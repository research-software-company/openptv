from cython cimport view
from cython.view cimport array as cvarray

from libc.stdlib cimport malloc #, free

import numpy as np

cdef extern from "optv/ray_tracing.h":
    void c_ray_tracing_v2  "ray_tracing_v2" (double x , double y , Exterior Ex , Interior I , Glass G , mm_np mm , double *Xb2 , double *Yb2 , double *Zb2 , double *a3 , double *b3 , double *c3)


cdef extern from "optv/calibration.h":
    int c_read_ori "read_ori" (Exterior *Ex, Interior *I, Glass *G, char *ori_file, ap_52 *addp, char *add_file, char *add_fallback)

cdef extern from "optv/parameters.h":
    control_par* c_read_control_par "read_control_par"(char *filename)
    void c_free_control_par "free_control_par"(control_par *cp)




class pExterior:
    def __init__(self):
        self.x0, self.y0, self.z0 = 0,0,0
        self.omega, self.phi, self.kappa = 0,0,0
        self.dm = np.zeros([3,3])

class pmm_np:
    def __init__(self):
        self.nlay = 0
        self.n1 = 0
        self.n2 = np.zeros(3)
        self.d = np.zeros(3)
        self.n3 = 0
        self.lut = 0

def get_dummy_mm_np():
    ret = pmm_np()
    ret.nlay = 3
    ret.n1 = 1.
    ret.n2 = np.array([1.49,0.,0])
    ret.d = np.array([5.,0.,0])
    ret.n3 = 1.33
    ret.lut = 1

    return ret

#cdef mm_np get_dummy_mm_np():
# def get_dummy_mm_np():
#     cdef:
#         mm_np ret
#         view.array n2view = <double[:3]> ret.n2 
#         view.array dview = <double[:3]> ret.d 

#     ret.nlay = 3
#     ret.n1 = 1.
#     n2view = np.array([1.49,0.,0])
#     dview = np.array([5.,0.,0])
#     # for i in range(3):
#     #     ret.n2[i] = n2view[i]
#     #     ret.d[i] = dview[i]
#     ret.n3 = 1.33
#     ret.lut = 1

#     return ret

    



def get_dummy_Interior():
    return {'xh':0. , 'yh':0. , 'cc':100.}


def get_dummy_Glass():
    return {'vec_x': 0.0001,'vec_y': 0.00001,'vec_z': 1.}


def get_dummy_Exterior():
    ret = pExterior()
    ret.z0 = 100
    ret.dm = np.array([[1.0, 0.2, -0.3], 
        [0.2, 1.0, 0.0],
        [-0.3, 0.0, 1.0]])

    return ret


        
# def ray_tracing(x, Ex, Interior I_, Glass G_,mm_np mm_):
#     cdef:
#        double x_, y_ , Xb2, Yb2, Zb2 , a3, b3, c3
#        Exterior Ex_
#     Ex_.x0,Ex_.y0,Ex_.z0 = Ex.x0,Ex.y0,Ex.z0
#     Ex_.omega,Ex_.phi,Ex_.kappa = Ex.omega,Ex.phi,Ex.kappa
#     for i in range(3):
#         for j in range(3):
#             Ex_.dm[i][j] = Ex.dm[i][j]
#     x_ = x[0]
#     y_ = x[1]   
#     c_ray_tracing_v2( x_, y_, Ex_, I_, G_, mm_, &Xb2 , &Yb2 , &Zb2 , &a3, &b3 , &c3 )
#     return (Xb2, Yb2, Zb2), (a3, b3, c3) #ALL OK!




def ray_tracing(x , Ex ,  I , G, mm):
    cdef:
       double x_, y_ , Xb2, Yb2, Zb2 , a3, b3, c3
       Exterior Ex_
       Interior I_ = I
       Glass G_ = G
       mm_np mm_

    Ex_.x0,Ex_.y0,Ex_.z0 = Ex.x0,Ex.y0,Ex.z0
    Ex_.omega,Ex_.phi,Ex_.kappa = Ex.omega,Ex.phi,Ex.kappa
#    cdef double [: , :] dm_view = Ex.dm

    for i in range(3):
        for j in range(3):
            Ex_.dm[i][j] = Ex.dm[i][j]

    
    mm_.nlay = mm.nlay
    mm_.n1 = mm.n1
    for i in range(3):
        mm_.n2[i] = mm.n2[i]
        mm_.d[i] = mm.d[i]
        
    mm_.n3 = mm.n3
    mm_.lut = mm.lut  
        

    x_ = x[0]
    y_ = x[1]   

    c_ray_tracing_v2( x_, y_, Ex_, I_, G_, mm_, &Xb2 , &Yb2 , &Zb2 , &a3, &b3 , &c3 )

    return (Xb2, Yb2, Zb2), (a3, b3, c3) #ALL OK!


cdef class Ray_tracing:
    
    def __init__(self):
        self._owns_par_control_par = 0

    def read_ori(self,ori_file):
        return c_read_ori( &self.par_calibration.ext_par
                           , &self.par_calibration.int_par
                           , &self.par_calibration.glass_par                           
                           , ori_file
                           , &self.par_calibration.added_par
                           , NULL
                           , NULL)

    def read_control_par(self,filename):
        cdef:
            control_par* ret 
        ret = c_read_control_par(filename)
        if ret != NULL:            
            self.par_control_par = ret
            self._owns_par_control_par = 1

    def trace(self,x):
        cdef:
            double x_ ,  y_ , Xb2 , Yb2 , Zb2 , a3, b3, c3
        
        x_ = x[0]
        y_ = x[1]   

        c_ray_tracing_v2( x_, y_
                          , self.par_calibration.ext_par
                          , self.par_calibration.int_par
                          , self.par_calibration.glass_par
                          , self.par_control_par[0].mm[0]
                          , &Xb2 , &Yb2 , &Zb2
                          , &a3, &b3 , &c3 )

        return (Xb2, Yb2, Zb2), (a3, b3, c3) #ALL OK!


    def force_set_exterior(self,Ex):
        self.par_calibration.ext_par.x0 \
        , self.par_calibration.ext_par.y0 \
        , self.par_calibration.ext_par.z0 = Ex.x0, Ex.y0, Ex.z0

        self.par_calibration.ext_par.omega \
        , self.par_calibration.ext_par.phi \
        , self.par_calibration.ext_par.kappa = Ex.omega,Ex.phi,Ex.kappa

        for i in range(3):
            for j in range(3):
                self.par_calibration.ext_par.dm[i][j] = Ex.dm[i][j]


    def force_set_interior(self,In):
        self.par_calibration.int_par = In

    def force_set_glass(self,Gl):
        self.par_calibration.glass_par = Gl

    def force_set_mm_np(self,mm):
        self.dealloc_par_control_par()
        self.par_control_par = <control_par *>malloc(sizeof(control_par))
        self.par_control_par[0].mm = <mm_np *> malloc(sizeof(mm_np))

        cdef:
            mm_np * mm_

        mm_ = self.par_control_par[0].mm
            
         
        mm_[0].nlay = mm.nlay
        mm_[0].n1 = mm.n1
        for i in range(3):
            mm_[0].n2[i] = mm.n2[i]
            mm_[0].d[i] = mm.d[i]
        
        mm_[0].n3 = mm.n3
        mm_[0].lut = mm.lut 

    def __dealloc__(self):
        self.dealloc_par_control_par()

    def dealloc_par_control_par(self):        
        if self._owns_par_control_par == 1:
            c_free_control_par(self.par_control_par)
            self._owns_par_control_par = 0
        
        

        
        

    #def set_parameters(exterior = None):
    #    pass




    
    
        
       
    

