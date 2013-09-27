from cython.view cimport array as cvarray
import numpy as np

cdef extern from "optv/ray_tracing.h":
    void c_ray_tracing_v2  "ray_tracing_v2" (double x , double y , Exterior Ex , Interior I , Glass G , mm_np mm , double *Xb2 , double *Yb2 , double *Zb2 , double *a3 , double *b3 , double *c3)

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

def get_dummy_Interior():
    return {'xh':0. , 'yh':0. , 'cc':100}


def get_dummy_Glass():
    return {'vec_x': 0.0001,'vec_y': 0.00001,'vec_z': 1.}


def get_dummy_Exterior():
    ret = pExterior()
    ret.zo = 100
    ret.dm = np.array([[1.0, 0.2, -0.3], 
        [0.2, 1.0, 0.0],
        [-0.3, 0.0, 1.0]])

    return ret

def get_dummy_mm_np():
    ret = pmm_np()
    ret.nlay = 3
    ret.n1 = 1.
    ret.d2 = np.array([1.49,0.,0])
    ret.d = np.array([5.,0.,0])
    ret.n3 = 1.33
    ret.lut = 1

    return ret

        
def ray_tracing(x , Ex , Interior I , Glass G, mm):
    cdef:
       double x_, y_ , Xb2, Yb2, Zb2 , a3, b3, c3
       Exterior Ex_
       #Interior I_
       #Glass G_
       mm_np mm_

    Ex_.x0,Ex_.y0,Ex_.z0 = Ex.x0,Ex.y0,Ex.z0
    Ex_.omega,Ex_.phi,Ex_.kappa = Ex.omega,Ex.phi,Ex.kappa
    cdef double [: , :] dm_view = Ex.dm

    for i in range(3):
        for j in range(3):
            Ex_.dm[i][j] = Ex.dm[i][j]

    print Ex_.omega, Ex_.z0
    print Ex_.dm[0][2]

    mm_.nlay = mm.nlay
    mm_.n1 = mm.n1
    for i in range(3):
        mm_.n2[i] = mm.n2[i]
        mm_.d[i] = mm.d[i]
        
    mm_.n3 = mm.n3
    mm_.lut = mm.lut  
    
        

    x_ = x[0]
    y_ = x[1]

    print y_


    c_ray_tracing_v2( x_, y_, Ex_, I, G, mm_, &Xb2 , &Yb2 , &Zb2 , &a3, &b3 , &c3 )

    print Xb2


    return (Xb2, Yb2, Zb2), (a3, b3, c3)


    
       
    

