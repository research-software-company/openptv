from cython cimport view
from cython.view cimport array as cvarray

from libc.stdlib cimport malloc , free

import numpy as np


# cdef extern from "optv/calibration.h":
#     int c_read_ori "read_ori" (Exterior *Ex, Interior *I, Glass *G, char *ori_file, ap_52 *addp, char *add_file, char *add_fallback)
#     Calibration* c_read_calibration "read_calibration" (char* ori_file , char* add_file , char* fallback_file)
# 
# cdef extern from "optv/parameters.h":
#     control_par* c_read_control_par "read_control_par"(char *filename)
#     void c_free_control_par "free_control_par"(control_par *cp)
#     volume_par* c_read_volume_par "read_volume_par" (char* filename)
	

cdef extern from "../liboptv/include/imgcoord.h":
    void c_img_coord "img_coord" (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y)
    void c_img_xy_mm_geo "img_xy_mm_geo" (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y)


        

cdef class ImgCoord:    
    def __init__(self):
        self._owns_par_control_par = 0
        self._owns_par_volume_par = 0

#     def read_ori(self,ori_file,add_file):
#         return c_read_ori( &self.par_calibration.ext_par
#                            , &self.par_calibration.int_par
#                            , &self.par_calibration.glass_par                           
#                            , ori_file
#                            , &self.par_calibration.added_par
#                            , add_file
#                            , NULL)
# 
#     def read_control_par(self,filename):
#         cdef:
#             control_par* ret 
#         ret = c_read_control_par(filename)
#         if ret != NULL:            
#             self.par_control_par = ret
#             self._owns_par_control_par = 1
#         else:
#             print "An error while reading", filename, "occurred."
# 
#     def read_volume_par(self,filename):
#         cdef:
#             volume_par* ret
#         ret = c_read_volume_par(filename)
#         if ret != NULL:
#             self.par_volume_par = ret
#             self._owns_par_volume_par = 1
#         else:
#             print "An error while reading", filename, "occurred."

       
    def img_coord(self, x, int_cam):
        cdef:
            double x_, y_, z_, xp, yp
            int int_cam_
            
        x_ = x[0]
        y_ = x[1]
        z_ = x[2]
        
        int_cam_ = int_cam
        
        c_img_coord(x_ 
                , y_
                , z_
                , &self.par_calibration
                , self.mm
                , int_cam_
                , &self.par_mmlut
                , &xp
                , &yp)
                
        return [xp, yp]
                
    
    def img_xy_mm_geo(self,x,int_cam):
        cdef:
            double x_, y_, z_, xp, yp
            int int_cam_
            
        x_ = x[0]
        y_ = x[1]
        z_ = x[2]
        
        int_cam_ = int_cam
        
        c_img_xy_mm_geo(x_ 
                , y_
                , z_
                , &self.par_calibration
                , self.mm
                , int_cam_
                , &self.par_mmlut
                , &xp
                , &yp)
                
        return [xp, yp]
        

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


    def force_set_interior(self,Interior In):
        self.par_calibration.int_par = In

    def force_set_glass(self,Glass Gl):
        self.par_calibration.glass_par = Gl
        
    def force_set_addpar(self, ap_52 addpar):
        self.par_calibration.added_par = addpar
    
    def force_set_calibration(self, Calibration cal):
        self.par_calibration.ext_par = cal.ext_par
        self.par_calibration.int_par = cal.int_par
        self.par_calibration.glass_par = cal.glass_par
        self.par_calibration.added_par = cal.added_par
     
#    def force_set_mmlut(self, mmlut mmLUT):
#        self.par_mmlut = mmLUT
        
    def force_set_only_mm_np(self,mm):
        self.par_mm = mm    
    

    def force_set_mm_np(self,mm):
        self.dealloc_par_control_par()
        self.par_control_par = <control_par *> malloc(sizeof(control_par))
        self._owns_par_control_par = 1 
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

#     def __dealloc__(self):
#         self.dealloc_par_control_par()
# 
#     def dealloc_par_control_par(self):        
#         if self._owns_par_control_par == 1:
#             c_free_control_par(self.par_control_par)
#             self._owns_par_control_par = 0
        

    #def set_parameters(exterior = None):
    #    pass


################## Classes used to force values of calibration or parameters instead of loading them from file

class pExterior:
    def __init__(self):
        self.x0, self.y0, self.z0 = 0,0,0
        self.omega, self.phi, self.kappa = 0,0,0
        self.dm = np.zeros([3,3])

class pAddpar:
    def __init__(self):
        self.k1, self.k2, self.k3 = 0,0,0
        self.p1, self.p2 = 0,0
        self.scx = 1	
        self.she = 0
        
class pOrigin:
    def __init__(self):
        self.x, self.y, self.z = 0,0,0

class pInterior:
    def __init__(self):
        self.xh, self.yh, self.cc = 0,0,0

class pGlass:
    def __init__(self):
        self.vec_x, self.vec_y, self.vec_z = 0,0,0
        
class pCalibration:
    def __init__(self):
        self.ext_par = self.pExterior()
        self.int_par = self.pInterior()
        self.glass_par = self.pGlass()
        self.added_par = self.pAddpar()
        

class pmm_np:
    def __init__(self):
        self.nlay = 0
        self.n1 = 0
        self.n2 = np.zeros(3)
        self.d = np.zeros(3)
        self.n3 = 0
        self.lut = 0
        
class pmmlut:
    def __init__(self):
        self.nr, self.nz, self.nw = 2,2,2
        self.origin = pOrigin()
        self.origin.x, self.origin.y, self.origin.z = 0,0,0
        self.data = np.zeros([self.nr*self.nz])
