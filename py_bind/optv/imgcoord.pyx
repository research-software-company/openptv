from __future__ import division

# from libc.stdlib cimport malloc , free

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "../liboptv/include/imgcoord.h":
    void img_coord (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y)
    void img_xy_mm_geo (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y)



# @cython.boundscheck(False)
# @cython.wraparound(False)
def py_imgcoord(np.ndarray[np.float_t, ndim=1] pos1, np.ndarray[np.float_t, ndim=1] vec1, np.ndarray[np.float_t, ndim=1] pos2, np.ndarray[np.float_t, ndim=1] vec2):
    cdef np.ndarray[np.float_t, ndim=1] X = np.empty(3,dtype=np.float)
    intersect_rt(&pos1[0], &vec1[0], &pos2[0], &vec2[0], &X[0])
    print X
    return X
    
    
def py_img_coord(X, cal, mm, i_cam, mmlut, x, y):
    cdef:
       double x_, y_ , z_, x, y
       Calibration* cal_
       mm_np mm_
       int i_cam_
       mmlut* mmlut_

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
        

    x_ = X[0]
    y_ = X[1]
    z_ = X[2]   

    img_coord( x_, y_, z_, cal_, mm_, i_cam_, mmlut_, &x, &y )

    return (x, y) 