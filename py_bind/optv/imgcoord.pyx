from __future__ import division

# from libc.stdlib cimport malloc , free

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "../liboptv/include/imgcoord.h":
    void img_coord (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, 
    mmlut *mmLUT, double *x, double *y);
    void img_xy_mm_geo (double X, double Y, double Z, Calibration *cal, mm_np mm, 
    int i_cam, mmlut *mmLUT, double *x, double *y);



# @cython.boundscheck(False)
# @cython.wraparound(False)
def py_imgcoord(np.ndarray[np.float_t, ndim=1] pos1, np.ndarray[np.float_t, ndim=1] vec1,\
np.ndarray[np.float_t, ndim=1] pos2, np.ndarray[np.float_t, ndim=1] vec2):
    cdef np.ndarray[np.float_t, ndim=1] X = np.empty(3,dtype=np.float)
    intersect_rt(&pos1[0], &vec1[0], &pos2[0], &vec2[0], &X[0])
    print X
    return X