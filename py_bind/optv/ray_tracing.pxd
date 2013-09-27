cdef extern from "optv/calibration.h":
    ctypedef double Dmatrix[3][3]

    ctypedef struct Exterior:
        double  x0, y0, z0
        double  omega, phi, kappa
        Dmatrix dm

    ctypedef struct Interior:
        double xh, yh
        double cc

    ctypedef struct Glass:
        double vec_x,vec_y,vec_z



cdef extern from "optv/parameters.h":
    ctypedef struct mm_np:
        int  	nlay 
        double  n1
        double  n2[3]
        double  d[3]
        double  n3
        int     lut




 

