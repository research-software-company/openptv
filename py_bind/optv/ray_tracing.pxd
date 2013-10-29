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

    ctypedef struct ap_52:
        double k1,k2,k3,p1,p2,scx,she

    ctypedef struct Calibration:
        Exterior ext_par
        Interior int_par
        Glass glass_par
        ap_52 added_par

    

cdef extern from "optv/parameters.h":   

    ctypedef struct mm_np:
        int  	nlay 
        double  n1
        double  *n2 #[3]
        double  *d #[3]
        double  n3
        int     lut

    ctypedef struct control_par:
        mm_np *mm

cdef class Ray_tracing:
    cdef:        
        Calibration par_calibration
        int _owns_par_control_par
        control_par* par_control_par 

        





 

