# Using a copy of Alessandro's ray_tracing.pxd


cdef extern from "../liboptv/include/calibration.h":
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

    

cdef extern from "../liboptv/include/parameters.h":   

    ctypedef struct mm_np:
        int  	nlay 
        double  n1
        double  *n2 #[3]
        double  *d #[3]
        double  n3
        int     lut

    ctypedef struct control_par:
        mm_np *mm

    ctypedef struct volume_par:
        pass
        

cdef extern from "../liboptv/include/multimed.h":
    ctypedef struct Origin:
        double x, y, z

    ctypedef struct mmlut:
        Origin origin
        int    nr, nz, rw
        double *data


cdef class ImgCoord:
    cdef:        
        Calibration par_calibration
        
        int _owns_par_control_par
        control_par* par_control_par

        volume_par* par_volume_par
        int _owns_par_volume_par
        
        mmlut par_mmlut
        int _owns_mmlut
        
        mm_np mm

