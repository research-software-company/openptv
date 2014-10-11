/****************************************************************************

Routine:	       	imgcoord.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
                        ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	22.4.88

Description:	       	computes x', y' from given Point and orientation
	       		(see: Kraus)

Routines contained:

****************************************************************************/
/*test-commit*/
#include "imgcoord.h"


/* Calculates projection from coordinates in
 * world space to pixelcoordinates in image space
 */
void img_coord (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y){
    double deno, r, dx, dy;
	Calibration *cal_t;
    Exterior Ex_t;
    double X_t,Y_t,Z_t,cross_p[3],cross_c[3];
	
    /*calculate tilted positions and copy them to X_t, Y_t and Z_t*/
	trans_Cam_Point(cal->ext_par,*mm,cal->glass_par,X,Y,Z,&Ex_t,&X_t,&Y_t,&Z_t,&cross_p,&cross_c);
    
    multimed_nlay (&Ex_t, mm, X_t,Y_t,Z_t, &X_t,&Y_t, i_cam, mmLUT);
    /*
     void  multimed_nlay (Exterior *ex
                   , mm_np *mm
                   , double X
                   , double Y
                   , double Z
                   , double *Xq
                   , double *Yq
                   , int i_cam
                   , mmlut *mmLUT);
    */
    back_trans_Point(X_t,Y_t,Z_t,*mm, cal->glass_par,cross_p,cross_c,&X,&Y,&Z);
	
    X -= cal->ext_par.x0;  Y -= cal->ext_par.y0;  Z -= cal->ext_par.z0;

    deno = cal->ext_par.dm[0][2] * X + cal->ext_par.dm[1][2] * Y + cal->ext_par.dm[2][2] * Z;
    *x = cal->int_par.xh - cal->int_par.cc * (cal->ext_par.dm[0][0]*X + cal->ext_par.dm[1][0]*Y + cal->ext_par.dm[2][0]*Z) / deno;
    *y = cal->int_par.yh - cal->int_par.cc * (cal->ext_par.dm[0][1]*X + cal->ext_par.dm[1][1]*Y + cal->ext_par.dm[2][1]*Z) / deno;

    r = sqrt (*x * *x + *y * *y);
	
    dx = (*x) * (cal->added_par.k1*r*r + cal->added_par.k2*r*r*r*r + cal->added_par.k3*r*r*r*r*r*r)
      + cal->added_par.p1 * (r*r + 2*(*x)*(*x)) + 2*cal->added_par.p2*(*x)*(*y);
    dy = (*y) * (cal->added_par.k1*r*r + cal->added_par.k2*r*r*r*r + cal->added_par.k3*r*r*r*r*r*r)
      + cal->added_par.p2 * (r*r + 2*(*y)*(*y)) + 2*cal->added_par.p1*(*x)*(*y);

    *x += dx;
    *y += dy;

    *x = cal->added_par.scx * (*x) - sin(cal->added_par.she) * (*y);
    *y = cos(cal->added_par.she) * (*y);

}

/* Todo: Add short description
 * 
 */
void img_xy (double X, double Y, double Z, Calibration *cal, double *x, double *y)
{
  double deno;

  deno = cal->ext_par.dm[0][2] * (X-cal->ext_par.x0)
    + cal->ext_par.dm[1][2] * (Y-cal->ext_par.y0)
    + cal->ext_par.dm[2][2] * (Z-cal->ext_par.z0);

  *x = cal->int_par.xh - cal->int_par.cc *  (cal->ext_par.dm[0][0] * (X-cal->ext_par.x0)
		       + cal->ext_par.dm[1][0] * (Y-cal->ext_par.y0)
		       + cal->ext_par.dm[2][0] * (Z-cal->ext_par.z0)) / deno;

  *y = cal->int_par.yh - cal->int_par.cc *  (cal->ext_par.dm[0][1] * (X-cal->ext_par.x0)
		       + cal->ext_par.dm[1][1] * (Y-cal->ext_par.y0)
		       + cal->ext_par.dm[2][1] * (Z-cal->ext_par.z0)) / deno;
}

/* Todo: Add short description
 * 
 */
void img_xy_mm_geo (double X, double Y, double Z, Calibration *cal, mm_np *mm, int i_cam, \
mmlut *mmLUT, double *x, double *y)
{

  double deno;
  Exterior Ex_t;
  double X_t,Y_t,Z_t,cross_p[3],cross_c[3],Xh,Yh,Zh;

  /*calculate tilted positions and copy them to X_t, Y_t and Z_t*/
  trans_Cam_Point(cal->ext_par,*mm,cal->glass_par,X,Y,Z,&Ex_t,&X_t,&Y_t,&Z_t,&cross_p,&cross_c);
  
  multimed_nlay (&Ex_t, mm, X_t,Y_t,Z_t, &X_t,&Y_t, i_cam, mmLUT);

  back_trans_Point(X_t,Y_t,Z_t,*mm, cal->glass_par,cross_p,cross_c,&X,&Y,&Z);

  deno = cal->ext_par.dm[0][2] * (X-cal->ext_par.x0)
    + cal->ext_par.dm[1][2] * (Y-cal->ext_par.y0)
    + cal->ext_par.dm[2][2] * (Z-cal->ext_par.z0);

  *x = - cal->int_par.cc *  (cal->ext_par.dm[0][0] * (X-cal->ext_par.x0)
		  + cal->ext_par.dm[1][0] * (Y-cal->ext_par.y0)
		  + cal->ext_par.dm[2][0] * (Z-cal->ext_par.z0)) / deno;

  *y = - cal->int_par.cc *  (cal->ext_par.dm[0][1] * (X-cal->ext_par.x0)
		  + cal->ext_par.dm[1][1] * (Y-cal->ext_par.y0)
		  + cal->ext_par.dm[2][1] * (Z-cal->ext_par.z0)) / deno;
}

