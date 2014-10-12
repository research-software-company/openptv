#include "epi.h"


 /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */



void epi_mm_2D ( double* xp           /*output coordinates*/
		, double* yp
		, double* zp
	        , double x1          /* input coordinates */
		, double y1          
		, Calibration* calib /*Ex1, I1, G1*/
		, mm_np* mm
		, volume_par* vpar
		)

{
 
  double a, b, c;
  double X1,Y1,Z1,X,Y,Z;
  
  double Zmin, Zmax;

  ray_tracing_v2 (x1, y1
		  , calib->ext_par
		  , calib->int_par
		  , calib->glass_par
		  , *mm
		  , &X1, &Y1, &Z1, &a, &b, &c);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = vpar->Zmin_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmin_lay[1] - vpar->Zmin_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);
  Zmax = vpar->Zmax_lay[0]
    + (X1 - vpar->X_lay[0]) * (vpar->Zmax_lay[1] - vpar->Zmax_lay[0]) /
    (vpar->X_lay[1] - vpar->X_lay[0]);

  Z = 0.5*(Zmin+Zmax);   
  X = X1 + (Z-Z1) * a/c;   
  Y = Y1 + (Z-Z1) * b/c;
  
  *xp=X; *yp=Y; *zp=Z;
  
  printf("x,y = %3.2f, %3.2f, X,Y,Z = %3.2f,%3.2f,%3.2f\n",x1,y1,X,Y,Z);
  printf("a,b,c = %3.2f, %3.2f, %3.2f, X1,Y1,Z1 = %3.2f, %3.2f, %3.2f\n", a,b,c, X1,Y1,Z1);
  

  
}
