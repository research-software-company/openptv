
#ifndef EPI_H
#define EPI_H

#include "ray_tracing.h"


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
		);



#endif
