#ifndef IMGCOORD_H
#define IMGCOORD_H

#include <math.h>
#include "parameters.h"
#include "calibration.h"
#include "tracking_frame_buf.h"
#include "lsqadj.h"
#include "ray_tracing.h"
#include "trafo.h"

void img_coord (double X, double Y, double Z, Exterior Ex, Interior I, Glas G, \
                ap_52 ap, mm_np mm, double *x, double *y);
void img_xy (double X, double Y, double Z, Exterior Ex, Interior I, Glass G, \
             double *x, double *y);
void img_xy_mm_geo (double X, double Y, double Z, Exterior Ex, Interior I, Glass G, \
mm_np mm, double *x, double *y);


#endif

