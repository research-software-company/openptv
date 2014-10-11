#ifndef IMGCOORD_H
#define IMGCOORD_H

#include <math.h>
#include "parameters.h"
#include "calibration.h"
#include "tracking_frame_buf.h"
#include "lsqadj.h"
#include "ray_tracing.h"
#include "trafo.h"
#include "multimed.h"
#include "epi.h"

void img_coord (double X, double Y, double Z, Calibration *cal, mm_np mm, int i_cam, mmlut *mmLUT, double *x, double *y);
void img_xy    (double X, double Y, double Z, Calibration *cal, double *x, double *y);
void img_xy_mm_geo (double X, double Y, double Z, Calibration *cal, mm_np *mm, int i_cam, mmlut *mmLUT, double *x, double *y);



#endif