/* Definitions for tracking routines. */

#ifndef TRACK_H
#define TRACK_H

#include "tracking_frame_buf.h"
#include "parameters.h"
#include "trafo.h"
#include "tracking_run.h"
#include "vec_utils.h"
#include "imgcoord.h"
#include "multimed.h"
#include "orientation.h"


typedef struct /* struct for what was found to corres */
{
 int ftnr, freq, whichcam[4];
}
foundpix;

int candsearch_in_pix(target  next[], int num_targets, double x, double y,
    double dl, double dr, double du, double dd, int p[4], control_par *cpar);
int candsearch_in_pixrest(target  next[], int num_targets, double x, double y,
    double dl, double dr, double du, double dd, int p[4], control_par *cpar);
void sortwhatfound (foundpix item[16], int *zaehler, int num_cams);
void searchquader(vec3d point, double xr[4], double xl[4], double yd[4], \
    double yu[4], track_par *tpar, control_par *cpar, Calibration *cal);
void predict(vec2d a, vec2d b, vec2d c);
void det_lsq_3d (Calibration *cals, mm_np mm, vec2d v[], double *Xp, \
    double *Yp, double *Zp, int num_cams);
void sort(int n, float a[], int b[]);

tracking_run* trackcorr_c_init();
void trackcorr_c_loop (tracking_run *run_info, int step, int display, \
    Calibration *cal);
void trackcorr_c_finish(tracking_run *run_info, int step, int display);
void trackback_c(tracking_run *run_info, int step, int display, \
    Calibration *cal);

#endif
