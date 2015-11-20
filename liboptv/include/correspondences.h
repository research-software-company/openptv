/* Utilities for handling sortgrid */

#ifndef CORRESPONDENCES_H
#define CORRESPONDENCES_H

#include "calibration.h"
#include "parameters.h"
#include "tracking_frame_buf.h"
#include "imgcoord.h"
#include "multimed.h"   
#include "trafo.h"
#include "vec_utils.h"
#include "sortgrid.h"
#include "epi.h"
#include "btree.h"
#include <stdio.h>

#define maxcand 200 
#define nmax 20240

typedef struct
{
  int    	p1;	       	    /* point number of master point */
  int    	n;	       	    /* # of candidates */
  int    	p2[maxcand];	/* point numbers of candidates */
  double	corr[maxcand];	/* feature based correlation coefficient */
  double	dist[maxcand];	/* distance perpendicular to epipolar line */
}
correspond;	       	        /* correspondence candidates */

typedef struct
{
  int     p[4];             /* 4 is probably the number of cameras, to be fixed to N */
  double  corr;
}
n_tupel;


int advance_nd_iterator(int nd, int *current, int *count);
n_tupel **correspondences(frame *frm, Calibration **calib, volume_par *vpar);

void quicksort_coord2d_x (coord_2d	*crd, int num);
void qs_coord2d_x (coord_2d	*crd, int left, int right);


#endif

