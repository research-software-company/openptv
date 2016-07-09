/* Utilities for handling simple 3d vectors implemented as an array of 3
doubles.
*/

#ifndef CORRESPONDENCES_H
#define CORRESPONDENCES_H

#include "parameters.h"
#include "epi.h"

#define maxcand 200


typedef struct
{
  int     p[4];
  double  corr;
}
n_tupel;

typedef struct
{
  int    	p1;	       	/* point number of master point */
  int    	n;	       	/* # of candidates */
  int    	p2[maxcand];	/* point numbers of candidates */
  double	corr[maxcand];	/* feature based correlation coefficient */
  double	dist[maxcand];	/* distance perpendicular to epipolar line */
}
correspond;	       	/* correspondence candidates */

void correspondences (volume_par *vpar, control_par *cpar);

#endif

