#ifndef TRACK_H
#define TRACK_H

#include "tracking_frame_buf.h"
#include "parameters.h"
#include "trafo.h"
#include "tracking_run.h"
#include "vec_utils.h"
#include "imgcoord.h"
#include "multimed.h"


typedef struct /* struct for what was found to corres */
{
 int ftnr, freq, whichcam[4];
}
foundpix;



#endif

