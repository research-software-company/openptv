/*******************************************************************

Routine:	      	track.c

Author/Copyright:     	Jochen Willneff

Address:	       	Institute of Geodesy and Photogrammetry
		       	    ETH - Hoenggerberg
	               	CH - 8093 Zurich

Creation Date:		Beginning: February '98
                        End: far away

Description:   	        Tracking of particles in image- and objectspace

Routines contained:    	trackcorr_c

*******************************************************************/

/* References:
[1] http://en.wikipedia.org/wiki/Gradian
*/

#include "tracking_run.h"
#include "track.h"



/* Global variables marked extern in 'globals.h' and not defined elsewhere: */
int intx0_tr[4][10000], inty0_tr[4][10000], intx1_tr[4][10000],\
    inty1_tr[4][10000], intx2_tr[4][10000], inty2_tr[4][10000], \
    pnr1_tr[4][10000], pnr2_tr[4][10000], m1_tr;
int seq_step_shake;
double pnr3_tr[4][10000];
double npart, nlinks;


/* The buffer space required for this algorithm:

Note that MAX_TARGETS is taken from the global M, but I want a separate
definition because the fb created here should be local, not used outside
this file.

MAX_CANDS is the max number of candidates sought in search volume for next
link.
*/
#define BUFSPACE 4
#define MAX_TARGETS 20000
#define MAX_CANDS 4

/* trackcorr_c_init - initializes the tracking frame buffer with the parameters
   from sequence, track, criteria and ptv.par files
   Arguments:
   Calibration *cal - pointer to the calibration data of all the cameras

   returns: tracking_run type of a structure
*/

tracking_run* trackcorr_c_init(Calibration *cal) {
    int step;
    tracking_run *ret;

    /* Remaining globals:
    see below for communication globals.
    */

    ret = (tracking_run *) malloc(sizeof(tracking_run));
    tr_init(ret, "parameters/sequence.par", "parameters/track.par",
        "parameters/criteria.par", "parameters/ptv.par");

    fb_init(ret->fb, 4, ret->cpar->num_cams, MAX_TARGETS,
        "res/rt_is", "res/ptv_is", "res/added", ret->seq_par->img_base_name);

    /* Prime the buffer with first frames */
    for (step = ret->seq_par->first; step < ret->seq_par->first + 3; step++) {
        fb_read_frame_at_end(ret->fb, step, 0);
        fb_next(ret->fb);
    }
    fb_prev(ret->fb);

    ret->lmax = norm((ret->tpar->dvxmin - ret->tpar->dvxmax), \
        (ret->tpar->dvymin - ret->tpar->dvymax), \
        (ret->tpar->dvzmin - ret->tpar->dvzmax));
    volumedimension (&(ret->vpar->X_lay[1]), &(ret->vpar->X_lay[0]), &(ret->ymax),
        &(ret->ymin), &(ret->vpar->Zmax_lay[1]), &(ret->vpar->Zmin_lay[0]),
        ret->vpar, ret->cpar, cal);

    // Denis - globals below are used in trackcorr_finish
    npart=0;
    nlinks=0;

    return ret;
}

/* reset_foundpix_array() sets default values for foundpix objects in an array.
 *
 * Arguments:
 * foundpix *arr - the array to reset
 * int arr_len - array length
 * int num_cams - number of places in the whichcam member of foundpix.
 */
void reset_foundpix_array(foundpix *arr, int arr_len, int num_cams) {
    int i, cam;
    for (i = 0; i < arr_len; i++) {
	    arr[i].ftnr = -1;
        arr[i].freq = 0;

        for(cam = 0; cam < num_cams; cam++) {
            arr[i].whichcam[cam] = 0;
	    }
    }
}

/* copy_foundpix_array() copies foundpix objects from one array to another.
 *
 * Arguments:
 * foundpix *dest, *src - src is the array to copy, dest receives it.
 * int arr_len - array length
 * int num_cams - number of places in the whichcam member of foundpix.
 */
void copy_foundpix_array(foundpix *dest, foundpix *src, int arr_len,
    int num_cams)
{
    int i, cam;
    for (i = 0; i < arr_len; i++) {
        dest[i].ftnr = src[i].ftnr;
        dest[i].freq = src[i].freq;
        for (cam = 0; cam < num_cams; cam++) {
            dest[i].whichcam[cam] = src[i].whichcam[cam];
        }
    }
}

/* register_closest_neighbs() finds candidates for continuing a particle's
 * path in the search volume, and registers their data in a foundpix array
 * that is later used by the tracking algorithm.
 * TODO: the search area can be in a better data structure.
 *
 * Arguments:
 * target *targets - the targets list to search.
 * int num_targets - target array length.
 * int cam - the index of the camera we're working on.
 * double cent_x, cent_y - image coordinates of search area, [pixel]
 * double dl, dr, du, dd - respectively the left, right, up, down distance to
 *   the search area borders from its center, [pixel]
 * foundpix *reg - an array of foundpix objects, one for each possible
 *   neighbour. Output array.
 */
void register_closest_neighbs(target *targets, int num_targets, int cam,
    double cent_x, double cent_y, double dl, double dr, double du, double dd,
    foundpix *reg, control_par *cpar)
{
    int cand, all_cands[MAX_CANDS];

    cand = candsearch_in_pix (targets, num_targets, cent_x, cent_y, dl, dr,
        du, dd, all_cands, cpar);

    for (cand = 0; cand < MAX_CANDS; cand++) {
        if(all_cands[cand] == -999) {
            reg[cand].ftnr = -1;
        } else {
            reg[cand].whichcam[cam] = 1;
            reg[cand].ftnr = targets[all_cands[cand]].tnr;
        }
    }
}

/* search_volume_center_moving() finds the position of the center of the search
 * volume for a moving particle using the velocity of last step.
 *
 * Arguments:
 * vec3d prev_pos - previous position
 * vec3d curr_pos - current position
 * vec3d *output - output variable, for the  calculated
 *   position.
 */
void search_volume_center_moving(vec3d prev_pos, vec3d curr_pos, vec3d output)
{
    vec_scalar_mul(curr_pos, 2, output);
    vec_subt(output, prev_pos, output);

}

/* predict is used in display loop (only) of track.c to predict the position of
 * a particle in the next frame, using the previous and current positions
 * Arguments:
 * vec2d prev_pos, curr_pos are 2d positions at previous/current frames
 * vec2d output - output of the 2D positions of the particle in the next frame.
*/
void predict (vec2d prev_pos, vec2d curr_pos, vec2d output)
{
  output[0] = 2*curr_pos[0] - prev_pos[0];
  output[1] = 2*curr_pos[1] - prev_pos[1];
}


/* pos3d_in_bounds() checks that all components of a pos3d are in their
   respective bounds taken from a track_par object.

   Arguments:
   vec3d pos - the 3-component array to check.
   track_par *bounds - the struct containing the bounds specification.

   Returns:
   True if all components in bounds, false otherwise.
 */
int pos3d_in_bounds(vec3d pos, track_par *bounds) {
    return (
        bounds->dvxmin < pos[0] && pos[0] < bounds->dvxmax &&
        bounds->dvymin < pos[1] && pos[1] < bounds->dvymax &&
        bounds->dvzmin < pos[2] && pos[2] < bounds->dvzmax );
}

/* angle_acc() calculates the angle between the (1st order) numerical velocity
   vectors to the predicted next position and to the candidate actual position.
   The angle is calculated in [gon], see [1].
   The predicted position is the position if the particle continued at current
   velocity.

   Arguments:
   vec3d start, pred, cand - the particle start position, predicted position,
      and possible actual position, respectively.
   double *angle - output variable, the angle between the two velocity
      vectors, [gon]
   double *acc - output variable, the 1st-order numerical acceleration embodied
      in the deviation from prediction.
 */
void angle_acc(vec3d start, vec3d pred, vec3d cand, double *angle, double *acc)
{
    vec3d v0, v1;

    vec_subt(pred, start, v0);
    vec_subt(cand, start, v1);

    *acc = vec_diff_norm(v0, v1);


    if ((v0[0] == -v1[0]) && (v0[1] == -v1[1]) && (v0[2] == -v1[2])) {
        *angle = 200;
    } else if ((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2])) {
        *angle = 0; // otherwise it returns NaN
    } else {
        *angle = (200./M_PI) * acos(vec_dot(v0, v1) / vec_norm(v0) \
            / vec_norm(v1));
    }
}

void trackcorr_c_loop (tracking_run *run_info, int step, int display, Calibration **cal)
{
   /* sequence loop */
    char  buf[256];
    int j, h, k, mm, kk, invol=0;
    int counter1, counter2, philf[4][MAX_CANDS];
    int count1=0, count2=0, count3=0, zusatz=0;
    int intx0, intx1, inty0, inty1;
    int intx2, inty2;
    int quali=0;
    vec3d diff_pos, X[7]; /* 7 reference points used in the algorithm, TODO: check if can reuse some */
    double angle, acc, angle0, acc0,  dl;
    double xr[4], xl[4], yd[4], yu[4], angle1, acc1;
    vec2d p[4], c[4], n[4]; // would replace xp,yp and xc,yc and xn,yn by p[i][0],p[i][1]
    vec2d v1[4], v2[4];     // replaced x1[4], y1[4], x2[4], y2[4] by v1[j][0], v1[j][1]
    double rr;
    int flag_m_tr=0;
    int    	zoom_x[4], zoom_y[4], zoom_f[4];  /* zoom parameters */


    /* Shortcuts to inside current frame */
    P *curr_path_inf, *ref_path_inf;
    corres *curr_corres, *ref_corres;
    target **curr_targets, **ref_targets;
    int _ix; /* For use in any of the complex index expressions below */
    int _frame_parts; /* number of particles in a frame */

    /* Shortcuts into the tracking_run struct */
    framebuf *fb;
    track_par *tpar;
    volume_par *vpar;
    control_par *cpar;

    /* Remaining globals:
    all those in trackcorr_c_init.
    calibration globals.
    */

    foundpix *w, *wn, p16[4*MAX_CANDS];
    sprintf (buf, "Time step: %d, seqnr: %d, Particle info:",
        step - run_info->seq_par->first, step);
    count1=0; zusatz=0;

    fb = run_info->fb;
    tpar = run_info->tpar;
    vpar = run_info->vpar;
    cpar = run_info->cpar;
    curr_targets = fb->buf[1]->targets;

    /* try to track correspondences from previous 0 - corp, variable h */
    for (h = 0; h < fb->buf[1]->num_parts; h++) {
        for (j = 0; j < 7; j++) vec_init(X[j]);

        curr_path_inf = &(fb->buf[1]->path_info[h]);
        curr_corres = &(fb->buf[1]->correspond[h]);

	    curr_path_inf->inlist = 0;
        reset_foundpix_array(p16, 16, fb->num_cams);

	    /* 3D-position */
	    vec_copy(X[1], curr_path_inf->x);

	    /* use information from previous to locate new search position
	       and to calculate values for search area */
	    if (curr_path_inf->prev >= 0) {
            ref_path_inf = &(fb->buf[0]->path_info[curr_path_inf->prev]);
	        vec_copy(X[0], ref_path_inf->x);
            search_volume_center_moving(ref_path_inf->x, curr_path_inf->x, X[2]);

	        for (j = 0; j < fb->num_cams; j++) {
                img_coord(X[2], cal[j], cpar->mm, &n[j][0], &n[j][1]);
		        metric_to_pixel(&v1[j][0], &v1[j][1], n[j][0], n[j][1], cpar);
	        }
	    } else {
            vec_copy(X[2], X[1]);
	        for (j=0; j < fb->num_cams; j++) {
	            if (curr_corres->p[j] == -1) {
                    img_coord(X[2], cal[j], cpar->mm, &n[j][0], &n[j][1]);
                    metric_to_pixel(&v1[j][0], &v1[j][1], n[j][0], n[j][1], cpar);
	            } else {
                    _ix = curr_corres->p[j];
                    v1[j][0] = curr_targets[j][_ix].x;
                    v1[j][1] = curr_targets[j][_ix].y;
                }
            }
	    }

	    /* calculate searchquader and reprojection in image space */
	    searchquader(X[2], xr, xl, yd, yu, tpar, cpar, cal);

	    /* search in pix for candidates in next time step */
	    for (j = 0; j < fb->num_cams; j++) {
            register_closest_neighbs(fb->buf[2]->targets[j],
                fb->buf[2]->num_targets[j], j, v1[j][0], v1[j][1],
                xl[j], xr[j], yu[j], yd[j], &p16[j*MAX_CANDS], cpar);
	    }

	    /* fill and sort candidate struct */
	    sortwhatfound(p16, &counter1, fb->num_cams);
	    w = (foundpix *) calloc (counter1, sizeof (foundpix));

	    if (counter1 > 0) count2++;
        copy_foundpix_array(w, p16, counter1, fb->num_cams);
	    /*end of candidate struct */

	    /* check for what was found */
	    for (mm=0; mm<counter1;mm++) { /* counter1-loop */
	        /* search for found corr of current the corr in next
		    with predicted location */

            reset_foundpix_array(p16, 16, fb->num_cams);

	        /* found 3D-position */
            ref_path_inf = &(fb->buf[2]->path_info[w[mm].ftnr]);
            vec_copy(X[3], ref_path_inf->x);

	        if (curr_path_inf->prev >= 0) {
                for (j = 0; j < 3; j++)
                    X[5][j] = 0.5*(5.0*X[3][j] - 4.0*X[1][j] + X[0][j]);
	        } else {
                search_volume_center_moving(X[1], X[3], X[5]);
            }
            searchquader(X[5], xr, xl, yd, yu, tpar, cpar, cal);

	        for (j = 0; j < fb->num_cams; j++) {
                img_coord(X[5], cal[j], cpar->mm, &n[j][0], &n[j][1]);
		        metric_to_pixel(&v2[j][0], &v2[j][1], n[j][0], n[j][1], cpar);
	        }

	        /* search for candidates in next time step */
	        for (j=0; j < fb->num_cams; j++) {
	            counter2 = candsearch_in_pix (fb->buf[3]->targets[j],
                    fb->buf[3]->num_targets[j], v1[j][0], v1[j][1],
					xl[j], xr[j], yu[j], yd[j], philf[j], cpar);

		        for(k = 0; k < 4; k++) {
				    if (philf[j][k] == -999) {
                        p16[j*4+k].ftnr=-1;
				    } else {
				        if (fb->buf[3]->targets[j][philf[j][k]].tnr != -1) {
                            _ix = philf[j][k];
                            p16[j*4+k].ftnr = fb->buf[3]->targets[j][_ix].tnr;
                            p16[j*4+k].whichcam[j] = 1;
					    }
				    }
		        }
		    }
	        /* end of search in pix */

	        /* fill and sort candidate struct */
	        sortwhatfound(p16, &counter2, fb->num_cams);
	        wn = (foundpix *) calloc (counter2, sizeof (foundpix));
	        if (counter2 > 0) count3++;
            copy_foundpix_array(wn, p16, counter2, fb->num_cams);

	        /*end of candidate struct */
	        /* ************************************************ */
	        for (kk=0; kk < counter2; kk++)  { /* counter2-loop */
                ref_path_inf = &(fb->buf[3]->path_info[wn[kk].ftnr]);
                vec_copy(X[4], ref_path_inf->x);

                vec_subt(X[4], X[3], diff_pos);
                if ( pos3d_in_bounds(diff_pos, tpar)) {
                    angle_acc(X[3], X[4], X[5], &angle1, &acc1);

                    if (curr_path_inf->prev >= 0) {
                        angle_acc(X[1], X[2], X[3], &angle0, &acc0);
		            } else {
                        acc0=acc1; angle0=angle1;
                    }

                    acc=(acc0+acc1)/2; angle=(angle0+angle1)/2;
                    quali=wn[kk].freq+w[mm].freq;

                    if ((acc < tpar->dacc && angle < tpar->dangle) || \
                        (acc < tpar->dacc/10))
                    {
                        dl = (vec_diff_norm(X[1], X[3]) +
                            vec_diff_norm(X[4], X[3]) )/2;
                        rr = (dl/run_info->lmax + acc/tpar->dacc + angle/tpar->dangle)/(quali);
                        register_link_candidate(curr_path_inf, rr, w[mm].ftnr);
                    }
		        }
	        }   /* end of counter2-loop */

	        /* creating new particle position */
	        /* *************************************************************** */
	        for (j = 0;j < fb->num_cams; j++) {
                img_coord(X[5], cal[j], cpar->mm, &n[j][0], &n[j][1]);
		        metric_to_pixel(&n[j][0], &n[j][1], n[j][0], n[j][1], cpar);
	        }

	        /* reset img coord because of num_cams < 4 */
	        for (j=0;j < fb->num_cams; j++) { v2[j][0]=-1e10; v2[j][1]=-1e10; }

	        /* search for unused candidates in next time step */
	        for (j = 0;j < fb->num_cams; j++) {
		        /* use fix distance to define xl, xr, yu, yd instead of searchquader */
		        xl[j]= xr[j]= yu[j]= yd[j] = 3.0;

	            counter2 = candsearch_in_pix (fb->buf[3]->targets[j],
                    fb->buf[3]->num_targets[j], n[j][0], n[j][1],
					xl[j], xr[j], yu[j], yd[j], philf[j], cpar);
		        if(counter2>0 ) {
                    _ix = philf[j][0];
		            v2[j][0] = fb->buf[3]->targets[j][_ix].x;
                    v2[j][1] = fb->buf[3]->targets[j][_ix].y;
		        }
		    }
	        quali=0;

	        for (j = 0; j < fb->num_cams; j++) {
		        if (v2[j][0] != -1e10 && v2[j][1] != -1e10) {
		        pixel_to_metric(&v2[j][0],&v2[j][1], v2[j][0],v2[j][1], cpar);
                quali++;
		        }
		    }

	        if ( quali >= 2) {
                vec_copy(X[4], X[5]);
		        invol=0;

		        det_lsq_3d (*cal, *(cpar->mm), v2, \
                            &(X[4][0]), &(X[4][1]), &(X[4][2]), fb->num_cams);

		        /* volume check */
                if ( vpar->X_lay[0] < X[4][0] && X[4][0] < vpar->X_lay[1] &&
		            run_info->ymin < X[4][1] && X[4][1] < run_info->ymax &&
		            vpar->Zmin_lay[0] < X[4][2] && X[4][2] < vpar->Zmax_lay[1]) {invol=1;}

                vec_subt(X[3], X[4], diff_pos);
                if ( invol == 1 && pos3d_in_bounds(diff_pos, tpar) ) {
                    angle_acc(X[3], X[4], X[5], &angle, &acc);

                    if ((acc < tpar->dacc && angle < tpar->dangle) || \
                        (acc < tpar->dacc/10))
                    {
                        dl=(vec_diff_norm(X[1], X[3]) +
                            vec_diff_norm(X[4], X[3]) )/2;
                        rr = (dl/run_info->lmax + acc/tpar->dacc + angle/tpar->dangle) /
                            (quali+w[mm].freq);
                        register_link_candidate(curr_path_inf, rr, w[mm].ftnr);

                        if (tpar->add) {
                            ref_path_inf = &(fb->buf[3]->path_info[
                                fb->buf[3]->num_parts]);
                            vec_copy(ref_path_inf->x, X[4]);
                            reset_links(ref_path_inf);

                            _frame_parts = fb->buf[3]->num_parts;
                            ref_corres = &(fb->buf[3]->correspond[_frame_parts]);
                            ref_targets = fb->buf[3]->targets;
			                for (j = 0; j < fb->num_cams; j++) {
				                ref_corres->p[j]=-1;

				                if(philf[j][0]!=-999) {
                                    _ix = philf[j][0];
                                    ref_targets[j][_ix].tnr = _frame_parts;
                                    ref_corres->p[j] = _ix;
                                    ref_corres->nr = _frame_parts;
				                }
			                }
			                fb->buf[3]->num_parts++;
                            zusatz++;
                        }
                    }
                }
		        //invol=0;
	        }
	        //quali=0;

	        /* end of creating new particle position */
	        /* *************************************************************** */

	        /* try to link if kk is not found/good enough and prev exist */
	        if ( curr_path_inf->inlist == 0 && curr_path_inf->prev >= 0 ) {
                vec_subt(X[3], X[1], diff_pos);

                if (pos3d_in_bounds(diff_pos, tpar)) {
                    angle_acc(X[1], X[2], X[3], &angle, &acc);

                    if ( (acc < tpar->dacc && angle < tpar->dangle) || \
                        (acc < tpar->dacc/10) )
                    {
                        quali = w[mm].freq;
                        dl = (vec_diff_norm(X[1], X[3]) +
                            vec_diff_norm(X[0], X[1]) )/2;
                        rr = (dl/run_info->lmax + acc/tpar->dacc + angle/tpar->dangle)/(quali);
                        register_link_candidate(curr_path_inf, rr, w[mm].ftnr);
			        }
		        }
	        }

	        free(wn);
        } /* end of counter1-loop */

	    /* begin of inlist still zero */
	    if (tpar->add) {
	        if ( curr_path_inf->inlist == 0 && curr_path_inf->prev >= 0 ) {
                for (j = 0; j < fb->num_cams; j++) {
                    img_coord(X[2], cal[j], cpar->mm, &n[j][0], &n[j][1]);
		            metric_to_pixel(&n[j][0], &n[j][1], n[j][0], n[j][1], cpar);
		            v2[j][0]=-1e10;
                    v2[j][1]=-1e10;
                }

    		    /* search for unused candidates in next time step */
		        for (j = 0; j < fb->num_cams; j++) {
		            /*use fix distance to define xl, xr, yu, yd instead of searchquader */
		            xl[j]= xr[j]= yu[j]= yd[j] = 3.0;
	                counter2 = candsearch_in_pix (fb->buf[2]->targets[j],
                        fb->buf[2]->num_targets[j], n[j][0], n[j][1],
					    xl[j], xr[j], yu[j], yd[j], philf[j], cpar);
		            if(counter2 > 0) {
                        _ix = philf[j][0];
	    	            v2[j][0] = fb->buf[2]->targets[j][_ix].x;
                        v2[j][1] = fb->buf[2]->targets[j][_ix].y;
		            }
		        }
		        quali=0;

		        for (j = 0; j < fb->num_cams; j++) {
		            if (v2[j][0] !=-1e10 && v2[j][1] != -1e10) {
		                pixel_to_metric(&v2[j][0], &v2[j][1], v2[j][0], v2[j][1], cpar);
                        quali++;
		            }
		        }

		        if (quali>=2) {
                    vec_copy(X[3], X[2]);
		            invol=0;

	    	        det_lsq_3d (*cal, *(cpar->mm),v2, \
                        &(X[3][0]), &(X[3][1]), &(X[3][2]), \
                        fb->num_cams);

		            /* in volume check */
		            if ( vpar->X_lay[0] < X[3][0] && X[3][0] < vpar->X_lay[1] &&
                        run_info->ymin < X[3][1] && X[3][1] < run_info->ymax &&
                        vpar->Zmin_lay[0] < X[3][2] &&
                        X[3][2] < vpar->Zmax_lay[1]) {invol = 1;}

                    vec_subt(X[2], X[3], diff_pos);
                    if ( invol == 1 && pos3d_in_bounds(diff_pos, tpar) ) {
                        angle_acc(X[1], X[2], X[3], &angle, &acc);

                        if ( (acc < tpar->dacc && angle < tpar->dangle) || \
                            (acc < tpar->dacc/10) )
                        {
                            dl = (vec_diff_norm(X[1], X[3]) +
                                vec_diff_norm(X[0], X[1]) )/2;
                            rr = (dl/run_info->lmax + acc/tpar->dacc + angle/tpar->dangle)/(quali);

                            ref_path_inf = &(fb->buf[2]->path_info[
                                fb->buf[2]->num_parts]);
                            vec_copy(ref_path_inf->x, X[3]);
                            reset_links(ref_path_inf);

                            _frame_parts = fb->buf[2]->num_parts;
                            register_link_candidate(curr_path_inf, rr, _frame_parts);

                            ref_corres = &(fb->buf[2]->correspond[_frame_parts]);
                            ref_targets = fb->buf[2]->targets;
			                for (j = 0;j < fb->num_cams; j++) {
                                ref_corres->p[j]=-1;

                                if(philf[j][0]!=-999) {
                                    _ix = philf[j][0];
                                    ref_targets[j][_ix].tnr = _frame_parts;
                                    ref_corres->p[j] = _ix;
                                    ref_corres->nr = _frame_parts;
                                }
                            }
                            fb->buf[2]->num_parts++;
                            zusatz++;
                        }
                    }
		            //invol=0;
		        } // if quali >= 2
            }
        }
	    /* end of inlist still zero */
	    /***********************************/

	    free(w);
	} /* end of h-loop */

    /* sort decis and give preliminary "finaldecis"  */
    for (h = 0; h < fb->buf[1]->num_parts; h++) {
        curr_path_inf = &(fb->buf[1]->path_info[h]);

	    if(curr_path_inf->inlist > 0 ) {
	        sort(curr_path_inf->inlist, (float *) curr_path_inf->decis,
                curr_path_inf->linkdecis);
      	    curr_path_inf->finaldecis = curr_path_inf->decis[0];
	        curr_path_inf->next = curr_path_inf->linkdecis[0];
	    }
	}

    /* create links with decision check */
    for (h = 0;h < fb->buf[1]->num_parts; h++) {
        curr_path_inf = &(fb->buf[1]->path_info[h]);

	    if(curr_path_inf->inlist > 0 ) {
            ref_path_inf = &(fb->buf[2]->path_info[curr_path_inf->next]);

	        if (ref_path_inf->prev == -1) {
	            /* best choice wasn't used yet, so link is created */
                ref_path_inf->prev = h;
            } else {
	            /* best choice was already used by mega[2][mega[1][h].next].prev */
	            /* check which is the better choice */
	            if ( fb->buf[1]->path_info[ref_path_inf->prev].finaldecis > \
                    curr_path_inf->finaldecis)
                {
		            /* remove link with prev */
		            fb->buf[1]->path_info[ref_path_inf->prev].next = NEXT_NONE;
                    ref_path_inf->prev = h;
		        } else {
		            curr_path_inf->next = NEXT_NONE;
	            }
	        }
        }
        if (curr_path_inf->next != -2 ) count1++;
    }
    /* end of creation of links with decision check */
    /* ******** Draw links now ******** */
    m1_tr = 0;

    if (display) {
        for (h = 0; h < fb->buf[1]->num_parts; h++) {
            curr_path_inf = &(fb->buf[1]->path_info[h]);
            curr_corres = &(fb->buf[1]->correspond[h]);
            ref_corres = &(fb->buf[2]->correspond[curr_path_inf->next]);

            if (curr_path_inf->next != -2 ) {

                // sprintf(buf ,"green"); all the display is in Python

                for (j = 0; j < fb->num_cams; j++) {
                    if (curr_corres->p[j] > 0 && ref_corres->p[j] > 0) {
                        flag_m_tr=1;
                        p[j][0] = curr_targets[j][curr_corres->p[j]].x;
               		    p[j][1] = curr_targets[j][curr_corres->p[j]].y;
               		    c[j][0] = fb->buf[2]->targets[j][ref_corres->p[j]].x;
               		    c[j][1] = fb->buf[2]->targets[j][ref_corres->p[j]].y;

               		    predict (p[j], c[j], n[j]);

               		    if ( ( fabs(p[j][0]-zoom_x[j]) < cpar->imx/(2*zoom_f[j]))
                            && (fabs(p[j][1]-zoom_y[j]) < cpar->imy/(2*zoom_f[j])))
                        {
                           	// sprintf(val ,"orange"); all the display is in Python

                        	intx0 = (int)(cpar->imx/2 + zoom_f[j]*(p[j][0] - zoom_x[j]));
	                        inty0 = (int)(cpar->imy/2 + zoom_f[j]*(p[j][1] - zoom_y[j]));
                            intx1 = (int)(cpar->imx/2 + zoom_f[j]*(c[j][0] - zoom_x[j]));
	                        inty1 = (int)(cpar->imy/2 + zoom_f[j]*(c[j][1] - zoom_y[j]));
	                        intx2 = (int)(cpar->imx/2 + zoom_f[j]*(n[j][0] - zoom_x[j]));
	                        inty2 = (int)(cpar->imy/2 + zoom_f[j]*(n[j][1] - zoom_y[j]));

	                        intx0_tr[j][m1_tr]=intx0;
	                        inty0_tr[j][m1_tr]=inty0;
	                        intx1_tr[j][m1_tr]=intx1;
	                        inty1_tr[j][m1_tr]=inty1;
	                        intx2_tr[j][m1_tr]=intx2;
	                        inty2_tr[j][m1_tr]=inty2;
	                        pnr1_tr[j][m1_tr]=-1;
	                        pnr2_tr[j][m1_tr]=-1;
	                        pnr3_tr[j][m1_tr]=-1;

	                        if (curr_path_inf->finaldecis > 0.2) {
	            	            pnr1_tr[j][m1_tr] = h;
		                        pnr2_tr[j][m1_tr] = curr_path_inf->next;
		                        pnr3_tr[j][m1_tr] = curr_path_inf->finaldecis;
		                    }
                        }
                    }

                    if (flag_m_tr==0)  {
                        intx0_tr[j][m1_tr]=0;
    	                inty0_tr[j][m1_tr]=0;
	                    intx1_tr[j][m1_tr]=0;
	                    inty1_tr[j][m1_tr]=0;
	                    intx2_tr[j][m1_tr]=0;
	                    inty2_tr[j][m1_tr]=0;
	                    pnr1_tr[j][m1_tr]=-1;
	                    pnr2_tr[j][m1_tr]=-1;
	                    pnr3_tr[j][m1_tr]=-1;
                    }
                    flag_m_tr=0;
                }
                m1_tr++;
            }
        }
    }
    /* ******** End of Draw links now ******** */
    sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
        step, fb->buf[1]->num_parts, fb->buf[2]->num_parts, count1,
        fb->buf[1]->num_parts - count1, zusatz);

    /* for the average of particles and links */
    npart = npart + fb->buf[1]->num_parts;
    nlinks = nlinks + count1;

    fb_next(fb);
    fb_write_frame_from_start(fb, step);
    if(step < run_info->seq_par->last - 2) {
        fb_read_frame_at_end(fb, step + 3, 0);
    }
} /* end of sequence loop */

void trackcorr_c_finish(tracking_run *run_info, int step, int display)
{
  int range = run_info->seq_par->last - run_info->seq_par->first;

  /* average of all steps */
  npart /= range;
  nlinks /= range;
  printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
	  npart, nlinks, npart-nlinks);

  fb_next(run_info->fb);
  fb_write_frame_from_start(run_info->fb, step);

  fb_free(run_info->fb);

  /* reset of display flag */
  //display = 1;
}

/*     track backwards */
void trackback_c (tracking_run *run_info, int step, int display, Calibration **cal)
{
    char  buf[256];
    int i, j, h, k, invol=0;
    int counter1, philf[4][MAX_CANDS];
    int count1=0, count2=0, zusatz=0;
    int quali=0;
    double  angle, acc, lmax, dl;
    double xr[4], xl[4], yd[4], yu[4];
    vec3d diff_pos, X[7]; /* 7 reference points used in the algorithm, TODO: check if can reuse some */
    // double xn[4], yn[4];
    vec2d n[4], v2[4]; // replaces xn,yn, x2[4], y2[4],
    double rr, Ymin=0, Ymax=0;
    double npart=0, nlinks=0;
    foundpix *w, p16[4*MAX_CANDS];

    sequence_par *seq_par;
    track_par *tpar;
    volume_par *vpar;
    control_par *cpar;
    framebuf *fb;

    /* Shortcuts to inside current frame */
    P *curr_path_inf, *ref_path_inf;
    corres *ref_corres;
    target **ref_targets;
    int _ix; /* For use in any of the complex index expressions below */
    int _frame_parts; /* number of particles in a frame */

    //display = 1;
    /* read data */
    seq_par = read_sequence_par("parameters/sequence.par", 4);
    tpar = read_track_par("parameters/track.par");
    vpar = read_volume_par("parameters/criteria.par");
    cpar = read_control_par("parameters/ptv.par");

    fb = (framebuf *) malloc(sizeof(framebuf));
    fb_init(fb, 4, cpar->num_cams, MAX_TARGETS,
        "res/rt_is", "res/ptv_is", "res/added", seq_par->img_base_name);

    /* Prime the buffer with first frames */
    for (step = seq_par->last; step > seq_par->last - 4; step--) {
        fb_read_frame_at_end(fb, step, 1);
        fb_next(fb);
    }
    fb_prev(fb);

    lmax = norm((tpar->dvxmin - tpar->dvxmax), (tpar->dvymin - tpar->dvymax),
	    (tpar->dvzmin - tpar->dvzmax));
    volumedimension (&(vpar->X_lay[1]), &(vpar->X_lay[0]), &Ymax,
        &Ymin, &(vpar->Zmax_lay[1]), &(vpar->Zmin_lay[0]), vpar, cpar, *cal);

    /* sequence loop */
    for (step = seq_par->last - 1; step > seq_par->first; step--) {
        sprintf (buf, "Time step: %d, seqnr: %d, Particle info:",
            step - seq_par->first, step);

        for (h = 0; h < fb->buf[1]->num_parts; h++) {
            curr_path_inf = &(fb->buf[1]->path_info[h]);

            /* We try to find link only if the forward search failed to. */
            if ((curr_path_inf->next < 0) || (curr_path_inf->prev != -1)) continue;

            for (j = 0; j < 7; j++) vec_init(X[j]);
            curr_path_inf->inlist = 0;
            reset_foundpix_array(p16, 16, fb->num_cams);

            /* 3D-position of current particle */
            vec_copy(X[1], curr_path_inf->x);

            /* use information from previous to locate new search position
            and to calculate values for search area */
            ref_path_inf = &(fb->buf[0]->path_info[curr_path_inf->next]);
	        vec_copy(X[0], ref_path_inf->x);
            search_volume_center_moving(ref_path_inf->x, curr_path_inf->x, X[2]);

            for (j=0; j < fb->num_cams; j++) {
                img_coord(X[2], cal[j], cpar->mm, &n[j][0], &n[j][1]);
                metric_to_pixel(&n[j][0], &n[j][1], n[j][0], n[j][1], cpar);
            }

            /* calculate searchquader and reprojection in image space */
            searchquader(X[2], xr, xl, yd, yu, tpar, cpar, cal);


            for (j = 0; j < fb->num_cams; j++) {
                counter1 = candsearch_in_pix (
                    fb->buf[2]->targets[j], fb->buf[2]->num_targets[j], n[j][0], n[j][1],
                    xl[j], xr[j], yu[j], yd[j], philf[j], cpar);

                for(k=0; k<4; k++) {
                    if( counter1>0) {
                        if (philf[j][k] == -999){
                            p16[j*4+k].ftnr=-1;
                        } else {
                            p16[j*4+k].ftnr=fb->buf[3]->targets[j][philf[j][k]].tnr;
                            p16[j*4+k].whichcam[j]=1;
                        }
                    }
                }
            }
            /* search in pix for candidates in next time step */
            //for (j = 0; j < fb->num_cams; j++) {
            //    register_closest_neighbs(fb->buf[2]->targets[j],
            //        fb->buf[2]->num_targets[j], j, n[j][0], n[j][1],
            //        xl[j], xr[j], yu[j], yd[j], &p16[j*MAX_CANDS]);
            //}

            /* fill and sort candidate struct */
            sortwhatfound(p16, &counter1, fb->num_cams);
            w = (foundpix *) calloc (counter1, sizeof (foundpix));

            /*end of candidate struct */
            if (counter1 > 0) count2++;
            copy_foundpix_array(w, p16, counter1, fb->num_cams);

            for (i = 0; i < counter1; i++) {
                ref_path_inf = &(fb->buf[2]->path_info[w[i].ftnr]);
                vec_copy(X[3], ref_path_inf->x);

                vec_subt(X[1], X[3], diff_pos);
                if (pos3d_in_bounds(diff_pos, tpar)) {
                    angle_acc(X[1], X[2], X[3], &angle, &acc);

                    /* *********************check link *****************************/
                    if ((acc < tpar->dacc && angle < tpar->dangle) || \
                        (acc < tpar->dacc/10))
                    {
                        dl = (vec_diff_norm(X[1], X[3]) +
                            vec_diff_norm(X[0], X[1]) )/2;
                        quali=w[i].freq;
                        rr = (dl/lmax + acc/tpar->dacc + angle/tpar->dangle)/quali;
                        register_link_candidate(curr_path_inf, rr, w[i].ftnr);
                    }
                }
            }

            free(w);
            /******************/
            quali=0;

            /* reset img coord because num_cams < 4 */
            for (j=0;j<4;j++) { v2[j][0]=-1e10; v2[j][1]=-1e10;}

            /* if old wasn't found try to create new particle position from rest */
            if (tpar->add) {
                if ( curr_path_inf->inlist == 0) {
                    for (j = 0; j < fb->num_cams; j++) {
                        /* use fix distance to define xl, xr, yu, yd instead of searchquader */
                        xl[j]= xr[j]= yu[j]= yd[j] = 3.0;

                        counter1 = candsearch_in_pix (fb->buf[2]->targets[j],
                            fb->buf[2]->num_targets[j], n[j][0], n[j][1],
                            xl[j], xr[j], yu[j], yd[j], philf[j], cpar);
                        if(counter1 > 0) {
                            _ix = philf[j][0];
                            v2[j][0] = fb->buf[2]->targets[j][_ix].x;
                            v2[j][1] = fb->buf[2]->targets[j][_ix].y;
                        }
                    }

                    for (j = 0; j < fb->num_cams; j++) {
                        if (v2[j][0] !=-1e10 && v2[j][1] != -1e10) {
                            pixel_to_metric(&v2[j][0], &v2[j][1], v2[j][0],v2[j][1], cpar);
                            quali++;
                        }
                    }

                    if (quali>=2) {
                        vec_copy(X[3], X[2]);
                        invol=0;

                        det_lsq_3d (*cal, *(cpar->mm), v2,
                            &(X[3][0]), &(X[3][1]), &(X[3][2]), \
                            fb->num_cams);

                        /* volume check */
                        if ( vpar->X_lay[0] < X[3][0] && X[3][0] < vpar->X_lay[1] &&
                            Ymin < X[3][1] && X[3][1] < Ymax &&
                            vpar->Zmin_lay[0] < X[3][2] && X[3][2] < vpar->Zmax_lay[1])
                                {invol = 1;}

                        vec_subt(X[1], X[3], diff_pos);
                        if (invol == 1 && pos3d_in_bounds(diff_pos, tpar)) {
                            angle_acc(X[1], X[2], X[3], &angle, &acc);

                            if ( (acc<tpar->dacc && angle<tpar->dangle) || \
                                (acc<tpar->dacc/10) )
                            {
                                dl = (vec_diff_norm(X[1], X[3]) +
                                    vec_diff_norm(X[0], X[1]) )/2;
                                rr =(dl/lmax+acc/tpar->dacc + angle/tpar->dangle)/(quali);

                                ref_path_inf = &(fb->buf[2]->path_info[
                                    fb->buf[2]->num_parts]);
                                vec_copy(ref_path_inf->x, X[3]);
                                reset_links(ref_path_inf);

                                _frame_parts = fb->buf[2]->num_parts;
                                register_link_candidate(curr_path_inf, rr,
                                    _frame_parts);

                                ref_corres = &(fb->buf[2]->correspond[_frame_parts]);
                                ref_targets = fb->buf[2]->targets;

                                for (j = 0;j < fb->num_cams; j++) {
                                    ref_corres->p[j]=-1;

                                    if(philf[j][0]!=-999) {
                                        _ix = philf[j][0];
                                        ref_targets[j][_ix].tnr = _frame_parts;
                                        ref_corres->p[j] = _ix;
                                        ref_corres->nr = _frame_parts;
                                    }
                                }
                                fb->buf[2]->num_parts++;
                            }
                        }
                        //invol=0;
                    }
                }
            } /* end of if old wasn't found try to create new particle position from rest */
        } /* end of h-loop */

        for (h = 0; h < fb->buf[1]->num_parts; h++) {
            curr_path_inf = &(fb->buf[1]->path_info[h]);

            if(curr_path_inf->inlist > 0 ) {
                sort(curr_path_inf->inlist, (float *)curr_path_inf->decis,
                    curr_path_inf->linkdecis);
            }
        }

        /* create links with decision check */
        count1=0; zusatz=0;
        for (h = 0; h < fb->buf[1]->num_parts; h++) {
            curr_path_inf = &(fb->buf[1]->path_info[h]);

            if (curr_path_inf->inlist > 0 ) {
                /* if old/new and unused prev == -1 and next == -2 link is created */
                ref_path_inf = &(fb->buf[2]->path_info[curr_path_inf->linkdecis[0]]);

                if ( ref_path_inf->prev == PREV_NONE && \
                    ref_path_inf->next == NEXT_NONE )
                {
                    curr_path_inf->finaldecis = curr_path_inf->decis[0];
                    curr_path_inf->prev = curr_path_inf->linkdecis[0];
                    fb->buf[2]->path_info[curr_path_inf->prev].next = h;
                    zusatz++;
                }

                /* old which link to prev has to be checked */
                if ((ref_path_inf->prev != PREV_NONE) && \
                    (ref_path_inf->next == NEXT_NONE) )
                {
                    vec_copy(X[0], fb->buf[0]->path_info[curr_path_inf->next].x);
                    vec_copy(X[1], curr_path_inf->x);
                    vec_copy(X[3], ref_path_inf->x);
                    vec_copy(X[4], fb->buf[3]->path_info[ref_path_inf->prev].x);
                    for (j = 0; j < 3; j++)
                        X[5][j] = 0.5*(5.0*X[3][j] - 4.0*X[1][j] + X[0][j]);

                    angle_acc(X[3], X[4], X[5], &angle, &acc);

                    if ( (acc<tpar->dacc && angle<tpar->dangle) ||  (acc<tpar->dacc/10) ) {
                        curr_path_inf->finaldecis = curr_path_inf->decis[0];
                        curr_path_inf->prev = curr_path_inf->linkdecis[0];
                        fb->buf[2]->path_info[curr_path_inf->prev].next = h;
                        zusatz++;
                    }
                }
            }

            if (curr_path_inf->prev != -1 ) count1++;
        } /* end of creation of links with decision check */

        sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
            step, fb->buf[1]->num_parts, fb->buf[2]->num_parts, count1,
            fb->buf[1]->num_parts - count1, zusatz);

        /* for the average of particles and links */
        npart = npart + fb->buf[1]->num_parts;
        nlinks = nlinks + count1;

        fb_next(fb);
        fb_write_frame_from_start(fb, step);
        if(step > seq_par->first + 2) { fb_read_frame_at_end(fb, step - 3, 1); }
    } /* end of sequence loop */

    /* average of all steps */
    npart /= (seq_par->last - seq_par->first - 1);
    nlinks /= (seq_par->last - seq_par->first - 1);

    printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
    npart, nlinks, npart-nlinks);

    fb_next(fb);
    fb_write_frame_from_start(fb, step);

    fb_free(fb);
    free(fb);
    free(tpar);

    /* reset of display flag */
    //display = 1;
}


/* candsearch_in_pix searches of four (4) near candidates in target list
 * Arguments:
 * target next[] array of targets (pointer, x,y, n, nx,ny, sumg, track ID),
 * has to be y sorted ?!! this is not tested in the function.
 * int num_targets - target array length.
 * double cent_x, cent_y - image coordinates of the position of a particle [pixel]
 * double dl, dr, du, dd - respectively the left, right, up, down distance to
 *   the search area borders from its center, [pixel]
 * int p[] array of integer pointers
 * control_par *cpar array of parameters (cpar->imx,imy are needed)
 * Returns the integer counter of the number of candidates, between 0 - 3
*/

int candsearch_in_pix (target next[], int num_targets, double cent_x, double cent_y,
double dl, double dr, double du, double dd, int p[4], control_par *cpar) {

  int  	  j, j0, dj;
  int  counter=0, p1, p2, p3, p4;
  double  d, dmin=1e20, xmin, xmax, ymin, ymax;
  double d1, d2, d3, d4;

  xmin = cent_x - dl;  xmax = cent_x + dr;  ymin = cent_y - du;  ymax = cent_y + dd;

  if(xmin<0.0) xmin=0.0;
  if(xmax > cpar->imx)
        xmax = cpar->imx;
  if(ymin<0.0) ymin=0.0;
  if(ymax > cpar->imy)
    ymax = cpar->imy;

  p1 = p2 = p3 = p4 = -999;
  d1 = d2 = d3 = d4 = dmin;

  if (cent_x >= 0.0 && cent_x <= cpar->imx ) {
        if (cent_y >= 0.0 && cent_y <= cpar->imy ) {

      /* binarized search for start point of candidate search */
      for (j0=num_targets/2, dj=num_targets/4; dj>1; dj/=2)
        {
          if (next[j0].y < ymin) j0 += dj;
          else j0 -= dj;
        }

      j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
      for (j=j0; j<num_targets; j++) {	       	        /* candidate search */
        if (next[j].tnr != -1 ) {
            if (next[j].y > ymax )  break;	    /* finish search */

            if (next[j].x > xmin && next[j].x < xmax \
            && next[j].y > ymin && next[j].y < ymax){
                d = sqrt ((cent_x-next[j].x)*(cent_x-next[j].x) + \
                          (cent_y-next[j].y)*(cent_y-next[j].y));

                if (d < dmin) {
                   dmin = d;
                }

                if ( d < d1 ) {
                   p4=p3; p3=p2; p2=p1; p1=j;
                   d4=d3; d3=d2; d2=d1; d1=d;
                  }
                else if ( d1 < d &&  d < d2 ){
                   p4=p3; p3=p2; p2=j;
                   d4=d3; d3=d2; d2=d;
                }
                else if ( d2 < d && d < d3 ){
                   p4=p3; p3=j;
                   d4=d3; d3=d;
                }
                else if ( d3 < d && d < d4 ){
                   p4=j;
                   d4=d;
                }
              }
        }
      }

      p[0]=p1;
      p[1]=p2;
      p[2]=p3;
      p[3]=p4;

      for (j=0; j<4; j++) if ( p[j] != -999 ) { counter++; }
      } /* if x is within the image boundaries */
    }   /* if y is within the image boundaries */
  return (counter);
}




/* searchquader defines the search region, using tracking parameters
 * dvxmin, ... dvzmax (but within the image boundaries), per camera
 * Arguments:
 * vec3d point position in physical space
 * track_par *tpar set of tracking parameters
 * control_par *cpar set of control parameters for the num_cams
 * Calibration *cal calibration per camera to find a projection of a 3D vertex
 * of a cuboid in the image space.
 * Returns the arrays xr,xl,yd,yu (right, left, down, up) per camera
 * for the search of a quader (cuboid).
*/

void searchquader(vec3d point, double xr[4], double xl[4], double yd[4], \
    double yu[4], track_par *tpar, control_par *cpar, Calibration **cal){
  int i, pt, dim;
  vec3d mins, maxes;
  double x, y, xz, yz;
  vec3d quader[8];



  vec_set(mins, tpar->dvxmin, tpar->dvymin, tpar->dvzmin);
  vec_set(maxes, tpar->dvxmax, tpar->dvymax, tpar->dvzmax);
  /* 3D positions of search volume - eight corners of a box */
  for (pt = 0; pt < 8; pt++) {
    vec_copy(quader[pt], point);
    for (dim = 0; dim < 3; dim++) {
        if (pt & 1<<dim) {
            quader[pt][dim] += maxes[dim];
        } else {
            quader[pt][dim] += mins[dim];
        }
        // printf("quader[%d][%d]=%f \n", pt,dim,quader[pt][dim]);
    }
  }



  /* calculation of search area in each camera */
  for (i = 0; i < cpar->num_cams; i++) {
      xr[i]=0;
      xl[i] = cpar->imx;
      yd[i]=0;
      yu[i] = cpar->imy;

      img_coord (point, cal[i], cpar->mm, &xz, &yz);
      metric_to_pixel (&xz, &yz, xz, yz, cpar);

      for (pt = 0; pt < 8; pt++) {
        img_coord (quader[pt], cal[i], cpar->mm, &x, &y);
	    metric_to_pixel (&x, &y, x, y, cpar);

	    if (x <xl[i] ) xl[i]=x;
	    if (y <yu[i] ) yu[i]=y;
	    if (x >xr[i] ) xr[i]=x;
	    if (y >yd[i] ) yd[i]=y;
	  }
      if (xl[i] < 0 ) xl[i]=0;
      if (yu[i] < 0 ) yu[i]=0;
      if (xr[i] > cpar->imx)
        xr[i] = cpar->imx;
      if (yd[i] > cpar->imy)
        yd[i] = cpar->imy;



      xr[i]=xr[i]-xz;
      xl[i]=xz-xl[i];
      yd[i]=yd[i]-yz;
      yu[i]=yz-yu[i];
  }
}



void sortwhatfound (foundpix item[16], int *counter, int num_cams)
{
  int i,j,m, different;
  foundpix temp;

  different=0;

  /* where what was found */
  for (i=0; i<16; i++)
    for (j=0; j<4; j++)
      for (m=0; m<4; m++)
	if(item[i].ftnr == item[4*j+m].ftnr)
	  {
	    item[i].whichcam[j]=1;
	  }

  /* how often was ftnr found */
  for (i=0; i<16; i++)
    for (j=0; j < num_cams; j++)
      if (item[i].whichcam[j] == 1 && item[i].ftnr !=-1) item[i].freq++;

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }

  for (i=0; i<16; i++)
    for (j=i+1; j<16; j++)
      {
	if (item[i].ftnr == item[j].ftnr || item[j].freq <2)
	  {
	    item[j].freq=0;
	    item[j].ftnr=-1;
	  }
      }

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }
  for (i=0; i<16; ++i) if(item[i].freq != 0) different++;
  *counter=different;

}

/* sorts a float array a and an integer array b both of length n
 * Arguments:
 *  float array a (returned sorted in the ascending order)
 *  integer array b (returned sorted according to float array a)
 *  int n (length of a)
*/
void sort(int n, float a[], int b[]){
  int flag = 0, i, itemp;
  float ftemp;

  do {
    flag = 0;
    for(i=0; i<(n-1); i++)
      if(a[i] > a[i+1]) {
	ftemp =  a[i];
	itemp =  b[i];
	a[i] = a[i+1];
	b[i] = b[i+1];
	a[i+1] = ftemp;
	b[i+1] = itemp;
        flag = 1;
      }
  }while(flag);
}

void det_lsq_3d (Calibration *cals, mm_np mm, vec2d v[], double *Xp, double *Yp, double *Zp, int num_cams) {
    int     i,count_inner=0,n,m, flag[4] = {0., 0., 0., 0.};
    double  d_inner=0.,x,y;
    double X[4][3], a[4][3];
    double dist,X_pos[6],Y_pos[6],Z_pos[6],XX,YY,ZZ,si0,sqX,sqY,sqZ;
    double rmsX, rmsY, rmsZ, mean_sigma0;
    vec3d res;
    int cam;
    
    rmsX = rmsY = rmsZ = mean_sigma0 = 0.0;



	    for (cam = 0; cam < num_cams; cam++){
          if(v[cam][0] > -999){
			 flag[cam]=1;
             dist_to_flat(v[cam][0], v[cam][1], &(cals[cam]),
                 &x, &y, 100000);
		     ray_tracing(x,y, &(cals[cam]), mm, X[cam], a[cam]);

		  }
		}


		count_inner=0;
		for (n = 0; n < num_cams; n++){
			for(m = n+1; m < num_cams; m++){
				if(flag[n]==1 && flag[m]==1){
                    dist = skew_midpoint(X[n], a[n], X[m], a[m], res);
                    d_inner += dist;
					X_pos[count_inner]=res[0];
                    Y_pos[count_inner]=res[1];
                    Z_pos[count_inner]=res[2];
					count_inner++;
				}
			}
		}
        //d_inner/=(double)count_inner;
		XX=0.;YY=0.;ZZ=0.;
		for(i=0;i<count_inner;i++){
           XX+=X_pos[i];
		   YY+=Y_pos[i];
		   ZZ+=Z_pos[i];
		}
		XX/=(double)count_inner;YY/=(double)count_inner;ZZ/=(double)count_inner;
		//end of new det_lsq
		*Xp=XX;
		*Yp=YY;
		*Zp=ZZ;

		//statistics
		si0=0.;sqX=0.;sqY=0.;sqZ=0.;
		for(i=0;i<count_inner;i++){
           si0+=pow(X_pos[i]-XX,2.)+pow(Y_pos[i]-YY,2.)+pow(Z_pos[i]-ZZ,2.);
           sqX+=pow(X_pos[i]-XX,2.);
		   sqY+=pow(Y_pos[i]-YY,2.);
		   sqZ+=pow(Z_pos[i]-ZZ,2.);
		}
		si0/=(double)count_inner;sqX/=(double)count_inner;sqY/=(double)count_inner;sqZ/=(double)count_inner;

		mean_sigma0 += pow(si0,0.5);
        rmsX += pow(sqX,0.5);
        rmsY += pow(sqY,0.5);
        rmsZ += pow(sqZ,0.5);
		//end of statistics

}
