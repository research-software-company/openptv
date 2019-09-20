/* main.c 
I need this file to start preparing some structure in my head. Alex
*/

#define MAXTARGETS 2048
#define BUFFER_LENGTH 4 // we do something very simple and studpid here
#include "main.h"

// These functions are part of the a test suite, see under /tests

void read_all_calibration(Calibration *calib[4], control_par *cpar)
{
    char ori_tmpl[] = "cal/cam%d.tif.ori";
    char added_name[] = "cal/cam1.tif.addpar";
    char ori_name[40];
    int cam;

    for (cam = 0; cam < cpar->num_cams; cam++)
    {
        sprintf(ori_name, ori_tmpl, cam + 1);
        calib[cam] = read_calibration(ori_name, added_name, NULL);
    }
}

/*  correct_frame() performs the transition from pixel to metric to flat 
    coordinates and x-sorting as required by the correspondence code.
    
    Arguments:
    frame *frm - target information for all cameras.
    control_par *cpar - parameters of image size, pixel size etc.
    tol - tolerance parameter for iterative flattening phase, see 
        trafo.h:correct_brown_affine_exact().
*/
coord_2d **correct_frame(frame *frm, Calibration *calib[], control_par *cpar,
                         double tol)
{
    coord_2d **corrected;
    int cam, part;

    corrected = (coord_2d **)malloc(cpar->num_cams * sizeof(coord_2d *));
    for (cam = 0; cam < cpar->num_cams; cam++)
    {
        corrected[cam] = (coord_2d *)malloc(
            frm->num_targets[cam] * sizeof(coord_2d));
        if (corrected[cam] == NULL)
        {
            /* roll back allocations and fail */
            for (cam -= 1; cam >= 0; cam--)
                free(corrected[cam]);
            free(corrected);
            return NULL;
        }

        for (part = 0; part < frm->num_targets[cam]; part++)
        {
            pixel_to_metric(&corrected[cam][part].x, &corrected[cam][part].y,
                            frm->targets[cam][part].x, frm->targets[cam][part].y,
                            cpar);

            dist_to_flat(corrected[cam][part].x, corrected[cam][part].y,
                         calib[cam], &corrected[cam][part].x, &corrected[cam][part].y,
                         tol);

            corrected[cam][part].pnr = frm->targets[cam][part].pnr;
        }

        /* This is expected by find_candidate() */
        quicksort_coord2d_x(corrected[cam], frm->num_targets[cam]);
    }
    return corrected;
}

int main(int argc, const char *argv[])
{
    // initialize variables

    int i, ntargets;
    // DIR *dirp;
    // struct dirent *dp;
    char file_name[256];
    int step, cam;
    unsigned char *img, *img_hp;
    target pix[MAXTARGETS], targ_t[MAXTARGETS];
    coord_2d **corrected;
    int match_counts[4];
    n_tupel *con;
    tracking_run *run;

    // read parameters from the working directory
    // for simplicity all names are default and hard coded (sorry)

    // 1. process inputs: directory, first frame, last frame

    printf("This program was called with \"%s\".\n", argv[0]);

    if (argc != 2 && argc != 4)
    {
        printf("Wrong number of inputs, expecting: \n");
        printf(" ./openptv test_cavity \n");
        printf(" or \n");
        printf(" ./openptv test_cavity 10000 10004 \n");
        return 0;
    }

    // change directory to the user-supplied working folder
    chdir(argv[1]);

    printf("changed directory to %s\n", argv[1]);

    // 2. read parameters and calibrations
    Calibration *calib[4]; // sorry only for 4 cameras now

    control_par *cpar = read_control_par("parameters/ptv.par");
    read_all_calibration(calib, cpar);
    free_control_par(cpar);
    printf("read calibrations\n");

    run = tr_new_legacy("parameters/sequence.par",
                        "parameters/track.par", "parameters/criteria.par",
                        "parameters/ptv.par", calib);

    if (argc == 4)
    {
        run->seq_par->first = atoi(argv[2]);
        run->seq_par->last = atoi(argv[3]);
    }
    printf("from frame %d to frame %d \n", run->seq_par->first, run->seq_par->last);

    // target_par *targ_read = read_target_par("parameters/targ_rec.par");

    // initialize memory buffers

    // for (step = 0; step < N_FRAMES_IN_DIRECTORY-BUFFER_LENGTH-1; step+BUFFER_LENGTH){
    // MAIN LOOP - see below we will just give inputs of 10000 10004 as a very simple approach

    // for each camera and for each time step the images are processed
    for (step = run->seq_par->first; step < run->seq_par->last + 1; step++)
    {
        for (cam = 0; cam < run->cpar->num_cams; cam++)
        {
            // we decided to focus just on the _targets, so we will read them from the
            // test directory test_cavity
            printf("reading targets from cam %d %s\n", cam, run->fb->target_file_base[cam]);
            printf("step is %d\n", step);

            /*
    target tbuf[2]; 
    target t1 = {0, 1127.0000, 796.0000, 13320, 111, 120, 828903, 1};
    target t2 = {1, 796.0000, 809.0000, 13108, 113, 116, 658928, 0};
    
    char *file_base = "testing_fodder/sample_";
    int frame_num = 42;
    int targets_read = 0;
    
    targets_read = read_targets(tbuf, file_base, frame_num);

    */

            run->fb->buf[step - run->seq_par->first]->num_targets[cam] = read_targets(
                run->fb->buf[step - run->seq_par->first]->targets[cam], run->fb->target_file_base[cam], step);
            printf("run->fb->buf[step]->num_targets[cam] is %d\n", run->fb->buf[step - run->seq_par->first]->num_targets[cam]);
            printf("and first x is %f\n", run->fb->buf[step - run->seq_par->first]->targets[cam][0].x);
        } // inner loop is per camera
        corrected = correct_frame(run->fb->buf[step - run->seq_par->first], calib, cpar, 0.0001);
        con = correspondences(run->fb->buf[step - run->seq_par->first], corrected, run->vpar, run->cpar, calib, match_counts);
        run->fb->buf[step - run->seq_par->first]->num_parts = match_counts[3]; // sum of all matches?
        printf("number of matched points is %d \n ", run->fb->buf[step - run->seq_par->first]->num_parts);
        // so here is missing frame into run->frame ?
        // WORK HERE

        // we need to push con into the run->fb somehow . 
        

    } // external loop is through frames

    // ok, theoretically we have now a buffer full of stuff from 4 frames
    // it's a good buffer on which we can just track stuff
    // and then we need to jump to a next chunk, remove all and start over.
    // the missing part is how to "chain the chunks" or make a smart use of
    // memory and buffers, it's beyond me now

    run->tpar->add = 0;
    track_forward_start(run);

    for (step = run->seq_par->first + 1; step < run->seq_par->last; step++)
    {
        trackcorr_c_loop(run, step);
    }
    trackcorr_c_finish(run, run->seq_par->last);

    // probably here we need to send something to plot it
    // or store in some memory for the next chunk?
    // basically we have a human in the loop here - his/her brain
    // will simply follow the velocity values in time like a movie
    // and we will store it to binary files. Later if someone wants to do
    // tracking, our simmple solution is not good enough. we kind of doing 3D-PIV here
    // of 4 frames and show the vectors. The quasi-vectors are not really connected. if we
    // will create nice animation - then the user will build trajectories him/herself.

    for (cam -= 1; cam >= 0; cam--) free(corrected[cam]);
    
    free(corrected);
    free(con);
    free(run->vpar);
    free_control_par(run->cpar);

    return 0;

} // should be end of main now