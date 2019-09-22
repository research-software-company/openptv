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

// int main(int argc, const char *argv[])
int main()
{
    // initialize variables

    int i, cam, ntargets, step, geo_id;
    coord_2d **corrected;
    int match_counts[4];
    n_tupel *corresp_buf;
    tracking_run *run;
    vec3d res;
    vec2d targ;
    framebuf *frm; 

    corres t_corres = { 3, {96, 66, 26, 26} };
    P t_path = {
        .x = {45.219, -20.269, 25.946},
        .prev = -1,
        .next = -2,
        .prio = 4,
        .finaldecis = 1000000.0,
        .inlist = 0.
    };

    // read parameters from the working directory
    // for simplicity all names are default and hard coded (sorry)

    // 1. process inputs: directory, first frame, last frame

    // printf("This program was called with \"%s\".\n", argv[0]);

    // argc = 4;

    // argv[1] = '../tests/testing_fodder/test_cavity/';
    // argv[2] = 10001;
    // argv[3] = 10004;

    // if (argc != 2 && argc != 4)
    // {
    //     printf("Wrong number of inputs, expecting: \n");
    //     printf(" ./main ../tests/testing_fodder/test_cavity/ 10001 10004 \n");
    //     return 0;
    // }

    // argv[1] = '../tests/testing_fodder/test_cavity/';
    // argv[2] = 10001;
    // argv[3] = 10004;

    // change directory to the user-supplied working folder
    // chdir(argv[1]);

    chdir("../tests/testing_fodder/test_cavity/");

    // printf("changed directory to %s\n", argv[1]);

    // 2. read parameters and calibrations
    Calibration *calib[4]; // sorry only for 4 cameras now

    control_par *cpar = read_control_par("parameters/ptv.par");
    read_all_calibration(calib, cpar);
    free_control_par(cpar);
    printf("read calibrations\n");

    run = tr_new_legacy("parameters/sequence.par",
                        "parameters/track.par", "parameters/criteria.par",
                        "parameters/ptv.par", calib);


    frm = *(run->fb);

    // if (argc == 4)
    // {
    //     run->seq_par->first = atoi(argv[2]);
    //     run->seq_par->last = atoi(argv[3]);
    // }
    run->seq_par->first = 10001;
    run->seq_par->last = 10004;

    printf("from frame %d to frame %d \n", run->seq_par->first, run->seq_par->last);

    // target_par *targ_read = read_target_par("parameters/targ_rec.par");

    // initialize memory buffers

    // for (step = 0; step < N_FRAMES_IN_DIRECTORY-BUFFER_LENGTH-1; step+BUFFER_LENGTH){
    // MAIN LOOP - see below we will camust give inputs of 10000 10004 as a very simple approach

    // for each camera and for each time step the images are processed
    for (step = run->seq_par->first; step < run->seq_par->last + 1; step++)
    {
        for (cam = 0; cam < run->cpar->num_cams; cam++)
        {
            // we decided to focus camust on the _targets, so we will read them from the
            // test directory test_cavity

            printf("reading targets from %s%d\n", run->fb->target_file_base[cam], step);

            run->fb->buf[step - run->seq_par->first]->num_targets[cam] = read_targets(
                run->fb->buf[step - run->seq_par->first]->targets[cam], run->fb->target_file_base[cam], step);

            quicksort_target_y(run->fb->buf[step - run->seq_par->first]->targets[cam], run->fb->buf[step - run->seq_par->first]->num_targets[cam]);

            for (i = 0; i < run->fb->buf[step - run->seq_par->first]->num_targets[cam]; i++)
                run->fb->buf[step - run->seq_par->first]->targets[cam][i].pnr = i;

            // debugging purposes print the status of targets - see below another print.
            printf("%d targets and the first is %f,%f \n ",
                   run->fb->buf[step - run->seq_par->first]->num_targets[cam],
                   run->fb->buf[step - run->seq_par->first]->targets[cam][0].x,
                   run->fb->buf[step - run->seq_par->first]->targets[cam][0].y);

        } // inner loop is per camera
        corrected = correct_frame(run->fb->buf[step - run->seq_par->first], calib, cpar, 0.0001);
        corresp_buf = correspondences(run->fb->buf[step - run->seq_par->first], corrected, run->vpar, run->cpar, calib, match_counts);
        run->fb->buf[step - run->seq_par->first]->num_parts = match_counts[run->cpar->num_cams - 1];
        printf("number of matched points is %d \n ", run->fb->buf[step - run->seq_par->first]->num_parts);

        // first we need to create 3d points after correspondences and fill it into the buffer
        // use point_position and loop through the num_parts
        // probably feed it directly into the buffer

        // so we split into two parts:
        // first i copy the code from correspondences.pyx
        // and create the same types of arrays in C
        // then we will convert those to 3D using similar approach to what is now in Python


        // shortcut
        int num_parts = run->fb->buf[step - run->seq_par->first]->num_parts;
        int p[4];
        float x[4],y[4],skew_dist; 


        for (i=0; i<num_parts; i++) {
            for (cam=0; cam < run->cpar->num_cams; cam++) {
                if (corresp_buf[i].p[cam] >= 0)  
                    p[cam] = corrected[cam][corresp_buf[i].p[cam]].pnr;
                else
                    p[cam] = -1;
            }


            // printf("corrected %d %f %f \n ",corrected[cam][corresp_buf[i].p[cam]].pnr, corrected[cam][corresp_buf[i].p[cam]].x,corrected[cam][corresp_buf[i].p[cam]].y);
            if (p[cam] > -1){
                pixel_to_metric(&targ[0], &targ[1], \
                    run->fb->buf[step - run->seq_par->first]->targets[cam][p[cam]].x, \
                    run->fb->buf[step - run->seq_par->first]->targets[cam][p[cam]].y, \
                    run->cpar);
            } else {
                targ[0] = 1e-10;
                targ[1] = 1e-10;
            }
            skew_dist = point_position(&targ, run->cpar->num_cams, run->cpar->mm, calib, res);

            t_corres.nr = i;
            for (cam=0; cam < run->cpar->num_cams; cam++){
                t_corres.p[cam] = run->fb->buf[step - run->seq_par->first]->targets[cam][p[cam]].pnr;
            }
            run->fb->buf[step - run->seq_par->first]->correspond[i] = t_corres;

            t_path.x[0] = res[0];
            t_path.x[1] = res[1];
            t_path.x[2] = res[2];

             run->fb->buf[step - run->seq_par->first]->path_info[i] = t_path;
        }

    } // external loop is through frames

    run->tpar->add = 0;
    // we do not need to read frames - it's all in memory now
    // track_forward_start(run); 

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
    // will create nice animation - then the user will build tracamectories him/herself.

    for (cam -= 1; cam >= 0; cam--)
        free(corrected[cam]);

    free(corrected);
    free(corresp_buf);
    free(run->vpar);
    free_control_par(run->cpar);

    return 0;

} // should be end of main now