//
//  openptv.c
//  OpenPTV
//
//  Created by Alex Liberzon on 11/05/2017.
//
//

#include "openptv.h"
void read_all_calibration(Calibration *calib[], int num_cams);

/* 
    Compile it using: 
       make 
    Run it using: 
       ./openptv test_cavity 10000 10004
*/
int main(int argc, char *argv[])
{
    

    int count;
    DIR *dirp;
    struct dirent *dp;
    char file_name[256];
    
    // 1. process inputs: directory, first frame, last frame
  
  printf ("This program was called with \"%s\".\n",argv[0]);
  
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
    
    // 2. read parameters and calibrations
    Calibration *calib[4];
    
    control_par *cpar = read_control_par("parameters/ptv.par");
    read_all_calibration(calib, cpar->num_cams);
    
    tracking_run *run = tr_new_legacy("parameters/sequence.par",
                                      "parameters/track.par", "parameters/criteria.par",
                                      "parameters/ptv.par", calib);
    
    if (argc == 4)
    {
        run->seq_par->first = atoi(argv[2]);
        run->seq_par->last = atoi(argv[3]);
    }
    printf("from frame %d to frame %d \n", run->seq_par->first, run->seq_par->last);
    
    target_par *targ_read = read_target_par("parameters/targ_rec.par");
    

    // 3. sequence (init, set images, loop)
    // to be completed
    
    // 4. tracking (init, loop, finish)
    track_forward_start(run);
    //trackcorr_c_loop(run, run->seq_par->first);
    
    int step;
    for (step = run->seq_par->first; step < run->seq_par->last; step++) {
        trackcorr_c_loop(run, step);
    }
    trackcorr_c_finish(run, run->seq_par->last);
    trackback_c(run, run->seq_par->last);
    return 0;
}

/* read_all_calibration function helps to read the num_cams calibrations from
   the ORI and ADDPAR files in the folder. 
   the name convention is cal/cam1.tif.ori and so on. 
*/
void read_all_calibration(Calibration *calib[], int num_cams) {
    char ori_tmpl[] = "cal/cam%d.tif.ori";
    char added_tmpl[] = "cal/cam%d.tif.addpar";
    char ori_name[256],added_name[256];
    int cam;
    
    for (cam = 0; cam < num_cams; cam++) {
        sprintf(ori_name, ori_tmpl, cam + 1);
        sprintf(added_name, added_tmpl, cam + 1);
        calib[cam] = read_calibration(ori_name, added_name, NULL);
    }
}

