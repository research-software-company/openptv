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

//  if (argc > 1)
//    {
//      for (count = 1; count < argc; count++)
//	    {
//	        printf("User input: argv[%d] = %s \n", count, argv[count]);
//	    }
//    }
    
    // 2. init_proc - irrelevant?
    // 3. start_proc - irrelevant?
    // 4. read parameters

    // change directory to the user-supplied working folder
    chdir(argv[1]);
    // dirp = opendir(argv[1]);
    
    control_par *cpar = read_control_par("parameters/ptv.par");
//    printf(" ------------------ \n");
//    printf("Control parameters \n");
//    printf(" ------------------ \n");
//    printf(" highpass flag = %d \n", cpar->hp_flag);
//    printf(" use all cameras flag = %d \n", cpar->allCam_flag);
//    printf(" TIFF flag = %d \n", cpar->tiff_flag);
//    printf(" image is  %d x %d \n", cpar->imx, cpar->imy);
//    printf(" pixel size  %4.3f x %4.3f \n", cpar->pix_x, cpar->pix_y);
//    printf(" chfield flag = %d \n", cpar->chfield);
//    printf(" Multimedia parameters: \n");
//    printf(" gas/air index of refraction = %3.2f \n", cpar->mm->n1);
//    printf(" glass/perspex index of refraction = %3.2f \n", cpar->mm->n2[0]);
//    printf(" water/liquid index of refraction = %3.2f \n", cpar->mm->n3);
//    printf(" glass thickness = %3.2f \n",cpar->mm->d[0]);
    
    
    volume_par *vpar = read_volume_par("parameters/criteria.par");
    track_par *tpar = read_track_par("parameters/track.par");
    target_par *targ_read = read_target_par("parameters/targ_rec.par");
    sequence_par *seqp = read_sequence_par("parameters/sequence.par", cpar->num_cams);
    
    if (argc == 4)
    {
        seqp->first = atoi(argv[2]);
        seqp->last = atoi(argv[3]);
    }
    printf("from frame %d to frame %d \n", seqp->first, seqp->last);
    

    // 5. sequence (init, set images, loop)
    // 6. tracking (init, loop, finish)
    
    tracking_run *run;
    Calibration *calib[4];
    int step;
    
    read_all_calibration(calib, cpar->num_cams);
    
    run = tr_new_legacy("parameters/sequence.par",
                        "parameters/track.par", "parameters/criteria.par",
                        "parameters/ptv.par", calib);
    track_forward_start(run);
    trackcorr_c_loop(run, run->seq_par->first);
    for (step = run->seq_par->first + 1; step < run->seq_par->last; step++) {
        trackcorr_c_loop(run, step);
    }
    trackcorr_c_finish(run, run->seq_par->last);
    
    return 0;
}

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

