//
//  openptv.c
//  OpenPTV
//
//  Created by Alex Liberzon on 11/05/2017.
//
//

#include "openptv.h"

/* 
    Compile it using: 
       make 
    Run it using: 
       ./openptv test_cavity 10000 10004
*/
int main(int argc, char *argv[])
{
    control_par *cpar;
    track_par *tpar;
    volume_par *vpar;

    int count;
    DIR *dirp;
    struct dirent *dp;
    char file_name[256];
    
    // 1. process inputs: directory, first frame, last frame
  
  printf ("This program was called with \"%s\".\n",argv[0]);
  
  if (argc != 4) 
    {
        printf("Wrong number of inputs, expecting: ./openptv test_cavity 10000 10004 \n");
        return 0;
    }

  if (argc > 1)
    {
      for (count = 1; count < argc; count++)
	    {
	        printf("argv[%d] = %s \n", count, argv[count]);
	    }
    }
    
    // 2. init_proc - irrelevant?
    // 3. start_proc - irrelevant?
    // 4. read parameters

    // change directory to the user-supplied working folder
    chdir(argv[1]);
    // dirp = opendir(argv[1]);
    
    cpar = read_control_par("parameters/ptv.par");
    printf(" ------------------ \n");
    printf("Control parameters \n");
    printf(" ------------------ \n");
    printf(" highpass flag = %d \n", cpar->hp_flag);
    printf(" use all cameras flag = %d \n", cpar->allCam_flag);
    printf(" TIFF flag = %d \n", cpar->tiff_flag);
    printf(" image is  %d x %d \n", cpar->imx, cpar->imy);
    printf(" pixel size  %4.3f x %4.3f \n", cpar->pix_x, cpar->pix_y);
    printf(" chfield flag = %d \n", cpar->chfield);
    printf(" Multimedia parameters: \n");
    printf(" gas/air index of refraction = %3.2f \n", cpar->mm->n1);
    printf(" glass/perspex index of refraction = %3.2f \n", cpar->mm->n2[0]);
    printf(" water/liquid index of refraction = %3.2f \n", cpar->mm->n3);
    printf(" glass thickness = %3.2f \n",cpar->mm->d[0]);
    
    
    vpar = read_volume_par("parameters/criteria.par");
    tpar = read_track_par("parameters/track.par");

    // 5. sequence (init, set images, loop)
    // 6. tracking (init, loop, finish)
        
    return 0;
}
