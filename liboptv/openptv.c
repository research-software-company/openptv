//
//  openptv.c
//  OpenPTV
//
//  Created by Alex Liberzon on 11/05/2017.
//
//

#include "openptv.h"
void read_all_calibration(Calibration *calib[], int num_cams);
void imread(unsigned char *img, char *filename);
coord_2d **correct_frame(frame *frm, Calibration *calib[], control_par *cpar,
                         double tol);

/* 
    Compile it using: 
       make 
    Run it using: 
       ./openptv test_cavity 10000 10004
*/
int main(int argc, char *argv[])
{
    int i, count, ntargets;
    DIR *dirp;
    struct dirent *dp;
    char file_name[256];
    int step;
    unsigned char *img, *img_hp;
    target pix[1024];
    coord_2d **corrected;
    int match_counts[4];
    frame *frm;
    int max_targets = 1024;
    
    
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
    free_control_par(cpar);
    
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
    
    // allocate memory for images and highpass(prepared) images
    img = (unsigned char *) malloc(run->cpar->imx*run->cpar->imy* \
                                                     sizeof(unsigned char));
    img_hp = (unsigned char *) malloc(run->cpar->imx*run->cpar->imy* \
                                                     sizeof(unsigned char));
                                                     
    
    // for each camera and for each time step the images are processed
    for (i = 1; i<run->cpar->num_cams+1; i++) {
        for (step = run->seq_par->first; step < run->seq_par->last+1; step++) {
            // a. read image
            sprintf(file_name, "img/cam%d.%d", i, step);
            imread(img, file_name);
            // b. highpass
            if (run->cpar->hp_flag)
            {
                prepare_image(img, img_hp, 1, 0, 0, run->cpar);
            } else {
                memcpy(img_hp, img, run->cpar->imx*run->cpar->imy);
            }
            // c. segmentation
            // detection
            //ntargets = peak_fit(img_hp, targ_read, 0, run->cpar->imx, 0, run->cpar->imy, run->cpar, 1, pix);
            ntargets = targ_rec(img_hp, targ_read, 0, run->cpar->imx, 0, run->cpar->imy, run->cpar, 1, pix);
            // write the target files, used later by tracking
            sprintf(file_name, "img/cam%d.", i, step);
            write_targets(pix, ntargets, file_name, step);
       }
    }
    
    // d. prepare the frame buffer, read targets and find correspondences
    for (step = run->seq_par->first; step < run->seq_par->last+1; step++){
            //frame_init(frm, run->cpar->num_cams, max_targets);
        for (i = 1; i<run->cpar->num_cams+1; i++) {
            sprintf(file_name, "img/cam%d.", i, step);
            ntargets = read_targets(pix, file_name, step);
            // printf(" read %d targets from %s \n", ntargets, file_name);
            //frm->num_targets[i] = ntargets;
            //frm->targets[i] = pix;
           }
        // corrected = correct_frame(frm, calib, run->cpar, 0.0001);
        // correspondences(frm, corrected, run->vpar, run->cpar, calib, match_counts);
    }
    
    // to be completed
    
    // 4. tracking (init, loop, finish)
    track_forward_start(run);

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

/* imread - reads greyscale TIFF file using libtiff library subroutines
   Arguments: 
   unsigned char *img - output of the image array
   char *filename - name of the TIFF file to read
*/

void imread(unsigned char *img, char *filename)
{
    TIFF* tif = TIFFOpen(filename, "r");
    if (tif) {
        uint32 imagelength;
        tsize_t scanline;
        tdata_t buf;
        uint32 row;
        uint32 col;
        
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
        scanline = TIFFScanlineSize(tif);
        // buf = _TIFFmalloc(scanline);
        for (row = 0; row < imagelength; row++)
        {
            TIFFReadScanline(tif, &img[row*scanline], row, 1);
//            for (col = 0; col < scanline; col++)
//            {
//                printf("%u", img[row*scanline+col]);
//            }
//            printf("\n");
        }
        // _TIFFfree(buf);
        TIFFClose(tif);
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
    
    corrected = (coord_2d **) malloc(cpar->num_cams * sizeof(coord_2d *));
    for (cam = 0; cam < cpar->num_cams; cam++) {
        corrected[cam] = (coord_2d *) malloc(
                                             frm->num_targets[cam] * sizeof(coord_2d));
        if (corrected[cam] == NULL) {
            /* roll back allocations and fail */
            for (cam -= 1; cam >= 0; cam--) free(corrected[cam]);
            free(corrected);
            return NULL;
        }
        
        for (part = 0; part < frm->num_targets[cam]; part++) {
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


