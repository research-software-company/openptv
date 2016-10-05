/* Unit tests for reading and writing parameter files. */

#include <check.h>
#include <stdio.h>
#include "segmentation.h"
#include "tracking_frame_buf.h"



/* create_empty_img creates a black (full of 0) image by 
 * allocating an array and filling with zeros.
 * Arguments: 
 * integer width, height of the image
 * Returns:
 * unsigned char* pointer to the image array
 */
unsigned char * create_empty_img(int width,int height){
    
    unsigned char* img = malloc(width*height*sizeof(unsigned char));
    
    for (int row=0;row<height;row++){
        for (int col=0;col<width;col++){
            img[row*width+col] = 0;
        }
    }
    return img;
}

/* create_blob_img creates a black (full of 0) image by
 * allocating an array and filling with zeros and the fills the center
 * by 255 (white pixels). Only a 1-pixel boundaries remain black.
 * Arguments:
 * integer width, height of the image
 * Returns:
 * unsigned char* pointer to the image array
 */
unsigned char * create_blob_img(int width,int height){
    
    unsigned char* img = create_empty_img(width,height);
    
    for (int row=1;row<height-1;row++){
        for (int col=1;col<width-1;col++){
            img[row*width+col] = 255;
        }
    }
    return img;
}

START_TEST(test_peak_fit)
{
    int ntargets; 
    
    control_par cpar = {
        .imx = 5,
        .imy = 5,
    };
    
    unsigned char* img = create_blob_img(cpar.imx,cpar.imy);

    
    
    target pix[1024];
    



    target_par targ_par= { 
        .gvthres = {250, 100, 20, 20}, 
        .discont = 5,
        .nnmin = 1, .nnmax = 10,
        .nxmin = 1, .nxmax = 10,
        .nymin = 1, .nymax = 10, 
        .sumg_min = 12, 
        .cr_sz = 13 };
    
                
   ntargets = peak_fit(img, &targ_par, 0, cpar.imx, 0, cpar.imy, &cpar, 0, pix);
   fail_unless(ntargets == 1);
   fail_unless(pix[0].n == 9);
   
   /* test the two objects */
    unsigned char* img1 = create_empty_img(cpar.imx,cpar.imy);
    img1[6] = 255;
    img1[18] = 251;
    
   ntargets = peak_fit(img1, &targ_par, 0, cpar.imx, 0, cpar.imy, &cpar, 1, 
        pix);
   fail_unless(ntargets == 2);
   
   targ_par.gvthres[1] = 252; 
   ntargets = peak_fit((unsigned char *)img1, &targ_par, 0, cpar.imx, 
        0, cpar.imy, &cpar, 1, pix);
   fail_unless(ntargets == 1);

}
END_TEST

START_TEST(test_targ_rec)
{
    int ntargets; 
   
    target pix[1024];
    
    control_par cpar = {
        .imx = 5,
        .imy = 5,
    }; 
    unsigned char* img = create_blob_img(cpar.imx,cpar.imy);

    target_par targ_par= { 
        .gvthres = {250, 100, 20, 20}, 
        .discont = 5,
        .nnmin = 1, .nnmax = 10,
        .nxmin = 1, .nxmax = 10,
        .nymin = 1, .nymax = 10, 
        .sumg_min = 12, 
        .cr_sz = 13 };
    
                
    ntargets = targ_rec (img, &targ_par, 0, cpar.imx, 0, cpar.imy, &cpar, 0, pix);
    fail_unless(ntargets == 1);
    fail_unless(pix[0].n == 9);
    fail_unless(pix[0].tnr == CORRES_NONE);
   
    /* test the two objects */
    unsigned char* img1 = create_empty_img(cpar.imx,cpar.imy);
    img1[6] = 255;
    img1[18] = 251;
    
    
    ntargets = targ_rec (img1, &targ_par, 0, cpar.imx, 0, cpar.imy, &cpar, 1, 
        pix);
    fail_unless(ntargets == 2);
   
    targ_par.gvthres[1] = 252; 
    ntargets = targ_rec (img1, &targ_par, 0, cpar.imx,
        0, cpar.imy, &cpar, 1, pix);
    fail_unless(ntargets == 1);

    /* Trip a segfault writing over the edge. */
    img1[4*cpar.imx+4] = 255;
    ntargets = targ_rec (img1, &targ_par, 0, cpar.imx, 0, cpar.imy, &cpar, 1,
        pix);
    /* If execution reached here, test passed. */
}
END_TEST



Suite* fb_suite(void) {
    Suite *s = suite_create ("Segmentation");
    
    TCase *tc = tcase_create ("Target recording");
    tcase_add_test(tc, test_targ_rec);
    suite_add_tcase (s, tc);

    tc = tcase_create ("Peak Fitting");
    tcase_add_test(tc, test_peak_fit);
    suite_add_tcase (s, tc);
    
    
    return s;
}

int main(void) {
    int number_failed;
    Suite *s = fb_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_ENV);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

