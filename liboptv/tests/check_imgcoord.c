/* Unit tests for ray tracing. */

#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "imgcoord.h"

#define EPS 1E-4



START_TEST(test_img_xy_mm_geo)
{

    /* input */
    double x = 100.0;
    double y = -50.0;
    double z =  10.0;
    
    double xa,ya;  
    
    int i, i_cam = 0;      
        
    /* parameters */
    Calibration *cal;

    char ori_file[] = "testing_fodder/cal/cam3.tif.ori";
    char add_file[] = "testing_fodder/cal/cam3.tif.addpar";    
    cal = read_calibration(ori_file, add_file, NULL);
    
        
    volume_par *vpar;
    vpar = read_volume_par("testing_fodder/parameters/criteria_2.par");
    
    
    control_par *cpar;
    char filename[] = "testing_fodder/parameters/ptv_2.par";
    cpar = read_control_par(filename);
    /* two default values which are not in the parameter file */
    cpar->mm->lut = 0;
    cpar->mm->nlay = 1;
        



     mmlut test_mmlut[4], correct_mmlut[4];  
     
     correct_mmlut[0].origin.x = 0.0;
     correct_mmlut[0].origin.y = 0.0;
     correct_mmlut[0].origin.z = -70.0;
     correct_mmlut[0].nr = 9;
     correct_mmlut[0].nz = 47;
     correct_mmlut[0].rw = 2;
             
     init_mmLUT (vpar
               , cpar
               , cal
               , test_mmlut);
                             
     for (i=0; i<cpar->num_cams; i++){
       ck_assert_msg( 
                    fabs(test_mmlut[i].origin.x - correct_mmlut[0].origin.x) < EPS && 
                    fabs(test_mmlut[i].origin.y - correct_mmlut[0].origin.y) < EPS && 
                    fabs(test_mmlut[i].origin.z - correct_mmlut[0].origin.z)  < EPS &&
                    test_mmlut[i].nr == correct_mmlut[i].nr &&
                    test_mmlut[i].nz == correct_mmlut[i].nz &&
                    test_mmlut[i].rw ==  correct_mmlut[i].rw &&
                    fabs(test_mmlut[i].data[35] - 1.000000) < EPS &&
                    fabs(test_mmlut[i].data[82] - 1.003130) < EPS,
         "\n Expected different correct_mmlut values \n  \
         but found %4.3f %4.3f %4.3f %d %d %d %4.3 %4.3f in camera %d\n", \
         test_mmlut[i].origin.x, test_mmlut[i].origin.y, test_mmlut[i].origin.z, \
         test_mmlut[i].nr, test_mmlut[i].nz, test_mmlut[i].rw, 1.0, 1.003130, i);
          
        }     
    
    
    
    img_xy_mm_geo (x, y, z, cal, *cpar->mm, i_cam, (mmlut *) test_mmlut, &xa, &ya); 
   
    ck_assert_msg( fabs(xa - 71.1275) < EPS &&
                   fabs(ya + 35.5637) < EPS,
         "\n Expected 71.1275,-35.5637  \n  \
         but found %6.4f %6.4f \n", xa, ya);	
 
}
END_TEST

START_TEST(test_img_coord)
{

    /* input */
    double x = 100.0;
    double y = -50.0;
    double z =  10.0;
    
    double xa,ya;  
    
    int i, i_cam = 0;      
        
    /* parameters */
    Calibration *cal;

    char ori_file[] = "testing_fodder/cal/cam3.tif.ori";
    char add_file[] = "testing_fodder/cal/cam3.tif.addpar";    
    cal = read_calibration(ori_file, add_file, NULL);
    
        
    volume_par *vpar;
    vpar = read_volume_par("testing_fodder/parameters/criteria_2.par");
    
    
    control_par *cpar;
    char filename[] = "testing_fodder/parameters/ptv_2.par";
    cpar = read_control_par(filename);
    /* two default values which are not in the parameter file */
    cpar->mm->lut = 0;
    cpar->mm->nlay = 1;
        



     mmlut test_mmlut[4], correct_mmlut[4];  
     
     correct_mmlut[0].origin.x = 0.0;
     correct_mmlut[0].origin.y = 0.0;
     correct_mmlut[0].origin.z = -70.0;
     correct_mmlut[0].nr = 9;
     correct_mmlut[0].nz = 47;
     correct_mmlut[0].rw = 2;
             
     init_mmLUT (vpar
               , cpar
               , cal
               , test_mmlut);
                             
     for (i=0; i<cpar->num_cams; i++){
       ck_assert_msg( 
                    fabs(test_mmlut[i].origin.x - correct_mmlut[0].origin.x) < EPS && 
                    fabs(test_mmlut[i].origin.y - correct_mmlut[0].origin.y) < EPS && 
                    fabs(test_mmlut[i].origin.z - correct_mmlut[0].origin.z)  < EPS &&
                    test_mmlut[i].nr == correct_mmlut[i].nr &&
                    test_mmlut[i].nz == correct_mmlut[i].nz &&
                    test_mmlut[i].rw ==  correct_mmlut[i].rw &&
                    fabs(test_mmlut[i].data[35] - 1.000000) < EPS &&
                    fabs(test_mmlut[i].data[82] - 1.003130) < EPS,
         "\n Expected different correct_mmlut values \n  \
         but found %4.3f %4.3f %4.3f %d %d %d %4.3 %4.3f in camera %d\n", \
         test_mmlut[i].origin.x, test_mmlut[i].origin.y, test_mmlut[i].origin.z, \
         test_mmlut[i].nr, test_mmlut[i].nz, test_mmlut[i].rw, 1.0, 1.003130, i);
          
        }     
    
    
    
    img_coord (x, y, z, cal, *cpar->mm, i_cam, (mmlut *) test_mmlut, &xa, &ya); 
   
    ck_assert_msg( fabs(xa - 69.5061) < EPS &&
                   fabs(ya + 35.5628) < EPS,
         "\n Expected 69.5061,-35.5628  \n  \
         but found %6.4f %6.4f \n", xa, ya);
         
    /* Note that with this precision, 0.001 is almost 0.0, will probably fail on 
    64 bit machines
    */     
    x = 0.001;
    y = 0.001;
    z = -10000.0; 	
    
    img_coord (x, y, z, cal, *cpar->mm, i_cam, (mmlut *) test_mmlut, &xa, &ya); 
   
    ck_assert_msg( fabs(xa - 0.0) < EPS &&
                   fabs(ya + 0.0) < EPS,
         "\n Expected 0.0, 0.0 at least on 32-bit machine  \n  \
         but found %f %f \n", xa, ya);

 
}
END_TEST


Suite* fb_suite(void) {
    Suite *s = suite_create ("imgcoord");
    TCase *tc = tcase_create ("imgcoord_test");
    tcase_add_test(tc, test_img_xy_mm_geo);
    tcase_add_test(tc, test_img_coord);    
    suite_add_tcase (s, tc);   
    return s;
}

int main(void) {
    int number_failed;
    Suite *s = fb_suite ();
    SRunner *sr = srunner_create (s);
    //srunner_run_all (sr, CK_ENV);
    //srunner_run_all (sr, CK_SUBUNIT);
    srunner_run_all (sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

