/*  Unit tests for the tracking. Uses the Check
    framework: http://check.sourceforge.net/

    To run it, type "make verify" when in the top C directory, src_c/
    If that doesn't run the tests, use the Check tutorial:
    http://check.sourceforge.net/doc/check_html/check_3.html
*/

/* Unit tests for ray tracing. */

#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "track.h"
#include "calibration.h"

#define EPS 1E-5

/* Tests of correspondence components and full process using dummy data */

void read_all_calibration(Calibration *calib[4], control_par *cpar) {
    char ori_tmpl[] = "testing_fodder/cal/sym_cam%d.tif.ori";
    char added_name[] = "testing_fodder/cal/cam1.tif.addpar";
    char ori_name[40];
    int cam;

    for (cam = 0; cam < cpar->num_cams; cam++) {
        sprintf(ori_name, ori_tmpl, cam + 1);
        calib[cam] = read_calibration(ori_name, added_name, NULL);
    }
}

Calibration test_cal(void) {
    Exterior correct_ext = {
        105.2632, 102.7458, 403.8822,
        -0.2383291, 0.2442810, 0.0552577,
        {{0.9688305, -0.0535899, 0.2418587},
        {-0.0033422, 0.9734041, 0.2290704},
        {-0.2477021, -0.2227387, 0.9428845}}};
    Interior correct_int = {-2.4742, 3.2567, 100.0000};
    Glass correct_glass = {0.0001, 0.00001, 150.0};
    ap_52 correct_addp = {0., 0., 0., 0., 0., 1., 0.};
    Calibration correct_cal = {correct_ext, correct_int, correct_glass,
        correct_addp};
    rotation_matrix(&correct_cal.ext_par);

    return correct_cal;
}


START_TEST(test_predict)
{
    vec2d prev_pos = {1.1, 0.6};
    vec2d curr_pos = {2.0, -0.8};
    vec2d result = {2.9, -2.2};

    vec2d c;

    predict(prev_pos,curr_pos,c);


    ck_assert_msg( fabs(c[0] - result[0]) < EPS,
             "Was expecting 2.9 but found %f \n", fabs(c[0]));

    ck_assert_msg( fabs(c[1] - result[1]) < EPS,
             "Was expecting -2.2 but found %f \n", fabs(c[1]));

}
END_TEST


START_TEST(test_search_volume_center_moving)
{
    vec3d prev_pos = {1.1, 0.6, 0.1};
    vec3d curr_pos = {2.0, -0.8, 0.2};
    vec3d result = {2.9, -2.2, 0.3};

    vec3d c;

    search_volume_center_moving(prev_pos, curr_pos, c);


    ck_assert_msg( fabs(c[0] - result[0]) < EPS,
             "Was expecting 2.9 but found %f \n", c[0]);

    ck_assert_msg( fabs(c[1] - result[1]) < EPS,
             "Was expecting -2.2 but found %f \n", c[1]);

    ck_assert_msg( fabs(c[2] - result[2]) < EPS,
                      "Was expecting 0.3 but found %f \n", c[2]);

}
END_TEST

START_TEST(test_pos3d_in_bounds)
{
    vec3d inside = {1.0,-1.0,0.0};
    vec3d outside = {2.0, -0.8, 2.1};

    track_par bounds[] = {
        {0.4, 120, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 0., 0., 0., 0., 1.}
    };

    int result;
    result = pos3d_in_bounds(inside, bounds);

    ck_assert_msg( result == 1,
             "Was expecting True but found %f \n", result);

    result = pos3d_in_bounds(outside, bounds);

    ck_assert_msg( result == 0,
            "Was expecting False but found %f \n", result);

}
END_TEST

START_TEST(test_angle_acc)
{
    vec3d start = {0.0, 0.0, 0.0};
    vec3d pred  = {1.0, 1.0, 1.0};
    vec3d cand  = {1.1, 1.0, 1.0};

    double angle, acc;

    angle_acc(start, pred, cand, &angle, &acc);
    ck_assert_msg( fabs(angle - 2.902234) < EPS,
             "Was expecting 2.902234 but found %f \n", angle);

    ck_assert_msg( fabs(acc - 0.1) < EPS,
             "Was expecting 0.1 but found %f \n", acc);


    angle_acc(start, pred, pred, &angle, &acc);
    ck_assert_msg( fabs(acc) < EPS,
                      "Was expecting 0.0 but found %f \n", acc);
    ck_assert_msg( fabs(angle) < EPS,
                 "Was expecting 0.0 but found %f \n", angle);

    vec_scalar_mul(pred,-1,cand);
    angle_acc(start, pred, cand, &angle, &acc);

     ck_assert_msg( fabs(angle - 200.0) < EPS,
                    "Was expecting 200.0 but found %f \n", angle);

}
END_TEST




START_TEST(test_candsearch_in_pix)
{
    double cent_x, cent_y, dl, dr, du, dd;
    int p[4], counter = 0;

    target test_pix[] = {
        {0, 0.0, -0.2, 5, 1, 2, 10, -999},
        {6, 0.2, 0.2, 10, 8, 1, 20, -999},
        {3, 0.2, 0.3, 10, 3, 3, 30, -999},
        {4, 0.2, 1.0, 10, 3, 3, 40, -999},
        {1, -0.7, 1.2, 10, 3, 3, 50, -999},
        {7, 1.2, 1.3, 10, 3, 3, 60, -999},
        {5, 10.4, 2.1, 10, 3, 3, 70, -999}
    };
    int num_targets = 7;


    /* prepare test control parameters, basically for pix_x  */
    int cam;
    char img_format[] = "cam%d";
    char cal_format[] = "cal/cam%d.tif";
    control_par *test_cpar;

    test_cpar = new_control_par(4);
    for (cam = 0; cam < 4; cam++) {
        sprintf(test_cpar->img_base_name[cam], img_format, cam + 1);
        sprintf(test_cpar->cal_img_base_name[cam], cal_format, cam + 1);
    }
    test_cpar->hp_flag = 1;
    test_cpar->allCam_flag = 0;
    test_cpar->tiff_flag = 1;
    test_cpar->imx = 1280;
    test_cpar->imy = 1024;
    test_cpar->pix_x = 0.02; /* 20 micron pixel */
    test_cpar->pix_y = 0.02;
    test_cpar->chfield = 0;
    test_cpar->mm->n1 = 1;
    test_cpar->mm->n2[0] = 1.49;
    test_cpar->mm->n3 = 1.33;
    test_cpar->mm->d[0] = 5;

    cent_x = cent_y = 0.2;
    dl = dr = du = dd = 0.1;

    counter = candsearch_in_pix (test_pix, num_targets, cent_x, cent_y, \
                                 dl, dr, du, dd, p, test_cpar);

    printf("counter %d \n",counter);
    printf("candidates: \n");
    for (int i=0;i<counter;i++){
        printf("%f,%f\n",test_pix[p[i]].x,test_pix[p[i]].y);
    }
    fail_unless(counter == 2);


    cent_x = 0.5;
    cent_y = 0.3;
    dl = dr = du = dd = 10.2;

    counter = candsearch_in_pix (test_pix, num_targets, cent_x, cent_y, \
                                 dl, dr, du, dd, p, test_cpar);
    printf("counter %d \n",counter);
    printf("candidates:\n");
    for (int i=0;i<counter;i++){
        printf("%f,%f\n",test_pix[p[i]].x,test_pix[p[i]].y);
    }

    fail_unless(counter == 4);

}
END_TEST

START_TEST(test_sort)
{
    float test_array[] = {1.0, 2200.2, 0.3, -0.8, 100.0};
    int ix_array[] = {0,5,13,2,124};
    int len_array = 5;

    sort(len_array,test_array,ix_array);

    ck_assert_msg( fabs(test_array[0] + 0.8) < EPS,
             "Was expecting -0.8 but found %f \n", fabs(test_array[0]));

    ck_assert_msg( ix_array[len_array-1] != 1,
             "Was expecting 1 but found %f \n", ix_array[len_array-1]);

    printf("Sorted array:\n");
    for (int i=0;i<len_array;i++){
        printf("test_array[%d]=%f\n",ix_array[i],test_array[i]);
    }

}
END_TEST

START_TEST(test_copy_foundpix_array)
{
    foundpix src[] = {  {1,1,{1,0}},
                        {2,5,{1,1}}
                    };
    foundpix *dest;
    int arr_len = 2;
    int num_cams = 2;

    dest = (foundpix *) calloc (arr_len, sizeof (foundpix));

    reset_foundpix_array(dest, arr_len, num_cams);
    ck_assert_msg( dest[1].ftnr == -1 ,
             "Was expecting dest[1].ftnr == -1 but found %d \n", dest[1].ftnr);
    ck_assert_msg( dest[0].freq == 0 ,
              "Was expecting dest.freq == 0 but found %d \n", dest[0].freq);
    ck_assert_msg( dest[1].whichcam[0] == 0 ,
                       "Was expecting 0 but found %d \n", dest[1].whichcam[0]);



    copy_foundpix_array(dest, src, 2, 2);

    ck_assert_msg( dest[1].ftnr == 2 ,
             "Was expecting dest[1].ftnr == 2 but found %d \n", dest[1].ftnr);

    printf(" destination foundpix array\n");
    for (int i=0; i<arr_len; i++){
        printf("ftnr = %d freq=%d whichcam = %d %d\n", dest[i].ftnr, dest[i].freq, \
        dest[i].whichcam[0],dest[i].whichcam[1]);
    }


}
END_TEST


START_TEST(test_searchquader)
{
    vec3d point = {0.0, 0.0, 0.0};
    double xr[4], xl[4], yd[4], yu[4];
    Calibration *calib[4];
    control_par *cpar;

    Calibration cal;

    cal = test_cal();

    fail_if((cpar = read_control_par("testing_fodder/parameters/ptv.par"))== 0);

    /* see check_correspondences for explanations */
    cpar->mm->n2[0] = 1.0001;
    cpar->mm->n3 = 1.0001;

    printf (" number of cameras in searchquader test: %d \n",cpar->num_cams);

    track_par tpar[] = { {0.4, 120, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 0., 0., 0., 0., 1.} };

    read_all_calibration(calib, cpar);

    searchquader(point, xr, xl, yd, yu, tpar, cpar, calib);

    // printf("searchquader returned:\n");
    // for (int i=0; i<cpar->num_cams;i++){
    //     printf("%f %f %f %f\n",xr[i],xl[i],yd[i],yu[i]);
    // }
    ck_assert_msg( fabs(xr[0] - 47.060202)<EPS ,
             "Was expecting 47.06 but found %d \n", xr[0]);
    ck_assert_msg( fabs(yu[3] - 33.512680)<EPS ,
                      "Was expecting 33.512680 but found %d \n", yu[3]);

}
END_TEST





Suite* fb_suite(void) {
    Suite *s = suite_create ("ttools");

    TCase *tc = tcase_create ("predict test");
    tcase_add_test(tc, test_predict);
    suite_add_tcase (s, tc);

    tc = tcase_create ("search_volume_center_moving");
    tcase_add_test(tc, test_search_volume_center_moving);
    suite_add_tcase (s, tc);

    tc = tcase_create ("pos3d_in_bounds");
    tcase_add_test(tc, test_pos3d_in_bounds);
    suite_add_tcase (s, tc);

    tc = tcase_create ("angle_acc");
    tcase_add_test(tc, test_angle_acc);
    suite_add_tcase (s, tc);

    tc = tcase_create ("candsearch_in_pix");
    tcase_add_test(tc, test_candsearch_in_pix);
    suite_add_tcase (s, tc);

    tc = tcase_create ("sort");
    tcase_add_test(tc, test_sort);
    suite_add_tcase (s, tc);

    tc = tcase_create ("copy_foundpix_array");
    tcase_add_test(tc, test_copy_foundpix_array);
    suite_add_tcase (s, tc);

    tc = tcase_create ("searchquader");
    tcase_add_test(tc, test_searchquader);
    suite_add_tcase (s, tc);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s = fb_suite ();
    SRunner *sr = srunner_create (s);
    // srunner_run_all (sr, CK_ENV);
    srunner_run_all (sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
