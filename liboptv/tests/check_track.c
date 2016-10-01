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

#define EPS 1E-5


START_TEST(test_predict)
{
    vec2d a = {1.1, 0.6};
    vec2d b = {2.0, -0.8};
    vec2d result = {2.9, -2.2};

    vec2d c;

    predict(a,b,c);

    ck_assert_msg( fabs(c[0] - result[0]) < EPS,
             "Was expecting 2.9 but found %f \n", fabs(c[0]));

    ck_assert_msg( fabs(c[1] - result[1]) < EPS,
             "Was expecting -2.2 but found %f \n", fabs(c[1]));

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



Suite* fb_suite(void) {
    Suite *s = suite_create ("ttools");

    TCase *tc = tcase_create ("predict test");
    tcase_add_test(tc, test_predict);
    suite_add_tcase (s, tc);

    tc = tcase_create ("candsearch_in_pix");
    tcase_add_test(tc, test_candsearch_in_pix);
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
