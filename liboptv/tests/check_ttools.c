/*  Unit tests for the vector utilities. Uses the Check
    framework: http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
    If that doesn't run the tests, use the Check tutorial:
    http://check.sourceforge.net/doc/check_html/check_3.html
*/

#include <check.h>
#include <stdlib.h>
#include <math.h>
#include "ttools.h"

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



Suite* fb_suite(void) {
    Suite *s = suite_create ("ttools");
 
    TCase *tc = tcase_create ("predict test");
    tcase_add_test(tc, test_predict);
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

