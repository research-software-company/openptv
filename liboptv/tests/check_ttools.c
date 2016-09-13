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


START_TEST(test_ttools)
{
    vec3d a = {1.0, 0.0, 0.0};
    vec3d b = {1.0, 0.0, 0.0};
    
    ck_assert_msg( fabs(a - a) < EPS,
             "Was expecting a-a to be 0.0 but found %f \n", fabs(a-a));

}
END_TEST



Suite* fb_suite(void) {
    Suite *s = suite_create ("ttools");
 
    TCase *tc = tcase_create ("dummy test");
    tcase_add_test(tc, test_ttools);
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

