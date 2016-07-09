/* Unit tests for ray tracing. */

#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "epi.h"
#include "multimed.h"
#include "correspondences.h"

#define EPS 1E-5

START_TEST(test_correspondences)
{

    
}
END_TEST


Suite* fb_suite(void) {
    Suite *s = suite_create ("correspondences");
    
    TCase *tc = tcase_create ("test_corrrespondences");
    tcase_add_test(tc, test_corrrespondences);    
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

