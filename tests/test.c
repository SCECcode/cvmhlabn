/**
 * @file test.c
 * @brief Bootstraps the test framework for the CVMHLABN library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the CVMHLABN library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cvmhlabn.h"

int compare_double(double f1, double f2) {
  double precision = 0.00001;
  if (((f1 - precision) < f2) && ((f1 + precision) > f2)) {
    return 1;
    } else {
      return 0;
  }
}

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	cvmhlabn_point_t pt;
	cvmhlabn_properties_t ret;

	// Initialize the model. 
        // try to use Use UCVM_INSTALL_PATH
        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
	   assert(cvmhlabn_init(envstr, "cvmhlabn") == 0);
           } else {
	     assert(cvmhlabn_init("..", "cvmhlabn") == 0);
        }

	printf("Loaded the model successfully.\n");

        // id = 99, don't care 
        cvmhlabn_setparam(99,CVMHLABN_PARAM_QUERY_MODE,CVMHLABN_COORD_GEO_DEPTH);

	// Query a point.
	pt.longitude = -118.0701;
	pt.latitude = 34.155;
	pt.depth = 2000;

	cvmhlabn_query(&pt, &ret, 1);

        //printf("vs : %lf\n",ret.vs);
        //printf("vp : %lf\n",ret.vp);
        //printf("rho: %lf\n",ret.rho);

	assert(ret.vs == -1);
	assert(ret.vp == -1);
	assert(ret.rho == -1);

 // next point

	// Query a point.
	pt.longitude = -118.1;
	pt.latitude = 34.0;
	pt.depth = 1500;

	cvmhlabn_query(&pt, &ret, 1);

        //printf("vs : %lf\n",ret.vs);
        //printf("vp : %lf\n",ret.vp);
        //printf("rho: %lf\n",ret.rho);
 
        assert(compare_double(ret.vs, 1569.190063));
        assert(compare_double(ret.vp, 3180.260498));
        assert(compare_double(ret.rho, 2261.115808));

	printf("Query was successful.\n");

	// Close the model.
	assert(cvmhlabn_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL CVMHLABN TESTS PASSED\n");

	return 0;
}
