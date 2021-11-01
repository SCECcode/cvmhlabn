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

	// Query a point.
	pt.longitude = -118.0701;
	pt.latitude = 34.155;
	pt.depth = 2000;

	cvmhlabn_query(&pt, &ret, 1);

        printf("vs : %lf\n",ret.vs);
        printf("vp : %lf\n",ret.vp);
        printf("rho: %lf\n",ret.rho);

// LABN as LR, LABN as LR&HR
	assert(ret.vs == -1);
	assert(ret.vp == -1);
	assert(ret.rho == -1);
// for CVM_LR.vo
	//assert(ret.vs == 2722.399);
	//assert(ret.vp == 4907.421);
	//assert(ret.rho == 2520.679);

	printf("Query was successful.\n");

	// Close the model.
	assert(cvmhlabn_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL CVMHLABN TESTS PASSED\n");

	return 0;
}
