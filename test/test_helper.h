#ifndef TEST_HELPER_H
#define TEST_HELPER_H

#include "vx_sub.h"
#include "cvmhlabn.h"

/* Constants */
#define MAX_TEST_POINTS 10
#define PLACEHOLDER -99999.0

/* modes of operation */
#define MODE_NONE 0
#define MODE_ELEVATION 2
#define MODE_DEPTH 4

/* Test data sets */
typedef enum { VX_TEST_DATASET_NOBKG = 0, 
	       VX_TEST_DATASET_BKG, 
	       VX_TEST_DATASET_NOGTL } vx_test_dataset_t;

int test_assert_file_exist(const char* filename);
double get_preset_ucvm_surface(double, double);

/* Retrieve basic test points */
int get_test_points(double *x, double *y, double *z, 
                    vx_coord_t *coord_types);

/* Retrieve expected surface elev at the test points */
int get_surf_values(double *surf_values);

/* Retrieve expected mat props at the test points */
int get_mat_props(float *vp, float *vs, double *rho, vx_test_dataset_t ds);

/* Save test points to file */
int save_depth_test_points(const char* filename);
int save_elevation_test_points(const char* filename);

/* run with cvmhlabn model api */
int runCVMHLABN(const char *bindir, const char *cvmdir, 
	  const char *infile, const char *outfile,
          int mode);

/* Execute vx_cvmhlabn as a child process */
int runVXCVMHLABN(const char *bindir, const char *cvmdir, 
	  const char *infile, const char *outfile,
          int mode);

/* Execute vx_lite_cvmhlabn as a child process */
int runVXLiteCVMHLABN(const char *bindir, const char *cvmdir, 
	      const char *infile, const char *outfile,
	      int mode);

#endif
