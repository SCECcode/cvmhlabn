/**  
   test_vx_cvmhlabn_exec.c

   invokes src/run_vx_cvmhlabn.sh/vx_cvmhlabn
     which uses cvmhlabn api,
       cvmhlabn_init, cvmhlabn_setparam, cvmhlabn_query, cvmhlabn_finalize
**/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <getopt.h>
#include "vx_sub.h"
#include "unittest_defs.h"
#include "test_helper.h"
#include "test_vx_cvmhlabn_exec.h"

int VX_TESTS=3;

int test_vx_cvhmlabn_points_elevation()
{
  char infile[1280];
  char outfile[1280];
  char reffile[1280];
  char currentdir[1000];

  printf("Test: vx_cvmhlabn executable with elevation option\n");

  /* Save current directory */
  getcwd(currentdir, 1000);

  sprintf(infile, "%s/%s", currentdir, "./inputs/test_elev.in");
  sprintf(outfile, "%s/%s", currentdir, 
	  "test-10-point-vx-cvmhlabn-extract-elev.out");
  sprintf(reffile, "%s/%s", currentdir, 
	  "./ref/test-10-point-vx-cvmhlabn-extract-elev.ref");

  if (test_assert_int(save_elevation_test_points(infile), 0) != 0) {
    return(1);
  }

  if (test_assert_int(runVXCVMHLABN(BIN_DIR, MODEL_DIR, infile, outfile, 
				MODE_ELEVATION), 0) != 0) {
    printf("vx_cvmhlabn failure\n");
    return(1);
  }

  /* Perform diff btw outfile and ref */
  if (test_assert_file(outfile, reffile) != 0) {
    return(1);
  }

//  unlink(infile);
//  unlink(outfile);

  printf("PASS\n");
  return(0);
}



int test_vx_cvhmlabn_points_depth()
{
  char infile[1280];
  char outfile[1280];
  char reffile[1280];
  char currentdir[1000];

  printf("Test: vx_cvmhlabn executable with depth option\n");

  /* Save current directory */
  getcwd(currentdir, 1000);

  sprintf(infile, "%s/%s", currentdir, "./inputs/test_depth.in");
  sprintf(outfile, "%s/%s", currentdir, 
	  "test-10-point-vx-cvmhlabn-extract-depth.out");
  sprintf(reffile, "%s/%s", currentdir, 
	  "./ref/test-10-point-vx-cvmhlabn-extract-depth.ref");

  if (test_assert_int(save_depth_test_points(infile), 0) != 0) {
    printf("save test point failure\n");
    return(1);
  }

  if (test_assert_int(runVXCVMHLABN(BIN_DIR, MODEL_DIR, infile, outfile, 
				MODE_DEPTH), 0) != 0) {
    printf("vx_cvmhlabn failure\n");
    return(1);
  }  

  /* Perform diff btw outfile and ref */
  if (test_assert_file(outfile, reffile) != 0) {
    printf("diff failure\n");
    return(1);
  }

//  unlink(infile);
//  unlink(outfile);

  printf("PASS\n");
  return(0);
}

int test_vx_cvhmlabn_points_offset()
{
  char infile[1280];
  char outfile[1280];
  char reffile[1280];
  char currentdir[1000];

  printf("Test: vx_cvmhlabn executable with offset(none) option\n");

  /* Save current directory */
  getcwd(currentdir, 1000);

  sprintf(infile, "%s/%s", currentdir, "./inputs/test_offset.in");
  sprintf(outfile, "%s/%s", currentdir, 
	  "test-10-point-vx-cvmhlabn-extract-offset.out");
  sprintf(reffile, "%s/%s", currentdir, 
	  "./ref/test-10-point-vx-cvmhlabn-extract-offset.ref");

  if (test_assert_int(save_elevation_test_points(infile), 0) != 0) {
    return(1);
  }

  if (test_assert_int(runVXCVMHLABN(BIN_DIR, MODEL_DIR, infile, outfile, 
				MODE_NONE), 0) != 0) {
    printf("vx_cvmhlabn failure\n");
    return(1);
  }

  /* Perform diff btw outfile and ref */
  if (test_assert_file(outfile, reffile) != 0) {
    return(1);
  }

//  unlink(infile);
//  unlink(outfile);

  printf("PASS\n");
  return(0);
}


int suite_vx_cvmhlabn_exec(const char *xmldir)
{
  suite_t suite;
  char logfile[1280];
  FILE *lf = NULL;

  /* Setup test suite */
  strcpy(suite.suite_name, "suite_vx_cvhmlabn_exec");

  suite.num_tests = VX_TESTS;
  suite.tests = malloc(suite.num_tests * sizeof(test_t));
  if (suite.tests == NULL) {
    fprintf(stderr, "Failed to alloc test structure\n");
    return(1);
  }
  test_get_time(&suite.exec_time);

  /* Setup test cases */
  strcpy(suite.tests[0].test_name, "test_vx_cvhmlabn_points_elevation");
  suite.tests[0].test_func = &test_vx_cvhmlabn_points_elevation;
  suite.tests[0].elapsed_time = 0.0;

  strcpy(suite.tests[1].test_name, "test_vx_cvhmlabn_points_depth");
  suite.tests[1].test_func = &test_vx_cvhmlabn_points_depth;
  suite.tests[1].elapsed_time = 0.0;

  strcpy(suite.tests[1].test_name, "test_vx_cvhmlabn_points_offset");
  suite.tests[2].test_func = &test_vx_cvhmlabn_points_offset;
  suite.tests[2].elapsed_time = 0.0;

  if (test_run_suite(&suite) != 0) {
    fprintf(stderr, "Failed to execute tests\n");
    return(1);
  }

  if (xmldir != NULL) {
    sprintf(logfile, "%s/%s.xml", xmldir, suite.suite_name);
    lf = init_log(logfile);
    if (lf == NULL) {
      fprintf(stderr, "Failed to initialize logfile\n");
      return(1);
    }
    
    if (write_log(lf, &suite) != 0) {
      fprintf(stderr, "Failed to write test log\n");
      return(1);
    }
    
    close_log(lf);
  }

  free(suite.tests);

  return 0;
}
