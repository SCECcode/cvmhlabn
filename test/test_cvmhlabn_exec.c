/**
   test_cvmhlabn_exec.c

   uses cvmhlabn's model api,
       model_init, model_setparam, model_query, model_finalize
**/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include "cvmhlabn.h"
#include "unittest_defs.h"
#include "test_helper.h"
#include "test_cvmhlabn_exec.h"


int test_setup()
{
  printf("Test: model_init and model_finalize\n");

  char *envstr=getenv("UCVM_INSTALL_PATH");
  if(envstr != NULL) {
    if (test_assert_int(model_init(envstr, "cvmhlabn"), 0) != 0) {
      return(1);
    }
  } else if (test_assert_int(model_init("..", "cvmhlabn"), 0) != 0) {
    return(1);
  }

  if (test_assert_int(model_finalize(), 0) != 0) {
    return(1);
  }

  printf("PASS\n");
  return(0);
}

int test_setparam()
{
  printf("Test: model_setparam() with depth\n");

// Initialize the model, try to use Use UCVM_INSTALL_PATH
  char *envstr=getenv("UCVM_INSTALL_PATH");
  if(envstr != NULL) {
    if (test_assert_int(model_init(envstr, "cvmhlabn"), 0) != 0) {
      return(1);
    }
  } else if (test_assert_int(model_init("..", "cvmhlabn"), 0) != 0) {
    return(1);
  }

  int zmode = CVMHLABN_COORD_GEO_DEPTH;
  if (test_assert_int(model_setparam(0, CVMHLABN_PARAM_QUERY_MODE, zmode), 0) != 0) {
      return(1);
  }

  // Close the model.
  assert(model_finalize() == 0);

  printf("PASS\n");
  return(0);
}


int test_query_by_depth()
{
  printf("Test: model_query() by depth\n");

  cvmhlabn_point_t pt;
  cvmhlabn_properties_t ret;

// Initialize the model, try to use Use UCVM_INSTALL_PATH
  char *envstr=getenv("UCVM_INSTALL_PATH");
  if(envstr != NULL) {
    if (test_assert_int(model_init(envstr, "cvmhlabn"), 0) != 0) {
      return(1);
    }
  } else if (test_assert_int(model_init("..", "cvmhlabn"), 0) != 0) {
    return(1);
  }

  int zmode = CVMHLABN_COORD_GEO_DEPTH;
  if (test_assert_int(model_setparam(0, CVMHLABN_PARAM_QUERY_MODE, zmode), 0) != 0) {
      return(1);
  }

  // Query a point.
  pt.longitude = -118.1;
  pt.latitude = 34.0;
  pt.depth = 1500;

  if (test_assert_int(model_query(&pt, &ret, 1), 0) != 0) {
      return(1);
  }

  if ((test_assert_double(ret.vs, 1569.190063), 0) != 0) {
      return(1);
  }
  if ((test_assert_double(ret.vp, 3180.260498), 0) != 0) {
      return(1);
  }
  if ((test_assert_double(ret.vp, 3180.260498), 0) != 0) {
      return(1);
  }
  // printf("Query was successful.\n");

  // Close the model.
  assert(model_finalize() == 0);

  printf("PASS\n");
  return(0);
}

int test_query_by_elevation()
{
  printf("Test: model_query() by elevation\n");

  cvmhlabn_point_t pt;
  cvmhlabn_properties_t ret;

// Initialize the model, try to use Use UCVM_INSTALL_PATH
  char *envstr=getenv("UCVM_INSTALL_PATH");
  if(envstr != NULL) {
    if (test_assert_int(model_init(envstr, "cvmhlabn"), 0) != 0) {
      return(1);
    }
  } else if (test_assert_int(model_init("..", "cvmhlabn"), 0) != 0) {
    return(1);
  }

  int zmode = CVMHLABN_COORD_GEO_ELEV;
  if (test_assert_int(model_setparam(0, CVMHLABN_PARAM_QUERY_MODE, zmode), 0) != 0) {
      return(1);
  }

  // Query a point.
  pt.longitude = -118.1;
  pt.latitude = 34.0;
  pt.depth = -1500;

  if (test_assert_int(model_query(&pt, &ret, 1), 0) != 0) {
      return(1);
  }

  if ((test_assert_double(ret.vs, 1569.190063), 0) != 0) {
      return(1);
  }
  if ((test_assert_double(ret.vp, 3180.260498), 0) != 0) {
      return(1);
  }
  if ((test_assert_double(ret.vp, 3180.260498), 0) != 0) {
      return(1);
  }
  //printf("Query was successful.\n");

  // Close the model.
  assert(model_finalize() == 0);

  printf("PASS\n");
  return(0);
}

int suite_cvmhlabn_exec(const char *xmldir)
{
  suite_t suite;
  char logfile[256];
  FILE *lf = NULL;

  /* Setup test suite */
  strcpy(suite.suite_name, "suite_vcmhlabn_exec");
  suite.num_tests = 4;
  suite.tests = malloc(suite.num_tests * sizeof(test_t));
  if (suite.tests == NULL) {
    fprintf(stderr, "Failed to alloc test structure\n");
    return(1);
  }
  test_get_time(&suite.exec_time);

  /* Setup test cases */
  strcpy(suite.tests[0].test_name, "test_setup()");
  suite.tests[0].test_func = &test_setup;
  suite.tests[0].elapsed_time = 0.0;

  strcpy(suite.tests[1].test_name, "test_separam()");
  suite.tests[1].test_func = &test_setparam;
  suite.tests[1].elapsed_time = 0.0;

  strcpy(suite.tests[2].test_name, "test_query_by_depth()");
  suite.tests[2].test_func = &test_query_by_depth;
  suite.tests[2].elapsed_time = 0.0;

  strcpy(suite.tests[3].test_name, "test_query_by_elevation()");
  suite.tests[3].test_func = &test_query_by_elevation;
  suite.tests[3].elapsed_time = 0.0;

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
