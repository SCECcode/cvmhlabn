#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "unittest_defs.h"
#include "test_helper.h"

int debug_mode=0;

int test_assert_file_exist(const char* filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fclose(fp);
    return(1);
  }
  return(0);
}

double get_preset_ucvm_surface(double y, double x) {
   if(y == -118.1 && x == 34.0) {
      return 55.827;
   }
   return 0;
}

/* Retrieve 10 test points */
int get_test_points(double *x, double *y, double *z, 
                    vx_coord_t *coord_types)
{
  x[0] = -125.0;
  y[0] = 35.0;
  z[0] = -7777.0;
  coord_types[0] = VX_COORD_GEO;

  x[1] = -118.56;
  y[1] = 32.55;
  z[1] = -2450.0;
  coord_types[1] = VX_COORD_GEO;

  x[2] = 360061.0;
  y[2] = 3750229.0;
  z[2] = -1400.0;
  coord_types[2] = VX_COORD_UTM;

  x[3] = -118.52;
  y[3] = 34.12;
  z[3] = -1400.0;
  coord_types[3] = VX_COORD_GEO;

  x[4] = -116.40;
  y[4] = 32.34;
  z[4] = -1000.0;
  coord_types[4] = VX_COORD_GEO;

  x[5] = 376592.0;
  y[5] = 3773379.0;
  z[5] = -1700.0;
  coord_types[5] = VX_COORD_UTM;

  x[6] = 376592.0;
  y[6] = 3773379.0;
  z[6] = -17700.0;
  coord_types[6] = VX_COORD_UTM;

  x[7] = 408669.0;
  y[7] = 3766189.0;
  z[7] = -3000.0;
  coord_types[7] = VX_COORD_UTM;

// mt wilson
  x[8] = -118.0642;
  y[8] = 34.2264;
  z[8] = -1000.0;
  coord_types[8] = VX_COORD_GEO;

// cvmhlabn
  x[9] = 421000.0;
  y[9] = 3712000.0;
  z[9] = -2000.0;
  coord_types[9] = VX_COORD_UTM;

/* depth ..
  pt.longitude = -118.1;
  pt.latitude = 34.0;
  pt.depth = 1500;

  (ret.vs, 1569.190063)
  (ret.vp, 3180.260498)
  (ret.rho, 2261.115808)
*/

  return(0);
}

/* Retrieve expected surface elev at the test points */
// 8th column
int get_surf_values(double *surf_values)
{
  surf_values[0] = PLACEHOLDER;
  surf_values[1] = -1114.907715;
  surf_values[2] = -56.939297;
  surf_values[3] = 491.459229;
  surf_values[4] = 780.433289;
  surf_values[5] = 99.369789;
  surf_values[6] = 99.369789;
  surf_values[7] = 93.892632;
  surf_values[8] = 93.892632;
  surf_values[9] = 93.892632;

  return(0);
}


int save_elevation_test_points(const char* filename)
{
  double x[MAX_TEST_POINTS], y[MAX_TEST_POINTS], z[MAX_TEST_POINTS];
  vx_coord_t coord_types[MAX_TEST_POINTS];
  FILE *fp;
  int i;
  char line[128]; 
  size_t retval;

  // generate if not there
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fclose(fp);
    return(0);
  }

  get_test_points(x, y, z, coord_types);

  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("FAIL: cannot open %s\n", filename);
    return(1);
  }

  for (i = 0; i < MAX_TEST_POINTS; i++) {
    sprintf(line, "%f %f %f\n", x[i], y[i], z[i]);
    retval = fwrite(line, strlen(line), 1, fp);
    if (retval != 1) {
      printf("FAIL: write failed\n");
      return(1);
    }
  }
  fclose(fp);
  return(0);
}

int save_depth_test_points(const char* filename)
{
  double x[MAX_TEST_POINTS], y[MAX_TEST_POINTS], z[MAX_TEST_POINTS];
  double surf[MAX_TEST_POINTS];
  vx_coord_t coord_types[MAX_TEST_POINTS];
  FILE *fp;
  int i;
  char line[128]; 
  size_t retval;

  // generate if not there
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fclose(fp);
    return(0);
  }

  get_test_points(x, y, z, coord_types);
  get_surf_values(surf);

  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("FAIL: cannot open %s\n", filename);
    return(1);
  }

  for (i = 0; i < MAX_TEST_POINTS; i++) {
    sprintf(line, "%f %f %f\n", x[i], y[i], z[i]+surf[i]);
    retval = fwrite(line, strlen(line), 1, fp);
    if (retval != 1) {
      printf("FAIL: write failed\n");
      return(1);
    }
  }
  fclose(fp);
  return(0);
}

/*************************************************************************/
int runCVMHLABN(const char *bindir, const char *cvmdir, 
	  const char *infile, const char *outfile, int mode)
{
  cvmhlabn_point_t pt;
  cvmhlabn_properties_t ret;

  FILE *infp, *outfp;
  char line[1000];

  char *envstr=getenv("UCVM_INSTALL_PATH");
  if(envstr != NULL) {
    if (test_assert_int(model_init(envstr, "cvmhlabn"), 0) != 0) {
      return(1);
    }
  } else if (test_assert_int(model_init("..", "cvmhlabn"), 0) != 0) {
    return(1);
  }

  int zmode = UCVM_COORD_GEO_ELEVOFF;
  switch (mode) {
    case MODE_ELEVATION:
      zmode = UCVM_COORD_GEO_ELEV;
      break;
    case MODE_DEPTH:
      zmode = UCVM_COORD_GEO_DEPTH;
      break;
    case MODE_NONE:
      break; // default
  }

  if (test_assert_int(model_setparam(0, UCVM_PARAM_QUERY_MODE, zmode), 0) != 0) {
      return(1);
  }

  /* open infile, outfile */
  infp = fopen(infile, "r");
  if (infp == NULL) {
    printf("FAIL: cannot open %s\n", infile);
    return(1);
  }
  outfp = fopen(outfile, "w");
  if (outfp == NULL) {
    printf("FAIL: cannot open %s\n", outfile);
    return(1);
  }

/* process one term at a time */
  while(fgets(line, 1000, infp) != NULL) {
    if(line[0] == '#') continue; // a comment 
    if (sscanf(line,"%lf %lf %lf",
         &pt.longitude,&pt.latitude,&pt.depth) == 3) {
      if (test_assert_int(model_query(&pt, &ret, 1), 0) == 0) {
         fprintf(outfp,"%lf %lf %lf\n",ret.vs, ret.vp, ret.rho);
      }
    }
  }
  fclose(infp);
  fclose(outfp);
                
  if (test_assert_int(model_finalize(),0) != 0) {
      return(1);
  }

  return(0);
}

int runVXCVMHLABN(const char *bindir, const char *cvmdir, 
	  const char *infile, const char *outfile, int mode)
{
  char currentdir[1280];
  char runpath[1280];
  char flags[1280]="";

  sprintf(runpath, "%s/run_vx_cvmhlabn.sh", bindir);

  switch (mode) {
    case MODE_ELEVATION:
      sprintf(flags, "-z elev ");
      break;
    case MODE_DEPTH:
      sprintf(flags, "-z dep ");
      break;
    case MODE_NONE:
      sprintf(flags, "-z off ");
      break;
  }

  if(debug_mode) { strcat(flags, "-g "); }

  /* Save current directory */
  getcwd(currentdir, 1280);
  
  /* Fork process */
  pid_t pid;
  pid = fork();
  if (pid == -1) {
    perror("fork");
    return(1);
  } else if (pid == 0) {
    /* Change dir to cvmdir */
    if (chdir(bindir) != 0) {
      printf("FAIL: Error changing dir in run_vx_cvmhlabn.sh\n");
      return(1);
    }

    if (strlen(flags) == 0) {
      execl(runpath, runpath, infile, outfile, (char *)0);
    } else {
      execl(runpath, runpath, flags, infile, outfile, (char *)0);
    }

    perror("execl"); /* shall never get to here */
    printf("FAIL: CVM exited abnormally\n");
    return(1);
  } else {
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) {
      return(0);
    } else {
      printf("FAIL: CVM exited abnormally\n");
      return(1);
    }
  }

  return(0);
}


int runVXLiteCVMHLABN(const char *bindir, const char *cvmdir, 
	      const char *infile, const char *outfile,
	      int mode)
{
  char currentdir[1280];
  char flags[1280]="";

  char runpath[1280];

  sprintf(runpath, "./run_vx_lite_cvmhlabn.sh");

  sprintf(flags, "-m %s ", cvmdir);

  switch (mode) {
     case MODE_ELEVATION:
       strcat(flags, "-z elev ");
       break;
     case MODE_DEPTH:
       strcat(flags, "-z dep ");
       break;
    case MODE_NONE:
       strcat(flags, "-z off ");
       break;
  }

  if(debug_mode) { strcat(flags, "-g "); }

  /* Save current directory */
  getcwd(currentdir, 1280);
  
  /* Fork process */
  pid_t pid;
  pid = fork();
  if (pid == -1) {
    perror("fork");
    printf("FAIL: unable to fork\n");
    return(1);
  } else if (pid == 0) {

    /* Change dir to bindir */
    if (chdir(bindir) != 0) {
      printf("FAIL: Error changing dir in run_vx_lite_cvmhlabn.sh\n");
      return(1);
    }

    if (strlen(flags) == 0) {
      execl(runpath, runpath, infile, outfile, (char *)0);
    } else {
      execl(runpath, runpath, flags, infile, outfile, (char *)0);
    }
    perror("execl"); /* shall never get to here */
    printf("FAIL: CVM exited abnormally\n");
    return(1);
  } else {
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) {
      return(0);
    } else {
      printf("FAIL: CVM exited abnormally\n");
      return(1);
    }
  }

  return(0);
}
