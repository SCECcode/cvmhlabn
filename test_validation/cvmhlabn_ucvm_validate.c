/*
 * @file cvmhlabn_ucvm_validate.c
 * @brief test with a full set of validation points in depth
 * @author - SCEC
 * @version 1.0
 *
 * Tests the CVMHLABN library by running with UCVM
 *
 *
 *  ./cvmhlabn_ucvm_validate -c ucvm.conf -f validate_api_good.txt
 *
 *  test 2 mode: query-by-elevation and query-by-depth
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include "ucvm.h"
#include "ucvm_utils.h"

#include "cvmhlabn.h"

#define NUM_POINTS 10

int validate_debug = 0;

/*********************************************/
typedef struct dat_entry_t 
{
  int x_idx;
  int y_idx;
  int z_idx;
  int vp_idx;
  int vs_idx;
  int rho_idx;
  int depth_idx;
} dat_entry_t;

typedef struct dat_data_t 
{
  double x;
  double y;
  double z;
  double vp;
  double vs;
  double rho;
  double depth;
} dat_data_t;

dat_entry_t dat_entry;

/*
X,Y,Z,depth,vp63_basin,vs63_basin
383000.000000,3744000.000000,-15000.000000,0.000000,-99999.000000,-99999.000000
384000.000000,3744000.000000,-15000.000000,0.000000,-99999.000000,-99999.000000
*/
FILE *_process_datfile(char *fname) {

  char dat_line[1028];
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr,"VALIDATE_UCVM: FAIL: Unable to open the validation data file %s\n", fname);
    exit(1);
  }
  /* read the title line */
  if (fgets(dat_line, 1028, fp) == NULL) {
    fprintf(stderr,"VALIDATE_UCVM: FAIL: Unable to extract validation data file %s\n", fname);
    fclose(fp);
    exit(1);
  }

  /* Strip terminating newline */
  int slen = strlen(dat_line);
  if ((slen > 0) && (dat_line[slen-1] == '\n')) {
    dat_line[slen-1] = '\0';
  }

  char delimiter[] = ",";
  char *p = strtok(dat_line, delimiter);
  int counter=0;

// X,Y,Z,depth,vp63_basin,vs63_basin
  while(p != NULL)
  {
    if(validate_debug) { printf("VALIDATE_UCVM:'%s'\n", p); }
    if(strcmp(p,"X")==0)
      dat_entry.x_idx=counter;
    else if(strcmp(p,"Y")==0)
      dat_entry.y_idx=counter;
    else if(strcmp(p,"Z")==0)
      dat_entry.z_idx=counter;
    else if(strcmp(p,"vp63_basin")==0)
      dat_entry.vp_idx=counter;
    else if(strcmp(p,"vs63_basin")==0)
      dat_entry.vs_idx=counter;
    else if(strcmp(p,"depth")==0)
      dat_entry.depth_idx=counter;
    p = strtok(NULL, delimiter);
    counter++;
  }
  return fp;
}

int _next_datfile(FILE *fp, dat_data_t *dat) {

  char dat_line[1028];
  if (fgets(dat_line, 1028, fp) == NULL) {
    return(1); 
  }

  char delimiter[] = ",";
  char *p = strtok(dat_line, delimiter);
  int counter=0;

//X,Y,Z,depth,vp63_basin,vs63_basin
//383000.000000,3744000.000000,-15000.000000,0.000000,-99999.000000,-99999.000000
  while(p != NULL) {
    double val = atof(p);
    if(counter == dat_entry.x_idx)
        dat->x=val;
      else if (counter == dat_entry.y_idx)
        dat->y=val;
      else if (counter == dat_entry.z_idx)
        dat->z=val;
      else if (counter == dat_entry.vs_idx)
        dat->vs=val;
      else if (counter == dat_entry.vp_idx)
        dat->vp=val;
      else if (counter == dat_entry.depth_idx)
        dat->depth=val;
    p = strtok(NULL, delimiter);
    counter++;
  }
  return(0);
}


/*********************************************/

int _compare_double(double f1, double f2) {
  double precision = 0.00001;
  if (((f1 - precision) < f2) && ((f1 + precision) > f2)) {
    return 0; // good
    } else {
      return 1; // bad
  }
}

/* Usage function */
void usage() {
  printf("     vx_cvmhlabn_valiate - (c) SCEC\n");
  printf("Extract velocities from a simple GOCAD voxet. Accepts\n");
  printf("geographic coordinates and UTM Zone 11, NAD27 coordinates in\n");
  printf("X Y Z depth columns. Z is expressed as elevation by default.\n\n");
  printf("\tusage: vx_cvmhlabn_validate [-d] -c ucvm.conf -f file.dat\n\n");
  printf("Flags:\n");
  printf("\t-c ucvm.conf\n\n");
  printf("\t-f point.dat\n\n");
  printf("\t-d enable debug/verbose mode\n\n");
  printf("Output format is:\n");
  printf("\tvp vs rho\n\n");
  exit (0);
}

extern char *optarg;
extern int optind, opterr, optopt;

/**
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, char* const argv[]) {
        int rc=0;
        int opt;

        char datfile[250]="";
        char configfile[250]="";
        ucvm_ctype_t cmode=UCVM_COORD_GEO_DEPTH;

        dat_data_t *dat;
        ucvm_point_t *pnts;
        ucvm_data_t *props;

        int idx=0;
        int tcount=0;  // total entry
        int mcount=0;  // real mismatch
        int mmcount=0; // fake mismatch -- no data
        int okcount=0;

        FILE *ofp= fopen("validate_ucvm_bad.txt","w");
        fprintf(ofp,"X,Y,Z,depth,vp63_basin,vs63_basin,vs,vp\n");

        /* Parse options */
        while ((opt = getopt(argc, argv, "df:c:h")) != -1) {
          switch (opt) {
          case 'c':
            strcpy(configfile, optarg);
            break;
          case 'f':
            strcpy(datfile, optarg);
            break;
          case 'd':
            validate_debug=1;
            break;
          case 'h':
            usage();
            exit(0);
            break;
          default: /* '?' */
            usage();
            exit(1);
          }
        }

      
        FILE *fp=_process_datfile(datfile);

        /* Initialize interface */
        if (ucvm_init(configfile) != UCVM_CODE_SUCCESS) {
          fprintf(stderr, "Failed to initialize UCVM API\n");
          return(1);
        }

        /* Add models */
        if (ucvm_add_model_list("cvmhlabn") != UCVM_CODE_SUCCESS) {
          fprintf(stderr, "Failed to enable model list: cvmhlabn\n");
          return(1);
        }

        /* Set z mode for depth*/
        if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, cmode) != UCVM_CODE_SUCCESS) {
          fprintf(stderr, "Failed to set z mode\n");
          return(1);
        }


        /* Allocate buffers */
        dat = malloc(NUM_POINTS * sizeof(dat_data_t));
        pnts = malloc(NUM_POINTS * sizeof(ucvm_point_t));
        props = malloc(NUM_POINTS * sizeof(ucvm_data_t));

        while(rc==0 && idx < NUM_POINTS && tcount < 3) {
              memset(&(pnts[idx]), 0, sizeof(ucvm_point_t));
              memset(&(dat[idx]), 0, sizeof(dat_data_t));
              rc=_next_datfile(fp, &dat[idx]); 
              tcount++;
              pnts[idx].coord[0]=dat[idx].x;
              pnts[idx].coord[1]=dat[idx].y;
              pnts[idx].coord[2]=dat[idx].depth;
              idx++;

              if(idx == NUM_POINTS) {
                /* Query the UCVM */
                if (ucvm_query(NUM_POINTS, pnts, props) != UCVM_CODE_SUCCESS) {
                  fprintf(stderr, "Query CVM Failed\n");
                  return(1);
                }
                // compare result
                idx=0;
                // is result matching ?
                for(int j=0; j<NUM_POINTS; j++) {
                  if(_compare_double(props[j].cmb.vs, dat[j].vs) || _compare_double(props[j].cmb.vp, dat[j].vp)) { 
                    fprintf(ofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].x,dat[j].y,dat[j].z,dat[j].depth,dat[j].vp,dat[j].vs,props[j].cmb.vs,props[j].cmb.vp);
                    } else {
                      okcount++;
                  }
                }
                idx=0;
              }
        }
        if(idx > 0) {
            /* Query the UCVM */
            if (ucvm_query(idx, pnts, props) != UCVM_CODE_SUCCESS) {
              fprintf(stderr, "Query CVM Failed\n");
              return(1);
            }
            // compare result
            idx=0;
            // is result matching ?
            for(int j=0; j<idx; j++) {
              if(_compare_double(props[j].cmb.vs, dat[j].vs) || _compare_double(props[j].cmb.vp, dat[j].vp)) {
                fprintf(ofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].x,dat[j].y,dat[j].z,dat[j].depth,dat[j].vp,dat[j].vs,props[j].cmb.vs,props[j].cmb.vp);
                } else {
                  okcount++;
              }
            }
        }

        fprintf(stderr,"VALIDATE_UCVM: %d mismatch out of %d \n", mcount, tcount);
        fprintf(stderr,"VALIDATE_UCVM: good with matching values(%d) mmcount(%d) \n",okcount, mmcount );
	assert(cvmhlabn_finalize() == 0);
	printf("VALIDATE_UCVM:Model closed successfully.\n");

        fclose(fp);
        fclose(ofp);

	return 0;
}



