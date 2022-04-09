/*
 * @file vx_cvmhlabn_lite_validate.c
 * @brief test with a full set of validation points
 * @author - SCEC
 * @version 1.0
 *
 * Tests the CVMHLABN as vx_lite_cvmhlabn would and recreate a
 * new input dataset that is depth based
 *
 *
 *  ./vx_cvmhlabn_lite_validate -m modeldir -f CVMHB-Los-Angeles-Basin.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "vx_sub_cvmhlabn.h"

int validate_debug = 0;
int local_query_fail_count=0;


/*********************************************/
typedef struct dat_entry_t 
{
  int x_idx;
  int y_idx;
  int z_idx;
  int vp_idx;
  int vs_idx;
  int rho_idx;
} dat_entry_t;

typedef struct dat_data_t 
{
  double x;
  double y;
  double z;
  double vp;
  double vs;
  double rho;
} dat_data_t;

dat_entry_t dat_entry;

/*
X,Y,Z,tag61_basin,vp63_basin,vs63_basin
383000.000000,3744000.000000,-15000.000000,-99999.000000,-99999.000000,-99999.000000
384000.000000,3744000.000000,-15000.000000,-99999.000000,-99999.000000,-99999.000000
*/
FILE *_process_datfile(char *fname) {

  char dat_line[1028];
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr,"VALIDATE_LITE: FAIL: Unable to open the validation data file %s\n", fname);
    exit(1);
  }
  /* read the title line */
  if (fgets(dat_line, 1028, fp) == NULL) {
    fprintf(stderr,"VALIDATE_LITE: FAIL: Unable to extract validation data file %s\n", fname);
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

// X,Y,Z,tag61_basin,vp63_basin,vs63_basin
  while(p != NULL)
  {
    if(validate_debug) { printf("VALIDATE_LITE:'%s'\n", p); }
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

//X,Y,Z,tag61_basin,vp63_basin,vs63_basin
//383000.000000,3744000.000000,-15000.000000,-99999.000000,-99999.000000,-99999.000000
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
  printf("     vx_cvmhlabn_lite_valiate - (c) SCEC\n");
  printf("Extract velocities from a simple GOCAD voxet. Accepts\n");
  printf("geographic coordinates and UTM Zone 11, NAD27 coordinates in\n");
  printf("X Y Z columns. Z is expressed as elevation offset by default.\n\n");
  printf("\tusage: vx_cvmhlabn_lite_validate [-d] -m modeldir -f file.dat\n\n");
  printf("Flags:\n");
  printf("\t-f point.dat\n\n");
  printf("\t-d enable debug/verbose mode\n\n");
  printf("\t-m path to cvmhlabn model directory \n\n");
  exit (0);
}

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc, char* const argv[]) {

	// Declare the structures.
        vx_entry_t entry;
        char modeldir[500];
        FILE *ofp = fopen("validate_lite_bad.out", "w");
        FILE *oofp = fopen("validate_lite_foo.out", "w");

        fprintf(ofp,"X,Y,Z,depth,vp63_basin,vs63_basin\n");
        fprintf(oofp,"X,Y,Z,depth,vp63_basin,vs63_basin\n");

        vx_zmode_t zmode=VX_ZMODE_ELEV;
        strcpy(modeldir,".");

        int rc;
        int opt;
        char datfile[100]="";
        dat_data_t dat;
        int tcount=0;  // total entry
        int mcount=0;  // real mismatch
        int mmcount=0; // fake mismatch -- no data
        int okcount=0;

        local_query_fail_count=0;

        /* Parse options */
        while ((opt = getopt(argc, argv, "df:m:h")) != -1) {
          switch (opt) {
          case 'f':
            strcpy(datfile, optarg);
            break;
          case 'm':
            strcpy(modeldir, optarg);
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

        /* Perform setup */
        if (vx_setup(modeldir) != 0) {
          fprintf(stderr, "Failed to init vx\n");
          exit(1);
        }

        /* Set zmode */
        vx_setzmode(zmode);

        rc=_next_datfile(fp, &dat);
        while(rc==0) {
              tcount++;
              entry.coor[0] = dat.x;
              entry.coor[1] = dat.y;
              entry.coor[2] = dat.z;

              /* In case we got anything like degrees */
              if ((entry.coor[0]<360.) && (fabs(entry.coor[1])<90)) {
                 entry.coor_type = VX_COORD_GEO;
                 } else {
                    entry.coor_type = VX_COORD_UTM;
               }

               vx_getcoord(&entry);

               if(validate_debug) fprintf(stderr, "VALIDATE_LITE:   with.. %lf %lf %lf\n", entry.coor[0],entry.coor[1],entry.coor[2]);

              if(1) { // sometimes rcc could be 1 because it can not find right surface or of area or of area
                if(validate_debug) {
                   fprintf(stderr,"VALIDATE_LITE:     vs:%lf vp:%lf rho:%lf\n\n",entry.vs, entry.vp, entry.rho);
                }

                // is result matching ?
                if(_compare_double(entry.vs, dat.vs) ||
                             _compare_double(entry.vp, dat.vp)) { 
/*** special case.. only in lr
VALIDATE_LITE:356000.000000,3754000.000000,-100.000114
VALIDATE_LITE: dat.vs(-99999.000000),dat.vp(1480.000000)
VALIDATE_LITE:   ret vs:(-1.000000) ret vp:(-1.000000)
**/
                     // okay if ( dat.vp == -99999, dat.vs== -99999 ) and (entry.vs == -1, entry.vp == -1) 
                     if (!_compare_double(entry.vs, -1.0) && !_compare_double(entry.vp, -1.0) &&
                              !_compare_double(dat.vs, -99999.0) && !_compare_double(dat.vp, -99999.0)) {
                       mmcount++;  // just -1 vs -99999
                       } else {
                         fprintf(stderr,"\nVALIDATE_LITE:Mismatching -\n");
                         fprintf(stderr,"VALIDATE_LITE:%lf,%lf,%lf\n",dat.x, dat.y, dat.z);
                         fprintf(stderr,"VALIDATE_LITE: dat.vs(%lf),dat.vp(%lf)\n",dat.vs, dat.vp);
                         fprintf(stderr,"VALIDATE_LITE:   entry vs:(%lf) entry vp:(%lf)\n",entry.vs, entry.vp);
                         mcount++;  // real mismatch
                      }
                      fprintf(ofp,"%lf,%lf,%lf,%lf,%lf,%lf\n",entry.coor[0],entry.coor[1],entry.coor[2],entry.depth,dat.vp,dat.vs);
                    } else {
                         okcount++;
                         fprintf(oofp,"%lf,%lf,%lf,%lf,%lf,%lf\n",entry.coor[0],entry.coor[1],entry.coor[2],entry.depth,dat.vp,dat.vs);
                  }
              }
          rc=_next_datfile(fp, &dat);
        }

        fprintf(stderr,"VALIDATE_LITE: %d mismatch out of %d (mmcount %d)\n", mcount, tcount,mmcount);
        fprintf(stderr,"VALIDATE_LITE: good with matching values(%d) \n",okcount );
 
       // vx_cleanup();

        fclose(fp);
        fclose(ofp);
        fclose(oofp);

	return 0;
}
