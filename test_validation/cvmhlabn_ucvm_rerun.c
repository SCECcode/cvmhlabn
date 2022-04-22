/*
 * @file cvmhlabn_ucvm_rerun.c
 * @brief rerun the failed case from cvmhlabn_ucvm_validate
 * @author - SCEC
 * @version 1.0
 *
 * Tests the CVMHLABN library by running with UCVM
 *
 *
 *  ./cvmhlabn_ucvm_rerun -e -c ucvm.conf -f validate_ucvm_debug.txt
 *
 *  test mode: query-by-depth with lat lon depth
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

#define NUM_POINTS 10000

int validate_debug = 0;
int has_negative_depth = 0;

// from gctpc
void gctp();

/*********************************************/
typedef struct dat_entry_t 
{
  int id_idx;
  int lat_idx;
  int lon_idx;
  int depth_idx;
  int vp_idx;
  int vs_idx;
} dat_entry_t;

typedef struct dat_data_t 
{
  long id;
  double lat;
  double lon;
  double depth;
  double vp;
  double vs;
} dat_data_t;

dat_entry_t dat_entry;

/*
lon,lat,depth,vp63_basin,vs63_basin
*/
FILE *_process_datfile(char *fname) {

  char dat_line[1028];
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr,"RERUN_UCVM: FAIL: Unable to open the validation data file %s\n", fname);
    exit(1);
  }

  int done=0;
  while(!done) {

  /* read the title line */
    if (fgets(dat_line, 1028, fp) == NULL) {
      fprintf(stderr,"RERUN_UCVM: FAIL: Unable to extract validation data file %s\n", fname);
      fclose(fp);
      exit(1);
    }

    if(dat_line[0]=='#') continue;
    done=1;

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
      if(validate_debug) { printf("RERUN_UCVM:'%s'\n", p); }
      if(strcmp(p,"lon")==0)
        dat_entry.lon_idx=counter;
      else if(strcmp(p,"id")==0)
        dat_entry.id_idx=counter;
      else if(strcmp(p,"lat")==0)
        dat_entry.lat_idx=counter;
      else if(strcmp(p,"depth")==0)
        dat_entry.depth_idx=counter;
      else if(strcmp(p,"vp63_basin")==0)
        dat_entry.vp_idx=counter;
      else if(strcmp(p,"vs63_basin")==0)
        dat_entry.vs_idx=counter;
      p = strtok(NULL, delimiter);
      counter++;
    }
  }
  return fp;
}
int _next_datfile(FILE *fp, dat_data_t *dat, int positive_depth) {

  char dat_line[1028];
  int done=0;

  while(!done) {
    if (fgets(dat_line, 1028, fp) == NULL) {
      return(1); 
    }
    if(dat_line[0]=='#') continue;

    char delimiter[] = ",";
    char *p = strtok(dat_line, delimiter);
    int counter=0;
  
  //lon,lat,depth,vp63_basin,vs63_basin
    while(p != NULL) {
      double val = atof(p);
      if(counter == dat_entry.lon_idx)
          dat->lon=val;
        else if (counter == dat_entry.id_idx)
          dat->id=(long) val;
        else if (counter == dat_entry.lat_idx)
          dat->lat=val;
        else if (counter == dat_entry.depth_idx)
          dat->depth=val;
        else if (counter == dat_entry.vs_idx)
          dat->vs=val;
        else if (counter == dat_entry.vp_idx)
          dat->vp=val;
      p = strtok(NULL, delimiter);
      counter++;
    }

    // only return selected line
    if(positive_depth) { // skip negative depth
      if(dat->depth < 0) { 
         has_negative_depth=1;
         continue; 
      } 
      } else {           // skip positive depth
        if(dat->depth >= 0) {
           continue;
        }
    }
    done=1;
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
  printf("     cvmhlabn_ucvm_rerun - (c) SCEC\n");
  printf("\tusage: cvmhlabn_ucvm_rerun [-d] -c ucvm.conf -f file.dat\n\n");
  printf("Flags:\n");
  printf("\t-c ucvm.conf\n\n");
  printf("\t-f point.dat\n\n");
  printf("\t-d enable debug/verbose mode\n\n");
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
        long tcount=0;  // total entry
        int mcount=0;  // real mismatch
        int mmcount=0; // fake mismatch -- no data
        int okcount=0;

        FILE *ofp= fopen("ucvm_rerun_bad.txt","w");
        fprintf(ofp,"id,lon,lat,depth,vp63_basin,vs63_basin,vp,vs\n");
        FILE *oofp= fopen("ucvm_rerun_good.txt","w");
        fprintf(oofp,"lon,lat,depth,vp63_basin,vs63_basin,vp,vs\n");

        /* Parse options */
        while ((opt = getopt(argc, argv, "edf:c:h")) != -1) {
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

      

// USING UCVM to query for points with positive depth
{
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

        if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, cmode) != UCVM_CODE_SUCCESS) {
          fprintf(stderr, "Failed to set z mode\n");
          return(1);
        }

        /* Allocate buffers */
        dat = malloc(NUM_POINTS * sizeof(dat_data_t));
        pnts = malloc(NUM_POINTS * sizeof(ucvm_point_t));
        props = malloc(NUM_POINTS * sizeof(ucvm_data_t));

        while(rc==0 && idx < NUM_POINTS) {
              memset(&(pnts[idx]), 0, sizeof(ucvm_point_t));
              memset(&(dat[idx]), 0, sizeof(dat_data_t));

              rc=_next_datfile(fp, &dat[idx], 1 /* just for positive depth */); 
              if(rc) continue;

if(validate_debug) {
fprintf(stderr,"lon %lf lat %lf dep %lf\n", dat[idx].lon, dat[idx].lat, dat[idx].depth);
}
              tcount++;

              pnts[idx].coord[0]=dat[idx].lon;
              pnts[idx].coord[1]=dat[idx].lat;
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

                     // okay if ( dat[j].vp == -99999, dat[j].vs== -99999 ) and (props[j].cmb.vs == 0, props[j].cmb.vp == 0)
                     if (!_compare_double(props[j].cmb.vs, 0) && !_compare_double(props[j].cmb.vp, 0) &&
                              !_compare_double(dat[j].vs, -99999.0) && !_compare_double(dat[j].vp, -99999.0)) {
                       mmcount++;  // just 0 vs -99999
                       fprintf(oofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].lon,dat[j].lat,dat[j].depth,
                          dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);
                       } else {
                          fprintf(ofp,"%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].id,
                          dat[j].lon,dat[j].lat,dat[j].depth,
                          dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);
                          mcount++;
                     }

                    } else {
                      okcount++;
                      fprintf(oofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].lon,dat[j].lat,dat[j].depth,
                          dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);
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
            // is result matching ?
            for(int j=0; j<idx; j++) {

if(validate_debug && j==0) {
fprintf(stderr," VV|| pnts  << coord0(%lf) coord1(%lf) coord2(%lf) \n", pnts[j].coord[0], pnts[j].coord[1], pnts[j].coord[2]);
fprintf(stderr," VV|| props << surf(%lf) vs30(%lf) depth(%lf) \n", props[j].surf, props[j].vs30, props[j].depth);
fprintf(stderr," Vv|| props << vp(%lf) vs(%lf) rho(%lf) \n", props[j].cmb.vp, props[j].cmb.vs, props[j].cmb.rho);
}

              if(_compare_double(props[j].cmb.vs, dat[j].vs) || _compare_double(props[j].cmb.vp, dat[j].vp)) {

                // okay if ( dat[j].vp == -99999, dat[j].vs== -99999 ) and (props[j].cmb.vs == 0, props[j].cmb.vp == 0)
                 if (!_compare_double(props[j].cmb.vs, 0) && !_compare_double(props[j].cmb.vp, 0) &&
                              !_compare_double(dat[j].vs, -99999.0) && !_compare_double(dat[j].vp, -99999.0)) {
                   mmcount++;  // just 0 vs -99999
                   fprintf(oofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                      dat[j].lon,dat[j].lat,dat[j].depth,
                      dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);
                   } else {
                      fprintf(ofp,"%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                      dat[j].id,dat[j].lon,dat[j].lat,dat[j].depth,
                      dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);
                      mcount++;
                 }
                } else {
                  okcount++;
                  fprintf(oofp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat[j].lon,dat[j].lat,dat[j].depth,
                          dat[j].vp,dat[j].vs,props[j].cmb.vp,props[j].cmb.vs);

              }
            }
        }
        assert(cvmhlabn_finalize() == 0);
        printf("RERUN_UCVM:Model closed successfully.\n");
        fclose(fp);
}

//USE cvmhlabn_query to process points with negative depth
if(has_negative_depth) {

        cvmhlabn_point_t pt;
        cvmhlabn_properties_t ret;

        FILE *fp=_process_datfile(datfile);

        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
          assert(model_init(envstr, "cvmhlabn")==0);
          } else {
            assert(model_init("..", "cvmhlabn")==0);
        }
      
        assert(model_setparam(0, UCVM_PARAM_QUERY_MODE, cmode)==0);
      
        int done=0;
      /* process one term at a time */
       
        while(!done) {

          dat_data_t dat;
          int rc=_next_datfile(fp, &dat, 0 /* just for negative depth */);
          if(rc) { done=1; continue; }
          tcount++;

          pt.longitude=dat.lon;
          pt.latitude=dat.lat;
          pt.depth=dat.depth;
          rc= model_query(&pt, &ret, 1);
          if(rc) {
            done=1;
            } else {
if(validate_debug) { fprintf(stderr,"%lf %lf %lf\n",ret.vs, ret.vp, ret.rho); }
              if(_compare_double(ret.vs, dat.vs) || _compare_double(ret.vp, dat.vp)) {

                 if (!_compare_double(ret.vs, 0) && !_compare_double(ret.vp, 0) &&
                              !_compare_double(dat.vs, -99999.0) && !_compare_double(dat.vp, -99999.0)) {
                   mmcount++;  // just 0 vs -99999
                   fprintf(oofp,"%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                      dat.id, dat.lon,dat.lat,dat.depth,
                      dat.vp,dat.vs,ret.vp,ret.vs);
                   } else {
                      fprintf(ofp,"%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                      dat.id, dat.lon,dat.lat,dat.depth,
                      dat.vp,dat.vs,ret.vp,ret.vs);
                      mcount++;
                 }
                } else {
                  okcount++;
                  fprintf(oofp,"%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                          dat.id, dat.lon,dat.lat,dat.depth,
                          dat.vp,dat.vs,ret.vp,ret.vs);

              }
            }
        }

          // Close the model.
        assert(model_finalize()==0);
        fclose(fp);
}


    fprintf(stderr,"RERUN_UCVM: %d mismatch out of %ld \n", mcount, tcount);
    fprintf(stderr,"RERUN_UCVM: good with matching values(%d) mmcount(%d) \n",okcount, mmcount );

    fclose(ofp);
    fclose(oofp);

     return 0;
}



