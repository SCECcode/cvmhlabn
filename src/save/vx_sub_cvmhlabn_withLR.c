/** vx_sub_cvmhlabn.c - Query interface to GoCAD voxet volumes and GTL. Supports
    queries for material properties and topography. Accepts Geographic Coordinates 
    or UTM Zone 11 coordinates.

01/2010: PES: Derived from original VX interface, vx.c. 
              Added Vs30 Derived GTL, 1D background, smoothing
07/2011: PES: Extracted io into separate module from vx_sub.c
**/

/*
b.dat
X,Y,Z,vp63_basin,vs63_basin
421000.000000,3712000.000000,-1906.000000,-1.0000,-1.0000
354000.000000,3756000.000000,-2000.000366,3326.486084,1695.246582
g.dat
X,Y,Z,vp63_basin,vs63_basin
359000.000000,3766000.000000,-100.000114,1710.692505,302.321136
s.dat
X,Y,Z,vp63_basin,vs63_basin
359000.000000,3766000.000000,-100.000114,1710.692505,302.321136
354000.000000,3756000.000000,-2000.000366,3326.486084,1695.246582
464000.000000,3662000.000000,-1900.000000,3438.379395,1791.699829

bad-in
354000.000000 3756000.000000-2000.000366
./vx_lite -z elev -m ../data/cvmhlabn  < bad-in

./vx_cvmhlabn_validate  -z elev -f b.dat
#./vx_cvmhlabn_validate  -z elev -f g.dat

#./vx_cvmhlabn_validate  -z elev -f s.dat
#./vx_cvmhlabn_validate  -z elev -f ss.dat
#./vx_cvmhlabn_validate -z elev -f /Users/mei/scec/UCVM_MODELS/Harvard/WORKSPACE/cvmhlabn/CVMHB-Los-Angeles-Basin.dat

interesting set..
355000.000000 3755000.000000 -100.000114
*/

int MEI=1 ;
int LABN=1; 

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "params.h"
#include "coor_para.h"
#include "voxet.h"
#include "proj.h"
#include "cproj.h"
#include "utils.h"
#include "vx_io.h"
#include "vx_sub_cvmhlabn.h"

/* Smoothing parameters for SCEC 1D */
#define SCEC_SMOOTH_DIST 50.0 // km

/* Topography filtering DEPRECATED */
#define ELEV_EPSILON 0.01
#define MAX_ITER_ELEV 4

float COOR_UTM2=-99999.0;
 

int cvmhlabn_debug=0;
int _debug=1;
int surface_nodata_count=0; // tracing how many NO_DATA_VALUE

float p0_NO_DATA_VALUE = -99999.0; // p0.NO_DATA_VALUE 
//float p0_NO_DATA_VALUE = 0.0; // p0.NO_DATA_VALUE 
int p0_ESIZE = 4; // p0.ESIZE

/* Function declarations */
void gctp();
int voxbytepos(int *, int* ,int);
double calc_rho(float vp, vx_src_t data_src);

/* User-defined background model function pointer */
int (*callback_bkg)(vx_entry_t *entry, vx_request_t req_type) = NULL;

/* Model state variables */
static int is_setup = False;
vx_zmode_t vx_zmode = VX_ZMODE_ELEV; // default
struct axis lr_a, mr_a, hr_a, to_a;

int pkey=0;
// ?? struct property p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13;
struct property p_vp_hr; // HR p2
struct property p_tag_hr; // HR p9
struct property p_vs_hr; // HR p12 
struct property p_topo_dem; // interface p4 
struct property p_moho; // interface p5
struct property p_base; // interface p6
struct property p_model_top; // interface p13
struct property p_vp_lr; // LR p0
struct property p_tag_lr; // LR p7
struct property p_vs_lr; // LR p11

float step_to[3], step_hr[3], step_lr[3];
float step0hr, step1hr, step2hr;

/* Model buffers */
static char *lrbuffer = NULL;
static char *hrbuffer = NULL;
static char *tobuffer = NULL;
static char *mobuffer = NULL;
static char *babuffer = NULL;
static char *mtopbuffer = NULL;
static char *lrtbuffer = NULL;
static char *hrtbuffer = NULL;
static char *lrvsbuffer = NULL;
static char *hrvsbuffer = NULL;

/* Data source labels */
char *VX_SRC_NAMES[7] = {"nr", "hr", "lr", "cm", "to", "bk", "gt"};


/* Setup function to be called prior to querying points */
int vx_setup(const char *data_dir)
{
  int NCells;
  int n;

  /* zero-out inparm for gctpc */
  for(n=0;n>15;n++) inparm[n]=0;

  /* Initialize buffer pointers to NULL */
  hrbuffer = NULL;
  tobuffer = mobuffer = babuffer = mtopbuffer = NULL;
  hrtbuffer = NULL;
  hrvsbuffer = NULL;

  char LR_PAR[CMLEN];
  sprintf(LR_PAR, "%s/CVM_LR.vo", data_dir);
  
  char HR_PAR[CMLEN];
if(LABN) {
sprintf(HR_PAR, "%s/CVMHB-Los-Angeles-Basin.vo", data_dir);
} else {
sprintf(HR_PAR, "%s/CVM_HR.vo", data_dir);
}
  
  char CM_PAR[CMLEN];
  sprintf(CM_PAR, "%s/CVM_CM.vo", data_dir);

  char TO_PAR[CMLEN];
  sprintf(TO_PAR, "%s/interfaces.vo", data_dir);

  /**** First we load the LowRes File****/
  if (vx_io_init(LR_PAR) != 0) {
    fprintf(stderr, "Failed to load LR param file %s. Check that the model path is correct.\n", LR_PAR);
    return(1);
  }

if(_debug) fprintf(stderr, "load LR param file %s. \n", LR_PAR);

  vx_io_getvec("AXIS_O",lr_a.O);
  vx_io_getvec("AXIS_U",lr_a.U);
  vx_io_getvec("AXIS_V",lr_a.V);
  vx_io_getvec("AXIS_W",lr_a.W);
  vx_io_getvec("AXIS_MIN",lr_a.MIN);
  vx_io_getvec("AXIS_MAX",lr_a.MAX);
  vx_io_getdim("AXIS_N",lr_a.N);

  /** AP: AXIS_MIN and AXIS_MAX are currently not used and need to be 0
      0 0 and 1 1 1, respectively, in the .vo file. The AXIS_UVW would
      need to be adjusted accordingly in the .vo file.
  **/

  NCells=lr_a.N[0]*lr_a.N[1]*lr_a.N[2];
  sprintf(p_vp_lr.NAME,"vint");
  vx_io_getpropname("PROP_FILE",1,p_vp_lr.FN);
  vx_io_getpropsize("PROP_ESIZE",1,&p_vp_lr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",1,&p_vp_lr.NO_DATA_VALUE);

if(_debug) { fprintf(stderr,"using LR VP file ..%s\n\n",p_vp_lr.FN); }

  lrbuffer=(char *)malloc(NCells*p_vp_lr.ESIZE);
  if (lrbuffer == NULL) {
    fprintf(stderr, "Failed to allocate LR Vp buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_vp_lr.FN,
		       p_vp_lr.ESIZE,NCells,lrbuffer) != 0) {
    fprintf(stderr, "Failed to load LR Vp volume\n");
    return(1);
  }

  // and the tags
  sprintf(p_tag_lr.NAME,"tag");
  vx_io_getpropname("PROP_FILE",2,p_tag_lr.FN);
  vx_io_getpropsize("PROP_ESIZE",2,&p_tag_lr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",2,&p_tag_lr.NO_DATA_VALUE);

if(_debug) { fprintf(stderr,"using LR TAG file ..%s\n\n",p_tag_lr.FN); }

  lrtbuffer=(char *)malloc(NCells*p_tag_lr.ESIZE);
  if (lrtbuffer == NULL) {
    fprintf(stderr, "Failed to allocate LR tag buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_tag_lr.FN,
		       p_tag_lr.ESIZE,NCells,lrtbuffer) != 0) {
    fprintf(stderr, "Failed to load LR tag volume\n");
    return(1);
  }

  // and vs
  sprintf(p_vs_lr.NAME,"vs");
  vx_io_getpropname("PROP_FILE",3,p_vs_lr.FN);
  vx_io_getpropsize("PROP_ESIZE",3,&p_vs_lr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",3,&p_vs_lr.NO_DATA_VALUE);

if(_debug) { fprintf(stderr,"using LR VS file ..%s\n\n",p_vs_lr.FN); }

  lrvsbuffer=(char *)malloc(NCells*p_vs_lr.ESIZE);
  if (lrvsbuffer == NULL) {
    fprintf(stderr, "Failed to allocate LR Vs buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_vs_lr.FN,
		       p_vs_lr.ESIZE,NCells,lrvsbuffer) != 0) {
    fprintf(stderr, "Failed to load LR Vs volume\n");
    return(1);
  }
  vx_io_finalize();


  /**** First we load the HighRes File *****/
  if (vx_io_init(HR_PAR) != 0) {
    fprintf(stderr, "Failed to load HR param file %s\n", HR_PAR);
    return(1);
  }

  vx_io_getvec("AXIS_O",hr_a.O);
  vx_io_getvec("AXIS_U",hr_a.U);
  vx_io_getvec("AXIS_V",hr_a.V);
  vx_io_getvec("AXIS_W",hr_a.W);
  vx_io_getvec("AXIS_MIN",hr_a.MIN);
  vx_io_getvec("AXIS_MAX",hr_a.MAX);
  vx_io_getdim("AXIS_N ",hr_a.N);

  if(_debug) {
    fprintf(stderr," From:\n");
    fprintf(stderr,"    hr_a.O %.0f %.0f %.0f\n", hr_a.O[0],hr_a.O[1],hr_a.O[2]);
    fprintf(stderr,"    hr_a.U %.0f %.0f %.0f\n", hr_a.U[0],hr_a.U[1],hr_a.U[2]);
    fprintf(stderr,"    hr_a.V %.0f %.0f %.0f\n", hr_a.V[0],hr_a.V[1],hr_a.V[2]);
    fprintf(stderr,"    hr_a.W %.0f %.0f %.0f\n", hr_a.W[0],hr_a.W[1],hr_a.W[2]);
    fprintf(stderr,"    hr_a.MAX %f %f %f\n", hr_a.MAX[0],hr_a.MAX[1],hr_a.MAX[2]);
    fprintf(stderr,"    hr_a.MIN %f %f %f\n", hr_a.MIN[0],hr_a.MIN[1],hr_a.MIN[2]);
    fprintf(stderr,"    hr_a.N %d %d %d\n", hr_a.N[0],hr_a.N[1],hr_a.N[2]);
  }

  if( MEI || (hr_a.MIN[0] != 0) || (hr_a.MAX[0] != 1)) {
/* origin */
      float origin0hr= hr_a.O[0] + hr_a.U[0] * hr_a.MIN[0];
      float origin1hr= hr_a.O[1] + hr_a.V[1] * hr_a.MIN[1];
      float origin2hr= hr_a.O[2] + hr_a.W[2] * hr_a.MIN[2];
/* Umax */
      float umax0hr=hr_a.O[0] + (hr_a.U[0] * hr_a.MAX[0]);
      float umax1hr=origin1hr;
      float umax2hr=origin2hr;
/* Vmax */
      float vmax0hr=origin0hr;
      float vmax1hr=hr_a.O[1] + (hr_a.V[1] * hr_a.MAX[1]);
      float vmax2hr=origin2hr;
/* Wmax */
      float wmax0hr=origin0hr;
      float wmax1hr=origin1hr;
      float wmax2hr=hr_a.O[2] + (hr_a.W[2] * hr_a.MAX[2]);
/* cell size */
      step0hr=((hr_a.MAX[0] - hr_a.MIN[0]) * hr_a.U[0]) / (hr_a.N[0]-1); 
      step1hr=((hr_a.MAX[1] - hr_a.MIN[1]) * hr_a.V[1]) / (hr_a.N[1]-1); 
      step2hr=((hr_a.MAX[2] - hr_a.MIN[2]) * hr_a.W[2]) / (hr_a.N[2]-1); 

      hr_a.O[0]=origin0hr; hr_a.O[1]=origin1hr; hr_a.O[2]=origin2hr;
      hr_a.U[0]=umax0hr; hr_a.U[1]=umax1hr; hr_a.U[2]=umax2hr;
      hr_a.V[0]=vmax0hr; hr_a.V[1]=vmax1hr; hr_a.V[2]=vmax2hr;
      hr_a.W[0]=wmax0hr; hr_a.W[1]=wmax1hr; hr_a.W[2]=wmax2hr;

      if(_debug) {
        fprintf(stderr,">>>Info: newly  calculated HR ------\n");
        fprintf(stderr," origin %.0f %.0f %.0f\n", hr_a.O[0],hr_a.O[1],hr_a.O[2]);
        fprintf(stderr," umax %.0f %.0f %.0f\n", hr_a.U[0],hr_a.U[1],hr_a.U[2]);
        fprintf(stderr," vmax %.0f %.0f %.0f\n", hr_a.V[0],hr_a.V[1],hr_a.V[2]);
        fprintf(stderr," wmax %.0f %.0f %.0f\n", hr_a.W[0],hr_a.W[1],hr_a.W[2]);
        fprintf(stderr," step %f %f %f\n\n", step0hr, step1hr, step2hr);
      }

      } else { // get the info out
        if(_debug) {
          step0hr=((hr_a.MAX[0] - hr_a.MIN[0]) * hr_a.U[0]) / (hr_a.N[0]-1);
          step1hr=((hr_a.MAX[1] - hr_a.MIN[1]) * hr_a.V[1]) / (hr_a.N[1]-1);
          step2hr=((hr_a.MAX[2] - hr_a.MIN[2]) * hr_a.W[2]) / (hr_a.N[2]-1);
          fprintf(stderr,">>>Info: read HR------\n");
          fprintf(stderr," origin %.0f %.0f %.0f\n", hr_a.O[0],hr_a.O[1],hr_a.O[2]);
          fprintf(stderr," umax %.0f %.0f %.0f\n", hr_a.U[0],hr_a.U[1],hr_a.U[2]);
          fprintf(stderr," vmax %.0f %.0f %.0f\n", hr_a.V[0],hr_a.V[1],hr_a.V[2]);
          fprintf(stderr," wmax %.0f %.0f %.0f\n", hr_a.W[0],hr_a.W[1],hr_a.W[2]);
          fprintf(stderr," step %f %f %f\n\n", step0hr, step1hr, step2hr);
      }

   }

  NCells=hr_a.N[0]*hr_a.N[1]*hr_a.N[2];
if(LABN) {
sprintf(p_vp_hr.NAME,"vp63_basin");
} else { 
sprintf(p_vp_hr.NAME,"vp63");
}
  pkey=vx_io_getpropkey(p_vp_hr.NAME);

  vx_io_getpropname("PROP_FILE",pkey,p_vp_hr.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_vp_hr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_vp_hr.NO_DATA_VALUE);

  if(_debug) { fprintf(stderr,"using HR VP file ..%s\n\n",p_vp_hr.FN); }

  hrbuffer=(char *)malloc(NCells*p_vp_hr.ESIZE);
  if (hrbuffer == NULL) {
    fprintf(stderr, "Failed to allocate HR Vp file\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_vp_hr.FN,
		       p_vp_hr.ESIZE,NCells,hrbuffer) != 0) {
    fprintf(stderr, "Failed to load HR Vp volume\n");
    return(1);
  }

  // and the tags
if(LABN) {
sprintf(p_tag_hr.NAME,"tag61_basin");
} else {
sprintf(p_tag_hr.NAME,"tag");
}

  pkey=vx_io_getpropkey(p_tag_hr.NAME);

  vx_io_getpropname("PROP_FILE",pkey,p_tag_hr.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_tag_hr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_tag_hr.NO_DATA_VALUE);

if(_debug) {fprintf(stderr,"using HR TAG file..%s\n\n",p_tag_hr.FN); }

  hrtbuffer=(char *)malloc(NCells*p_tag_hr.ESIZE);
  if (hrtbuffer == NULL) {
    fprintf(stderr, "Failed to allocate HR tag file\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_tag_hr.FN,
		       p_tag_hr.ESIZE,NCells,hrtbuffer) != 0) {
    fprintf(stderr, "Failed to load HR tag volume\n");
    return(1);
  }

  // and vs
if(LABN) {
sprintf(p_vs_hr.NAME,"vs63_basin");
} else {
sprintf(p_vs_hr.NAME,"vs63");
}
  pkey=vx_io_getpropkey(p_vs_hr.NAME);

  vx_io_getpropname("PROP_FILE",pkey,p_vs_hr.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_vs_hr.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_vs_hr.NO_DATA_VALUE);

if(_debug) {fprintf(stderr,"using HR VS file..%s\n\n",p_vs_hr.FN); }

  hrvsbuffer=(char *)malloc(NCells*p_vs_hr.ESIZE);
  if (hrvsbuffer == NULL) {
    fprintf(stderr, "Failed to allocate HR Vs file\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_vs_hr.FN,
		       p_vs_hr.ESIZE,NCells,hrvsbuffer) != 0) {
    fprintf(stderr, "Failed to load HR Vs volume\n");
    return(1);
  }

  vx_io_finalize();

  /**** Now we load the topo, moho, base, model top File *****/
  if (vx_io_init(TO_PAR) != 0) {
    fprintf(stderr, "Failed to load topo param file %s\n", TO_PAR);
    return(1);
  }

  vx_io_getvec("AXIS_O",to_a.O);
  vx_io_getvec("AXIS_U",to_a.U);
  vx_io_getvec("AXIS_V",to_a.V);
  vx_io_getvec("AXIS_W",to_a.W);
  vx_io_getvec("AXIS_MIN",to_a.MIN);
  vx_io_getvec("AXIS_MAX",to_a.MAX);
  vx_io_getdim("AXIS_N ",to_a.N);

  NCells=to_a.N[0]*to_a.N[1]*to_a.N[2];

  // topo
  sprintf(p_topo_dem.NAME,"topo_dem");
  pkey=vx_io_getpropkey(p_topo_dem.NAME);

  vx_io_getpropname("PROP_FILE",pkey,p_topo_dem.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_topo_dem.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_topo_dem.NO_DATA_VALUE);

  tobuffer=(char *)malloc(NCells*p_topo_dem.ESIZE);
  if (tobuffer == NULL) {
    fprintf(stderr, "Failed to allocate topo dem buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_topo_dem.FN, 
		       p_topo_dem.ESIZE, NCells, tobuffer) != 0) {
    fprintf(stderr, "Failed to load topo volume\n");
    return(1);
  }

  // moho
  sprintf(p_moho.NAME,"moho");
  pkey=vx_io_getpropkey(p_moho.NAME);
  vx_io_getpropname("PROP_FILE",pkey,p_moho.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_moho.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_moho.NO_DATA_VALUE);

  mobuffer=(char *)malloc(NCells*p_moho.ESIZE);
  if (mobuffer == NULL) {
    fprintf(stderr, "Failed to allocate topo moho buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_moho.FN, 
		       p_moho.ESIZE, NCells, mobuffer) != 0) {
    fprintf(stderr, "Failed to load moho volume\n");
    return(1);
  }

  // basement
  sprintf(p_base.NAME,"base");
  pkey=vx_io_getpropkey(p_base.NAME);
  vx_io_getpropname("PROP_FILE",pkey,p_base.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_base.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_base.NO_DATA_VALUE);

  babuffer=(char *)malloc(NCells*p_base.ESIZE);
  if (babuffer == NULL) {
    fprintf(stderr, "Failed to allocate topo basement buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_base.FN, 
		       p_base.ESIZE, NCells, babuffer) != 0) {
    fprintf(stderr, "Failed to load basement volume\n");
    return(1);
  }

  // top elevation of model
  sprintf(p_model_top.NAME,"model_top");
  pkey=vx_io_getpropkey(p_model_top.NAME);

  vx_io_getpropname("PROP_FILE",pkey,p_model_top.FN);
  vx_io_getpropsize("PROP_ESIZE",pkey,&p_model_top.ESIZE);
  vx_io_getpropval("PROP_NO_DATA_VALUE",pkey,&p_model_top.NO_DATA_VALUE);

  mtopbuffer=(char *)malloc(NCells*p_model_top.ESIZE);
  if (mtopbuffer == NULL) {
    fprintf(stderr, "Failed to allocate topo model_top buffer\n");
    return(1);
  }
  if (vx_io_loadvolume(data_dir, p_model_top.FN, 
		       p_model_top.ESIZE, NCells, mtopbuffer) != 0) {
    fprintf(stderr, "Failed to load model_top volume\n");
    return(1);
  }

  vx_io_finalize();

  // compute steps
  step_to[0]=to_a.U[0]/(to_a.N[0]-1);
  step_to[1]=to_a.V[1]/(to_a.N[1]-1);
  step_to[2]=0.0;

  step_lr[0]=lr_a.U[0]/(lr_a.N[0]-1);
  step_lr[1]=lr_a.V[1]/(lr_a.N[1]-1);
  step_lr[2]=lr_a.W[2]/(lr_a.N[2]-1);

  if( MEI || (hr_a.MIN[0] != 0) || (hr_a.MAX[0] != 1)) {
      step_hr[0]=step0hr;
      step_hr[1]=step1hr;
      step_hr[2]=step2hr;
      } else {
          step_hr[0]=hr_a.U[0]/(hr_a.N[0]-1);
          step_hr[1]=hr_a.V[1]/(hr_a.N[1]-1);
          step_hr[2]=hr_a.W[2]/(hr_a.N[2]-1);
  }
  
  is_setup = True;

  return(0);
}


/* Cleanup function to free resources and restore state */
int vx_cleanup()
{
  if (!is_setup) {
    return(1);
  }

  free(hrbuffer);
  free(lrbuffer);
  free(tobuffer);
  free(mobuffer);
  free(babuffer);
  free(mtopbuffer);

  free(hrtbuffer);
  free(hrvsbuffer);
  free(lrtbuffer);
  free(lrvsbuffer);

  vx_zmode = VX_ZMODE_ELEV;

  return(0);
}


/* Set query mode: elevation, elevation offset, depth */
int vx_setzmode(vx_zmode_t m) {
  vx_zmode = m;
  return(0);
}


/* Query material properties and topography at desired point. 
   Coordinates may be Geo or UTM */
int vx_getcoord(vx_entry_t *entry) {
  return(vx_getcoord_private(entry, True));
}


/* Private query function for material properties. Allows caller to 
   disable advanced features like depth/offset query modes.
*/ 
int vx_getcoord_private(vx_entry_t *entry, int enhanced) {
if(cvmhlabn_debug) fprintf(stderr,"CALLING --- vx_getcoord_private (enhanced %d)\n",enhanced);
  int j=0;
  double SP[2],SPUTM[2];
  int gcoor[3];
  // fall into bkg
  int do_bkg = False;
  float surface;
  double elev, depth, zt, topo_gap;
  double incoor[3];

  /* Initialize variables */
  elev = 0.0;
  depth = 0.0;
  zt = 0.0;
  topo_gap = 0.0;
  gcoor[0]=0;
  gcoor[1]=0;
  gcoor[2]=0;


  /* Proceed only if setup has been performed */
  if ((entry == NULL) || (is_setup != True)) {
    return(1);
  }

  /* Make copy of original input coordinates */
  memcpy(incoor, entry->coor, sizeof(double) * 3);

  /* Initialize entry structure */
  vx_init_entry(entry);

  /* Generate UTM coords */
  switch (entry->coor_type) {
  case VX_COORD_GEO:
    
    SP[0]=entry->coor[0];
    SP[1]=entry->coor[1];
    
    gctp(SP,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,efile,
	 SPUTM,&outsys,&outzone,inparm,&outunit,&outdatum,
	 file27, file83,&iflg);
    
    entry->coor_utm[0]=SPUTM[0];
    entry->coor_utm[1]=SPUTM[1];
    entry->coor_utm[2]=entry->coor[2];
    break;
  case VX_COORD_UTM:
    entry->coor_utm[0]=entry->coor[0];
    entry->coor_utm[1]=entry->coor[1];
    entry->coor_utm[2]=entry->coor[2];
    break;
  default:
    return(1);
    break;
  }

  if(enhanced) {  // park the Z utm..
     COOR_UTM2=entry->coor_utm[2];
     } else {
     COOR_UTM2=-99999.0;
  }
  if(_debug) { fprintf(stderr," UTM -> %lf %lf %lf\n", entry->coor_utm[0], entry->coor_utm[1], entry->coor_utm[2]);
}

  /* Now we have UTM Zone 11 */
  /*** Prevent all to obvious bad coordinates from being displayed */
  if (entry->coor_utm[1] < 10000000) {
     
    // we start with the elevations; the voxet does not have a vertical 
    // dimension
    gcoor[0]=round((entry->coor_utm[0]-to_a.O[0])/step_to[0]);
    gcoor[1]=round((entry->coor_utm[1]-to_a.O[1])/step_to[1]);
    gcoor[2]=0;

if(_debug){ fprintf(stderr,"CHECKING with elevation part %d %d %d\n", gcoor[0], gcoor[1], gcoor[2]); }
    
    //check if inside
    if(gcoor[0]>=0&&gcoor[1]>=0&&
       gcoor[0]<to_a.N[0]&&gcoor[1]<to_a.N[1]) {	      
      entry->elev_cell[0]= to_a.O[0]+gcoor[0]*step_to[0];
      entry->elev_cell[1]= to_a.O[1]+gcoor[1]*step_to[1];
      j=voxbytepos(gcoor,to_a.N,p_topo_dem.ESIZE);
      memcpy(&(entry->topo), &tobuffer[j], p_topo_dem.ESIZE);
      memcpy(&(entry->mtop), &mtopbuffer[j], p_topo_dem.ESIZE);
      memcpy(&(entry->base), &babuffer[j], p_topo_dem.ESIZE);
      memcpy(&(entry->moho), &mobuffer[j], p_topo_dem.ESIZE);
      if (((entry->topo - p0_NO_DATA_VALUE < 0.1) || 
	   (entry->mtop - p0_NO_DATA_VALUE < 0.1))) {
	do_bkg = True;
      }
    } else {
      do_bkg = True;
    }

    /* Convert depth/offset Z coordinate to elevation */
    if (enhanced == True) {
      elev = entry->coor_utm[2];
      vx_getsurface(entry->coor, entry->coor_type, &surface);
      if(_debug) { fprintf(stderr," cvmh surface -- %lf\n", surface); }
      if (surface < -90000.0) {
        surface_nodata_count++;
	return(1);
      }
      switch (vx_zmode) {
      case VX_ZMODE_ELEV:
	break;
      case VX_ZMODE_DEPTH:
	entry->coor[2] = surface - elev;
	entry->coor_utm[2] = entry->coor[2];
	break;
      case VX_ZMODE_ELEVOFF:
	entry->coor[2] = surface + elev;
	entry->coor_utm[2] = entry->coor[2];
	break;
      default:
	return(1);
	break;
      }
      depth = surface - entry->coor_utm[2];
    }

if(cvmhlabn_debug) fprintf(stderr,"\n  READY  \n");
if(cvmhlabn_debug) { fprintf(stderr,"Looking into (LAYER)>>>>>> entry->coor(%lf %lf %lf)\n",
                                                        entry->coor[0], entry->coor[1], entry->coor[2]); }
    if ((do_bkg == False) || (enhanced == False)) {
      /* AP: this calculates the cell numbers from the coordinates and 
	 the grid spacing. The -1 is necessary to do the counting 
	 correctly. The rounding is necessary because the data are cell 
	 centered, eg. they are valid half a cell width away from the 
	 data point */

      /* Extract vp/vs */      
      gcoor[0]=round((entry->coor_utm[0]-hr_a.O[0])/step_hr[0]);
      gcoor[1]=round((entry->coor_utm[1]-hr_a.O[1])/step_hr[1]);
      gcoor[2]=round((entry->coor_utm[2]-hr_a.O[2])/step_hr[2]);
      
if(cvmhlabn_debug) fprintf(stderr,"is this in HR? %d %d %d utm(%f %f %f) \n", 
gcoor[0], gcoor[1], gcoor[2] ,entry->coor_utm[0], entry->coor_utm[1], entry->coor_utm[2]);

      if(gcoor[0]>=0&&gcoor[1]>=0&&gcoor[2]>=0&&
	 gcoor[0]<hr_a.N[0]&&gcoor[1]<hr_a.N[1]&&gcoor[2]<hr_a.N[2]) {

if(cvmhlabn_debug) fprintf(stderr,"CHECKING into HR with %d %d %d\n", gcoor[0], gcoor[1], gcoor[2]); 

	/* AP: And here are the cell centers*/
	entry->vel_cell[0]= hr_a.O[0]+gcoor[0]*step_hr[0];
	entry->vel_cell[1]= hr_a.O[1]+gcoor[1]*step_hr[1];
	entry->vel_cell[2]= hr_a.O[2]+gcoor[2]*step_hr[2];

if(_debug) { fprintf(stderr,"  >with entry_vel_cell, %f %f %f\n", entry->vel_cell[0], entry->vel_cell[1], entry->vel_cell[2]); }
	j=voxbytepos(gcoor,hr_a.N,p_vp_hr.ESIZE);
	memcpy(&(entry->provenance), &hrtbuffer[j], p0_ESIZE);
	memcpy(&(entry->vp), &hrbuffer[j], p_vp_hr.ESIZE);
	memcpy(&(entry->vs), &hrvsbuffer[j], p_vp_hr.ESIZE);

        if(entry->vs == -99999.0 &&  entry->vp == -99999.0) { 
           if(cvmhlabn_debug) { fprintf(stderr,"  >FOUND IN HR but NODATA>>>>>> j(%d) gcoor(%d %d %d) vp(%f) vs(%f)\n",j, gcoor[0], gcoor[1], gcoor[2], entry->vp, entry->vs); }
//MEI, still set it ??
	  entry->data_src = VX_SRC_HR;

          } else {
	  entry->data_src = VX_SRC_HR;
if(cvmhlabn_debug) { fprintf(stderr,"  >DONE(In HR)>>>>>> j(%d) gcoor(%d %d %d) vp(%f) vs(%f)\n",j, gcoor[0], gcoor[1], gcoor[2], entry->vp, entry->vs); }
        }
      }

      if( entry->data_src != VX_SRC_HR ) { // try LR region	  

        gcoor[0]=round((entry->coor_utm[0]-lr_a.O[0])/step_lr[0]);
        gcoor[1]=round((entry->coor_utm[1]-lr_a.O[1])/step_lr[1]);
        gcoor[2]=round((entry->coor_utm[2]-lr_a.O[2])/step_lr[2]);

if(cvmhlabn_debug) fprintf(stderr,"is this in LR? %d %d %d utm(%f %f %f) \n", 
gcoor[0], gcoor[1], gcoor[2] ,entry->coor_utm[0], entry->coor_utm[1], entry->coor_utm[2]);

        if(gcoor[0]>=0&&gcoor[1]>=0&&gcoor[2]>=0&&
           gcoor[0]<lr_a.N[0]&&gcoor[1]<lr_a.N[1]&&gcoor[2]<lr_a.N[2]) {
          /* AP: And here are the cell centers*/
if(cvmhlabn_debug) fprintf(stderr,"CHECKING with LR with %d %d %d\n", gcoor[0], gcoor[1], gcoor[2]); 
          entry->vel_cell[0]= lr_a.O[0]+gcoor[0]*step_lr[0];
          entry->vel_cell[1]= lr_a.O[1]+gcoor[1]*step_lr[1];
          entry->vel_cell[2]= lr_a.O[2]+gcoor[2]*step_lr[2];
          j=voxbytepos(gcoor,lr_a.N,p_vp_lr.ESIZE);
          memcpy(&(entry->provenance), &lrtbuffer[j], p_vp_lr.ESIZE);
          memcpy(&(entry->vp), &lrbuffer[j], p_vp_lr.ESIZE);
          memcpy(&(entry->vs), &lrvsbuffer[j], p_vp_lr.ESIZE);
          entry->data_src = VX_SRC_LR;
if(cvmhlabn_debug) { fprintf(stderr,"  >Looking DONE(In LR)>>>>>> j(%d) gcoor(%d %d %d) vp(%f) vs(%f)\n",j, gcoor[0], gcoor[1], gcoor[2], entry->vp, entry->vs); }
        }
      }
      if( entry->data_src != VX_SRC_HR && entry->data_src != VX_SRC_LR) { // must be in bkg
        do_bkg = True;
      }
    }

//--  if it turns out to be do_bkg, that means return NO DATA 
    if ((enhanced == True) && (do_bkg == True)) {
      memcpy(entry->coor, incoor, sizeof(double) * 3);
      return(1);
    } else {
      /* Compute rho */
      entry->rho = calc_rho(entry->vp, entry->data_src);
    }
  }

  if(_debug) { fprintf(stderr,"   DONE(ALL)>>>>>> j(%d) gcoor(%d %d %d) vp(%f) vs(%f) rho(%f)\n",j, gcoor[0], gcoor[1], gcoor[2], entry->vp, entry->vs, entry->rho); }

  if(0) {
    fprintf(stderr,"KEEP warnings down: topo_gap(%lf) zt(%lf) depth(%lf)\n",topo_gap,zt,depth);
  }

  /* Restore original input coords */
  memcpy(entry->coor, incoor, sizeof(double) * 3);
if(cvmhlabn_debug) fprintf(stderr,"DONE --- vx_getcoord_private\n");
  return(0);
}

/* Get raw voxel information at the supplied voxel volume coordinates */
void vx_getvoxel(vx_voxel_t *voxel) {
  int gcoor[3];
  int j;

  /* Proceed only if setup has been performed */
  if ((voxel == NULL) || (is_setup != True)) {
    return;
  }

  // Initialize entry structure
  vx_init_voxel(voxel);

  gcoor[0] = voxel->coor[0];
  gcoor[1] = voxel->coor[1];
  gcoor[2] = voxel->coor[2];
	  
  switch (voxel->data_src) {
  case VX_SRC_TO:
    if(gcoor[0]>=0&&gcoor[1]>=0&&
       gcoor[0]<to_a.N[0]&&gcoor[1]<to_a.N[1])
      {
	voxel->elev_cell[0]= to_a.O[0]+gcoor[0]*step_to[0];
        voxel->elev_cell[1]= to_a.O[1]+gcoor[1]*step_to[1];
	j=voxbytepos(gcoor,to_a.N,p_topo_dem.ESIZE);
	memcpy(&(voxel->topo), &tobuffer[j], p_topo_dem.ESIZE);
	memcpy(&(voxel->mtop), &mtopbuffer[j], p_topo_dem.ESIZE);
	memcpy(&(voxel->base), &babuffer[j], p_topo_dem.ESIZE);
	memcpy(&(voxel->moho), &mobuffer[j], p_topo_dem.ESIZE);
      } 
    break;
  case VX_SRC_HR:

    if(gcoor[0]>=0&&gcoor[1]>=0&&gcoor[2]>=0&&
       gcoor[0]<hr_a.N[0]&&gcoor[1]<hr_a.N[1]&&gcoor[2]<hr_a.N[2])
      {
	voxel->vel_cell[0]= hr_a.O[0]+gcoor[0]*step_hr[0];
	voxel->vel_cell[1]= hr_a.O[1]+gcoor[1]*step_hr[1];
	voxel->vel_cell[2]= hr_a.O[2]+gcoor[2]*step_hr[2];
	j=voxbytepos(gcoor,hr_a.N,p_vp_hr.ESIZE);
	memcpy(&(voxel->provenance), &hrtbuffer[j], p0_ESIZE);
	memcpy(&(voxel->vp), &hrbuffer[j], p_vp_hr.ESIZE);
	memcpy(&(voxel->vs), &hrvsbuffer[j], p_vp_hr.ESIZE);	
      }

    break;

  case VX_SRC_LR:

    if(gcoor[0]>=0&&gcoor[1]>=0&&gcoor[2]>=0&&
       gcoor[0]<lr_a.N[0]&&gcoor[1]<lr_a.N[1]&&gcoor[2]<lr_a.N[2])
      {
        voxel->vel_cell[0]= lr_a.O[0]+gcoor[0]*step_lr[0];
        voxel->vel_cell[1]= lr_a.O[1]+gcoor[1]*step_lr[1];
        voxel->vel_cell[2]= lr_a.O[2]+gcoor[2]*step_lr[2];
        j=voxbytepos(gcoor,lr_a.N,p_vp_lr.ESIZE);
        memcpy(&(voxel->provenance), &lrtbuffer[j], p_vp_lr.ESIZE);
        memcpy(&(voxel->vp), &lrbuffer[j], p_vp_lr.ESIZE);
        memcpy(&(voxel->vs), &lrvsbuffer[j], p_vp_lr.ESIZE);
      }

    break;
  default:
    break;
  }

  return;
}


/* Query elevation of free surface at point 'coor' */
void vx_getsurface(double *coor, vx_coord_t coor_type, float *surface)
{
  vx_getsurface_private(coor, coor_type, surface);
  return;
}


/* Private function for querying elevation of free surface at point 'coor'.  */
int vx_getsurface_private(double *coor, vx_coord_t coor_type, float *surface)
{
if(cvmhlabn_debug) fprintf(stderr,"CALLING -- vx_getsurface_private\n");
  int gcoor[3];
  double SP[2],SPUTM[2];
  int j;
  vx_entry_t entry;
  int do_bkg = False;

  *surface = p0_NO_DATA_VALUE;

  entry.coor[0] = coor[0];
  entry.coor[1] = coor[1];
  entry.coor[2] = 0.0;
  entry.coor_type = coor_type;

  // Initialize entry structure
  vx_init_entry(&entry);

  switch (entry.coor_type) {
  case VX_COORD_GEO:
    
    SP[0]=entry.coor[0];
    SP[1]=entry.coor[1];
    
    gctp(SP,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,efile,
	 SPUTM,&outsys,&outzone,inparm,&outunit,&outdatum,
	 file27, file83,&iflg);
    
    entry.coor_utm[0]=SPUTM[0];
    entry.coor_utm[1]=SPUTM[1];
    entry.coor_utm[2]=entry.coor[2];
    break;
  case VX_COORD_UTM:
    entry.coor_utm[0]=entry.coor[0];
    entry.coor_utm[1]=entry.coor[1];
    entry.coor_utm[2]=entry.coor[2];
    break;
  default:
    return(1);
    break;
  }


  gcoor[0]=round((entry.coor_utm[0]-to_a.O[0])/step_to[0]);
  gcoor[1]=round((entry.coor_utm[1]-to_a.O[1])/step_to[1]);
  gcoor[2]=0;
    
  /* check if inside topo volume */
  if(gcoor[0]>=0&&gcoor[1]>=0&&
     gcoor[0]<to_a.N[0]&&gcoor[1]<to_a.N[1]) {	      

if(cvmhlabn_debug) fprintf(stderr," SURFACE: in topo\n");

    j=voxbytepos(gcoor,to_a.N,p_topo_dem.ESIZE);
    memcpy(&(entry.topo), &tobuffer[j], p_topo_dem.ESIZE);
    memcpy(&(entry.mtop), &mtopbuffer[j], p_topo_dem.ESIZE);

    // Test points for topo - mtop gap
    // -118.75 36.8 0.0
    // -118.5 36.8 0.0
    // 345500.000000  4059000.0 0.0

    { // this is from where (vx_use_gtl != TRUE)
      /* check for valid topo values */
      if ((entry.topo - p0_NO_DATA_VALUE > 0.1) && 
	  (entry.mtop - p0_NO_DATA_VALUE > 0.1)) {
	
	if (entry.topo > entry.mtop) {
	  *surface = entry.mtop - ELEV_EPSILON;
	} else {
	  *surface = entry.topo - ELEV_EPSILON;
	}

        int MAX_ITER=abs(round((COOR_UTM2 - *surface)/step_hr[2]));
	
	int flag = 0;
	int num_iter = 0;
if(cvmhlabn_debug) fprintf(stderr," SURFACE: LOOPING in here start %lf(%lf)\n", entry.coor[2], *surface);
	entry.coor[2] = *surface;
	while (!flag) {
//	  if (num_iter > MAX_ITER_ELEV) {
	  if (num_iter > MAX_ITER) {
	    *surface = p0_NO_DATA_VALUE;
	    flag = 1;
if(cvmhlabn_debug) fprintf(stderr," SURFACE: Out of MAX\n");
	  }
	  num_iter = num_iter + 1;
	  vx_getcoord_private(&entry, False);
	  if ((entry.vp < 0.0) || (entry.vs < 0.0)) {
	    switch (entry.data_src) {
	    case VX_SRC_HR:
	      entry.coor[2] -= fabs(step_hr[2]);
if(cvmhlabn_debug) fprintf(stderr, " SURFACE: HR offset %lf -> %lf\n", step_hr[2], entry.coor[2]);
	      break;
            case VX_SRC_LR:
              entry.coor[2] -= fabs(step_lr[2]);
if(cvmhlabn_debug) fprintf(stderr, " SURFACE: LR offset %lf -> %lf\n", step_hr[2], entry.coor[2]);
              break;
	    default:
	      do_bkg = True;
	      flag = 1;
	      break;
	    }
	  } else {
	    *surface = entry.coor[2];
	    flag = 1;
	  }
	}
      } else {
	do_bkg = True;
      }
    }
      
  } else {
    do_bkg = True;
  }

  if (do_bkg) {
    *surface = p0_NO_DATA_VALUE;
  }

if(cvmhlabn_debug)  fprintf(stderr,"DONE -- vx_getsurface_private surface found=%lf\n", *surface);
  return(0);
}


/* Return mtop at coordinates 'coor' in 'surface' */
void vx_model_top(double *coor, vx_coord_t coor_type, float *surface)
{
  int gcoor[3];
  double SP[2],SPUTM[2];
  int j;
  vx_entry_t entry;
  int do_bkg = False;

  *surface = p0_NO_DATA_VALUE;

  entry.coor[0] = coor[0];
  entry.coor[1] = coor[1];
  entry.coor[2] = 0.0;
  entry.coor_type = coor_type;

  // Initialize entry structure
  vx_init_entry(&entry);

  switch (entry.coor_type) {
  case VX_COORD_GEO:
    
    SP[0]=entry.coor[0];
    SP[1]=entry.coor[1];
    
    gctp(SP,&insys,&inzone,inparm,&inunit,&indatum,&ipr,efile,&jpr,efile,
	 SPUTM,&outsys,&outzone,inparm,&outunit,&outdatum,
	 file27, file83,&iflg);
    
    entry.coor_utm[0]=SPUTM[0];
    entry.coor_utm[1]=SPUTM[1];
    entry.coor_utm[2]=entry.coor[2];
    break;
  case VX_COORD_UTM:
    entry.coor_utm[0]=entry.coor[0];
    entry.coor_utm[1]=entry.coor[1];
    entry.coor_utm[2]=entry.coor[2];
    break;
  default:
    return;
    break;
  }

  gcoor[0]=round((entry.coor_utm[0]-to_a.O[0])/step_to[0]);
  gcoor[1]=round((entry.coor_utm[1]-to_a.O[1])/step_to[1]);
  gcoor[2]=0;
    
  /* check if inside topo volume */
  if(gcoor[0]>=0&&gcoor[1]>=0&&
     gcoor[0]<to_a.N[0]&&gcoor[1]<to_a.N[1]) {	      
    j=voxbytepos(gcoor,to_a.N,p_topo_dem.ESIZE);
    memcpy(&(entry.topo), &tobuffer[j], p_topo_dem.ESIZE);
    memcpy(&(entry.mtop), &mtopbuffer[j], p_topo_dem.ESIZE);
  } else {
    do_bkg = True;
  }

  if (!do_bkg) {

    /* check for valid topo values */
    if ((entry.topo - p0_NO_DATA_VALUE > 0.1) && 
	(entry.mtop - p0_NO_DATA_VALUE > 0.1)) {
      if (entry.topo > entry.mtop) {
	*surface = entry.mtop - ELEV_EPSILON;
      } else {
	*surface = entry.topo - ELEV_EPSILON;
      }
      
      int flag = 0;
      int num_iter = 0;
      entry.coor[2] = *surface;
      while (!flag) {
	if (num_iter > MAX_ITER_ELEV) {
	  *surface = p0_NO_DATA_VALUE;
	  flag = 1;
	}
	num_iter = num_iter + 1;
	vx_getcoord_private(&entry, False);
	if ((entry.vp < 0.0) || (entry.vs < 0.0)) {
	  switch (entry.data_src) {
	  case VX_SRC_HR:
	    entry.coor[2] -= fabs(step_hr[2]);
	    break;
          case VX_SRC_LR:
            entry.coor[2] -= fabs(step_lr[2]);
            break;
	  default:
	    do_bkg = True;
	    flag = 1;
	    break;
	  }
	} else {
	  *surface = entry.coor[2];
	  flag = 1;
	}
      }
      
    } else {
      do_bkg = True;
    }
  }

  if (do_bkg) {
    *surface = p0_NO_DATA_VALUE;
  }

  return;
}

/* Return the voxel 'voxel' that lies precisely at the coord/model
   specified in 'entry'. If point lies outside of volume, the returned
   voxel will be empty. */
void vx_voxel_at_coord(vx_entry_t *entry, vx_voxel_t *voxel)
{
  int j;
  int model_coor[3]; // x,y,z of closest voxel in volume
  int model_max[3]; // max size x,y,z of volume
  double gcoor[3]; // coord of point wrt volume
  float gcoor_min[3]; // UTM coord of volume origin
  float step[3];
  int esize=p0_ESIZE;

  vx_init_voxel(voxel);

  /* find min utm/max coord/step size for specified model */
  for (j = 0; j < 3; j++) {
    switch (entry->data_src) {
    case VX_SRC_TO:
      gcoor_min[j] = to_a.O[j];
      model_max[j] = to_a.N[j];
      step[j] = step_to[j];
      esize = p_topo_dem.ESIZE;
    default:
      return;
      break;
    }
  }

  /* Find coord of point wrt specified volume */
  gcoor[0]=(entry->coor_utm[0]-gcoor_min[0])/step[0];
  gcoor[1]=(entry->coor_utm[1]-gcoor_min[1])/step[1];
  if (entry->data_src != VX_SRC_TO) {
    gcoor[2]=(entry->coor_utm[2]-gcoor_min[2])/step[2];
  } else {
    gcoor[2] = 0.0;
  }

  /* Find voxel */
  for (j = 0; j < 3; j++) {
    if ((gcoor[j] < 0) || (gcoor[j] > (model_max[j] - 1))) {
      return;
    } else {
      model_coor[j] = round(gcoor[j]);
    }
  }

  /* Calc index byte offset in volume */
  j = voxbytepos(model_coor, model_max, esize);

  /* Get vp/vs for closest voxel */
  switch (entry->data_src) {
  case VX_SRC_TO:
    memcpy(&(voxel->topo), &tobuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->mtop), &mtopbuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->base), &babuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->moho), &mobuffer[j], p_topo_dem.ESIZE);
    break;
  default:
    return;
    break;
  }

  voxel->coor[0] = model_coor[0];
  voxel->coor[1] = model_coor[1];
  voxel->coor[2] = model_coor[2];
  voxel->data_src = entry->data_src;

  return;
}


/* Return closest voxel to coord/model specified in 'entry' */
void vx_closest_voxel_to_coord(vx_entry_t *entry, vx_voxel_t *voxel)
{
  int j;
  int model_coor[3]; // x,y,z of closest voxel in volume
  int model_max[3]; // max size x,y,z of volume
  double gcoor[3]; // coord of point wrt volume
  float gcoor_min[3]; // UTM coord of volume origin
  float step[3];
  int esize;
  //float testval;

  vx_init_voxel(voxel);

  /* find min utm/max coord/step size for specified model */
  for (j = 0; j < 3; j++) {
    switch (entry->data_src) {
    case VX_SRC_TO:
      gcoor_min[j] = to_a.O[j];
      model_max[j] = to_a.N[j];
      step[j] = step_to[j];
      esize = p_topo_dem.ESIZE;
      break;
    default:
      return;
      break;
    }
  }

  /* Find coord of point wrt specified volume */
  gcoor[0]=(entry->coor_utm[0]-gcoor_min[0])/step[0];
  gcoor[1]=(entry->coor_utm[1]-gcoor_min[1])/step[1];
  if (entry->data_src != VX_SRC_TO) {
    gcoor[2]=(entry->coor_utm[2]-gcoor_min[2])/step[2];
  } else {
    gcoor[2]=0.0;
  }

  /* Find closest voxel */
  for (j = 0; j < 3; j++) {
    if (gcoor[j] < 0) {
      model_coor[j] = 0;
    } else if (gcoor[j] > (model_max[j] - 1)) {
      model_coor[j] = model_max[j] - 1;
    } else {
      model_coor[j] = round(gcoor[j]);
    }
  }

  /* Calc index byte offset in volume */
  j = voxbytepos(model_coor, model_max, esize);

  /* Get vp/vs for closest voxel */
  switch (entry->data_src) {
  case VX_SRC_TO:
    memcpy(&(voxel->topo), &tobuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->mtop), &mtopbuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->base), &babuffer[j], p_topo_dem.ESIZE);
    memcpy(&(voxel->moho), &mobuffer[j], p_topo_dem.ESIZE);
    break;
  default:
    return;
    break;
  }

  voxel->coor[0] = model_coor[0];
  voxel->coor[1] = model_coor[1];
  voxel->coor[2] = model_coor[2];
  voxel->data_src = entry->data_src;

  return;
}


/* Calculate 2D/3D distance from point 'entry' to voxel 'voxel' */
void vx_dist_point_to_voxel(vx_entry_t *entry, vx_voxel_t *voxel, 
			    float *dist_2d, float *dist_3d)
{
  int j;
  int model_max[3]; // max size x,y,z of volume
  double gcoor[3]; // coord of point wrt volume
  float gcoor_min[3]; // UTM coord of volume origin
  float step[3];
  int esize;
  double dxyz[3];

  /* find min utm/max coord/step size for specified model */
  for (j = 0; j < 3; j++) {
    switch (voxel->data_src) {
    case VX_SRC_TO:
      gcoor_min[j] = to_a.O[j];
      model_max[j] = to_a.N[j];
      step[j] = step_to[j];
      esize = p_topo_dem.ESIZE;
      break;
    default:
      return;
      break;
    }
  }

  /* Find coord of point wrt closest volume */
  gcoor[0]=(entry->coor_utm[0]-gcoor_min[0])/step[0];
  gcoor[1]=(entry->coor_utm[1]-gcoor_min[1])/step[1];
  if (entry->data_src != VX_SRC_TO) {
    gcoor[2]=(entry->coor_utm[2]-gcoor_min[2])/step[2];
  } else {
    gcoor[2]=0.0;
  }

  for (j = 0; j < 3; j++) {
    dxyz[j] = fabs(gcoor[j] - (double)voxel->coor[j]);
  }

  if(0) {
      fprintf(stderr,"KEEP warnings down: esize(%d) model_max(%d,%d,%d)\n",esize,model_max[0],model_max[1],model_max[2]);
  }

  /* Calculate min distance from selected model */
  /* LR cell size is 1000x1000x100, CM is 10000x10000x1000 */
  *dist_2d = sqrt(pow(dxyz[0] * step[0], 2.0) + 
		  pow(dxyz[1] * step[1], 2.0));
  *dist_3d = sqrt(pow(dxyz[0] * step[0], 2.0) + 
		  pow(dxyz[1] * step[1], 2.0) + 
		  pow(dxyz[2] * step[2], 2.0));
  
  return;
}


/* Get voxel byte offset position by the index values 'ic'
   and datatype size 'esize' */
int voxbytepos(int *ic,int *gs,int esize) {
  int pos;

  pos=(gs[0]*gs[1]*(ic[2])+gs[0]*(ic[1])+ic[0])*esize;
  return pos;
}


/* Calculate density (rho) from vp and the source model */
double calc_rho(float vp, vx_src_t data_src)
{
  double rho = 0.0;
  float fl;

  /* Compute rho */
  switch (data_src) {
  case VX_SRC_HR:
  case VX_SRC_LR:
    /*** Density should be at least 1000 ***/
    if (vp!=1480) {
      if (vp>744.) {
	fl = vp/1000.0;
	rho = 1000.0*(fl*(1.6612 + 
			  fl*(-0.4721 + fl*(0.0671 + 
					    fl*(-0.0043 + 
						fl*0.000106)))));
      } else
	rho = 1000.0;
    } else
      rho = 1000.0;
    break;
  case VX_SRC_CM:
    fl = vp/1000;
    rho = 1000.0*(fl*(1.6612 + 
		      fl*(-0.4721 + fl*(0.0671 + 
					fl*(-0.0043 + 
					    fl*0.000106)))));
    break;
  default:
    rho = p0_NO_DATA_VALUE;
    break;
  }

  return(rho);  
}


/* Initialize contents of entry structure */
void vx_init_entry(vx_entry_t *entry) {
  int j;

  for(j = 0; j < 2; j++) {
    entry->coor_utm[j] = p0_NO_DATA_VALUE;
    entry->elev_cell[j] = p0_NO_DATA_VALUE;
    entry->vel_cell[j] = p0_NO_DATA_VALUE;
  }
  entry->vel_cell[2] = p0_NO_DATA_VALUE;
  entry->coor_utm[2] = p0_NO_DATA_VALUE;

  entry->topo = entry->mtop = entry->base = entry->moho = p0_NO_DATA_VALUE;
  entry->provenance = p0_NO_DATA_VALUE;
  entry->vp = entry->vs = entry->rho = p0_NO_DATA_VALUE;
  entry->data_src = VX_SRC_NR;

  return;
}


/* Initialize contents of voxel structure */
void vx_init_voxel(vx_voxel_t *voxel) {
  int j;

  // Initially set to no data
  for(j = 0; j < 2; j++) {
    voxel->elev_cell[j] = p0_NO_DATA_VALUE;
    voxel->vel_cell[j] = p0_NO_DATA_VALUE;
  }
  voxel->vel_cell[2] = p0_NO_DATA_VALUE;
  voxel->topo = voxel->mtop = voxel->base = voxel->moho = p0_NO_DATA_VALUE;
  voxel->provenance = p0_NO_DATA_VALUE;
  voxel->vp = voxel->vs = voxel->rho = p0_NO_DATA_VALUE;

  return;
}

