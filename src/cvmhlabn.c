/**
 * @file cvmhlabn.c
 * @brief Main file for CVMHLABN library.
 * @author - SCEC 
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers CVMH Los Angeles Basin Velocity Model
 *
 */

#include "cvmhlabn.h"

int debug=0;
int cvmhlabn_force_depth = 0;
vx_zmode_t cvmhlabn_zmode = VX_ZMODE_ELEVOFF;

/**
 * Initializes the CVMHLABN plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int cvmhlabn_init(const char *dir, const char *label) {
	char configbuf[512];

	// Initialize variables.
	cvmhlabn_configuration = calloc(1, sizeof(cvmhlabn_configuration_t));
	cvmhlabn_velocity_model = calloc(1, sizeof(cvmhlabn_model_t));

	// Configuration file location.
	sprintf(configbuf, "%s/model/%s/data/config", dir, label);

	// Read the configuration file.
	if (cvmhlabn_read_configuration(configbuf, cvmhlabn_configuration) != SUCCESS) {
                print_error("No configuration file was found to read from.");
                return FAIL;
        }

	// Set up the data directory.
	sprintf(cvmhlabn_data_directory, "%s/model/%s/data/%s", dir, label, cvmhlabn_configuration->model_dir);

        /* Init vx */
        if (vx_setup(cvmhlabn_data_directory) != 0) {
          return FAIL;
        }

	// Let everyone know that we are initialized and ready for business.
	cvmhlabn_is_initialized = 1;

	return SUCCESS;
}

/**  
  * 
**/

/* Setparam CVM-H */
int cvmhlabn_setparam(int id, int param, ...)
{
  va_list ap;
  int cvmhlabn_cmode;

  va_start(ap, param);

  switch (param) {
    case CVMHLABN_MODEL_PARAM_FORCE_DEPTH_ABOVE_SURF:
      cvmhlabn_force_depth = va_arg(ap, int);
      break;
    default:
      fprintf(stderr,"Unknown param for cvmhlabn_setparam()\n");
      break;
  }
  va_end(ap);
  return SUCCESS;
}


/**
 * Queries CVMHLABN at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @param cmode The cmode.
 * @return SUCCESS or FAIL.
 */
int cvmhlabn_query(cvmhlabn_point_t *points, cvmhlabn_properties_t *data, int numpoints, int cmode) {
  // setup >> points -> entry
  // retrieve >> entry -> data

  switch(cmode) {
      case CVMHLABN_COORD_GEO_DEPTH:
          cvmhlabn_zmode = VX_ZMODE_DEPTH;
          break;
      case CVMHLABN_COORD_GEO_ELEV:
          cvmhlabn_zmode = VX_ZMODE_ELEV;
          break;
      default:
          fprintf(stderr, "Unsupported coord type\n");
          return FAIL;
          break;
  }

  vx_setzmode(cvmhlabn_zmode);

  for(int i=0; i<numpoints; i++) {
      vx_entry_t entry;
      float vx_surf;

    /*
       By the time here, Conditions:

       Following condition must have met,
         1) Point data has not been filled in by previous model
         2) Point falls in crust or interpolation zone
         3) Point falls within the configured model region
     */

      /* Force depth mode if directed and point is above surface */
      if ((cvmhlabn_force_depth) && (cvmhlabn_zmode == VX_ZMODE_ELEV) &&
          (points[i].depth < 0.0)) {
        /* Setup point to query */
        entry.coor[0]=points[i].longitude;
        entry.coor[1]=points[i].latitude;
        vx_getsurface(&(entry.coor[0]), entry.coor_type, &vx_surf);
        if (vx_surf - VX_NO_DATA < 0.01) {
          /* Fallback to using UCVM topo */
          entry.coor[2]=points[i].depth;
        } else {
          entry.coor[2]=vx_surf - points[i].depth;
        }
      } else {
        /* Setup with direct point to query */
        entry.coor[0]=points[i].longitude;
        entry.coor[1]=points[i].latitude;
        entry.coor[2]=points[i].depth;
      }

      /* In case we got anything like degrees */
      if ((entry.coor[0]<360.) && (fabs(entry.coor[1])<90)) {
        entry.coor_type = VX_COORD_GEO;
      } else {
        entry.coor_type = VX_COORD_UTM;
      }

      /* Query the point */
      int rc=vx_getcoord(&entry);

      if(debug) {
        printf("%14.6f %15.6f %9.2f \n",
               entry.coor[0], entry.coor[1], entry.coor[2]);
        /* AP: Let's provide the computed UTM coordinates as well */
        printf("%10.2f %11.2f ", entry.coor_utm[0], entry.coor_utm[1]);

        printf("%10.2f %11.2f ", entry.elev_cell[0], entry.elev_cell[1]);
        printf("%9.2f ", entry.topo);
        printf("%9.2f ", entry.mtop);
        printf("%9.2f ", entry.base);
        printf("%9.2f ", entry.moho);
        printf("%s %10.2f %11.2f %9.2f ", VX_SRC_NAMES[entry.data_src],
               entry.vel_cell[0], entry.vel_cell[1], entry.vel_cell[2]);
        printf("%9.2f %9.2f %9.2f ", entry.provenance, entry.vp, entry.vs);
        printf("%9.2f\n", entry.rho);
      }

      if(debug) {
        fprintf(stderr,">>> a point..rc(%d)->",rc);
        switch(entry.data_src) {
          case VX_SRC_NR: {fprintf(stderr,"GOT VX_SRC_NR\n"); break; }
          case VX_SRC_HR: {fprintf(stderr,"GOT VX_SRC_HR\n"); break; }
          case VX_SRC_LR: {fprintf(stderr,"GOT VX_SRC_LR\n"); break; }
          case VX_SRC_CM: {fprintf(stderr,"GOT VX_SRC_CM\n"); break; }
          case VX_SRC_TO: {fprintf(stderr,"GOT VX_SRC_TO\n"); break; }
          case VX_SRC_BK: {fprintf(stderr,"GOT VX_SRC_BK\n"); break; }
          case VX_SRC_GT: {fprintf(stderr,"GOT VX_SRC_GT\n"); break; }
          default: {fprintf(stderr,"???\n"); break; }
        }
      }

      // 1 is bad, 0 is good and anything not in HR region/ie cvmhlabn 
      if(rc || entry.data_src != VX_SRC_HR) { 
        data[i].vp=-1;
        data[i].vs=-1;
        data[i].rho=-1;
        } else {
          data[i].vp=entry.vp;
          data[i].vs=entry.vs;
          data[i].rho=entry.rho;
      }

  }
  return SUCCESS;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int cvmhlabn_finalize() {
        vx_cleanup();
	cvmhlabn_is_initialized = 0;
	return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cvmhlabn_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(cvmhlabn_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, cvmhlabn_version_string, verlen);
  return 0;
}

/**
 * Reads the configuration file describing the various properties of CVMHLABNand populates
 * the configuration struct. This assumes configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The configuration file location on disk to read.
 * @param config The configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int cvmhlabn_read_configuration(char *file, cvmhlabn_configuration_t *config) {
	FILE *fp = fopen(file, "r");
	char key[40];
	char value[80];
	char line_holder[128];

	// If our file pointer is null, an error has occurred. Return fail.
	if (fp == NULL) {
		print_error("Could not open the configuration file.");
		return FAIL;
	}

	// Read the lines in the configuration file.
	while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
		if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
			sscanf(line_holder, "%s = %s", key, value);

			// Which variable are we editing?
			if (strcmp(key, "utm_zone") == 0)
  				config->utm_zone = atoi(value);
			if (strcmp(key, "model_dir") == 0)
				sprintf(config->model_dir, "%s", value);

		}
	}

	// Have we set up all configuration parameters?
	if (config->utm_zone == 0 || config->model_dir[0] == '\0' ) {
		print_error("One configuration parameter not specified. Please check your configuration file.");
		return FAIL;
	}

	fclose(fp);

	return SUCCESS;
}

/*
 * @param err The error string to print out to stderr.
 */
void print_error(char *err) {
	fprintf(stderr, "An error has occurred while executing CVMHLABN. The error was:\n\n");
	fprintf(stderr, "%s", err);
	fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
	fprintf(stderr, "about the computer you are running CVMHLABN on (Linux, Mac, etc.).\n");
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls cvmhlabn_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
	return cvmhlabn_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls cvmhlabn_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(cvmhlabn_point_t *points, cvmhlabn_properties_t *data, int numpoints) {
	return cvmhlabn_query(points, data, numpoints);
}

/**
 * Setparam function loaded and called by the UCVM library. Calls cvmhlabn_setparam.
 *
 * @param id  don'care
 * @param param 
 * @param val, it is actually just 1 int
 * @return Success or fail.
 */
int model_setparam(int id, int param, int val) {
	return cvmhlabn_setparam(id, param, val);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cvmhlabn_finalize.
 *
 * @return Success
 */
int model_finalize() {
	return cvmhlabn_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cvmhlabn_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
	return cvmhlabn_version(ver, len);
}

int (*get_model_init())(const char *, const char *) {
        return &cvmhlabn_init;
}
int (*get_model_query())(cvmhlabn_point_t *, cvmhlabn_properties_t *, int) {
         return &cvmhlabn_query;
}
int (*get_model_finalize())() {
         return &cvmhlabn_finalize;
}
int (*get_model_version())(char *, int) {
         return &cvmhlabn_version;
}
int (*get_model_setparam())(int, int, ...) {
         return &cvmhlabn_setparam;
}




#endif
