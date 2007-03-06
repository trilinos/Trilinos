/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef ZOLTAN_DRUM

/* Implementation of DRUM interface with Zoltan */
/* Compile with ZOLTAN_DRUM and set USE_DRUM parameter to 1 to use */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "params_const.h"
#include "zz_const.h"
#include <string.h>

/* These belong in drum.h, but we want to avoid the dependency on
   drum.h for Zoltan_Struct */
extern int DRUM_hierCreateCallbacks(DRUM_machineModel *dmm, 
				    struct Zoltan_Struct *zz);
extern void DRUM_hierSetCallbacks(struct Zoltan_Struct *zz);

/**********  parameters structure for DRUM-related methods **********/
static PARAM_VARS Drum_params[] = {
        { "USE_DRUM", NULL, "INT", 0 },
        { "DRUM_HIER", NULL, "INT", 0 },
        { "ZOLTAN_BUILD_DRUM_TREE", NULL, "INT", 0 },
        { "ZOLTAN_START_DRUM_MONITORS", NULL, "INT", 0 },
	{ "DRUM_MONITORING_FREQUENCY", NULL, "INT", 1 },
	{ "DRUM_DEBUG_LEVEL", NULL, "INT", 0 },
        { "DRUM_POWER_FILE_LOG", NULL, "STRING", 0 },
        { "DRUM_USE_NETWORK_POWERS", NULL, "INT", 0 },
	{ "DRUM_FIXED_NETWORK_WEIGHT", NULL, "FLOAT", 0 },
        { "DRUM_USE_NWS", NULL, "INT", 0 },
        { NULL, NULL, NULL, 0 } };

int Zoltan_Drum_Set_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, Drum_params, &result, &index);
    return(status);
}

int Zoltan_Drum_Init_Struct(struct Zoltan_Drum_Struct *zds) {

  zds->dmm = NULL;

  return ZOLTAN_OK;
}

int Zoltan_Drum_Init(ZZ *zz) {

  /* bind DRUM-related Zoltan parameters */
  Zoltan_Bind_Param(Drum_params, "USE_DRUM", (void *) &zz->Drum.use_drum);
  Zoltan_Bind_Param(Drum_params, "DRUM_HIER", (void *) &zz->Drum.drum_hier);
  Zoltan_Bind_Param(Drum_params, "ZOLTAN_BUILD_DRUM_TREE", 
		    (void *) &zz->Drum.build_tree);
  Zoltan_Bind_Param(Drum_params, "ZOLTAN_START_DRUM_MONITORS",
                    (void *) &zz->Drum.start_monitors);
  Zoltan_Bind_Param(Drum_params, "DRUM_MONITORING_FREQUENCY",
		    (void *) &zz->Drum.monitoring_frequency);
  Zoltan_Bind_Param(Drum_params, "DRUM_DEBUG_LEVEL",
		    (void *) &zz->Drum.debug_level);
  Zoltan_Bind_Param(Drum_params, "DRUM_POWER_FILE_LOG",
		    (void *) zz->Drum.power_filename);
  Zoltan_Bind_Param(Drum_params, "DRUM_USE_NETWORK_POWERS",
		    (void *) &zz->Drum.use_network_powers);
  Zoltan_Bind_Param(Drum_params, "DRUM_FIXED_NETWORK_WEIGHT",
		    (void *) &zz->Drum.fixed_network_weight);
  Zoltan_Bind_Param(Drum_params, "DRUM_USE_NWS",
		    (void *) &zz->Drum.use_nws);

  /* set default values */
  /* can't do this - this is called on each LB_Balance invocation */
  /*zz->Drum.dmm = NULL;*/
  zz->Drum.use_drum = 0;
  zz->Drum.drum_hier = 0;
  zz->Drum.build_tree = 1;
  zz->Drum.start_monitors = 1;
  zz->Drum.monitoring_frequency = 1;
  zz->Drum.debug_level = 0;
  zz->Drum.power_filename[0] = '\0';
  zz->Drum.use_network_powers = 0;
  zz->Drum.fixed_network_weight = 0.0;
  zz->Drum.use_nws = 0;

  Zoltan_Assign_Param_Vals(zz->Params, Drum_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);

  return ZOLTAN_OK;
}

int Zoltan_Drum_Create_Model(ZZ *zz) {
  char *yo = "Zoltan_Drum_Create_Model";
  int ierr;
  FILE *fp;
  char buf[80];

  /* check params */
  Zoltan_Drum_Init(zz);

  if (zz->Drum.use_drum && !zz->Drum.dmm) {
    if (zz->Drum.build_tree) {
      zz->Drum.dmm = DRUM_createMachineModel(zz->Communicator,
					     zz->Drum.debug_level);
      if (!zz->Drum.dmm) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo,
			   "Unable to create DRUM machine model");
	return ZOLTAN_FATAL;
      }

      ierr = DRUM_initMachineModel(zz->Drum.dmm);
      if (ierr == DRUM_FATAL || ierr == DRUM_MEMERR) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo,
			   "Unable to initialize DRUM machine model");
	return (ierr == DRUM_FATAL ? ZOLTAN_FATAL : ZOLTAN_MEMERR);
      }

      if( zz->Drum.drum_hier){
	Zoltan_Set_Param(zz, "LB_METHOD", "HIER");
	ierr = DRUM_hierCreateCallbacks(zz->Drum.dmm, zz);
	if(ierr != DRUM_OK){
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo,
			     "DRUM_hier_create_callbacks returned an error");
	  return (ierr == DRUM_FATAL ? ZOLTAN_FATAL : ZOLTAN_MEMERR);
	}
	DRUM_hierSetCallbacks(zz);
      }

      /* print the "power file" if it was requested */
      if (zz->Proc == 0 && strcmp(zz->Drum.power_filename,"")) {
	fp = fopen(zz->Drum.power_filename, "w");
	if (fp) {
	  DRUM_printMachineModel(zz->Drum.dmm, fp);
	  fclose(fp);
	}
	else {
	  ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open power file");
	}
      }
      
      DRUM_setMonitoringFrequency(zz->Drum.dmm, zz->Drum.monitoring_frequency);
      sprintf(buf, "%d", zz->Drum.use_network_powers);
      DRUM_setParam(zz->Drum.dmm, "USE_NETWORK_POWERS", buf);
      sprintf(buf, "%f", zz->Drum.fixed_network_weight);
      DRUM_setParam(zz->Drum.dmm, "FIXED_NETWORK_WEIGHT", buf);
      sprintf(buf, "%d", zz->Drum.use_nws);
      DRUM_setParam(zz->Drum.dmm, "USE_NWS", buf);


    }
    
    if (zz->Drum.start_monitors) {
      ierr = DRUM_startMonitoring(zz->Drum.dmm);
      if (ierr == DRUM_FATAL || ierr == DRUM_MEMERR) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo,
			   "Unable to start DRUM monitors");
	return (ierr == DRUM_FATAL ? ZOLTAN_FATAL : ZOLTAN_MEMERR);
      }
    }
  }

  return ZOLTAN_OK;
}

/****************************************************************************/
void Zoltan_Drum_Free_Structure(ZZ *zz) {

  if (zz->Drum.dmm) {
    DRUM_deleteMachineModel(zz->Drum.dmm);
  }
  zz->Drum.dmm = NULL;
}

/****************************************************************************/
void Zoltan_Drum_Copy_Struct(struct Zoltan_Drum_Struct *to,
			     struct Zoltan_Drum_Struct const *from) {

  to->dmm = from->dmm;
  to->use_drum = from->use_drum;
  to->build_tree = from->build_tree;
  to->start_monitors = from->start_monitors;
  to->monitoring_frequency = from->monitoring_frequency;
  to->debug_level = from->debug_level;
  strncpy(to->power_filename, from->power_filename, 256);
  to->use_network_powers = from->use_network_powers;
  to->fixed_network_weight = from->fixed_network_weight;
  to->use_nws = from->use_nws;
}

/****************************************************************************/

int Zoltan_Drum_Start_Monitors(ZZ *zz) {
  int ierr;
  char *yo = "Zoltan_Drum_Start_Monitors";

  if ((zz->Drum.use_drum == 0) || (zz->Drum.start_monitors == 0))
    return ZOLTAN_OK;

  if (!zz->Drum.dmm) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "DRUM not initialized");
    return ZOLTAN_FATAL;
  }

  ierr = DRUM_startMonitoring(zz->Drum.dmm);
  if (ierr == DRUM_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to start DRUM monitors");
    return ZOLTAN_FATAL;
  }
  return ZOLTAN_OK;
}

int Zoltan_Drum_Stop_Monitors(ZZ *zz) {
  int ierr;
  char *yo = "Zoltan_Drum_Stop_Monitors";
  FILE *fp;

  if ((zz->Drum.use_drum == 0) || (zz->Drum.start_monitors == 0))
    return ZOLTAN_OK;

  if (!zz->Drum.dmm) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "DRUM not initialized");
    return ZOLTAN_FATAL;
  }

  ierr = DRUM_stopMonitoring(zz->Drum.dmm);
  if (ierr == DRUM_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to stop DRUM monitors");
    return ZOLTAN_FATAL;
  }

  ierr = DRUM_computePowers(zz->Drum.dmm);
  if (ierr == DRUM_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to compute DRUM powers");
    return ZOLTAN_FATAL;
  }

  /* print the "power file" if it was requested */
  if (zz->Proc == 0 && strcmp(zz->Drum.power_filename,"")) {
    fp = fopen(zz->Drum.power_filename, "a");
    if (fp) {
      DRUM_printMachineModel(zz->Drum.dmm, fp);
      fclose(fp);
    }
    else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open power file");
      return ZOLTAN_WARN;
    }
  }
  else {
    if (zz->Proc == 0) {
      printf("Skipping power file output\n"); fflush(stdout);
    }
  }

  return ZOLTAN_OK;
}

int Zoltan_Drum_Set_Part_Sizes(ZZ *zz) {
  /* for now this work for a single weight dimension and exactly
     one partition per process.  It may make sense to enhance it, 
     particularly in relation to hierarchical balancing.

     Does it make sense just to use the same partition sizes for each
     weight in the context of multiple weights specified?
  */
  float part_size;
  char *yo = "Zoltan_Drum_Set_Part_Sizes";
  int ierr;
  int wgt_idx;
  int part_id;

  if (zz->Drum.use_drum == 0) return ZOLTAN_OK;

  part_size = DRUM_getLocalPartSize(zz->Drum.dmm, &ierr);
  if (ierr == DRUM_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Cannot get DRUM partition sizes");
    return ZOLTAN_FATAL;
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS && zz->Proc == 0) {
    printf("Part Sizes: "); fflush(stdout);
  }
  if (zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS) {
    printf("proc %d: %.2f ", zz->Proc, part_size); fflush(stdout);
  }
  if (zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS && zz->Proc == 0) {
    printf("\n"); fflush(stdout);
  }
  wgt_idx = 0;
  part_id = 0;
  /* clear old part size info, first */
  Zoltan_LB_Set_Part_Sizes(zz, 0, -1, NULL, NULL, NULL);
  /* now set the part size computed by DRUM */
  Zoltan_LB_Set_Part_Sizes(zz, 0, 1, &part_id, &wgt_idx, &part_size);

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* ZOLTAN_DRUM */
