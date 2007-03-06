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

#ifndef __HA_DRUM_H
#define __HA_DRUM_H

/* prototypes and structure definitions for Zoltan's DRUM interface */

/* include drum.h for DRUM_machineModel type definition and DRUM prototypes */
#ifdef ZOLTAN_DRUM
#include "drum.h"
#endif

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Information needed for Zoltan/DRUM interface */
struct Zoltan_Drum_Struct {
#ifdef ZOLTAN_DRUM
  DRUM_machineModel *dmm;            /* pointer to DRUM machine model
					structure */
#endif
  int use_drum;                      /* is DRUM to be used? */
  int drum_hier;                     /* will DRUM store hier info? */
  int build_tree;                    /* should Zoltan build the tree? */
  int start_monitors;                /* should Zoltan start/stop monitors? */
  int monitoring_frequency;          /* DRUM monitoring frequency */
  int debug_level;                   /* debug level to pass to DRUM */
  char power_filename[256];          /* file name to use to print power file */
  int use_network_powers;            /* include network powers? */
  float fixed_network_weight;        /* what weight for network if using it? */
  int use_nws;                       /* use NWS for monitoring? */
};
typedef struct Zoltan_Drum_Struct Zoltan_Drum;

/* prototype for set_param function needed by params/set_param.c */
extern int Zoltan_Drum_Set_Param(char *name, char *val);

extern int Zoltan_Drum_Init(struct Zoltan_Struct *zz);
extern int Zoltan_Drum_Init_Struct(struct Zoltan_Drum_Struct *zds);
extern int Zoltan_Drum_Create_Model(struct Zoltan_Struct *zz);
extern void Zoltan_Drum_Copy_Struct(struct Zoltan_Drum_Struct *to,
				    struct Zoltan_Drum_Struct const *from);
extern void Zoltan_Drum_Free_Structure(struct Zoltan_Struct *zz);
extern int Zoltan_Drum_Start_Monitors(struct Zoltan_Struct *zz);
extern int Zoltan_Drum_Stop_Monitors(struct Zoltan_Struct *zz);
extern int Zoltan_Drum_Set_Part_Sizes(struct Zoltan_Struct *zz);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif

#endif /* ZOLTAN_DRUM */
