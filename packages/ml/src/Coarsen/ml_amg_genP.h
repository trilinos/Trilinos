/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* functions for setting up AMG                                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : October, 2000                                        */
/* ******************************************************************** */

#ifndef __MLAMGGENP__
#define __MLAMGGENP__

#include "ml_common.h"
#include "ml_amg.h"
#include "ml_operator.h"

/* ******************************************************************** */
/* functions defined here                                               */
/* ******************************************************************** */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* ******************************************************************** */
/* functions called by users                                            */
/* -------------------------------------------------------------------- */

extern int ML_Gen_MGHierarchy_UsingAMG(ML *, int start, 
                       int increment_or_decrement, ML_AMG *);

/* ******************************************************************** */
/* internal functions called by developers                              */
/* -------------------------------------------------------------------- */

extern int ML_AMG_Gen_MGHierarchy(ML *, int fine_level,
               int (*next_level)(ML *, int, ML_Operator *, ML_AMG *),
               int (*user_gen_prolongator)(ML *,int,int,void *,ML_AMG*),
               void *data, int internal_or_external, ML_AMG *);
extern int ML_AMG_Gen_Prolongator(ML*,int ,int,void *data,ML_AMG*);
extern int ML_AMG_Increment_Level(ML *,int level,ML_Operator *Amat,ML_AMG*);
extern int ML_AMG_Decrement_Level(ML *,int level,ML_Operator *Amat,ML_AMG*);
extern int ML_AMG_Identity_Getrows(ML_Operator *data, int N_requested_rows, 
               int requested_rows[], int allocated_space, int columns[], 
               double values[], int row_lengths[]);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

