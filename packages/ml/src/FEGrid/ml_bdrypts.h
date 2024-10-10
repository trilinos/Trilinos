/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_BdryPts structure                              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLBDRYPTSH__
#define __MLBDRYPTSH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"

/* ******************************************************************** */
/* data definition for the ML_BdryPts Class                             */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* This data structure stores two arrays : one integer array for a list */
/* of boundary conditions in the grid space, and another integer array  */
/* for a list of boundary conditions in the equation space.  In this    */
/* implementation, only Dirichlet boundary conditions are stored.       */
/* -------------------------------------------------------------------- */

typedef struct ML_BdryPts_Struct ML_BdryPts;

struct ML_BdryPts_Struct {
   int    ML_id;
   int    Dirichlet_grid_CreateOrDup;
   int    Dirichlet_grid_length;
   int    *Dirichlet_grid_list;
   int    Dirichlet_eqn_CreateOrDup;
   int    Dirichlet_eqn_length;
   int    *Dirichlet_eqn_list;
};

/* ******************************************************************** */
/* function for accessing the ML_BdryPts Class                          */
/* ******************************************************************** */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int ML_BdryPts_Create(ML_BdryPts **);
extern int ML_BdryPts_Init(ML_BdryPts *);
extern int ML_BdryPts_Destroy(ML_BdryPts **);
extern int ML_BdryPts_Clean(ML_BdryPts *);
extern int ML_BdryPts_Check_Dirichlet_Grid(ML_BdryPts *);
extern int ML_BdryPts_Check_Dirichlet_Eqn(ML_BdryPts *);
extern int ML_BdryPts_Get_Dirichlet_Grid_Info(ML_BdryPts *, int *, int **);
extern int ML_BdryPts_Get_Dirichlet_Eqn_Info(ML_BdryPts *, int *, int **);
extern int ML_BdryPts_Load_Dirichlet_Grid(ML_BdryPts *,int,int*);
extern int ML_BdryPts_Load_Dirichlet_Eqn(ML_BdryPts *,int,int*);
extern int ML_BdryPts_Copy_Dirichlet_GridToEqn(ML_BdryPts *);
extern int ML_BdryPts_ApplyZero_Dirichlet_Grid(ML_BdryPts *, double *);
extern int ML_BdryPts_ApplyZero_Dirichlet_Eqn(ML_BdryPts *, double *);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
