/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_1Level structure                               */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __ML1LEVEL__
#define __ML1LEVEL__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* ******************************************************************** */
/* data structure type definition                                       */
/* ******************************************************************** */

typedef struct ML_1Level_Struct ML_1Level;

/* ******************************************************************** */
/* include files                                                        */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_bdrypts.h"
#include "ml_mapper.h"
#include "ml_grid.h"
#include "ml_comm.h"
#include "ml_comminfoop.h"
#include "ml_operator.h"
#include "ml_smoother.h"
#include "ml_csolve.h"
#include "ml_vec.h"

/* ******************************************************************** */
/* data definition for the ML_1Level Class                              */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* This data structure defines the components of a grid level in a      */
/* multilevel environment.                                              */
/* -------------------------------------------------------------------- */

struct ML_1Level_Struct
{
   int          id, levelnum;
   ML_Operator  *Amat, *Rmat, *Pmat;
   ML_BdryPts   *BCs;
   ML_Mapper    *eqn2grid;
   ML_Mapper    *grid2eqn;
   ML_Grid      *Grid;
   ML_DVector   *Amat_Normalization;
   ML_Smoother  *pre_smoother;
   ML_Smoother  *post_smoother;
   ML_CSolve    *csolve;
   ML_Comm      *comm;
};

#endif
