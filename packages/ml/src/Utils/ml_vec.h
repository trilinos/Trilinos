/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* local (to ML) data structure to hold vector information              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : August, 1997                                         */
/* ******************************************************************** */

#ifndef __MLVEC__
#define __MLVEC__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include <stdio.h>
/* #include <stdlib.h> */
#include "ml_common.h"
#include "ml_memory.h"
#include "ml_comm.h"
#include "ml_defs.h"

/* ******************************************************************** */
/* definition of the grid structure                                     */
/*  ML_id             : identification for a vector                     */
/*  VecLength         : length of vector                                */
/*  SetOrLoad         : a flag to see if storage is allocated.          */
/*  VecData           : holder for data                                 */
/* -------------------------------------------------------------------- */

typedef struct ML_DVector_Struct
{
  int     ML_id;
  ML_Comm *comm;
  int     VecLength;
  int     SetOrLoad;
  double  *VecData;

} ML_DVector;

/* ******************************************************************** */
/* functions to manipulate the vector structure                         */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int  ML_DVector_Create( ML_DVector **, ML_Comm *com );
extern int  ML_DVector_Init(ML_DVector *vec);
extern int  ML_DVector_Destroy( ML_DVector ** );
extern int  ML_DVector_Clean( ML_DVector *vec );
extern int  ML_DVector_LoadData( ML_DVector *, int, double * );
extern int  ML_DVector_SetData( ML_DVector *, int, double * );
extern int  ML_DVector_GetLength( ML_DVector * );
extern int  ML_DVector_GetData( ML_DVector *, int *, double * );
extern int  ML_DVector_GetDataPtr( ML_DVector *, double ** );
extern int  ML_DVector_Check( ML_DVector * );
extern int  ML_DVector_Scale( double, ML_DVector * );
extern int  ML_DVector_Copy( ML_DVector *, ML_DVector * );
extern int  ML_DVector_Axpy( double, ML_DVector *, ML_DVector * );
extern int  ML_DVector_Aypx( double, ML_DVector *, ML_DVector * );
extern int ML_DVector_Print(int length, double *data, char *label, ML_Comm *comm);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
