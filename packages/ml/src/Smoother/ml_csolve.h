/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_CSolve structure                               */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLCSOLVE__
#define __MLCSOLVE__

/* ******************************************************************** */
/* data structure type definition                                       */
/* ******************************************************************** */

typedef struct ML_CSolveFunc_Struct ML_CSolveFunc;
typedef struct ML_CSolve_Struct ML_CSolve;

/* ******************************************************************** */
/* local include files                                                  */
/* ******************************************************************** */

#include "ml_defs.h"
#include "ml_memory.h"
#include "ml_1level.h"

/* ******************************************************************** */
/* data definition for the ML_CSolve Class                              */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* These data structures define the coarse solver object.               */
/* -------------------------------------------------------------------- */

struct ML_CSolveFunc_Struct 
{
   int ML_id;
   int (*internal)(void *, int, double *, int, double *);
   int (*external)(void *, int, double *, int, double *);
};

struct ML_CSolve_Struct 
{
   int                     ML_id;
   struct ML_1Level_Struct *my_level;
   int                     ntimes;
   double                  tol;
   ML_CSolveFunc           *func;
   void                    *data;
   void                    (*data_destroy)(void *);
   double                  build_time, apply_time;
   char                    *label;
};

/* ******************************************************************** */
/* ******************************************************************** */
/*      User Interface Proto-types                                      */
/* ******************************************************************** */
/* ******************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif
extern int ML_CSolve_Create(ML_CSolve **);
extern int ML_CSolve_Set_Label(ML_CSolve *, char *label);
extern int ML_CSolve_Init(ML_CSolve *);
extern int ML_CSolve_Destroy(ML_CSolve **);
extern int ML_CSolve_Clean(ML_CSolve *);
extern int ML_CSolve_Check(ML_CSolve *);
extern int ML_CSolve_Set_1Level(ML_CSolve *, ML_1Level *);
extern int ML_CSolve_Apply(ML_CSolve *, int, double *, int, double *);

#ifdef __cplusplus
}
#endif

#endif
