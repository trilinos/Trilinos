/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the Grid and its access functions data structure      */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLGRID__
#define __MLGRID__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"
#include "ml_gridfunc.h"

typedef struct ML_Grid_Struct ML_Grid;

/* ******************************************************************** */
/* definition of the data structure for Grid                            */
/* -------------------------------------------------------------------- */

struct ML_Grid_Struct
{
   int           ML_id;
   void          *Grid;        /* user grid data structure        */
   ML_GridFunc   *gridfcn;     /* a set of grid access functions  */
   int           gf_SetOrLoad; /* see if gridfcn is created locally */
};

/* ******************************************************************** */
/* definition of the functions                                          */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int ML_Grid_Create( ML_Grid ** );
extern int ML_Grid_Init( ML_Grid * );
extern int ML_Grid_Destroy( ML_Grid ** );
extern int ML_Grid_Clean( ML_Grid * );
extern int ML_Grid_Set_Grid( ML_Grid *, void * );
extern int ML_Grid_Set_GridFunc( ML_Grid *, ML_GridFunc * );
extern int ML_Grid_Get_GridFunc( ML_Grid *, ML_GridFunc ** );
extern int ML_Grid_Create_GridFunc( ML_Grid * );

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
