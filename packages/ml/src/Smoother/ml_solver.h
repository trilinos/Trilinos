/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML structure                                      */
/* ******************************************************************** */
/* Author : Charles Tong (LLNL) and Raymond Tuminaro (SNL)              */
/* Date   : February, 1999                                              */
/* ******************************************************************** */

#ifndef _MLSOLVER_
#define _MLSOLVER_

#include <stdlib.h>

typedef struct ML_Solver_Struct 
{
   int          ML_id;               /* ID for the Solver structure     */
   int          reuse_flag;          /* flag for internal use           */
   void         (*func)();           /* function to perform the solve   */
   void         *Mat1;               /* primary matrix for the solver   */
   void         *Mat2;               /* L matrix (for direct solver)    */
   void         *Mat3;               /* U matrix (for direct solver)    */
   int          int_data;            /* integer data for the solver     */
   int          int_data2;           /* integer data for the solver     */
   double       dble_data;           /* double data for the solver      */
   int          *int_params1;        /* integer array for the solver    */
   int          *int_params2;        /* integer array for the solver    */
   double       *dble_params1;       /* double array for the solver     */
   double       *dble_params2;       /* double array for the solver     */
   void         *void_params1;       /* other data for the solver       */
   void         *void_params2;       /* other data for the solver       */
   void         *LUspl;              /* for direct solver               */
   void         *PERMspl;            /* for direct solver               */
   void         *grid;
   void         *ML_subcomm;

} ML_Solver;

#ifdef __cpluscplus
extern "C" {
#endif

int ML_Solver_Create( ML_Solver **sol );
int ML_Solver_Destroy( ML_Solver **sol );

#ifdef __cpluscplus
}
#endif
#endif

