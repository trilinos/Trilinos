/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* Functions for the ML_Solver structure                                     */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)            */
/* Date          : February, 1999                                            */
/* ************************************************************************* */

#include "ml_solver.h"
#include "ml_memory.h"
/* ************************************************************************* */
/* constructor                                                               */
/* ------------------------------------------------------------------------- */

int ML_Solver_Create( ML_Solver **sol )
{
   ML_Solver *solver;
  
   ML_memory_alloc((void**) sol, sizeof( ML_Solver ), "SO1" );
   solver = (*sol);
   solver->ML_id = ML_ID_SOLVER;
   solver->reuse_flag = 0;
   solver->func = NULL;
   solver->Mat1 = NULL;
   solver->Mat2 = NULL;
   solver->Mat3 = NULL;
   solver->int_data = 0;
   solver->int_data2 = 0;
   solver->dble_data = 0.0;
   solver->int_params1 = NULL;
   solver->int_params2 = NULL;
   solver->dble_params1 = NULL;
   solver->dble_params2 = NULL;
   solver->void_params1 = NULL;
   solver->void_params2 = NULL;
   solver->LUspl = NULL;
   solver->PERMspl = NULL;
   solver->gridtiles = NULL;

   return 0;
}

/* ************************************************************************* */
/* destructor                                                                */
/* ------------------------------------------------------------------------- */

int ML_Solver_Destroy( ML_Solver **sol )
{
   ML_Solver *solver;
  
   solver = (*sol);
   solver->ML_id = -1;
   solver->func = NULL;
   solver->Mat1 = NULL;
   solver->Mat2 = NULL;
   solver->Mat3 = NULL;
   if ( solver->int_params1 != NULL ) 
      ML_memory_free( (void**) &(solver->int_params1) );
   if ( solver->int_params2 != NULL ) 
      ML_memory_free( (void**) &(solver->int_params2) );
   if ( solver->dble_params1 != NULL ) 
      ML_memory_free( (void**) &(solver->dble_params1) );
   if ( solver->dble_params2 != NULL ) 
      ML_memory_free( (void**) &(solver->dble_params2) );
   if ( solver->LUspl != NULL )
      ML_memory_free( (void**) &(solver->LUspl) );
   if ( solver->PERMspl != NULL )
      ML_memory_free( (void**) &(solver->PERMspl) );
   ML_memory_free( (void**) sol );
   (*sol) = NULL;
   return 0;
}

/* ************************************************************************* */
/* check validity of the passed object                                       */
/* ------------------------------------------------------------------------- */

int ML_Solver_Check( ML_Solver *sol )
{
   if ( sol->ML_id != ML_ID_SOLVER )
   {
      printf("ML_Solver_Check ERROR : object unknown %d.\n",sol->ML_id);
      exit(1);
   }
   return 1;
}

