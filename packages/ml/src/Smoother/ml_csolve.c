/* *********************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact person, */
/* and disclaimer.                                                         */
/* *********************************************************************** */

/* *********************************************************************** */
/* Functions for the ML_CSolve structure                                   */
/* *********************************************************************** */
/* Author        : Charles Tong (LLNL)                                     */
/* Date          : March, 1999                                             */
/* *********************************************************************** */

#include <stdlib.h>
#include "ml_csolve.h"
#include <string.h>
#ifdef SUPERLU
#include "dsp_defs.h"
#endif
#include "ml_superlu.h"

/* *********************************************************************** */
/* Constructor                                                             */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Create(ML_CSolve **cs)
{
   ML_CSolve *ml_cs;

   ML_memory_alloc((void**) cs, sizeof(ML_CSolve), "CS1" );
   ml_cs = (*cs);
   ml_cs->ML_id = ML_ID_CSOLVE; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   ml_cs->data = NULL;
   ML_memory_alloc((void**) &(ml_cs->func),sizeof(ML_CSolveFunc),"CF1" );
   ml_cs->func->ML_id = ML_EMPTY;
   ml_cs->func->func_ptr = NULL;
   ml_cs->build_time = 0.0;
   ml_cs->apply_time = 0.0;
   ml_cs->label      = NULL;
   ml_cs->data_destroy = NULL;
   return 0;
} 

/* *********************************************************************** */
/* Initialize                                                              */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Init(ML_CSolve *ml_cs)
{
   ml_cs->ML_id = ML_ID_CSOLVE; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   ml_cs->data = NULL;
   ML_memory_alloc((void**) &(ml_cs->func),sizeof(ML_CSolveFunc),"CF2" );
   ml_cs->func->ML_id = ML_EMPTY;
   ml_cs->func->func_ptr = NULL;
   ml_cs->build_time = 0.0;
   ml_cs->apply_time = 0.0;
   ml_cs->label      = NULL;
   ml_cs->data_destroy = NULL;
   return 0;
} 

/* *********************************************************************** */
/* Destructor                                                              */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Destroy(ML_CSolve **cs)
{
   ML_CSolve_Clean(*cs);
   ML_memory_free( (void**) cs );
   (*cs) = NULL; 
   return 0;
}

/* *********************************************************************** */
/* partial destructor                                                      */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Clean(ML_CSolve *ml_cs)
{
#ifdef ML_TIMING_DETAILED
   double t1;
#endif

   if ( ml_cs->ML_id != ML_ID_CSOLVE ) 
   {
      printf("ML_CSolve_Clean error : wrong object.\n");
      exit(-1);
   }
#ifdef ML_TIMING_DETAILED
   if (ml_cs->label != NULL) {
      t1 = ML_gsum_double(ml_cs->build_time, global_comm);
      t1 = t1/((double) global_comm->ML_nprocs);
      if ( (global_comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Build time for %s\t= %e\n",ml_cs->label,t1);
   }

   if  (ml_cs->label != NULL) {
      t1 = ML_gsum_double(ml_cs->apply_time, global_comm);
      t1 = t1/((double) global_comm->ML_nprocs);
      if ( (global_comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s\t= %e\n",ml_cs->label,t1);
    }
#endif

   ml_cs->ML_id = -1; 
   ml_cs->my_level = NULL;
   ml_cs->ntimes = 0;
   ml_cs->tol = 0;
   if (ml_cs->data_destroy != NULL)
   {
      ml_cs->data_destroy( ml_cs->data );
      ml_cs->data = NULL;
   }
   if ( ml_cs->func->func_ptr == ML_SuperLU_Solve && ml_cs->data != NULL )
       ML_Clean_CSolveSuperLU( ml_cs->data, ml_cs->func );
   if (ml_cs->func->func_ptr == ML_CSolve_Aggr)
       ML_CSolve_Clean_Aggr( ml_cs->data, ml_cs->func );
   ML_memory_free( (void**) &(ml_cs->func) );
   ml_cs->data = NULL;
   ml_cs->func = NULL;
   ml_cs->data_destroy = NULL;
   if (ml_cs->label != NULL) { ML_free(ml_cs->label); ml_cs->label = NULL; }
   return 0;
}

/* *********************************************************************** */
/* Check to see if csolve function exists                                  */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Check(ML_CSolve *ml_cs)
{
   if ( ml_cs->ML_id != ML_ID_CSOLVE )
   {
      printf("ML_CSolve_Check : wrong object.\n");
      exit(1);
   }
   if ( ml_cs->func->func_ptr == NULL ) return 0;
   else                                  return 1;
}

/* *********************************************************************** */
/* set the 1Level parameter                                                */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Set_1Level(ML_CSolve *ml_cs, ML_1Level *mylevel)
{
   if ( ml_cs->ML_id != ML_ID_CSOLVE )
   {
      printf("ML_CSolve_Set_1Level error : wrong object.\n");
      exit(-1);
   }
   ml_cs->my_level = mylevel;
   return 0;
}

/* *********************************************************************** */
/* perform solve                                                           */
/* ----------------------------------------------------------------------- */

int ML_CSolve_Apply(ML_CSolve *csolve, int inlen, double din[], 
                    int outlen, double dout[])
{
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
   t0 = GetClock();
#endif
   if (csolve->func->func_ptr == NULL) 
      pr_error("ML_CSolve_Apply error : coarse solver not defined\n");

   csolve->func->func_ptr((ML_Solver *)csolve->data, inlen, din, outlen, dout);
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   csolve->apply_time += (GetClock() - t0);
#endif
   return 0;
}

/* ******************************************************************** */
/* give a label to this object                                          */
/* ******************************************************************** */

int ML_CSolve_Set_Label( ML_CSolve *csolve, char *label)
{
  int size;

   if (csolve->label != NULL) { ML_free(csolve->label); csolve->label = NULL; }
   size = strlen(label) + 1;
   csolve->label = (char *) ML_allocate(size*sizeof(char));
   if (csolve->label == NULL) pr_error("Not enough space in ML_CSolve_Set_Label\n");
   strncpy(csolve->label,label,size);
   return(1);
}

/* ************************************************************************* */
/* This subroutine calls the SuperLU subroutine to perform LU                */
/* factorization of a given matrix                                           */
/* ------------------------------------------------------------------------- */

int ML_CSolve_Aggr(ML_Solver *vsolver,int ilen,double *x,int olen,double *rhs)
{
   int            i, n, N_local, offset;
   double         *local_x, *local_rhs;
   ML_Comm        *comm;
   ML_Solver      *solver;
   ML             *ml_ptr;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen )
   {
      printf("ML_CSolve_Aggr ERROR : lengths not matched.\n");
      exit(1);
   }
   solver = (ML_Solver *) vsolver;

   ml_ptr   = (ML *) solver->void_params1;
   comm     = (ML_Comm *) solver->void_params2;
   N_local  = (int) solver->dble_params1[0];
   offset   = (int) solver->dble_params1[1];
   n        = (int) solver->dble_params1[2];

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   ML_memory_alloc((void**) &local_x,   n*sizeof(double),"LU2" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) comm, local_rhs, &i, n);
   for ( i = 0; i < n; i++ ) local_x[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* solve                                                         */
   /* ------------------------------------------------------------- */

   ML_Solve_AMGV( ml_ptr, local_rhs, local_x );

   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_x[i+offset];

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free( (void **) &local_x );
   ML_memory_free( (void **) &local_rhs );
   solver->reuse_flag = 1;
   return 0;
}

/* ************************************************************************* */
/* destructor for Aggregation based coarse solver                            */
/* ------------------------------------------------------------------------- */

int ML_CSolve_Clean_Aggr( void *vsolver, ML_CSolveFunc *func)
{
   ML_Solver *solver;

   solver = (ML_Solver *) vsolver;
   if ( solver != NULL )
   {
      if ( solver->dble_params1 != NULL )
         ML_memory_free( (void**) &(solver->dble_params1) );
      solver->dble_params1 = NULL;
      if ( solver->Mat1 != NULL )
      {
         ML_Matrix_DCSR_Destroy( (ML_Matrix_DCSR *) (solver->Mat1) );
         ML_memory_free( &(solver->Mat1) );
      }
      solver->Mat1 = NULL;
      if ( solver->void_params1 != NULL )
         ML_Destroy((ML**)&(solver->void_params1));
      solver->void_params1 = NULL;
      ML_memory_free( (void**) &solver );
   }
   else ML_avoid_unused_param( (void *) func);
   return 0;
}




#ifdef WKC
/* WKC  EPETRA STUFF TO FOLLOW */

/* *********************************************************************** */
/* perform solve                                                           */
/* ----------------------------------------------------------------------- */
#include <iostream>

extern int ML_SuperLU_Solve_WKC(void *vsolver,int ilen,double *x,int olen,
                            double *rhs);


int ML_CSolve_Apply(ML_CSolve *csolve, int inlen, Epetra_MultiVector &ep_din, 
                    int outlen, Epetra_MultiVector &ep_dout )
{

   double ** pp_din;
   double ** pp_dout;
   ep_din.ExtractView ( &pp_din );
   ep_dout.ExtractView ( &pp_dout );

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double t0;
   t0 = GetClock();
#endif
   if (csolve->func->func_ptr == NULL) 
      pr_error("ML_CSolve_Apply error : coarse solver not defined\n");

   if ( (void *) csolve->func->func_ptr == (void *)ML_SuperLU_Solve )
     ML_SuperLU_Solve_WKC ((ML_Solver *)csolve->data, inlen, (double *) &ep_din, outlen, 
                         (double *) &ep_dout); 
   else {
       for ( int KK = 0 ; KK != ep_din.NumVectors() ; KK++ ) {
         double *din = pp_din[KK];
         double *dout = pp_dout[KK];

         csolve->func->func_ptr ((ML_Solver *)csolve->data, inlen, din, outlen, dout); 
       }
   }
   

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   csolve->apply_time += (GetClock() - t0);
#endif
   return 0;
}

#endif
