#define MatrixProductHiptmair
/*
#define NoDampingFactor
#define MatrixFreeHiptmair
*/

/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* Functions for the ML_Smoother structure                                   */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)            */
/* Date          : April, 2000                                               */
/* ************************************************************************* */
/* ML_Smoother_Create                                                        */
/* ML_Smoother_Init                                                          */
/* ML_Smoother_Destroy                                                       */
/* ML_Smoother_Clean                                                         */
/* ML_Smoother_Set_Label                                                     */
/* ML_Smoother_Apply                                                         */
/* ML_Smoother_Set                                                           */
/* ML_Smoother_Jacobi                                                        */
/* ML_Smoother_GaussSeidel                                                   */
/* ML_Smoother_SGS                                                           */
/* ML_Smoother_BlockGS                                                       */
/* ML_Smoother_VBlockJacobi                                                  */
/* ML_Smoother_VBlockSGS                                                     */
/* ML_Smoother_VBlockSGSSequential                                           */
/* ML_Smoother_VBlockKrylovJacobi                                            */
/* ML_Smoother_OverlappedILUT                                                */
/* ML_Smoother_VBlockAdditiveSchwarz                                         */
/* ML_Smoother_VBlockMultiplicativeSchwarz                                   */
/* ML_Smoother_ParaSails                                                     */
/* ************************************************************************* */
/*  Methods available                                                        */
/*   - weighted Jacobi                                                       */
/*   - Gauss Seidel                                                          */
/*   - symmetric Gauss Seidel                                                */
/*   - symmetric Gauss Seidel (sequential)                                   */
/*   - block Gauss Seidel                                                    */
/*   - variable block Jacobi                                                 */
/*   - variable block Symmetric Gauss Seidel                                 */
/*   - variable block Symmetric Gauss Seidel (sequential)                    */
/*   - variable block Jacobi with Krylov                                     */
/*   - overlapped ILUT                                                       */
/*   - variable block additive Schwarz                                       */
/*   - variable block multiplicative Schwarz                                 */
/*   - sparse approximate inverse                                            */
/* ************************************************************************* */

#include <stdlib.h>
#include "ml_smoother.h"
#include "ml_lapack.h"
#include "ml_utils.h"

#ifdef out
/* ************************************************************************* */
/* include files for SuperLU and MPI                                         */
/* ------------------------------------------------------------------------- */

#ifdef SUPERLU
#include "dsp_defs.h"
#include "util.h"
#else
#ifdef DSUPERLU
#include "mpi.h"
#include "superlu_ddefs.h"
#endif
#endif
#endif

/* ************************************************************************* */
/* Constructor                                                               */
/* ************************************************************************* */

int ML_Smoother_Create(ML_Smoother **sm, ML_1Level *mylevel)
{
   ML_Smoother *ml_sm;
   ML_memory_alloc((void**) sm, sizeof(ML_Smoother), "SM1" );
   ml_sm = (*sm);
   ML_Smoother_Init(ml_sm, mylevel);
   return 0;
} 

/* ************************************************************************* */
/* Initialize                                                                */
/* ************************************************************************* */

int ML_Smoother_Init(ML_Smoother *ml_sm, ML_1Level *mylevel)
{
   ml_sm->ML_id = ML_ID_SMOOTHER; 
   ml_sm->my_level = mylevel;
   ml_sm->ntimes = 0;
   ml_sm->omega = 0;
   ml_sm->init_guess = ML_NONZERO;
   ml_sm->tol = 0;
   ML_memory_alloc((void**)&(ml_sm->smoother),sizeof(ML_SmootherFunc),"SF2");
   ml_sm->smoother->ML_id = ML_EMPTY; 
   ml_sm->smoother->internal = NULL; 
   ml_sm->smoother->external = NULL; 
   ml_sm->smoother->data = NULL; 
   ml_sm->data_destroy = NULL;
   ml_sm->build_time = 0.0;
   ml_sm->apply_time = 0.0;
   ml_sm->label      = NULL;
   ml_sm->pre_or_post = 0;
   ml_sm->envelope = NULL;
   return 0;
} 

/* ************************************************************************* */
/* Destructor                                                                */
/* ************************************************************************* */

int ML_Smoother_Destroy(ML_Smoother **sm)
{
   ML_Smoother *ml_sm;
   ml_sm = (*sm);
   ML_Smoother_Clean(ml_sm);
   ML_memory_free( (void**) sm );
   (*sm) = NULL; 
   return 0;
}

/* ************************************************************************* */
/* Cleaner                                                                   */
/* ************************************************************************* */

int ML_Smoother_Clean(ML_Smoother *ml_sm)
{
#ifdef ML_TIMING_DETAILED
   int    nprocs, mypid;
   double t1;
#endif
   /* ML_Sm_BGS_Data     *ml_data; */
   ML_Sm_ILUT_Data    *ilut_data;
#ifdef out
   ML_Sm_Schwarz_Data *schwarz_data;
#endif

#ifdef ML_TIMING_DETAILED
   if (ml_sm->my_level != NULL) {
     mypid  = ml_sm->my_level->comm->ML_mypid;
     nprocs = ml_sm->my_level->comm->ML_nprocs;
     t1 = ML_gsum_double(ml_sm->build_time, ml_sm->my_level->comm);
     t1 = t1/((double) nprocs);
     if ( (ml_sm->label != NULL) && ( t1 != 0.0) && mypid == 0)
       printf(" Build time for %s (average) \t= %e\n",ml_sm->label,t1);
     t1 = ML_gmax_double(ml_sm->build_time, ml_sm->my_level->comm);
     if ( (ml_sm->label != NULL) && ( t1 != 0.0) && mypid == 0)
       printf(" Build time for %s (maximum) \t= %e\n",ml_sm->label,t1);
     t1 = - ml_sm->build_time;
     t1 = ML_gmax_double(t1, ml_sm->my_level->comm);
     t1 = - t1;
     if ( (ml_sm->label != NULL) && ( t1 != 0.0) && mypid == 0)
       printf(" Build time for %s (minimum) \t= %e\n",ml_sm->label,t1);
   }
#endif

#ifdef ML_TIMING_DETAILED
   if (ml_sm->label != NULL) 
   {
      t1 = ML_gsum_double(ml_sm->apply_time, ml_sm->my_level->comm);
      t1 = t1/((double) nprocs);
      if ( (mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (average) \t= %e\n",ml_sm->label,t1);
      t1 = ML_gmax_double(ml_sm->apply_time, ml_sm->my_level->comm);
      if ( (mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (maximum) \t= %e\n",ml_sm->label,t1);
      t1 = - ml_sm->apply_time;
      t1 = ML_gmax_double(t1, ml_sm->my_level->comm);
      t1 = - t1;
      if ( (mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (minimum) \t= %e\n",ml_sm->label,t1);
   }
#endif

   ml_sm->ML_id = -1; 
   ml_sm->my_level = NULL;
   ml_sm->ntimes = 0;
   ml_sm->omega = 0;
   ml_sm->pre_or_post = 0;
   ml_sm->init_guess = ML_NONZERO;
   ml_sm->tol = 0;
   if ((ml_sm->data_destroy != NULL) && (ml_sm->smoother->data != NULL)) {
      ml_sm->data_destroy( ml_sm->smoother->data );
      ml_sm->smoother->data = NULL;
   }
   ml_sm->data_destroy = NULL;
   /* ********************************************************* */
   /* ********************************************************* */
   /* These have already been freed with the data destroy above */
   /* for VBlockSGS and VBlockJacobi!!!!!!!!!!!!!!!!!!!!!!!!!!  */
   /* Charles, we should probably talk about this --------------*/
   /* NOTE: the code below will not work if the pre and post    */
   /* smoother are sharing the same data pointer.               */
   /* ********************************************************* */
   /* ********************************************************* */
   /* ********************************************************* 
   if ( ml_sm->smoother->internal == ML_Smoother_VBlockSGS ||
        ml_sm->smoother->internal == ML_Smoother_VBlockJacobi )
   {
      ml_data = ml_sm->smoother->data;
      if ( ml_data != NULL )
      {
         ML_Smoother_Destroy_BGS_Data((ml_data));
         ml_sm->smoother->data = NULL;
      }
   }
   */
   if ( ml_sm->smoother->internal == ML_Smoother_OverlappedILUT &&
        ml_sm->smoother->data != NULL )
   {
     printf("This should be switched to using the data_destroy field\n");
     printf("Charles ... can we talk about fixing this -RST\n");
#ifdef ML_CPP
      ilut_data = (ML_Sm_ILUT_Data *)ml_sm->smoother->data;
#else
      ilut_data = ml_sm->smoother->data;
#endif 
      ML_Smoother_Destroy_ILUT_Data(ilut_data);
      ml_sm->smoother->data = NULL;
   }
   if ( ml_sm->smoother->internal == ML_Smoother_SGS &&
        ml_sm->smoother->data != NULL )
   {
     printf("This should be switched to using the data_destroy field\n");
     printf("Charles ... can we talk about fixing this -RST\n");
#ifdef ML_FAST
      ilut_data = ml_sm->smoother->data;
      ML_Smoother_Destroy_ILUT_Data(ilut_data);
      ml_sm->smoother->data = NULL;
#endif
   }
   if (((ml_sm->smoother->internal == ML_Smoother_VBlockAdditiveSchwarz) ||
        (ml_sm->smoother->internal == ML_Smoother_VBlockMultiplicativeSchwarz)) &&
        ml_sm->smoother->data != NULL )
   {
     printf("This should be switched to using the data_destroy field\n");
     printf("Charles ... can we talk about fixing this -RST\n");
#ifdef out
      schwarz_data = ml_sm->smoother->data;
      ML_Smoother_Destroy_Schwarz_Data(&(schwarz_data));
#endif
      ml_sm->smoother->data = NULL;
   }

   ML_memory_free((void**)&(ml_sm->smoother));
   if (ml_sm->label != NULL) { ML_free(ml_sm->label); ml_sm->label = NULL; }

   return 0;
}

/* ************************************************************************* */
/* set a label                                                               */
/* ************************************************************************* */

int ML_Smoother_Set_Label( ML_Smoother *smoo, char *label)
{
  int size;

   if (smoo->label != NULL) { ML_free(smoo->label); smoo->label = NULL; }
   size = strlen(label) + 1;
   smoo->label = (char *) ML_allocate(size*sizeof(char));
   if (smoo->label == NULL) 
      pr_error("Not enough space in ML_Smoother_Set_Label\n");
   strncpy(smoo->label,label,size);
   return(1);
}

/* ************************************************************************* */
/* smoothing operation                                                       */
/* ************************************************************************* */

int ML_Smoother_Apply(ML_Smoother *pre, int inlen, double sol[], 
                      int outlen, double rhs[], int init_guess)
{
   int         i, n;
   double      temp, *res, tol;
   ML_Operator *Amat;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double      t0;
   t0 = GetClock();
#endif

   if (pre->smoother->ML_id == ML_EMPTY) return 1;
   pre->init_guess = init_guess;

   if (pre->smoother->ML_id == ML_EXTERNAL)
        pre->smoother->external(pre->smoother->data,inlen,sol,outlen,rhs);
   else 
   {
      if (pre->ntimes == ML_CONVERGE) 
      {
         Amat = pre->my_level->Amat;
         n    = Amat->outvec_leng; 
         res  = (double *) ML_allocate( (n+1)*sizeof(double) );
         temp = sqrt(ML_gdot(n, rhs, rhs, pre->my_level->comm));
         tol  = temp*pre->tol;
         pre->ntimes = 100;
         while ( temp > tol ) 
         {
            pre->smoother->internal(pre,n,sol,n, rhs);
            ML_Operator_Apply(Amat, n, sol, n, res);
            for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];
            temp = sqrt(ML_gdot(n, res, res, pre->my_level->comm));
         }
         pre->ntimes = ML_CONVERGE;
         ML_free(res);
      }
      else pre->smoother->internal(pre,inlen,sol,outlen,rhs);
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   pre->apply_time += (GetClock() - t0);
#endif
   return 1;
}

/* ************************************************************************* */
/* set smoother                                                              */
/* ************************************************************************* */

int ML_Smoother_Set(ML_Smoother *smoo,int internal_or_external,void *data,
                    int (*internal)(void*,int,double*,int,double *),
                    int (*external)(void*,int,double*,int,double *), 
                    int ntimes, double omega, char *str)
{
   if (internal_or_external == ML_EXTERNAL) 
   {
      smoo->smoother->external = external;
      smoo->smoother->ML_id = ML_EXTERNAL;
   } 
   else 
   {
      smoo->smoother->internal = internal;
      smoo->smoother->ML_id = ML_INTERNAL;
   }
   smoo->smoother->data= data;
   smoo->ntimes = ntimes;
   smoo->omega = omega;
   if (str != NULL) ML_Smoother_Set_Label( smoo, str);
   return 0;
}

/* ************************************************************************* */
/* ML-supplied smoothers                                                     */
/* ************************************************************************* */
/* ************************************************************************* */

/* ************************************************************************* */
/* weighted Jacobi smoother                                                  */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Jacobi(void *sm,int inlen,double x[],int outlen,double rhs[])
{
   int i, j, n, *cols, allocated_space;
   ML_Operator *Amat;
   double *res,omega, *diagonal, *vals, *tdiag, *res2 = NULL;
   double r_z_dot, p_ap_dot;
   ML_Smoother  *smooth_ptr;

#ifdef ML_DEBUG_SMOOTHER
   ML_Comm      *comm;
   double res_norm;
#endif

   smooth_ptr = (ML_Smoother *) sm;

#ifdef ML_DEBUG_SMOOTHER
   comm = smooth_ptr->my_level->comm;
#endif

   omega = smooth_ptr->omega;
   Amat = smooth_ptr->my_level->Amat;
   if (Amat->matvec->ML_id == ML_EMPTY) {
         pr_error("Error(ML_Jacobi): Need matvec\n");
	 ML_avoid_unused_param((void *) &inlen);
	 ML_avoid_unused_param((void *) &outlen);
   }

   /* ----------------------------------------------------------------- */
   /* extract diagonal using getrow function if not found               */
   /* ----------------------------------------------------------------- */

   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(ML_Jacobi): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&n) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols); 
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < n; j++) 
               if (cols[j] == i) tdiag[i] = vals[j];
         }
         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);

   n     = Amat->outvec_leng;
   res   = (double *) ML_allocate(n*sizeof(double));
   if (smooth_ptr->omega == ML_ONE_STEP_CG) {
      res2  = (double *) ML_allocate(n*sizeof(double));
   }

#ifdef ML_DEBUG_SMOOTHER
   if (res2 == NULL) res2  = (double *) ML_allocate(n*sizeof(double));
   if (comm->ML_mypid == 0)
     printf("    ML_Jacobi, omega = %e\n", omega);
   res_norm = sqrt(ML_gdot(n, rhs, rhs, comm));
   if (comm->ML_mypid == 0)
     printf("       entering:  |rhs| = %e\n", res_norm);
   res_norm = sqrt(ML_gdot(n, diagonal, diagonal, comm));
   if (comm->ML_mypid == 0)
     printf("       entering: |diag| = %e\n", res_norm);
#endif
   for (j = 0; j < smooth_ptr->ntimes; j++) {
      ML_Operator_Apply(Amat, n, x, n, res);
      for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];
#ifdef ML_DEBUG_SMOOTHER
      res_norm = sqrt(ML_gdot(n, res, res, comm));
      if (comm->ML_mypid == 0)
        printf("               %d: |res| = %e\n", j, res_norm);
#endif
      for (i = 0; i < n; i++) res[i] /= diagonal[i];

      if (smooth_ptr->omega == ML_ONE_STEP_CG) {
         /* Compute damping parameter that corresonds to one step of CG. */
         r_z_dot = 0.;
         for (i = 0; i < n; i++) r_z_dot += res[i]*res[i]*diagonal[i];
         r_z_dot = ML_gsum_double(r_z_dot, smooth_ptr->my_level->comm);
         ML_Operator_Apply(Amat, n, res, n, res2);
         p_ap_dot = ML_gdot(n, res, res2, smooth_ptr->my_level->comm);
         if (p_ap_dot != 0.0) omega = r_z_dot/p_ap_dot;
         else omega = 1.;
      }
      for (i = 0; i < n; i++) x[i] += omega*res[i];

#ifdef ML_DEBUG_SMOOTHER
      ML_Operator_Apply(Amat, n, x, n, res2);
      for ( i = 0; i < n; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(n, res2, res2, comm));
      if (comm->ML_mypid == 0)
        printf("              %d: rnorm = %e\n", j, res_norm);
#endif
   }
   if (res2 != NULL) ML_free(res2);
   ML_free(res);
   return 0;
}

/* ************************************************************************* */
/* Gauss-Seidel smoother                                                     */
/* ------------------------------------------------------------------------- */

int ML_Smoother_GaussSeidel(void *sm, int inlen, double x[], int outlen, 
                            double rhs[])
{
   int iter, i, j, length, allocated_space, *cols, col;
   double dtemp, diag_term, *vals;
#ifdef ML_DEBUG_SMOOTHER
   double *res2, res_norm;
#endif
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2, omega;
   ML_Smoother  *smooth_ptr;
   smooth_ptr = (ML_Smoother *) sm;

   /* ----------------------------------------------------------------- */
   /* set up                                                            */
   /* ----------------------------------------------------------------- */

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;
   omega = smooth_ptr->omega;

   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_GaussSeidel): Need getrow() for GS smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_GaussSeidel(): Not enough space\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for GS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
				   *sizeof(double));
      if (x2 == NULL) 
      {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;

#ifdef ML_DEBUG_SMOOTHER
   res2 = (double*) ML_allocate(Nrows * sizeof(double));
#endif

   /* ----------------------------------------------------------------- */
   /* perform smoothing                                                 */
   /* ----------------------------------------------------------------- */

   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

      for (i = 0; i < Nrows; i++) 
      {
         dtemp = 0.0;
         diag_term = 0.0;
         ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
         for (j = 0; j < length; j++) 
         {
            col = cols[j];
            if (col == i) diag_term = vals[j];
            dtemp += vals[j]*x2[col];
         }
         if (diag_term == 0.0)
            pr_error("Error: GS() can not be used with a zero diagonal\n");

         x2[i] += omega * (rhs[i] - dtemp)/diag_term;
      }
#ifdef ML_DEBUG_SMOOTHER
      ML_Operator_Apply(Amat, Nrows, x2, Nrows, res2);
      for ( i = 0; i < Nrows; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(Nrows, res2, res2, comm));
      printf("      GS : iter = %2d, rnorm = %e\n", iter, res_norm);
#endif
   }
   if (getrow_comm != NULL) 
   {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+2) 
      Amat->max_nz_per_row = allocated_space;

#ifdef ML_DEBUG_SMOOTHER
   ML_free(res2);
#endif
   ML_free(vals); ML_free(cols);

   return 0;
}

/* ************************************************************************* */
/* Symmetric Gauss-Seidel smoother                                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_SGS(void *sm,int inlen,double x[],int outlen, double rhs[])
{
   int iter, i, j, length, allocated_space, *cols, col;
   double dtemp, diag_term, *vals, omega;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
#ifdef ML_DEBUG_SMOOTHER
   double *res2, res_norm, init_norm;
#endif
   ML_Smoother  *smooth_ptr;
#ifdef ML_FAST
   int *rptr, nnz;
   ML_Sm_ILUT_Data *ilut_data;
#endif

   /* ----------------------------------------------------------------- */
   /* fetch data                                                        */
   /* ----------------------------------------------------------------- */

   smooth_ptr = (ML_Smoother *) sm;
   omega = smooth_ptr->omega;
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

#ifdef ML_FAST
   ilut_data  = smooth_ptr->smoother->data;
#endif

   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_SGS): Need getrow() for SGS smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   /* ----------------------------------------------------------------- */
   /* if matrix not found, get it (more efficient implementation, but   */
   /* take up much more memory)                                         */
   /* ----------------------------------------------------------------- */

#ifdef ML_FAST
   ilut_data  = smooth_ptr->smoother->data;
   if ( ilut_data != NULL )
   {
      rptr = ilut_data->mat_ia;
      cols = ilut_data->mat_ja;
      vals = ilut_data->mat_aa;
   }
   else
   {
      allocated_space = (Amat->max_nz_per_row+2) * Nrows + 1;
      rptr = (int    *) ML_allocate(allocated_space*sizeof(int   ));
      cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
      vals = (double *) ML_allocate(allocated_space*sizeof(double));
      if (vals == NULL) pr_error("ML_Smoother_SGS ERROR: Not enough space\n");
      nnz = 0;
      rptr[0] = nnz;
      for (i = 0; i < Nrows; i++)
      {
         ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,&length,nnz);
         nnz += length;
         rptr[i+1] = nnz;
      }
      if ( nnz > allocated_space ) printf("ERROR in ML_Smoother_SGS memory\n");
      ML_memory_alloc((void**)&ilut_data, sizeof(ML_Sm_ILUT_Data), "ILU");
      ilut_data->Nrows = Nrows;
      ilut_data->mat_ia = rptr;
      ilut_data->mat_ja = cols;
      ilut_data->mat_aa = vals;
      ilut_data->getrow_comm = NULL;
      ilut_data->fillin = 0;
      smooth_ptr->smoother->data = ilut_data;
   }
#else
   allocated_space = Amat->max_nz_per_row + 2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_SymGaussSeidel: Not enough space\n");
#endif

   /* ----------------------------------------------------------------- */
   /* get an extended buffer for ghost data                             */
   /* ----------------------------------------------------------------- */

   getrow_comm = Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
				   *sizeof(double));
      if (x2 == NULL) 
      {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;

   /* ----------------------------------------------------------------- */
   /* perform smoothing (more memory efficient version)                 */
   /* ----------------------------------------------------------------- */

#ifdef ML_DEBUG_SMOOTHER
   res2 = (double*) ML_allocate(Nrows * sizeof(double));
#endif

#ifndef ML_FAST
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,
                          ML_OVERWRITE,smooth_ptr->envelope);

      for (i = 0; i < Nrows; i++) 
      {
         dtemp = 0.0;
         diag_term = 0.0;
         ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
         for (j = 0; j < length; j++) 
         {
            col = cols[j];
            dtemp += vals[j]*x2[col];
            if (col == i) diag_term = vals[j];
         }
         if (diag_term != 0.0)
            x2[i] += omega*(rhs[i] - dtemp)/diag_term;
      }
#ifdef ML_DEBUG_SMOOTHER
      ML_Operator_Apply(Amat, Nrows, x2, Nrows, res2);
      for ( i = 0; i < Nrows; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(Nrows, res2, res2, comm));
      printf("      SGS (for ) : iter = %2d, rnorm = %15.10e\n", iter,res_norm);
      res_norm = sqrt(ML_gdot(Nrows, x2, x2, comm));
      printf("                               ||x|| = %15.10e\n", iter, res_norm);
#endif

#ifdef ML_OBSOLETE
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,sm->envelope);
#endif

      for (i = Nrows-1; i >= 0; i--) 
      {
         dtemp = 0.0;
         diag_term = 0.0;
         ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
         for (j = 0; j < length; j++) 
         {
            col = cols[j];
            dtemp += vals[j]*x2[col];
            if (col == i) diag_term = vals[j];
         }
         if (diag_term != 0.0)
            x2[i] += omega*(rhs[i] - dtemp)/diag_term;
      }

#ifdef ML_DEBUG_SMOOTHER
      ML_Operator_Apply(Amat, Nrows, x2, Nrows, res2);
      for ( i = 0; i < Nrows; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(Nrows, res2, res2, comm));
      printf("      SGS (back) : iter = %2d, rnorm = %15.10e\n",
             iter, res_norm);
      res_norm = sqrt(ML_gdot(Nrows, x2, x2, comm));
      printf("                               ||x|| = %15.10e\n",
             iter, res_norm);
#endif
   }
#else

   /* ----------------------------------------------------------------- */
   /* perform smoothing (more CPU efficient version)                    */
   /* ----------------------------------------------------------------- */

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,sm->envelope);

      for (i = 0; i < Nrows; i++)
      {
         dtemp = rhs[i];
         diag_term = 0.0;
         for (j = rptr[i]; j < rptr[i+1]; j++)
         {
            col = cols[j];
            dtemp -= (vals[j]*x2[col]);
            if (col == i) diag_term = vals[j];
         }
         if (diag_term != 0.0)
            x2[i] += (omega * dtemp / diag_term);
      }

      for (i = Nrows-1; i >= 0; i--)
      {
         dtemp = rhs[i];
         diag_term = 0.0;
         for (j = rptr[i]; j < rptr[i+1]; j++)
         {
            col = cols[j];
            dtemp -= (vals[j]*x2[col]);
            if (col == i) diag_term = vals[j];
         }
         if (diag_term != 0.0)
            x2[i] += (omega * dtemp / diag_term);
      }
   }
#endif

   /* ----------------------------------------------------------------- */
   /* clean up                                                          */
   /* ----------------------------------------------------------------- */

#ifdef ML_DEBUG_SMOOTHER
   ML_free(res2);
#endif

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
#ifndef ML_FAST
   if (allocated_space != Amat->max_nz_per_row+2) {
      Amat->max_nz_per_row = allocated_space;
   }
   ML_free(vals); ML_free(cols);
#endif

   return 0;
#ifdef ML_DEBUG_SMOOTHER
#undef ML_DEBUG_SMOOTHER
#endif
}

/* ************************************************************************* */
/* Hiptmair smoother with matrix products                                    */
/* ------------------------------------------------------------------------- */
#ifdef MatrixProductHiptmair
/* Block Hiptmair Smoother */
int ML_Smoother_BlockHiptmair(void *sm, int inlen, double x[], int outlen, 
                            double rhs[])
{
  int iter, kk, Nrows, init_guess; 
   ML_Operator *Tmat, *Tmat_trans, *TtATmat, *Ke_mat;
   ML_Smoother  *smooth_ptr;
   ML_Sm_BlockHiptmair_Data *dataptr;
   ML_Comm_Envelope *envelope;
   int reduced_smoother_flag;
#ifdef ML_DEBUG_SMOOTHER
   int i,j;
   double *res2, res_norm;
   double t1,t2,t3,t4;
#endif
   double *res_edge;
   double *res_edge1, *res_edge2;
   double *rhs_nodal1, *rhs_nodal2;
   double *x_nodal1, *x_nodal2;
   double *edge_update1, *edge_update2;

   smooth_ptr = (ML_Smoother *) sm;


   Ke_mat = smooth_ptr->my_level->Amat;
   Nrows = Ke_mat->getrow->Nrows;

   /* pointer to private smoother data */
   dataptr = (ML_Sm_BlockHiptmair_Data *) smooth_ptr->smoother->data;

   Tmat = (ML_Operator *) dataptr->Tmat;
   Tmat_trans  = (ML_Operator *) dataptr->Tmat_trans;
   TtATmat = (ML_Operator *) dataptr->TtATmat;
   /*
      If true, do cheaper smoothing by smoothing (edges/nodes) on the pre
      smooth and smoothing (nodes/edges) on the post smooth.
   */
   reduced_smoother_flag = (int) dataptr->reduced_smoother;

   res_edge = (double *) dataptr->res_edge;
   res_edge1 = (double *) dataptr->res_edge1;
   res_edge2 = (double *) dataptr->res_edge2;
   edge_update1 = (double *) dataptr->edge_update1;
   edge_update2 = (double *) dataptr->edge_update2;
   rhs_nodal1 = (double *) dataptr->rhs_nodal1;
   rhs_nodal2 = (double *) dataptr->rhs_nodal2;
   x_nodal1 = (double *) dataptr->x_nodal1;
   x_nodal2 = (double *) dataptr->x_nodal2;

   if (Ke_mat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Hiptmair): Need getrow() for Hiptmair smoother\n");

#ifdef ML_DEBUG_SMOOTHER
   t1 = sqrt((ML_gdot(Nrows, x, x, Tmat_trans->comm)));
   t2 = sqrt((ML_gdot(Nrows, rhs, rhs, Tmat_trans->comm)));
   if (Tmat_trans->comm->ML_mypid == 0) {
     printf("\n--------------------------------\n");
     printf("Coming into matrix Hiptmair\n");
     printf("\t||x|| = %15.10e\n",  t1);
     printf("\t||rhs|| = %15.10e\n", t2);
     fflush(stdout);
   }
#endif

   ML_Comm_Envelope_Create(&envelope);
   smooth_ptr->envelope = envelope;

   /*
   if (smooth_ptr->pre_or_post != ML_TAG_PRESM)
      smooth_ptr->pre_or_post = ML_TAG_PRESM;
   else
      smooth_ptr->pre_or_post = ML_TAG_POSTSM;
   */
   
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
   printf("Coming into Hiptmair: pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif

   ML_Comm_Envelope_Set_Tag(envelope, smooth_ptr->my_level->levelnum,
                            smooth_ptr->pre_or_post); 

   init_guess = smooth_ptr->init_guess;
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (reduced_smoother_flag)
      {
         if (smooth_ptr->pre_or_post == ML_TAG_PRESM)
         {
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
            printf("Hiptmair: pre edge smoothing\n");
#endif
#endif
            ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			                  inlen, x, outlen, rhs, init_guess);
         }
      }
      else
      {
         ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, init_guess);
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
         printf("\t(1) smooth_ptr->pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif
      }
      init_guess = ML_NONZERO;

      for (kk = 0; kk < TtATmat->invec_leng; kk++) {
	x_nodal1[kk] = 0.;
	x_nodal2[kk] = 0.;
      }

      ML_Comm_Envelope_Increment_Tag(envelope);
   
      /* calculate initial residual */ 
      ML_Operator_Apply(Ke_mat, Ke_mat->invec_leng,
                        x, Ke_mat->outvec_leng,res_edge);
      /*
#ifdef ML_DEBUG_SMOOTHER
      ML_DVector_Print(Nrows, res_edge, "Atimesx", Ke_mat->comm);
      exit(1);
#endif
      printf("\t%d: ||A*x|| = %20.15e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Nrows, res_edge, res_edge, Tmat_trans->comm)));
      */
      for (kk = 0; kk < Nrows; kk++) res_edge[kk] = rhs[kk] - res_edge[kk];


      /* split residual for separate nodal calculations */

      for (kk = 0; kk < Nrows/2; kk++) res_edge1[kk] = -res_edge[kk];
      for (kk = Nrows/2; kk < Nrows; kk++) res_edge2[kk-Nrows/2] = res_edge[kk];
   
#ifdef ML_DEBUG_SMOOTHER
      t1 = sqrt(ML_gdot(Nrows, x, x, Tmat_trans->comm));
      t2 = sqrt(ML_gdot(Nrows,res_edge,res_edge,Tmat_trans->comm));
      t3 = sqrt(ML_gdot(Nrows/2,res_edge1,res_edge1,Tmat_trans->comm));
      t4 = sqrt(ML_gdot(Nrows/2,res_edge2,res_edge2,Tmat_trans->comm));
      if (Tmat_trans->comm->ML_mypid == 0) {
	printf("After SGS on edges\n");
	printf("\t ||x|| = %20.15e\n", t1);
	printf("\t ||res|| = %15.10e\n", t2);
	printf("\t ||res1|| = %15.10e\n", t3);
	printf("\t ||res2| = %15.10e\n", t4);
      }
#endif
   
      /****************************
      * Symmetric sweep on nodes. *
      ****************************/
      ML_Comm_Envelope_Increment_Tag(envelope);

      ML_Operator_Apply(Tmat_trans, Tmat_trans->invec_leng,
                        res_edge1, Tmat_trans->outvec_leng,rhs_nodal1);
      ML_Operator_Apply(Tmat_trans, Tmat_trans->invec_leng,
                        res_edge2, Tmat_trans->outvec_leng,rhs_nodal2);

#ifdef ML_DEBUG_SMOOTHER
      t1 = sqrt(ML_gdot(Tmat_trans->outvec_leng, x_nodal1,x_nodal1,
			Tmat_trans->comm));
      t2 =  sqrt(ML_gdot(Tmat_trans->outvec_leng,x_nodal2,x_nodal2,
			 Tmat_trans->comm));
      t3 =  sqrt(ML_gdot(Tmat_trans->outvec_leng,rhs_nodal1,rhs_nodal1,
			 Tmat_trans->comm));
      t4 =  sqrt(ML_gdot(Tmat_trans->outvec_leng,rhs_nodal2,rhs_nodal2,
			 Tmat_trans->comm));
      if (Tmat_trans->comm->ML_mypid == 0) {
	printf("Before SGS on nodes\n");
	printf("\t ||x_nodal1|| = %15.10e\n", t1);
	printf("\t ||x_nodal2|| = %15.10e\n", t2);
	printf("\t ||rhs_nodal1|| = %15.10e\n", t3);
	printf("\t ||rhs_nodal2|| = %15.10e\n", t4);
      }
#endif

      /* two solves on nodes */

      ML_Smoother_Apply(&(dataptr->ml_nodal->pre_smoother[0]),
			TtATmat->invec_leng, x_nodal1,
			TtATmat->outvec_leng, rhs_nodal2,ML_ZERO);
      ML_Smoother_Apply(&(dataptr->ml_nodal->pre_smoother[0]),
			TtATmat->invec_leng, x_nodal2,
			TtATmat->outvec_leng, rhs_nodal1,ML_ZERO);

#ifdef ML_DEBUG_SMOOTHER

      ML_Operator_Apply(TtATmat, TtATmat->invec_leng,
                        x_nodal1, TtATmat->outvec_leng,res_edge1);
      for (kk = 0; kk < Nrows/2; kk++) 
	res_edge1[kk] = rhs_nodal2[kk] - res_edge1[kk];
      ML_Operator_Apply(TtATmat, TtATmat->invec_leng,
                        x_nodal2, TtATmat->outvec_leng,res_edge2);
      for (kk = 0; kk < Nrows/2; kk++) 
	res_edge2[kk] = rhs_nodal1[kk] - res_edge2[kk];

      t1 = sqrt(ML_gdot(Tmat_trans->outvec_leng,x_nodal1,x_nodal1,
			Tmat_trans->comm));
      t2 = sqrt(ML_gdot(Tmat_trans->outvec_leng,x_nodal2,x_nodal2,
			Tmat_trans->comm));
      t3 = sqrt(ML_gdot(Tmat_trans->outvec_leng,res_edge1,res_edge1,
			Tmat_trans->comm));
      t4 = sqrt(ML_gdot(Tmat_trans->outvec_leng,res_edge2,res_edge2,
			Tmat_trans->comm));

      if (Tmat_trans->comm->ML_mypid == 0) {
	printf("After SGS on nodes\n");
	printf("\t ||x_nodal1|| = %15.10e\n", t1);
	printf("\t ||x_nodal2|| = %15.10e\n", t2);
	printf("\t ||rhs_nodal1|| = %15.10e\n", t3);
	printf("\t ||rhs_nodal2|| = %15.10e\n", t4);
      }
#endif
   
      /************************
      * Update edge solution. *
      ************************/
      ML_Comm_Envelope_Increment_Tag(envelope);

      ML_Operator_Apply(Tmat, Tmat->invec_leng,
                        x_nodal1, Tmat->outvec_leng,edge_update1);
      ML_Operator_Apply(Tmat, Tmat->invec_leng,
                        x_nodal2, Tmat->outvec_leng,edge_update2);
   
      for (kk=0; kk < Nrows/2; kk++) x[kk] += edge_update1[kk];
      for (kk=Nrows/2; kk < Nrows; kk++) x[kk] += edge_update2[kk-Nrows/2];
      
#ifdef ML_DEBUG_SMOOTHER
      t1 = sqrt(ML_gdot(Nrows,x,x,Tmat_trans->comm));
      if (Tmat_trans->comm->ML_mypid == 0) {
	printf("After edge correction\n");
	printf("\t ||x_corrected|| = %15.10e\n", t1);
      }
#endif

      if (reduced_smoother_flag)
      {
         if (smooth_ptr->pre_or_post == ML_TAG_POSTSM)
         {
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
            printf("Hiptmair: post edge smoothing\n");
#endif
#endif
            ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, ML_NONZERO);
         }
      }
      else 
      {
         ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, ML_NONZERO);
      
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
         printf("\t(2) smooth_ptr->pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif
      }
	 
#ifdef ML_DEBUG_SMOOTHER
      t1 = sqrt((ML_gdot(Nrows,x,x,Tmat_trans->comm)));
      if (Tmat_trans->comm->ML_mypid == 0) {
	printf("After updating edge solution\n");
	printf("\t ||x|| = %15.10e\n", t1);
	printf("--------------------------------\n");
      }
#endif

      ML_Comm_Envelope_Increment_Tag(envelope);
   
   } /*for (iter = 0; ...*/

   ML_Comm_Envelope_Destroy(envelope);
   return 0;
}
#ifdef ML_DEBUG_SMOOTHER
#undef ML_DEBUG_SMOOTHER
#endif

/******************************************************************************/

int ML_Smoother_Hiptmair(void *sm, int inlen, double x[], int outlen, 
                            double rhs[])
{
  int iter, kk, Nrows, init_guess;
   ML_Operator *Tmat, *Tmat_trans, *TtATmat, *Ke_mat;
   double *res_edge, *edge_update,
          *rhs_nodal, *x_nodal;
   ML_Smoother  *smooth_ptr;
   ML_Sm_Hiptmair_Data *dataptr;
   ML_Comm_Envelope *envelope;
   int reduced_smoother_flag;
#ifdef ML_DEBUG_SMOOTHER
   int i,j;
   double *res2, res_norm;
#endif

   smooth_ptr = (ML_Smoother *) sm;

   Ke_mat = smooth_ptr->my_level->Amat;
   Nrows = Ke_mat->getrow->Nrows;

   /* pointer to private smoother data */
   dataptr = (ML_Sm_Hiptmair_Data *) smooth_ptr->smoother->data;

   Tmat = (ML_Operator *) dataptr->Tmat;
   Tmat_trans  = (ML_Operator *) dataptr->Tmat_trans;
   TtATmat = (ML_Operator *) dataptr->TtATmat;
   /*
      If true, do cheaper smoothing by smoothing (edges/nodes) on the pre
      smooth and smoothing (nodes/edges) on the post smooth.
   */
   reduced_smoother_flag = (int) dataptr->reduced_smoother;

   res_edge = (double *) dataptr->res_edge;
   edge_update = (double *) dataptr->edge_update;
   rhs_nodal = (double *) dataptr->rhs_nodal;
   x_nodal = (double *) dataptr->x_nodal;

   if (Ke_mat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Hiptmair): Need getrow() for Hiptmair smoother\n");
#ifdef ML_DEBUG_SMOOTHER
   printf("\n--------------------------------\n");
   printf("Coming into matrix Hiptmair\n");
   printf("\t%d: ||x|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
           sqrt((ML_gdot(Nrows, x, x, Tmat_trans->comm))));
   printf("\t%d: ||rhs|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
          sqrt((ML_gdot(Nrows, rhs, rhs, Tmat_trans->comm))));
   fflush(stdout);
#endif

   ML_Comm_Envelope_Create(&envelope);
   smooth_ptr->envelope = envelope;

   /*
   if (smooth_ptr->pre_or_post != ML_TAG_PRESM)
      smooth_ptr->pre_or_post = ML_TAG_PRESM;
   else
      smooth_ptr->pre_or_post = ML_TAG_POSTSM;
   */

#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
   printf("Coming into Hiptmair: pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif

   ML_Comm_Envelope_Set_Tag(envelope, smooth_ptr->my_level->levelnum,
                            smooth_ptr->pre_or_post); 

   init_guess = smooth_ptr->init_guess;
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (reduced_smoother_flag)
      {
         if (smooth_ptr->pre_or_post == ML_TAG_PRESM)
         {
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
            printf("Hiptmair: pre edge smoothing\n");
#endif
#endif
            ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			                  inlen, x, outlen, rhs, init_guess);
         }
      }
      else
      {
         ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, init_guess);
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
         printf("\t(1) smooth_ptr->pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif
      }
      init_guess = ML_NONZERO;

      for (kk = 0; kk < TtATmat->invec_leng; kk++) x_nodal[kk] = 0.;

      ML_Comm_Envelope_Increment_Tag(envelope);
   
      /* calculate initial residual */ 
      ML_Operator_Apply(Ke_mat, Ke_mat->invec_leng,
                        x, Ke_mat->outvec_leng,res_edge);
      for (kk = 0; kk < Nrows; kk++) res_edge[kk] = rhs[kk] - res_edge[kk];
   
#ifdef ML_DEBUG_SMOOTHER
      printf("After SGS on edges\n");
      printf("\t%d: ||x|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Nrows, x, x, Tmat_trans->comm)));
      printf("\t%d: ||res|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Nrows,res_edge,res_edge,Tmat_trans->comm)));
#endif
   
      /****************************
      * Symmetric sweep on nodes. *
      ****************************/
      ML_Comm_Envelope_Increment_Tag(envelope);

      ML_Operator_Apply(Tmat_trans, Tmat_trans->invec_leng,
                        res_edge, Tmat_trans->outvec_leng,rhs_nodal);

#ifdef ML_DEBUG_SMOOTHER
      printf("Before SGS on nodes\n");

      printf("\t%d: ||x_nodal|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Tmat_trans->outvec_leng,
                          x_nodal,x_nodal,Tmat_trans->comm)));
      printf("\t%d: ||rhs_nodal|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Tmat_trans->outvec_leng,
                          rhs_nodal,rhs_nodal,Tmat_trans->comm)));
#endif

      ML_Smoother_Apply(&(dataptr->ml_nodal->pre_smoother[0]),
			TtATmat->invec_leng, x_nodal,
			TtATmat->outvec_leng, rhs_nodal,ML_ZERO);

#ifdef ML_DEBUG_SMOOTHER
      printf("After SGS on nodes\n");
      printf("\t%d: ||x_nodal|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt(ML_gdot(Tmat_trans->outvec_leng,
                  x_nodal,x_nodal,Tmat_trans->comm)));
      printf("\t%d: ||rhs_nodal|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
            sqrt(ML_gdot(Tmat_trans->outvec_leng,rhs_nodal,
                 rhs_nodal,Tmat_trans->comm)));
#endif
   
      /************************
      * Update edge solution. *
      ************************/
      ML_Comm_Envelope_Increment_Tag(envelope);
      ML_Operator_Apply(Tmat, Tmat->invec_leng,
                        x_nodal, Tmat->outvec_leng,edge_update);
   
      for (kk=0; kk < Nrows; kk++) x[kk] += edge_update[kk];

      if (reduced_smoother_flag)
      {
         if (smooth_ptr->pre_or_post == ML_TAG_POSTSM)
         {
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
            printf("Hiptmair: post edge smoothing\n");
#endif
#endif
            ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, ML_NONZERO);
         }
      }
      else 
      {
         ML_Smoother_Apply(&(dataptr->ml_edge->pre_smoother[0]),
			               inlen, x, outlen, rhs, ML_NONZERO);
#ifdef ML_DEBUG_SMOOTHER
#ifdef PRINTITNOW
         printf("\t(2) smooth_ptr->pre_or_post = %d\n",smooth_ptr->pre_or_post);
#endif
#endif
      }

#ifdef ML_DEBUG_SMOOTHER
      printf("After updating edge solution\n");
      printf("\t%d: ||x|| = %15.10e\n", Tmat_trans->comm->ML_mypid,
             sqrt((ML_gdot(Nrows,x,x,Tmat_trans->comm))));
      printf("--------------------------------\n");
#endif

      ML_Comm_Envelope_Increment_Tag(envelope);
   
   } /*for (iter = 0; ...*/

   ML_Comm_Envelope_Destroy(envelope);
   return 0;
}
#ifdef ML_DEBUG_SMOOTHER
#undef ML_DEBUG_SMOOTHER
#endif
#endif /* ifdef MatrixProductHiptmair */

#ifdef MatrixFreeHiptmair

/******************************************************************************/

int ML_Smoother_Hiptmair(void *sm, int inlen, double x[], int outlen, 
                         double rhs[])
{
   int i, j, kk, length, allocated_space = 0, *cols = NULL, Nrows;
   int iter;
   double *vals = NULL, numer, denom, alpha, *g, *res;
   double *local_x;
   ML_Operator *Amat, *Tmat_trans, *ATmat_trans;
   double *TtAT_diag;
   ML_Comm *comm;
   ML_Smoother  *smooth_ptr;
   ML_Sm_Hiptmair_Data *dataptr;
   ML_CommInfoOP *getrow_comm;
   int totalsize;
   struct ML_CSR_MSRdata *temp_csrdata;
   int *Tmat_trans_rowptr, *Tmat_trans_bindx;
   double *Tmat_trans_val;
   int *ATmat_trans_rowptr, *ATmat_trans_bindx;
   double *ATmat_trans_val;
   int ntimes;
#ifdef ML_DEBUG_SMOOTHER
   double *res2, res_norm, dtemp, *tmpvec1, *tmpvec2;
#endif

   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   /* Retrieve smoother specific data. */
   dataptr = (ML_Sm_Hiptmair_Data *) smooth_ptr->smoother->data;
   Tmat_trans  = (ML_Operator *) dataptr->Tmat_trans;
   ATmat_trans = (ML_Operator *) dataptr->ATmat_trans;
   TtAT_diag = (double *) dataptr->TtAT_diag;
#define SPECIAL
#ifdef SPECIAL
   temp_csrdata = (struct ML_CSR_MSRdata *) Tmat_trans->data;
   Tmat_trans_rowptr = temp_csrdata->rowptr;
   Tmat_trans_bindx = temp_csrdata->columns;
   Tmat_trans_val = temp_csrdata->values;
   temp_csrdata = (struct ML_CSR_MSRdata *) ATmat_trans->data;
   ATmat_trans_rowptr = temp_csrdata->rowptr;
   ATmat_trans_bindx = temp_csrdata->columns;
   ATmat_trans_val = temp_csrdata->values;
#endif

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Hiptmair): Need getrow() for Hiptmair smoother\n");

   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {

#ifdef ML_DEBUG_SMOOTHER
   printf("\n\n\t%d: (before SGS on edges) norm of global x = %15.10e\n",
          ATmat_trans->comm->ML_mypid, sqrt(ML_gdot(Nrows, x, x, comm)));
   printf("\t%d: (before SGS on edges) norm of rhs = %15.10e\n",
          ATmat_trans->comm->ML_mypid, sqrt((ML_gdot(Nrows, rhs, rhs, comm))));
   fflush(stdout);
#endif

   /************************
   * Symmetric GS on edges *
   ************************/

   /* Only do one smooth over the edges.  Without this, the following call
      will do N sweeps, where N is the number of Hiptmair sweeps. */
   ntimes = smooth_ptr->ntimes;
   smooth_ptr->ntimes = 1;

   /* fun = ML_Smoother_MSR_SGS; */
#ifdef SPECIAL
   ML_Smoother_MSR_SGSnodamping(sm, inlen, x, outlen, rhs);
#else
   ML_Smoother_SGS(sm, inlen, x, outlen, rhs); 
   /*ML_Smoother_Jacobi(sm, inlen, x, outlen, rhs);*/
#endif
   smooth_ptr->ntimes = ntimes;

#ifdef ML_DEBUG_SMOOTHER
   printf("\t%d:After SGS on edges, norm of x = %15.10e\n",
          ATmat_trans->comm->ML_mypid,
          sqrt(ML_gdot(Nrows, x, x, comm)));
   res2 = (double *) ML_allocate( Nrows * sizeof(double) );
   ML_Operator_Apply(Amat, Nrows, x, Nrows, res2);
   for (kk = 0; kk < outlen; kk++) res2[kk] = rhs[kk] - res2[kk];
   printf("\t%d:After SGS on edges, norm of res = %15.10e\n",
          ATmat_trans->comm->ML_mypid,
          sqrt(ML_gdot(outlen, res2, res2, comm)));
   fflush(stdout);
   ML_free(res2);
#endif

   getrow_comm= Tmat_trans->getrow->pre_comm;
   if (getrow_comm != NULL) 
      totalsize = inlen+getrow_comm->total_rcv_length;
   else
      totalsize = inlen;

   /* Local portion (local+ghost columns) of residual vector. */
   res = (double *) ML_allocate( (totalsize)*sizeof(double));
   for (kk = 0; kk < totalsize; kk++) res[kk] = 0.;

   /* Local portion (local+ghost columns) of solution vector. */
   local_x = (double *) ML_allocate( totalsize * sizeof(double));

   /* Copy global x to local_x. */
   for (kk = 0; kk < inlen; kk++) local_x[kk] = x[kk];
   for (kk = inlen; kk < totalsize; kk++) local_x[kk] = 0.;

   /* Calculate initial global residual. */
   ML_Operator_Apply(Amat, Nrows, x, Nrows , res);
   for (kk = 0; kk < outlen; kk++)
      res[kk] = rhs[kk] - res[kk];

   if (getrow_comm != NULL)
      ML_exchange_bdry(res,getrow_comm, inlen, comm, ML_OVERWRITE,NULL);

#ifdef ML_DEBUG_SMOOTHER
#ifdef vecdump
   if (ATmat_trans->comm->ML_mypid == 1)
   for (kk = 0; kk < ATmat_trans->outvec_leng; kk++)
   {
      printf("%d: TtAT_diag[%d] = %e\n",
	         ATmat_trans->comm->ML_mypid,
			 kk,
			 TtAT_diag[kk]);
      fflush(stdout);
   }
#endif
#undef vecdump

#ifdef vecdump
   if (ATmat_trans->comm->ML_mypid == 1)
   for (kk=0; kk < totalsize; kk++)
   {
      printf("%d: local_x[%d] = %d\n",
	         ATmat_trans->comm->ML_mypid,
			 kk,
			 (int) local_x[kk]);
      fflush(stdout);
   }
#endif
#endif

#ifdef ML_DEBUG_SMOOTHER
   tmpvec1 = (double *) ML_allocate( Tmat_trans->outvec_leng * sizeof(double) );
   ML_Operator_Apply(Tmat_trans, Tmat_trans->invec_leng, res,
                     Tmat_trans->outvec_leng,tmpvec1);
   printf("Before SGS on nodes\n");
   printf("\t||rhs|| = %15.10e\n",
          sqrt(ML_gdot(Tmat_trans->outvec_leng, tmpvec1,tmpvec1,comm))); 
   ML_free(tmpvec1);
#endif


   /*******************************
   * Forward GS sweep over nodes. *
   *******************************/
   for (i = 0; i < ATmat_trans->outvec_leng; i++)
   {
	  /* Grab row i of Tmat_trans. */
#ifdef SPECIAL
      length = Tmat_trans_rowptr[i+1]- Tmat_trans_rowptr[i];
      cols = &(Tmat_trans_bindx[Tmat_trans_rowptr[i]]);
      vals = &(Tmat_trans_val[Tmat_trans_rowptr[i]]);
#else
      ML_get_matrix_row(Tmat_trans, 1, &i , &allocated_space , &cols, &vals,
                        &length, 0);
#endif

	  /* Explicit dot product Tmat_trans(i,:)*res  */
	  numer = 0.;
	  for (kk = 0; kk < length; kk++) numer += vals[kk] * res[cols[kk]];
	  denom = TtAT_diag[i];
	  alpha = numer / denom;

#ifdef printcoeff
	  if (Amat->outvec_leng == 2265)
	  {
	     printf("%d (%d):  numer = %e   denom = %e    alpha = %e\n",i,
	            ATmat_trans->comm->ML_mypid,numer,denom,alpha);
	     for (kk = 0; kk < length; kk++)
	        printf("\tvals = %12.9e   res[%d] = %12.9e\n",vals[kk],cols[kk],
			       res[cols[kk]]);
      }
#endif

#ifdef ML_DEBUG_SMOOTHER
#define PIDNUM 0
#ifdef dump
      if (ATmat_trans->comm->ML_mypid == PIDNUM)
	  printf("%d (%d):  numer = %e   denom = %e    alpha = %e\n",i,
	         ATmat_trans->comm->ML_mypid,numer,denom,alpha);
      fflush(stdout);
#endif
#endif

	  /* Update edge-based local solution. */

      for (j=0; j < length; j++)
	  {
#ifdef ML_DEBUG_SMOOTHER
#ifdef dump
         if ((ATmat_trans->comm->ML_mypid == 1 && cols[j] == 45)
		     || (ATmat_trans->comm->ML_mypid == 0 && cols[j] == 270 ))
         /*
         if (cols[j] == 45)
		 */
		     printf("%d (%d): (before update) x[%d] = %e\n",
	                i,ATmat_trans->comm->ML_mypid,cols[j],local_x[cols[j]]);
#endif
#endif
	     local_x[cols[j]] += alpha * vals[j];
#ifdef ML_DEBUG_SMOOTHER
#ifdef dump
	         printf("%d (%d):  numer = %e   denom = %e    alpha = %e\n",i,
	                ATmat_trans->comm->ML_mypid,numer,denom,alpha);
/*
		     printf("%d (%d): (after update) x[%d] = %e, vals[%d] = %e\n",
	                i,ATmat_trans->comm->ML_mypid,cols[j],local_x[cols[j]],
					cols[j], vals[j]);
*/
             fflush(stdout);
#endif
#undef dump
#endif
      }
	     
      /* Get a row of transpose of matrix product AT (column of AT). */
#ifdef SPECIAL
length = ATmat_trans_rowptr[i+1]- ATmat_trans_rowptr[i];
cols = &(ATmat_trans_bindx[ATmat_trans_rowptr[i]]);
vals = &(ATmat_trans_val[ATmat_trans_rowptr[i]]);
#else
      ML_get_matrix_row(ATmat_trans, 1, &i , &allocated_space , &cols, &vals,
                        &length, 0);
#endif

	  /* Update edge-based local residual. */
      for (j=0; j < length; j++) res[cols[j]] -= alpha * vals[j];
   }

#ifdef ML_DEBUG_SMOOTHER
#ifdef vecdump
   if (ATmat_trans->comm->ML_mypid == PIDNUM)
   for (kk=0; kk < totalsize; kk++)
   {
      printf("%d before exchange: x[%d] = %e\n",
	         ATmat_trans->comm->ML_mypid,
			 kk,
			 local_x[kk]);
      fflush(stdout);
   }
#endif
#endif

#ifdef ML_DEBUG_SMOOTHER
#ifdef vecdump
   if (ATmat_trans->comm->ML_mypid == PIDNUM)
   for (kk=0; kk < totalsize; kk++)
   {
      printf("%d after exchange: x[%d] = %e\n",
	         ATmat_trans->comm->ML_mypid,
			 kk,
			 local_x[kk]);
      fflush(stdout);
   }
#endif
#endif

#ifdef ML_DEBUG_SMOOTHER
   printf("\t||r|| after forward sweep on nodes %15.10e\n",
	  sqrt(ML_gdot(Nrows,res,res,comm)));
   printf("\t||x_local|| after forward sweep on nodes %15.10e\n",
	  sqrt(ML_gdot(totalsize,local_x,local_x,comm)));
   printf("\t||x_glob|| after forward sweep on nodes %15.10e\n",
	  sqrt(ML_gdot(totalsize,x,x,comm)));
   fflush(stdout);
#endif

   /***********************************
   *  Backwards GS sweep over nodes.  *
   ***********************************/
   for (i = Tmat_trans->outvec_leng-1; i >= 0; i--)
   {

	  /* Grab row i of Tmat_trans. */
#ifdef SPECIAL
      length = Tmat_trans_rowptr[i+1]- Tmat_trans_rowptr[i];
      cols = &(Tmat_trans_bindx[Tmat_trans_rowptr[i]]);
      vals = &(Tmat_trans_val[Tmat_trans_rowptr[i]]);
#else
      ML_get_matrix_row(Tmat_trans, 1, &i , &allocated_space , &cols, &vals,
                        &length, 0);
#endif

	  numer = 0.;
	  for (kk = 0; kk < length; kk++)
	  {
	     numer += vals[kk] * res[cols[kk]];
      }
	  denom = TtAT_diag[i];
	  alpha = numer / denom;

	  /*
	  printf("      numer = %e   denom = %e    alpha = %e\n",numer,denom,alpha);
	  */

	  /* Update edge-based solution. */
      for (j=0; j < length; j++)
	     local_x[cols[j]] += alpha * vals[j];

      /* Get a row of transpose of matrix product AT. */
#ifdef SPECIAL
      length = ATmat_trans_rowptr[i+1]- ATmat_trans_rowptr[i];
      cols = &(ATmat_trans_bindx[ATmat_trans_rowptr[i]]);
      vals = &(ATmat_trans_val[ATmat_trans_rowptr[i]]);
#else
      ML_get_matrix_row(ATmat_trans, 1, &i , &allocated_space , &cols, &vals,
                        &length, 0);
#endif

	  /* Update edge-based residual. */
       for (kk=0; kk < length; kk++) res[cols[kk]] -= alpha * vals[kk]; 
   }
   /***************************************
   * end of backwards GS sweep over nodes *
   ***************************************/

   /* Exchange local_x information here. */
   ML_transposed_exchange_bdry(local_x,getrow_comm, inlen, comm, ML_ADD);

   /* copy global local_x to x here */
   for (kk = 0; kk < inlen; kk++) x[kk] = local_x[kk];

#ifdef ML_DEBUG_SMOOTHER

   dtemp = sqrt(ML_gdot(Nrows, res, res, comm));
   printf("\tnorm of residual  after backward sweep on nodes %15.10e\n",dtemp);
   printf("\tnorm of local solution after backward sweep on nodes %15.10e\n",
            sqrt(ML_gdot(totalsize, local_x,local_x, comm)));
   printf("\tnorm of global solution after forward sweep on nodes %15.10e\n\n",
	  sqrt(ML_gdot(totalsize,x,x,comm)));
   fflush(stdout);
#endif
   ML_free(res); ML_free(local_x);
   res = NULL; local_x = NULL;

   } /* for iter = 1... */

#ifndef SPECIAL
   ML_free(cols);
   ML_free(vals);
#endif

   return 0;
#ifdef ML_DEBUG_SMOOTHER
#undef ML_DEBUG_SMOOTHER
#endif
}
#endif /*ifdef MatrixFreeHiptmair*/

/* ************************************************************************* */
/* sequential Symmetric Gauss-Seidel smoother                                */
/* ------------------------------------------------------------------------- */

int ML_Smoother_SGSSequential(void *sm,int inlen,double x[],int outlen,
                              double rhs[])
{
   int           i, j, iter, Nrows, length, allocated_space, *cols, col;
   int           token, mypid, nprocs;
   double        dtemp, diag_term, *vals, omega, *x2;
   ML_Operator   *Amat;
   ML_Comm       *comm;
   ML_CommInfoOP *getrow_comm;
   ML_Smoother   *smooth_ptr;

   smooth_ptr = (ML_Smoother *) sm;
   omega      = smooth_ptr->omega;
   Amat       = smooth_ptr->my_level->Amat;
   comm       = smooth_ptr->my_level->comm;
   nprocs     = comm->ML_nprocs;
   mypid      = comm->ML_mypid;
   Nrows      = Amat->getrow->Nrows;

   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_SGSSequential): Need getrow() for SGS smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_SymGaussSeidel: Not enough space\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL)
   {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL)
      {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {
      token = 0;

      while ( token < nprocs )
      {
         if (getrow_comm != NULL)
            ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

         if ( token == mypid )
         {
            for (i = 0; i < Nrows; i++)
            {
               dtemp = 0.0;
               diag_term = 0.0;
               ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                 &length, 0);
               for (j = 0; j < length; j++)
               {
                  col = cols[j];
                  dtemp += vals[j]*x2[col];
                  if (col == i) diag_term = vals[j];
               }
               if (diag_term != 0.0)
                  x2[i] += omega*(rhs[i] - dtemp)/diag_term;
            }
         }
         token++;
         token = ML_gmax_int( token, comm);
      }
   }
   for (iter = smooth_ptr->ntimes; iter > 0; iter--)
   {
      token = nprocs - 1;

      while ( token >= 0 )
      {
         if (getrow_comm != NULL)
            ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

         if ( token == mypid )
         {
            for (i = Nrows-1; i >= 0; i--)
            {
               dtemp = 0.0;
               diag_term = 0.0;
               ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                 &length, 0);
               for (j = 0; j < length; j++)
               {
                  col = cols[j];
                  dtemp += vals[j]*x2[col];
                  if (col == i) diag_term = vals[j];
               }
               x2[i] += omega*(rhs[i] - dtemp)/diag_term;
            }
         }
         token--;
         token = ML_gmax_int( token, comm);
      }
   }
   if (getrow_comm != NULL)
   {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+2)
      Amat->max_nz_per_row = allocated_space;

   ML_free(vals); ML_free(cols);

   return 0;
}

/* ************************************************************************* */
/* Block Gauss-Seidel smoother                                               */
/* ------------------------------------------------------------------------- */

int ML_Smoother_BlockGS(void *sm,int inlen,double x[],int outlen,
                        double rhs[])
{
   int            iter, i, j, k, length, allocated_space, *cols, col, one;
   double         *vals, omega;
   ML_Operator    *Amat;
   ML_Comm        *comm;
   ML_CommInfoOP  *getrow_comm;
   int Nrows,     **perms, blocksize, Nblocks, row, info;
   double *x2,    **blockdata, *Atimesx, *correc;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp=0;
	 
   smooth_ptr = (ML_Smoother *) sm;

   omega = smooth_ptr->omega;
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;
   dataptr=(ML_Sm_BGS_Data *)smooth_ptr->smoother->data;
   perms = dataptr->perms;
   blockdata = dataptr->blockfacts;
   one=1; /* fortran needs pointers to everything.  I want to be able to
	     pass it a 1, indicating that I want solve matrix-vector systems.
	     I also need to pass fortran an "N": */
   strcpy(N,"N");

   blocksize=dataptr->blocksize;
   Nblocks=Nrows/blocksize;

   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_blockGaussSeidel): Need getrow() for smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   /* this is very space inefficient, but since lapack doesn't do sparsity,
      I'm not sure how to do better without rewriting their solver routine */
   Atimesx = (double *) ML_allocate(blocksize*sizeof(double));
   correc = (double *) ML_allocate(blocksize*sizeof(double));
   if (correc == NULL) pr_error("Error in ML_BlockGaussSeidel:Not enough space\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for BGS smoother\n");

   /* right now, I'm not sure how this will work - need to look up the
      matvec stuff.  It needs to be able to mutliply a set of rows of A
      by x.  Since matvec presumably won't do that, I'm not sure what to
      do about that */

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
				   *sizeof(double));
      if (x2 == NULL) 
      {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

      for (i = 0; i < Nblocks; i++) 
      {
	 for (k = 0; k < blocksize; k++) Atimesx[k]=0.0;
	 for (k = 0; k < blocksize; k++) 
         {
	    row=i*blocksize+k;
	    ML_get_matrix_row(Amat, 1, &row , &allocated_space , &cols, &vals,
	                                       &length, 0);
	    for (j = 0; j < length; j++) 
            {
               col = cols[j];
               Atimesx[k] += vals[j]*x2[col];
	    }
	    correc[k]=rhs[row]-Atimesx[k];
	 }
				
	 MLFORTRAN(dgetrs)(N, &blocksize, &one, blockdata[i], &blocksize, perms[i],
			   correc, &blocksize, &info, itmp);
	 for (k = 0; k < blocksize; k++)
	    x2[k+i*blocksize] += omega*correc[k];
      }
/* symmetrize it
      for (i = Nblocks-1; i >= 0; i--) {
         for (k = 0; k < blocksize; k++) Atimesx[k]=0.0;
         for (k = 0; k < blocksize; k++) {
            row=i*blocksize+k;
            ML_get_matrix_row(Amat, 1, &row , &allocated_space , &cols, &vals,
                                               &length, 0);
            for (j = 0; j < length; j++) {
               col = cols[j];
               Atimesx[k] += vals[j]*x2[col];
            }
            correc[k]=rhs[row]-Atimesx[k];
         }

         MLFORTRAN(dgetrs)(N, &blocksize, &one, blockdata[i], &blocksize, perms[i],
                           correc, &blocksize, &info, itmp);
         for (k = 0; k < blocksize; k++)
            x2[k+i*blocksize] += omega*correc[k];
      }
*/
   }
	 
   if (getrow_comm != NULL) 
   {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+2) {
      Amat->max_nz_per_row = allocated_space;
   }
	 
   ML_free(vals); ML_free(cols); ML_free(Atimesx); ML_free(correc);
	 
   return 0;
}

int ML_BlockDinv(ML_Sm_BGS_Data *BGS_Data, int inlen, double out[]) {
  unsigned int   itmp=0;
  int info, one = 1;
  int nblocks, **perms, blocksize, i, *blocklengths, index;
  double **blockdata, *dtemp;
  char N[2];
  int k,length;
  int maxBlocksize, *aggr_offset, *aggr_group, Nrows;
  int *block_indices;

  nblocks    = BGS_Data->Nblocks;
  blockdata  = BGS_Data->blockfacts;
  perms      = BGS_Data->perms;
  blocklengths =  BGS_Data->blocklengths;
  blocksize = BGS_Data->blocksize;
  block_indices = BGS_Data->blockmap;

  strcpy(N,"N");

  index = 0;
  if (blocksize > 0) { 
    /* constant blocksize with consecutive indices in blocks */

    for (i = 0; i < nblocks; i++) {
      length = blocksize;
      MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
			perms[i], &(out[index]), 
			&blocksize, &info, itmp);
      index += blocksize;
      if ( info != 0 ) {
	printf("dgetrs returns with %d at block %d\n",info,i); 
	exit(1);
      }
    }
  }
  else {
   maxBlocksize = 0;
   for ( i = 0; i < nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) ML_allocate( (nblocks+1) * sizeof(int) );
   aggr_offset[0] = 0;
   for (i = 1; i <= nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 
   Nrows = aggr_offset[nblocks];

   aggr_group  = (int *) ML_allocate( Nrows   * sizeof(int) );
   for (i = 0; i < Nrows; i++) 
      aggr_group[aggr_offset[block_indices[i]]++] = i; 
   aggr_offset[0] = 0;
   for (i = 1; i < nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 

   if ( maxBlocksize > 0 )    {
      dtemp   = (double *) ML_allocate( maxBlocksize * sizeof(double) );
   }


    for (i = 0; i < nblocks; i++) {
      blocksize = blocklengths[i];
      length = blocksize;

      for (k = 0; k < blocksize; k++) {
	dtemp[k] = out[aggr_group[aggr_offset[i]+k]];
      }
      MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
			perms[i], dtemp,
			&blocksize, &info, itmp);
      if ( info != 0 ) {
	printf("dgetrs returns with %d at block %d\n",info,i); 
	exit(1);
      }

      for (k = 0; k < blocksize; k++) {
	out[aggr_group[aggr_offset[i]+k]] = dtemp[k];
      }

    }
   ML_free(aggr_group);
   ML_free(aggr_offset);

  }
  return 0;
}
#include "ml_mls.h"

int ML_BlockScaledApply(ML_Operator *Amat, int inlen, double in[],
			int outlen, double out[]) {

  struct MLSthing *widget;


  widget     = (struct MLSthing *) Amat->data;

  ML_Operator_Apply(widget->unscaled_matrix, inlen, in, outlen, out);
  ML_BlockDinv(widget->block_scaling, outlen, out);

  return 0;
}


/* ************************************************************************* */
/* Variable size Block Jacobi smoother                                       */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockJacobi(void *sm, int inlen, double x[], int outlen, 
                             double rhs[])
{
   int            i, j, k, iter, blocksize, one=1, *aggr_offset, *block_indices;
   int            Nrows, **perms, Nblocks, info, *blocklengths, row, length;
   int            maxBlocksize, *aggr_group, *cols, allocated_space, col;
   int            *do_update = NULL;
   ML_Comm        *comm;
   ML_CommInfoOP  *getrow_comm;
   ML_Operator    *Amat;
   double         *x_ext, **blockdata, *vals, omega;
   double         *unprec_r = NULL, *Mr = NULL, *dtemp = NULL;
   double         r_z_dot, p_ap_dot;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp=0;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_VBlockJacobi): Need getrow() for smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for VBJacobi smoother\n");
   Nrows         = Amat->getrow->Nrows;
   dataptr       = (ML_Sm_BGS_Data *)smooth_ptr->smoother->data;
   perms         = dataptr->perms;
   blockdata     = dataptr->blockfacts;
   Nblocks       = dataptr->Nblocks;
   blocklengths  = dataptr->blocklengths;
   block_indices = dataptr->blockmap;

   /* ----------------------------------------------------- */
   /* set up for fetching the matrix                        */
   /* ----------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row+1000;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));

   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) ML_allocate( Nblocks * sizeof(int) );
   aggr_group  = (int *) ML_allocate( Nrows   * sizeof(int) );
   aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 
   for (i = 0; i < Nrows; i++) 
      aggr_group[aggr_offset[block_indices[i]]++] = i; 
   aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 

   /* ----------------------------------------------------- */
   /* create a long buffer for incoming data                */
   /* ----------------------------------------------------- */

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x_ext = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      dtemp   = (double *) ML_allocate( maxBlocksize * sizeof(double) );
   }
   if ( Nblocks > 0 )
   {
     unprec_r = (double *) ML_allocate(inlen*sizeof(double));

      do_update = (int *) ML_allocate( Nblocks * sizeof(int) );
      if ( do_update == NULL ) 
      {
         printf("ERROR : memory allocation.\n");
         exit(1);
      }
   }
   if (smooth_ptr->omega == ML_ONE_STEP_CG) {
     Mr = (double *) ML_allocate(inlen*sizeof(double));
     if (Mr == NULL) pr_error("ML_Smoother_VBlockJacobi: Out of space\n");
   }
   else Mr = unprec_r;

   /* ----------------------------------------------------- */
   /* iterate                                               */
   /* ----------------------------------------------------- */

   strcpy(N,"N");
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
     /* do we need this communication (as it is in the matvec ?) */
      if (getrow_comm != NULL)
         ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

      /* compute the residual */

      ML_Operator_Apply(Amat, inlen, x_ext, inlen, unprec_r);
      for (i = 0; i < inlen; i++) unprec_r[i] = rhs[i] - unprec_r[i];

      /* compute do_update */

      for (i = 0; i < Nblocks; i++) 
      {
         do_update[i] = 0;
         blocksize = blocklengths[i];

         for (k = 0; k < blocksize; k++) 
         {
            row = aggr_group[aggr_offset[i]+k];
            ML_get_matrix_row(Amat,1,&row,&allocated_space,&cols,&vals,
                                                &length,0);
            for (j = 0; j < length; j++) 
            {
               col = cols[j];
               if ( col == row ) do_update[i]++;
            }
         }
      }
 
      for (i = 0; i < Nblocks; i++) 
      {
         blocksize = blocklengths[i];
         if ( do_update[i] == blocksize && blocksize != 0 )
         {
	   for (k = 0; k < blocksize; k++) {
	     dtemp[k] = unprec_r[aggr_group[aggr_offset[i]+k]];
           }
            MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                              perms[i], dtemp, &blocksize, &info, itmp);
            if ( info != 0 ) 
            {
               printf("dgetrs returns with %d at block %d\n",info,i); 
               exit(1);
            }
            for (k = 0; k < blocksize; k++) {
	      /*
               x_ext[aggr_group[aggr_offset[i]+k]] += (omega * dtemp[k]);
	      */
               Mr[aggr_group[aggr_offset[i]+k]] = dtemp[k];
            }
         }
      }
      if (smooth_ptr->omega == ML_ONE_STEP_CG) {
         /* Compute damping parameter that corresonds to one step of CG. */
         r_z_dot = 0.;
         r_z_dot = ML_gdot(inlen, Mr, unprec_r, smooth_ptr->my_level->comm);
         ML_Operator_Apply(Amat, inlen, Mr, inlen, unprec_r);
         p_ap_dot = ML_gdot(inlen, Mr, unprec_r, smooth_ptr->my_level->comm);
         if (p_ap_dot != 0.0) omega = r_z_dot/p_ap_dot;
         else omega = 1.;
      }


      /* watch out for 'do_update' */
      for (i = 0; i < inlen; i++)  x_ext[i] += (omega * Mr[i]);
   }
	 
   /* ----------------------------------------------------- */
   /* copy data to output buffer                            */
   /* ----------------------------------------------------- */

   if (getrow_comm != NULL) 
   {
      for (i = 0; i < inlen; i++) x[i] = x_ext[i];
      ML_free(x_ext);
   }

   ML_free( vals ); 
   ML_free( cols ); 
   ML_free( aggr_offset );
   ML_free( aggr_group );
   if ( unprec_r != NULL) ML_free(unprec_r);
   if ((smooth_ptr->omega == ML_ONE_STEP_CG) && (Mr != NULL)) ML_free(Mr);
   if ( dtemp    != NULL) ML_free(dtemp);
   if ( Nblocks > 0 ) ML_free( do_update );
	 
   return 0;
}

/* ************************************************************************* */
/*  Variable size Block symmetric Gauss Seidel smoother                      */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockSGS(void *sm, int inlen, double x[], 
                          int outlen, double rhs[])
{
   int            i, j, k, iter, blocksize, one=1, *aggr_offset, *block_indices;
   int            Nrows, **perms, Nblocks, info, *blocklengths, row, length;
   int            maxBlocksize, *aggr_group, *cols, allocated_space, col;
   int            do_update;
   ML_Comm        *comm;
   ML_CommInfoOP  *getrow_comm;
   ML_Operator    *Amat;
   double         *x_ext, **blockdata, *Ax = NULL, *res=NULL, *vals, omega;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp=0;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_VBlockSymGS): Need getrow() for smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for VBSymGS smoother\n");
   Nrows         = Amat->getrow->Nrows;
   dataptr       = (ML_Sm_BGS_Data *)smooth_ptr->smoother->data;
   perms         = dataptr->perms;
   blockdata     = dataptr->blockfacts;
   Nblocks       = dataptr->Nblocks;
   blocklengths  = dataptr->blocklengths;
   block_indices = dataptr->blockmap;

   /* ----------------------------------------------------- */
   /* set up for fetching the matrix                        */
   /* ----------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row+1000;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) ML_allocate( Nblocks * sizeof(int) );
   aggr_group  = (int *) ML_allocate( Nrows   * sizeof(int) );
   if ( Nblocks > 0 ) aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 
   for (i = 0; i < Nrows; i++) 
      aggr_group[aggr_offset[block_indices[i]]++] = i; 
   if ( Nblocks > 0 ) aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 

   /* ----------------------------------------------------- */
   /* create a long buffer for incoming data                */
   /* ----------------------------------------------------- */

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x_ext = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      Ax  = (double *) ML_allocate( maxBlocksize * sizeof(double) );
      res = (double *) ML_allocate( maxBlocksize * sizeof(double) );
   }

   /* ----------------------------------------------------- */
   /* iterate                                               */
   /* ----------------------------------------------------- */

   strcpy(N,"N");
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

      for (i = 0; i < Nblocks; i++) 
      {
         do_update = 0;
         blocksize = blocklengths[i];
         for (k = 0; k < blocksize; k++) 
         {
            Ax[k] = 0.0;
            row = aggr_group[aggr_offset[i]+k];
            ML_get_matrix_row(Amat,1,&row,&allocated_space,&cols,&vals,
                                                &length,0);
            for (j = 0; j < length; j++) 
            {
               col = cols[j];
               if ( col == row ) do_update++;
               Ax[k] += vals[j]*x_ext[col];
            }
            res[k] = rhs[row] - Ax[k];
         }
				
         if ( do_update == blocksize && blocksize != 0 )
         {
            MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                              perms[i], res, &blocksize, &info, itmp);
            if ( info != 0 ) 
            {
               printf("dgetrs returns with %d at block %d(%d)\n",info,i,Nblocks); 
               exit(1);
            }
            for (k = 0; k < blocksize; k++)
            {
               x_ext[aggr_group[aggr_offset[i]+k]] += (omega * res[k]);
            }
         }
      }
      for (i = 0; i < Nrows; i++)
      {
         if ( block_indices[i] == -1 )
         {
            ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                                &length,0);
            Ax[0] = 0.0;
            for (j = 0; j < length; j++) Ax[0] += (vals[j]*x_ext[cols[j]]);
            x_ext[i] += (omega * (rhs[i] - Ax[0]));
         }
      }
      for (i = Nrows-1; i >= 0; i--)
      {
         if ( block_indices[i] == -1 )
         {
            ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                                &length,0);
            Ax[0] = 0.0;
            for (j = 0; j < length; j++) Ax[0] += (vals[j]*x_ext[cols[j]]);
            x_ext[i] += (omega * (rhs[i] - Ax[0]));
         }
      }
      for (i = Nblocks-1; i >= 0; i--) 
      {
         blocksize = blocklengths[i];
         do_update = 0;
         for (k = 0; k < blocksize; k++) 
         {
            Ax[k] = 0.0;
            row = aggr_group[aggr_offset[i]+k];
            ML_get_matrix_row(Amat,1,&row,&allocated_space,&cols,&vals,
                                                &length,0);
            for (j = 0; j < length; j++) 
            {
               col = cols[j];
               if ( col == row ) do_update++;
               Ax[k] += vals[j]*x_ext[col];
            }
            res[k] = rhs[row] - Ax[k];
         }
				
         if ( do_update == blocksize )
         {
            MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                              perms[i], res, &blocksize, &info, itmp);
            if ( info != 0 ) 
            {
               printf("dgetrs returns with %d at block %d(%d)\n",info,i,blocksize); 
               exit(1);
            }
            for (k = 0; k < blocksize; k++)
            {
               x_ext[aggr_group[aggr_offset[i]+k]] += (omega * res[k]);
            }
         }
      }
   }
	 
   /* ----------------------------------------------------- */
   /* copy data to output buffer                            */
   /* ----------------------------------------------------- */

   if (getrow_comm != NULL) 
   {
      for (i = 0; i < inlen; i++) x[i] = x_ext[i];
      ML_free(x_ext);
   }
   if ( maxBlocksize > 0 )
   {
      ML_free( Ax );
      ML_free( res );
   }
   ML_free( vals ); 
   ML_free( cols ); 
   if ( Nblocks > 0 ) ML_free( aggr_offset );
   if ( Nrows > 0 ) ML_free( aggr_group );
	 
   return 0;
}

/* ************************************************************************* */
/* Variable size Block Gauss Seidel smoother (sequential)                    */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockSGSSequential(void *sm, int inlen, double x[], 
                                    int outlen, double rhs[])
{
   int            i, j, k, iter, blocksize, one=1, *aggr_offset, *block_indices;
   int            Nrows, **perms, Nblocks, info, *blocklengths, row, length;
   int            maxBlocksize, *aggr_group, *cols, allocated_space, col;
   int            do_update, nprocs, mypid, token;
   ML_Comm        *comm;
   ML_CommInfoOP  *getrow_comm;
   ML_Operator    *Amat;
   double         *x_ext, **blockdata, *Ax = NULL, *res = NULL, *vals, omega;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp=0;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) {
      pr_error("Error(ML_VBlockSymGSSeq): Need getrow() for smoother\n");
      ML_avoid_unused_param((void *) &outlen);
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for VBSGSSeq smoother\n");
   Nrows         = Amat->getrow->Nrows;
   dataptr       = (ML_Sm_BGS_Data *)smooth_ptr->smoother->data;
   perms         = dataptr->perms;
   blockdata     = dataptr->blockfacts;
   Nblocks       = dataptr->Nblocks;
   blocklengths  = dataptr->blocklengths;
   block_indices = dataptr->blockmap;
   nprocs        = comm->ML_nprocs;
   mypid         = comm->ML_mypid;

   /* ----------------------------------------------------- */
   /* set up for fetching the matrix                        */
   /* ----------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row+1000;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) ML_allocate( Nblocks * sizeof(int) );
   aggr_group  = (int *) ML_allocate( Nrows   * sizeof(int) );
   if ( Nblocks > 0 ) aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 
   for (i = 0; i < Nrows; i++) 
      aggr_group[aggr_offset[block_indices[i]]++] = i; 
   if ( Nblocks > 0 ) aggr_offset[0] = 0;
   for (i = 1; i < Nblocks; i++) 
      aggr_offset[i] = aggr_offset[i-1] + blocklengths[i-1]; 

   /* ----------------------------------------------------- */
   /* create a long buffer for incoming data                */
   /* ----------------------------------------------------- */

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      x_ext = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      Ax  = (double *) ML_allocate( maxBlocksize * sizeof(double) );
      res = (double *) ML_allocate( maxBlocksize * sizeof(double) );
   }

   /* ----------------------------------------------------- */
   /* iterate                                               */
   /* ----------------------------------------------------- */

   strcpy(N,"N");
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      token = 0;

      while ( token < nprocs )
      {
         if (getrow_comm != NULL)
            ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

         if ( token == mypid )
         {
            for (i = 0; i < Nblocks; i++) 
            {
               do_update = 0;
               blocksize = blocklengths[i];
               for (k = 0; k < blocksize; k++) 
               {
                  Ax[k] = 0.0;
                  row = aggr_group[aggr_offset[i]+k];
                  ML_get_matrix_row(Amat,1,&row,&allocated_space,&cols,&vals,
                                                   &length,0);
                  for (j = 0; j < length; j++) {
                     col = cols[j];
                     if ( col == row ) do_update++;
                     Ax[k] += vals[j]*x_ext[col];
                  }
                  res[k] = rhs[row] - Ax[k];
               }
               if ( do_update == blocksize && blocksize != 0 )
               {
                  MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                                    perms[i], res, &blocksize, &info, itmp);
                  if ( info != 0 ) 
                  {
                     printf("dgetrs returns %d at blk %d(%d)\n",info,i,Nblocks); 
                     exit(1);
                  }
                  for (k = 0; k < blocksize; k++)
                  {
                     x_ext[aggr_group[aggr_offset[i]+k]] += (omega * res[k]);
                  }
               }
            }
            for (i = 0; i < Nrows; i++)
            {
               if ( block_indices[i] == -1 )
               {
                  ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                                      &length,0);
                  Ax[0] = 0.0;
                  for (j = 0; j < length; j++) Ax[0] += (vals[j]*x_ext[cols[j]]);
                  x_ext[i] += (omega * (rhs[i] - Ax[0]));
               }
            }
         }
         token++; 
         token = ML_gmax_int( token, comm);
      }
   }
   for (iter = smooth_ptr->ntimes; iter > 0; iter--) 
   {
      token = nprocs - 1;

      while ( token >= 0 )
      {
         if (getrow_comm != NULL)
            ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

         if ( token == mypid )
         {
            for (i = Nrows-1; i >= 0; i++)
            {
               if ( block_indices[i] == -1 )
               {
                  ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,
                                                      &length,0);
                  Ax[0] = 0.0;
                  for (j = 0; j < length; j++) Ax[0] += (vals[j]*x_ext[cols[j]]);
                  x_ext[i] += (omega * (rhs[i] - Ax[0]));
               }
            }
            for (i = Nblocks-1; i >= 0; i--) 
            {
               blocksize = blocklengths[i];
               do_update = 0;
               for (k = 0; k < blocksize; k++) 
               {
                  Ax[k] = 0.0;
                  row = aggr_group[aggr_offset[i]+k];
                  ML_get_matrix_row(Amat,1,&row,&allocated_space,&cols,&vals,
                                                      &length,0);
                  for (j = 0; j < length; j++) {
                     col = cols[j];
                     if ( col == row ) do_update++;
                     Ax[k] += vals[j]*x_ext[col];
                  }
                  res[k] = rhs[row] - Ax[k];
               }
				
               if ( do_update == blocksize )
               {
                  MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                                    perms[i], res, &blocksize, &info, itmp);
                  if ( info != 0 ) 
                  {
                     printf("dgetrs returns %d at blk %d(%d)\n",info,i,Nblocks); 
                     exit(1);
                  }
                  for (k = 0; k < blocksize; k++)
                  {
                     x_ext[aggr_group[aggr_offset[i]+k]] += (omega * res[k]);
                  }
               }
            }
         }
         token--;
         token = ML_gmax_int( token, comm);
      }
   }
	 
   /* ----------------------------------------------------- */
   /* copy data to output buffer                            */
   /* ----------------------------------------------------- */

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x_ext[i];
      ML_free(x_ext);
   }

   if ( maxBlocksize > 0 )
   {
      ML_free( Ax );
      ML_free( res );
   }
   ML_free( vals ); 
   ML_free( cols ); 
   if ( Nblocks > 0 ) ML_free( aggr_offset );
   if ( Nrows > 0 ) ML_free( aggr_group );
	 
   return 0;
}

/* ************************************************************************* */
/* Variable size Block Jacobi smoother with Krylov                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockKrylovJacobi(void *sm,int inlen,double x[],int outlen,
                                   double rhs[])
{
   ML_Comm        *comm;
   ML_Operator    *Amat;
   ML_Smoother    *smooth_ptr;
   ML_Krylov      *ml_kry;

   if (outlen == - 47) ML_avoid_unused_param((void *) &outlen);
   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;

   ml_kry = ML_Krylov_Create(comm);
   ML_Krylov_Set_Method(ml_kry, 0);
   ML_Krylov_Set_Amatrix(ml_kry, Amat);
   ML_Krylov_Set_MaxIterations(ml_kry, 3);
   ML_Krylov_Set_Precon(ml_kry, sm);
   ML_Krylov_Set_PrintFreq(ml_kry, 1000);
   ML_Krylov_Set_PreconFunc(ml_kry, ML_Smoother_VBlockJacobi);
   ML_Krylov_Solve(ml_kry, inlen, rhs, x);
   ML_Krylov_Destroy(&ml_kry);
   return 0;
}  

/* ************************************************************************* */
/* overlapped Domain decomposition                                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_OverlappedILUT(void *sm,int inlen,double x[],int outlen,
                               double rhs[])
{
   int             i, j, column, *idiag, *mat_ia, *mat_ja, extNrows;
   ML_Comm         *comm;
   ML_Operator     *Amat;
   ML_Smoother     *smooth_ptr;
   ML_Sm_ILUT_Data *dataptr;
   double          *dbuffer, ddata, *mat_aa;
   ML_CommInfoOP   *getrow_comm;

   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;
   dataptr    = (ML_Sm_ILUT_Data *) smooth_ptr->smoother->data;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_OverlappedILUT): Need getrow()\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for ML_OverlappedILUT\n");
   if ( dataptr == NULL ) 
      pr_error("Error(ML_OverlappedILUT): Need dataptr\n");

   extNrows = dataptr->Nrows;
   mat_ia   = dataptr->mat_ia;
   mat_ja   = dataptr->mat_ja;
   mat_aa   = dataptr->mat_aa;

   dbuffer = (double *) ML_allocate(extNrows * sizeof(double));
   idiag   = (int *)    ML_allocate(extNrows * sizeof(int));
   for ( i = 0; i < inlen; i++ ) dbuffer[i] = rhs[i];

   if ( extNrows > outlen )
   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(ML_OverlappedILUT): Need getrow()\n");
      if (Amat->getrow->post_comm != NULL)
         pr_error("Post communication not implemented for ML_OverlappedILUT\n");
      if ( dataptr == NULL ) 
         pr_error("Error(ML_OverlappedILUT): Need dataptr\n");

      getrow_comm= Amat->getrow->pre_comm;
      if (getrow_comm != NULL)
         ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);
   }

   for ( i = 0; i < extNrows; i++ )
   {
      ddata = 0.0;
      for ( j = mat_ia[i]; j < mat_ia[i+1]; j++ )
      {
         column = mat_ja[j];
         if ( column == i ) { idiag[i] = j; break;}
         ddata += mat_aa[j] * dbuffer[column];
      }
      dbuffer[i] -= ddata;
   }
   for ( i = extNrows-1; i >= 0; i-- )
   {
      ddata = 0.0;
      for ( j = idiag[i]+1; j < mat_ia[i+1]; j++ )
      {
         column = mat_ja[j];
         ddata += mat_aa[j] * dbuffer[column];
      }
      dbuffer[i] -= ddata;
      dbuffer[i] /= mat_aa[idiag[i]];
   }
   for ( i = 0; i < inlen; i++ ) x[i] = dbuffer[i];
   ML_free(dbuffer);
   ML_free(idiag);
      
   return 0;
}

#ifdef out
/* ************************************************************************* */
/* variable block additive Schwarz                                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockAdditiveSchwarz(void *sm, int inlen, double x[],
                                      int outlen, double rhs[])
{
#ifdef SUPERLU
   int                i, j, m, k, extNrows, nblocks, length, *indptr, ntimes;
   int                *blk_size, **blk_indices, max_blk_size;
   int                **aux_bmat_ia, **aux_bmat_ja;
   double             *dbuffer, **aux_bmat_aa, *rhsbuf, *solbuf, *xbuffer = NULL;
   ML_Comm            *comm;
   ML_Operator        *Amat;
   ML_Smoother        *smooth_ptr;
   ML_Sm_Schwarz_Data *dataptr;
   ML_CommInfoOP      *getrow_comm;
   int                info, *perm_r, *perm_c, *etree, panel_size, lwork;
   double             *R, *C, *ferr, *berr, rpg, rcond, dtemp;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* --------------------------------------------------------- */
   /* fetch parameters and check                                */
   /* --------------------------------------------------------- */

   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;
   dataptr    = (ML_Sm_Schwarz_Data *) smooth_ptr->smoother->data;
   ntimes     = smooth_ptr->ntimes;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Smoother_AdditiveSchwarz): Need getrow()\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for AdditiveSchwarz\n");
   if ( dataptr == NULL ) 
      pr_error("Error(AdditiveSchwarz): Need dataptr\n");

   getrow_comm = Amat->getrow->pre_comm;
   extNrows    = dataptr->Nrows;
   nblocks     = dataptr->nblocks;
   blk_indices = dataptr->blk_indices;
   blk_size    = dataptr->blk_size;
   aux_bmat_ia = dataptr->aux_bmat_ia;
   aux_bmat_ja = dataptr->aux_bmat_ja;
   aux_bmat_aa = dataptr->aux_bmat_aa;
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
      max_blk_size = (blk_size[i] > max_blk_size) ? blk_size[i] : max_blk_size;

   /* --------------------------------------------------------- */
   /* communicate the rhs and put into dbuffer                  */
   /* --------------------------------------------------------- */

   dbuffer = (double *) ML_allocate(extNrows * sizeof(double));
   for ( i = 0; i < outlen; i++ ) dbuffer[i] = rhs[i];
   for ( i = 0; i < inlen;  i++ ) x[i] = 0.0;

   if (extNrows > outlen && getrow_comm != NULL)
      ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

   /* --------------------------------------------------------- */
   /* set up for SuperLU solves                                 */
   /* --------------------------------------------------------- */

   rhsbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   solbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   *fact                    = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );

   /* --------------------------------------------------------- */
   /* the first pass                                            */
   /* --------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ )
   {
      indptr = blk_indices[i];
      length = blk_size[i];
      for ( j = 0; j < length; j++ ) rhsbuf[j] = dbuffer[indptr[j]];
      A = dataptr->slu_Amat[i];
      L = dataptr->slu_Lmat[i];
      U = dataptr->slu_Umat[i];
      perm_c = dataptr->perm_c[i];
      perm_r = dataptr->perm_r[i];
      dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
      dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
      for ( j = 0; j < length; j++ ) 
         /*if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   if (ntimes > 1) xbuffer = (double *) ML_allocate(extNrows * sizeof(double));

   for ( m = 1; m < ntimes; m++ )
   {
      for ( i = 0; i < inlen; i++ ) xbuffer[i] = x[i];
      if (extNrows > outlen && getrow_comm != NULL)
         ML_exchange_bdry(xbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

      for ( i = 0; i < nblocks; i++ )
      {
         indptr = blk_indices[i];
         length = blk_size[i];
         for ( j = 0; j < length; j++ )
         {
            dtemp = dbuffer[indptr[j]];
            for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
               dtemp -= (aux_bmat_aa[i][k] * xbuffer[aux_bmat_ja[i][k]]); 
            rhsbuf[j] = dtemp;
         }
         A = dataptr->slu_Amat[i];
         L = dataptr->slu_Lmat[i];
         U = dataptr->slu_Umat[i];
         perm_c = dataptr->perm_c[i];
         perm_r = dataptr->perm_r[i];
         dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
         dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
         dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
                equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
                ferr, berr, &mem_usage, &info);
         for ( j = 0; j < length; j++ ) 
            /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
         Destroy_SuperMatrix_Store(&B);
         Destroy_SuperMatrix_Store(&X);
      }
   }

   /* --------------------------------------------------------- */
   /* clean up                                                  */
   /* --------------------------------------------------------- */

   if (ntimes > 1) ML_free(xbuffer);
   ML_free( rhsbuf );
   ML_free( solbuf );
   ML_free( dbuffer );
   ML_free( etree );
   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   return 0;
#else
   printf("ML_Smoother_VBlockAdditiveSchwarz : not available.\n");
   exit(1);
   return 1;
#endif
}

/* ************************************************************************* */
/* variable block multiplicative Schwarz                                     */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockMultiplicativeSchwarz(void *sm, int inlen, double x[],
                                      int outlen, double rhs[])
{
#ifdef SUPERLU
   int                i, j, k, m, extNrows, nblocks, length, *indptr, ntimes;
   int                index, *blk_size, **blk_indices, max_blk_size;
   int                **aux_bmat_ia, **aux_bmat_ja;
   double             *dbuffer, **aux_bmat_aa, *rhsbuf, *solbuf, *xbuffer=NULL;
   ML_Comm            *comm;
   ML_Operator        *Amat;
   ML_Smoother        *smooth_ptr;
   ML_Sm_Schwarz_Data *dataptr;
   ML_CommInfoOP      *getrow_comm;
   int                info, *perm_r, *perm_c, *etree, panel_size, lwork;
   double             *R, *C, *ferr, *berr, rpg, rcond, dtemp;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* --------------------------------------------------------- */
   /* fetch parameters and check                                */
   /* --------------------------------------------------------- */

   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;
   dataptr    = (ML_Sm_Schwarz_Data *) smooth_ptr->smoother->data;
   ntimes     = smooth_ptr->ntimes;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Smoother_MultiplicativeSchwarz): Need getrow()\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for MultiplicativeSchwarz\n");
   if ( dataptr == NULL ) 
      pr_error("Error(MultiplicativeSchwarz): Need dataptr\n");
   getrow_comm= Amat->getrow->pre_comm;

   extNrows    = dataptr->Nrows;
   nblocks     = dataptr->nblocks;
   blk_indices = dataptr->blk_indices;
   blk_size    = dataptr->blk_size;
   aux_bmat_ia = dataptr->aux_bmat_ia;
   aux_bmat_ja = dataptr->aux_bmat_ja;
   aux_bmat_aa = dataptr->aux_bmat_aa;
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
      max_blk_size = (blk_size[i] > max_blk_size) ? blk_size[i] : max_blk_size;

   /* --------------------------------------------------------- */
   /* communicate the rhs and put into dbuffer                  */
   /* --------------------------------------------------------- */

   dbuffer = (double *) ML_allocate(extNrows * sizeof(double));
   for ( i = 0; i < outlen; i++ ) dbuffer[i] = rhs[i];
   for ( i = outlen; i < extNrows; i++ ) dbuffer[i] = 0.0;
   for ( i = 0; i < inlen;  i++ ) x[i] = 0.0;

   if (extNrows > outlen && getrow_comm != NULL)
      ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

   /* --------------------------------------------------------- */
   /* set up for SuperLU solves                                 */
   /* --------------------------------------------------------- */

   rhsbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   solbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   *fact                    = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );

   /* --------------------------------------------------------- */
   /* the first pass                                            */
   /* --------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ )
   {
      indptr = blk_indices[i];
      length = blk_size[i];
      for ( j = 0; j < length; j++ ) rhsbuf[j] = dbuffer[indptr[j]];
      A = dataptr->slu_Amat[i];
      L = dataptr->slu_Lmat[i];
      U = dataptr->slu_Umat[i];
      perm_c = dataptr->perm_c[i];
      perm_r = dataptr->perm_r[i];
      dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
      dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
      for ( j = 0; j < length; j++ ) 
         /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j]; */
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   if (ntimes > 1) xbuffer = (double *) ML_allocate(extNrows * sizeof(double));

   for ( m = 1; m < ntimes; m++ )
   {
      for ( i = 0; i < inlen; i++ ) xbuffer[i] = x[i];
      if (extNrows > outlen && getrow_comm != NULL)
         ML_exchange_bdry(xbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

      for ( i = 0; i < nblocks; i++ )
      {
         indptr = blk_indices[i];
         length = blk_size[i];
         for ( j = 0; j < length; j++ )
         {
            dtemp = dbuffer[indptr[j]];
            for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
            {
               index = aux_bmat_ja[i][k];
               if (index < inlen) dtemp -= (aux_bmat_aa[i][k] * x[index]); 
               else               dtemp -= (aux_bmat_aa[i][k] * xbuffer[index]); 
            }
            rhsbuf[j] = dtemp;
         }
         A = dataptr->slu_Amat[i];
         L = dataptr->slu_Lmat[i];
         U = dataptr->slu_Umat[i];
         perm_c = dataptr->perm_c[i];
         perm_r = dataptr->perm_r[i];
         dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
         dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
         dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
                equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
                ferr, berr, &mem_usage, &info);
         for ( j = 0; j < length; j++ ) 
            /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
         Destroy_SuperMatrix_Store(&B);
         Destroy_SuperMatrix_Store(&X);
      }
   }

   /* --------------------------------------------------------- */
   /* clean up                                                  */
   /* --------------------------------------------------------- */

   if (ntimes > 1) ML_free(xbuffer);
   ML_free( rhsbuf );
   ML_free( solbuf );
   ML_free( dbuffer );
   ML_free( etree );
   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   return 0;
#else
   printf("ML_Smoother_VBlockMultiplicativeSchwarz : not available.\n");
   exit(1);
   return 1;
#endif
}

#endif
/* ******************************************************************** */
/* ******************************************************************** */
/* setup routines for various smoothers                                 */
/* ******************************************************************** */
/* ******************************************************************** */
/* Constructor for Sm_Hiptmair_Data                                     */
/* ******************************************************************** */
 
int ML_Smoother_Create_Hiptmair_Data(ML_Sm_Hiptmair_Data **data)
{
   ML_Sm_Hiptmair_Data *ml_data;
 
   ML_memory_alloc((void**) data, sizeof(ML_Sm_Hiptmair_Data), "Hiptmair" );
   ml_data = (*data);
   ml_data->Tmat = NULL;
   ml_data->Tmat_trans = NULL;
   ml_data->ATmat_trans = NULL;
   ml_data->TtAT_diag = NULL;
   ml_data->TtATmat = NULL;
   ml_data->sm_nodal = NULL;
   ml_data->res_edge = NULL;
   ml_data->rhs_nodal = NULL;
   ml_data->x_nodal = NULL;
   ml_data->edge_update = NULL;
   ml_data->max_eig = 0.0;
   ml_data->omega = 1.0;
   ml_data->output_level = 2;
   ml_data->ml_nodal = NULL;
   ml_data->ml_edge = NULL;
   ml_data->reduced_smoother = 0;
   return(0);
}
 
int ML_Smoother_Create_BlockHiptmair_Data(ML_Sm_BlockHiptmair_Data **data)
{
   ML_Sm_BlockHiptmair_Data *ml_data;
 
   ML_memory_alloc((void**) data, sizeof(ML_Sm_BlockHiptmair_Data), "Hiptmair" );

   ml_data = (*data);
   ml_data->Tmat = NULL;
   ml_data->Tmat_trans = NULL;
   ml_data->ATmat_trans = NULL;
   ml_data->TtAT_diag = NULL;
   ml_data->TtATmat = NULL;
   ml_data->sm_nodal = NULL;
   ml_data->res_edge = NULL;
   ml_data->res_edge1 = NULL;
   ml_data->res_edge2 = NULL;
   ml_data->rhs_nodal1 = NULL;
   ml_data->rhs_nodal2 = NULL;
   ml_data->x_nodal1 = NULL;
   ml_data->x_nodal2 = NULL;
   ml_data->edge_update1 = NULL;
   ml_data->edge_update2 = NULL;
   ml_data->max_eig = 0.0;
   ml_data->omega = 1.0;
   ml_data->output_level = 2;
   ml_data->ml_nodal = NULL;
   ml_data->ml_edge = NULL;
   ml_data->reduced_smoother = 0;
   return(0);
}
 
/* ******************************************************************** */
/* Destructor for Block Sm_Hiptmair_Data                                */
/* ******************************************************************** */
 
void ML_Smoother_Destroy_BlockHiptmair_Data(void *data)
{
   ML_Sm_BlockHiptmair_Data *ml_data;
 
   ml_data = (ML_Sm_BlockHiptmair_Data *) data;
 
   if ( ml_data->ATmat_trans != NULL )
      ML_Operator_Destroy(&(ml_data->ATmat_trans));
 
   if ( ml_data->TtAT_diag != NULL )
      ML_free(ml_data->TtAT_diag);
 
   if ( ml_data->TtATmat != NULL )
      ML_Operator_Destroy(&(ml_data->TtATmat));

   if ( ml_data->res_edge != NULL )
      ML_free(ml_data->res_edge);

   if ( ml_data->res_edge1 != NULL )
      ML_free(ml_data->res_edge1);

   if ( ml_data->res_edge2 != NULL )
      ML_free(ml_data->res_edge2);

   if ( ml_data->rhs_nodal1 != NULL )
      ML_free(ml_data->rhs_nodal1);

   if ( ml_data->rhs_nodal2 != NULL )
      ML_free(ml_data->rhs_nodal2);

   if ( ml_data->x_nodal1 != NULL )
      ML_free(ml_data->x_nodal1);

   if ( ml_data->x_nodal2 != NULL )
      ML_free(ml_data->x_nodal2);

   if ( ml_data->edge_update1 != NULL )
      ML_free(ml_data->edge_update1);

   if ( ml_data->edge_update2 != NULL )
      ML_free(ml_data->edge_update2);

   if ( (ml_data->sm_nodal != NULL) && (ml_data->sm_nodal->my_level != NULL) )
   {
      ML_free(ml_data->sm_nodal->my_level);
   }

   if ( ml_data->sm_nodal != NULL )
      ML_Smoother_Destroy(&(ml_data->sm_nodal));

   if ( ml_data->ml_nodal != NULL ){
     ML_Destroy(&(ml_data->ml_nodal));
   }

   if ( ml_data->ml_edge != NULL ) {
     ML_Destroy(&(ml_data->ml_edge));
   }


   ML_memory_free((void**) &ml_data);
}                                                                               
 
/* ******************************************************************** */
/* Destructor for Sm_Hiptmair_Data                                          */
/* ******************************************************************** */
 
void ML_Smoother_Destroy_Hiptmair_Data(void *data)
{
   ML_Sm_Hiptmair_Data *ml_data;
 
   ml_data = (ML_Sm_Hiptmair_Data *) data;
 
   if ( ml_data->ATmat_trans != NULL )
      ML_Operator_Destroy(&(ml_data->ATmat_trans));
 
   if ( ml_data->TtAT_diag != NULL )
      ML_free(ml_data->TtAT_diag);
 
   if ( ml_data->TtATmat != NULL )
      ML_Operator_Destroy(&(ml_data->TtATmat));

   if ( ml_data->res_edge != NULL )
      ML_free(ml_data->res_edge);

   if ( ml_data->rhs_nodal != NULL )
      ML_free(ml_data->rhs_nodal);

   if ( ml_data->x_nodal != NULL )
      ML_free(ml_data->x_nodal);

   if ( ml_data->edge_update != NULL )
      ML_free(ml_data->edge_update);

   if ( (ml_data->sm_nodal != NULL) && (ml_data->sm_nodal->my_level != NULL) )
   {
      ML_free(ml_data->sm_nodal->my_level);
   }

   if ( ml_data->sm_nodal != NULL )
      ML_Smoother_Destroy(&(ml_data->sm_nodal));

   if ( ml_data->ml_nodal != NULL ){
     ML_Destroy(&(ml_data->ml_nodal));
   }

   if ( ml_data->ml_edge != NULL ) {
     ML_Destroy(&(ml_data->ml_edge));
   }


   ML_memory_free((void**) &ml_data);
}                                                                               
/* ******************************************************************** */
/* Constructor for Sm_BGS_Data                                          */
/* ******************************************************************** */

int ML_Smoother_Create_BGS_Data(ML_Sm_BGS_Data **data)
{
   ML_Sm_BGS_Data *ml_data;

   ML_memory_alloc((void**) data, sizeof(ML_Sm_BGS_Data), "BGS" );
   ml_data = (*data);
   ml_data->blockfacts = NULL;
   ml_data->perms = NULL;
   ml_data->blocksize = -1;
   ml_data->blocklengths = NULL;
   ml_data->blockmap = NULL;
   return 0;
}

/* ******************************************************************** */
/* Destructor for Sm_BGS_Data                                          */
/* ******************************************************************** */

void ML_Smoother_Destroy_BGS_Data(void *data)
{
   int i;
   ML_Sm_BGS_Data *ml_data;

   ml_data = (ML_Sm_BGS_Data *) data;
   if ( ml_data->blockfacts != NULL )
   {
      for ( i = 0; i < ml_data->Nblocks; i++ )
         if ( ml_data->blockfacts[i] != NULL ) ML_free(ml_data->blockfacts[i]);
      ML_free( ml_data->blockfacts );
   }
   if ( ml_data->perms != NULL )
   {
      for ( i = 0; i < ml_data->Nblocks; i++ )
         if ( ml_data->perms[i] != NULL ) ML_free(ml_data->perms[i]);
      ML_free( ml_data->perms );
   }
   if ( ml_data->blocklengths != NULL )
      ML_free( ml_data->blocklengths );
   ML_free(ml_data->blockmap);
   ML_memory_free((void**) &ml_data);
}

/* ************************************************************************* */
/* clean up the BGS data structure (naming unconventional)                   */
/* ************************************************************************* */

void ML_Smoother_Clean_BGS_Data(void *data)
{
   int            Nblocks, i, **perms;
   double         **blockfacts;
   ML_Sm_BGS_Data *dataptr;

   dataptr = (ML_Sm_BGS_Data *) data;

   Nblocks = dataptr->Nblocks;
   perms = dataptr->perms;
   blockfacts = dataptr->blockfacts;

   for (i=0; i < Nblocks; i++) {
      ML_free(perms[i]);
      ML_free(blockfacts[i]);
   }

   ML_free(perms);
   ML_free(blockfacts);
   ML_memory_free((void **) &dataptr);
}

/* ************************************************************************* */
/* Constructor for ML_Sm_ILUT_Data                                           */
/* ************************************************************************* */

int ML_Smoother_Create_ILUT_Data(ML_Sm_ILUT_Data **data)
{
   ML_Sm_ILUT_Data *ml_data;

   ML_memory_alloc((void**) data, sizeof(ML_Sm_ILUT_Data), "SMI");
   ml_data = (ML_Sm_ILUT_Data *) (*data);
   ml_data->mat_ia      = NULL;
   ml_data->mat_ja      = NULL;
   ml_data->mat_aa      = NULL;
   ml_data->getrow_comm = NULL;
   ml_data->Nrows       = 0;
   ml_data->fillin      = 0;
   ml_data->threshold   = 0;
   return 0;
}

/* ************************************************************************* */
/* Destructor for ML_Sm_ILUT_Data                                            */
/* ************************************************************************* */

void ML_Smoother_Destroy_ILUT_Data(void *data)
{
   ML_Sm_ILUT_Data *ml_data;

   ml_data = (ML_Sm_ILUT_Data *) data;
   if ( ml_data->mat_ia != NULL ) ML_free(ml_data->mat_ia);
   if ( ml_data->mat_ja != NULL ) ML_free(ml_data->mat_ja);
   if ( ml_data->mat_aa != NULL ) ML_free(ml_data->mat_aa);
   ML_memory_free( (void **) &ml_data);
}

#ifdef out
/* ************************************************************************* */
/* Constructor for ML_Sm_Schwarz_Data                                        */
/* ************************************************************************* */

int ML_Smoother_Create_Schwarz_Data(ML_Sm_Schwarz_Data **data)
{
   ML_Sm_Schwarz_Data *ml_data;

   ML_memory_alloc((void**) data, sizeof(ML_Sm_Schwarz_Data), "SMI");
   ml_data = (ML_Sm_Schwarz_Data *) (*data);
   ml_data->bmat_ia      = NULL;
   ml_data->bmat_ja      = NULL;
   ml_data->bmat_aa      = NULL;
   ml_data->aux_bmat_ia  = NULL;
   ml_data->aux_bmat_ja  = NULL;
   ml_data->aux_bmat_aa  = NULL;
   ml_data->getrow_comm  = NULL;
   ml_data->Nrows        = 0;
   ml_data->blk_size     = NULL;
   ml_data->blk_info     = NULL;
   ml_data->blk_indices  = NULL;
#ifdef SUPERLU
   ml_data->slu_Amat     = NULL;
   ml_data->slu_Lmat     = NULL;
   ml_data->slu_Umat     = NULL;
#endif 
   ml_data->perm_c       = NULL;
   ml_data->perm_r       = NULL;
   return 0;
}

/* ************************************************************************* */
/* Destructor for ML_Sm_Schwarz_Data                                         */
/* ************************************************************************* */

void ML_Smoother_Destroy_Schwarz_Data(void *data)
{
   int                i;
   ML_Sm_Schwarz_Data *ml_data;
#ifdef SUPERLU
   SuperMatrix        *A, *L, *U;
#endif

   ml_data = (ML_Sm_Schwarz_Data *) data;
   if ( ml_data->bmat_ia  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_ia[i]);
      ML_free(ml_data->bmat_ia);
   }
   if ( ml_data->bmat_ja  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_ja[i]);
      ML_free(ml_data->bmat_ja);
   }
   if ( ml_data->bmat_aa  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_aa[i]);
      ML_free(ml_data->bmat_aa);
   }
   if ( ml_data->aux_bmat_ia  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_ia[i]);
      ML_free(ml_data->aux_bmat_ia);
   }
   if ( ml_data->aux_bmat_ja  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_ja[i]);
      ML_free(ml_data->aux_bmat_ja);
   }
   if ( ml_data->aux_bmat_aa  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_aa[i]);
      ML_free(ml_data->aux_bmat_aa);
   }
   if ( ml_data->blk_size != NULL ) ML_free(ml_data->blk_size);
   if ( ml_data->blk_info != NULL ) ML_free(ml_data->blk_info);
   if ( ml_data->blk_indices != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->blk_indices[i] != NULL ) 
            ML_free( ml_data->blk_indices[i] );
   }
#ifdef SUPERLU
   if ( ml_data->slu_Amat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         A = ml_data->slu_Amat[i];
         if ( A != NULL )
         {
            SUPERLU_FREE( A->Store );
            ML_free(A);
         }
      }
      ML_free( ml_data->slu_Amat );
   }
   if ( ml_data->slu_Lmat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         L = ml_data->slu_Lmat[i];
         if ( L != NULL )
         {
            Destroy_SuperNode_Matrix(L);
            ML_free(L);
         }
      }
      ML_free( ml_data->slu_Lmat );
   }
   if ( ml_data->slu_Umat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         U = ml_data->slu_Umat[i];
         if ( U != NULL )
         {
            SUPERLU_FREE( ((NRformat *) U->Store)->colind);
            SUPERLU_FREE( ((NRformat *) U->Store)->rowptr);
            SUPERLU_FREE( ((NRformat *) U->Store)->nzval);
            SUPERLU_FREE( U->Store );
            ML_free(U);
         }
      }
      ML_free( ml_data->slu_Umat );
   }
   if ( ml_data->perm_c != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->perm_c[i] ) ML_free(ml_data->perm_c[i]);
      ML_free( ml_data->perm_c );
   }
   if ( ml_data->perm_r != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->perm_r[i] ) ML_free(ml_data->perm_r[i]);
      ML_free( ml_data->perm_r );
   }
#endif
   ML_memory_free( (void **) &ml_data);
}

#endif
/* ************************************************************************* */
/* Function to generate the matrix products needed in the Hiptmair smoother  */
/* on one level.                                                             */
/* ************************************************************************ */
#include "ml_operator_blockmat.h"

#ifdef MatrixProductHiptmair

/* Block Hiptmair stuff */
int ML_Smoother_Gen_BlockHiptmair_Data(ML_Sm_BlockHiptmair_Data **data,
		                 ML_Operator *Amat,
                                 ML_Operator *Tmat, ML_Operator *Tmat_trans,
                                 ML_Operator *Tmat_bc, int BClength,
				  int *BCindices, 
void *edge_smoother, void **edge_args, void *nodal_smoother, void **nodal_args)
{

   ML_Sm_BlockHiptmair_Data *dataptr;
   ML_Operator *tmpmat, *tmpmat2, *Mmat;
   ML_1Level *mylevel;
   struct ML_CSR_MSRdata *matdata;
   int *row_ptr, i, j, k;
   double *val_ptr;
   double *dbl_arg1, *diagonal;
   struct ML_Operator_blockmat_data *mat_data;
#ifdef GREG
   struct MLSthing *widget;
#endif
#ifdef ML_TIMING_DETAILED
   double t0;

   t0 = GetClock();
#endif

   mat_data = (struct ML_Operator_blockmat_data *) Amat->data;
   Mmat = mat_data->M_mat;
   dataptr = *data;
   dataptr->Tmat_trans = Tmat_trans;
   dataptr->Tmat = Tmat;
   dataptr->output_level = 2.0;
   dataptr->omega = 1.0;

   /* Get maximum eigenvalue for damping parameter. */

  if ( (edge_smoother == (void *) ML_Gen_Smoother_Jacobi) ||
       (edge_smoother == (void *) ML_Gen_Smoother_GaussSeidel) ||
       (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) ||
       (edge_smoother == (void *) ML_Gen_Smoother_VBlockJacobi) ||
       (edge_smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel) )
  {
    dbl_arg1 = (double *) ML_Smoother_Arglist_Get(edge_args, 1);
    dataptr->omega = *dbl_arg1;

    if ((( (int) dbl_arg1[0]) == ML_DEFAULT) && (Amat->comm->ML_nprocs != 1))
    {
      dataptr->max_eig = ML_Operator_GetMaxEig(Amat);
      dataptr->omega = 1.0 / dataptr->max_eig;
      if (Amat->comm->ML_mypid == 0
          && dataptr->output_level < ML_Get_PrintLevel())
        printf("E:Calculated max eigenvalue of %f.\n",dataptr->max_eig);
    }
    if (Amat->comm->ML_mypid == 0
        && dataptr->output_level < ML_Get_PrintLevel())
    {
	     printf("Ke: Total nonzeros = %d (Nrows = %d)\n",Amat->N_nonzeros,
		    Amat->invec_leng);
         printf("E:Using Hiptmair damping factor of %f.\n",dataptr->omega);
         fflush(stdout);
    }  
  }
   /* Check matrix dimensions. */

   if (Tmat_trans->invec_leng != Mmat->outvec_leng)
   {
     printf("In ML_Smoother_Gen_BlockHiptmair_Data: Tmat_trans and Mmat\n");
     printf("\tdimensions do not agree:\n");
     printf("\tTmat_trans->invec_leng = %d, Mmat->outvec_leng = %d\n",
			 Tmat_trans->invec_leng, Mmat->outvec_leng);
      exit(1);
   }
   if ( dataptr->Tmat_trans->invec_leng != Mmat->outvec_leng )
   {
     printf("In ML_Smoother_Gen_BlockHiptmair_Data: Tmat_trans and Mmat\n");
     printf("\tdimensions do not agree:\n");
     printf("\tATmat_trans->invec_leng = %d, Mmat->outvec_leng = %d\n",
			 dataptr->Tmat_trans->invec_leng, Mmat->outvec_leng);
      exit(1);
   }
   if ( Mmat->invec_leng != Tmat->outvec_leng )
   {
     printf("In ML_Smoother_Gen_BlockHiptmair_Data: Mmat and Tmat\n");
     printf("\tdimensions do not agree:\n");
     printf("\tMmat->invec_leng = %d, Tmat->outvec_leng = %d\n",
			 Mmat->invec_leng, Tmat->outvec_leng);
      exit(1);
   }
   ML_Smoother_HiptmairSubsmoother_Create(&(dataptr->ml_edge),Amat,
					  edge_smoother, edge_args, 
					  dataptr->omega);

   /* For the MLS ... we really want to switch to a complex version */
#ifdef GREG
   if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {

     dataptr->ml_edge->pre_smoother[0].smoother->internal = ML_complex_Cheby;
#ifdef IWANTONESTEP
     dataptr->ml_edge->pre_smoother[0].smoother->internal = ML_DiagScaled_1stepKrylov;
#endif
     widget = (struct MLSthing *) dataptr->ml_edge->pre_smoother[0].smoother->data;
     if (Amat->lambda_max_img != 0.0) {
       widget->beta_real = Amat->lambda_max;
       widget->beta_img  = Amat->lambda_max_img;
     }
     else {
       widget->beta_real = 3.265411492391633/1.1;
       widget->beta_img  = .02605229219923984/1.1;
     }
   }
#endif

   /* Triple matrix product T^{*}AT. */


   tmpmat = ML_Operator_Create(Mmat->comm);
   if (Tmat_bc != NULL)
   {
      tmpmat2 = ML_Operator_Create(Mmat->comm);
      /* Calculate matrix product Ke * T.  Postprocess to get (almost)
         the same result matrix as bc(Ke) * T, where bc(Ke) is Ke with
         Dirichlet b.c. applied to both rows and columns of Ke.  The
         only differences will be the b.c.  rows themselves, which
         will be zero. */
      ML_2matmult(Mmat,Tmat_bc,tmpmat2, ML_CSR_MATRIX);
      matdata = (struct ML_CSR_MSRdata *) (tmpmat2->data);
      row_ptr = matdata->rowptr;
      val_ptr = matdata->values;

      for (i=0; i < BClength; i++)
      {
         j = BCindices[i];
         /* Zero out corresponding row in product Ke * T. */
         for (k = row_ptr[j]; k < row_ptr[j+1]; k++)
            val_ptr[k] = 0.0;
      }
      ML_2matmult(Tmat_trans,tmpmat2,tmpmat, ML_CSR_MATRIX);
      ML_Operator_Destroy(&tmpmat2);
   }
   else
   {
      ML_rap(Tmat_trans, Mmat, Tmat, tmpmat, ML_MSR_MATRIX);
      /* Some garbage code to fix up the case when sigma is */
      /* very small and so tmpmat is very small. Probably   */
      /* something better should be put in here!!!!         */

      matdata = (struct ML_CSR_MSRdata *) (tmpmat->data);
      if (tmpmat->diagonal != NULL) {
	ML_DVector_GetDataPtr( tmpmat->diagonal, &diagonal);
	for (i = 0; i < tmpmat->outvec_leng; i++) {
	  if ( fabs(diagonal[i]) < 1.0e-10)  {
	    matdata->values[i] = 1.;
	    diagonal[i] = 1.;
	  }
	}
      }
   }
/*
   kdata = ML_Krylov_Create( tmpmat->comm );
   ML_Krylov_Set_ComputeEigenvalues( kdata );
   ML_Krylov_Set_PrintFreq( kdata, 0 );
   ML_Krylov_Set_Amatrix(kdata, tmpmat);
   ML_Krylov_Solve(kdata, tmpmat->outvec_leng, NULL, NULL);
   if (tmpmat->comm->ML_mypid == 0)
   {
      printf("N:Calculated max eigenvalue of %lf.\n",dataptr->max_eig);
      printf("N:Using Hiptmair damping factor of %lf.\n",dataptr->omega);
      fflush(stdout);
   }  
*/
   /* Create ML_Smoother data structure for nodes.
      This is used in symmetric GS sweep over nodes. */
   mylevel = (ML_1Level *) ML_allocate(sizeof(ML_1Level));
   ML_Smoother_Create(&(dataptr->sm_nodal), mylevel);
   dataptr->sm_nodal->ntimes = 1;
   dataptr->sm_nodal->omega = 1.0;

  if ( (nodal_smoother == (void *) ML_Gen_Smoother_Jacobi) ||
       (nodal_smoother == (void *) ML_Gen_Smoother_GaussSeidel) ||
       (nodal_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) ||
       (nodal_smoother == (void *) ML_Gen_Smoother_VBlockJacobi) ||
       (nodal_smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel) )
  {
    dbl_arg1 = (double *) ML_Smoother_Arglist_Get(nodal_args, 1);
    dataptr->omega = *dbl_arg1;

    if ( ((int) dbl_arg1[0]) == ML_DEFAULT) dataptr->omega = 1.0;
    else dataptr->omega = *dbl_arg1;
  }

/*
   dataptr->sm_nodal->omega = 1.0 / ML_Krylov_Get_MaxEigenvalue(kdata);
   ML_Krylov_Destroy(&kdata);
*/

   dataptr->sm_nodal->my_level->Amat = tmpmat;
   dataptr->sm_nodal->my_level->comm = tmpmat->comm;
   dataptr->TtATmat = tmpmat;

   ML_Smoother_HiptmairSubsmoother_Create(&(dataptr->ml_nodal),tmpmat,
					  nodal_smoother, nodal_args, 
					  dataptr->omega);

   /* Allocate some work vectors that are needed in the smoother. */

   dataptr->res_edge = (double*) ML_allocate( Amat->invec_leng*2*sizeof(double));
   dataptr->res_edge1 = (double*) ML_allocate( Amat->invec_leng*sizeof(double));
   dataptr->res_edge2 = (double*) ML_allocate( Amat->invec_leng*sizeof(double));
   dataptr->rhs_nodal1 = (double *)
                        ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->rhs_nodal2 = (double *)
                        ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->x_nodal1 = (double *)
                      ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->x_nodal2 = (double *)
                      ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->edge_update1 = (double * )
                          ML_allocate(Amat->invec_leng * sizeof(double));
   dataptr->edge_update2 = (double * )
                          ML_allocate(Amat->invec_leng * sizeof(double));
   /*
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
#endif       
   */
   return 0;
}

/******************************************************************************/

int ML_Smoother_Gen_Hiptmair_Data(ML_Sm_Hiptmair_Data **data, ML_Operator *Amat,
                                 ML_Operator *Tmat, ML_Operator *Tmat_trans,
                                 ML_Operator *Tmat_bc, int BClength,
				  int *BCindices, 
void *edge_smoother, void **edge_args, void *nodal_smoother, void **nodal_args)
{

   ML_Sm_Hiptmair_Data *dataptr;
   ML_Operator *tmpmat, *tmpmat2;
   ML_1Level *mylevel;
   struct ML_CSR_MSRdata *matdata;
   int *row_ptr, i, j, k;
   double *val_ptr;
   double *dbl_arg1, *diagonal;
#ifdef ML_TIMING_DETAILED
   double t0;

   t0 = GetClock();
#endif

   dataptr = *data;
   dataptr->Tmat_trans = Tmat_trans;
   dataptr->Tmat = Tmat;
   dataptr->output_level = 2.0;
   dataptr->omega = 1.0;

   /* Get maximum eigenvalue for damping parameter. */

  if ( (edge_smoother == (void *) ML_Gen_Smoother_Jacobi) ||
       (edge_smoother == (void *) ML_Gen_Smoother_GaussSeidel) ||
       (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) ||
       (edge_smoother == (void *) ML_Gen_Smoother_VBlockJacobi) ||
       (edge_smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel) )
  {
    dbl_arg1 = (double *) ML_Smoother_Arglist_Get(edge_args, 1);
    if ((( (int) dbl_arg1[0]) == ML_DEFAULT) && (Amat->comm->ML_nprocs != 1))
    {
      dataptr->max_eig = ML_Operator_GetMaxEig(Amat);
      dataptr->omega = 1.0 / dataptr->max_eig;
      if (Amat->comm->ML_mypid == 0
          && dataptr->output_level < ML_Get_PrintLevel())
        printf("E:Calculated max eigenvalue of %f.\n",dataptr->max_eig);
    }
    if (Amat->comm->ML_mypid == 0
        && dataptr->output_level < ML_Get_PrintLevel())
    {
	     printf("Ke: Total nonzeros = %d (Nrows = %d)\n",Amat->N_nonzeros,
		    Amat->invec_leng);
         printf("E:Using Hiptmair damping factor of %f.\n",dataptr->omega);
         fflush(stdout);
    }  
  }

   /* Check matrix dimensions. */

   if (Tmat_trans->invec_leng != Amat->outvec_leng)
   {
      printf("In ML_Smoother_Gen_Hiptmair_Data: Tmat_trans and Amat\n");
      printf("\tdimensions do not agree:\n");
      printf("\tTmat_trans->invec_leng = %d, Amat->outvec_leng = %d\n",
			 Tmat_trans->invec_leng, Amat->outvec_leng);
      exit(1);
   }
   if ( dataptr->Tmat_trans->invec_leng != Amat->outvec_leng )
   {
     printf("In ML_Smoother_Gen_Hiptmair_Data: Tmat_trans and Amat\n");
     printf("\tdimensions do not agree:\n");
     printf("\tATmat_trans->invec_leng = %d, Amat->outvec_leng = %d\n",
			 dataptr->Tmat_trans->invec_leng, Amat->outvec_leng);
      exit(1);
   }
   if ( Amat->invec_leng != Tmat->outvec_leng )
   {
     printf("In ML_Smoother_Gen_Hiptmair_Data: Amat and Tmat\n");
     printf("\tdimensions do not agree:\n");
     printf("\tAmat->invec_leng = %d, Tmat->outvec_leng = %d\n",
			 Amat->invec_leng, Tmat->outvec_leng);
      exit(1);
   }
   ML_Smoother_HiptmairSubsmoother_Create(&(dataptr->ml_edge),Amat,
					  edge_smoother, edge_args, 
					  dataptr->omega);


   /* Triple matrix product T^{*}AT. */


   tmpmat = ML_Operator_Create(Amat->comm);
   if (Tmat_bc != NULL)
   {
      tmpmat2 = ML_Operator_Create(Amat->comm);
      /* Calculate matrix product Ke * T.  Postprocess to get (almost)
         the same result matrix as bc(Ke) * T, where bc(Ke) is Ke with
         Dirichlet b.c. applied to both rows and columns of Ke.  The
         only differences will be the b.c.  rows themselves, which
         will be zero. */
      ML_2matmult(Amat,Tmat_bc,tmpmat2, ML_CSR_MATRIX);
      matdata = (struct ML_CSR_MSRdata *) (tmpmat2->data);
      row_ptr = matdata->rowptr;
      val_ptr = matdata->values;

      for (i=0; i < BClength; i++)
      {
         j = BCindices[i];
         /* Zero out corresponding row in product Ke * T. */
         for (k = row_ptr[j]; k < row_ptr[j+1]; k++)
            val_ptr[k] = 0.0;
      }
      ML_2matmult(Tmat_trans,tmpmat2,tmpmat, ML_CSR_MATRIX);
      ML_Operator_Destroy(&tmpmat2);
   }
   else
   {
      ML_rap(Tmat_trans, Amat, Tmat, tmpmat, ML_MSR_MATRIX);

      /* Some garbage code to fix up the case when sigma is */
      /* very small and so tmpmat is very small. Probably   */
      /* something better should be put in here!!!!         */

      matdata = (struct ML_CSR_MSRdata *) (tmpmat->data);
      if (tmpmat->diagonal != NULL) {
	ML_DVector_GetDataPtr( tmpmat->diagonal, &diagonal);
	for (i = 0; i < tmpmat->outvec_leng; i++) {
	  if ( fabs(diagonal[i]) < 1.0e-10)  {
	    matdata->values[i] = 1.;
	    diagonal[i] = 1.;
	  }
	}
      }
   }
/*
   kdata = ML_Krylov_Create( tmpmat->comm );
   ML_Krylov_Set_ComputeEigenvalues( kdata );
   ML_Krylov_Set_PrintFreq( kdata, 0 );
   ML_Krylov_Set_Amatrix(kdata, tmpmat);
   ML_Krylov_Solve(kdata, tmpmat->outvec_leng, NULL, NULL);
   if (tmpmat->comm->ML_mypid == 0)
   {
      printf("N:Calculated max eigenvalue of %lf.\n",dataptr->max_eig);
      printf("N:Using Hiptmair damping factor of %lf.\n",dataptr->omega);
      fflush(stdout);
   }  
*/
   /* Create ML_Smoother data structure for nodes.
      This is used in symmetric GS sweep over nodes. */
   mylevel = (ML_1Level *) ML_allocate(sizeof(ML_1Level));
   ML_Smoother_Create(&(dataptr->sm_nodal), mylevel);
   dataptr->sm_nodal->ntimes = 1;
   dataptr->sm_nodal->omega = 1.0;

/*
   dataptr->sm_nodal->omega = 1.0 / ML_Krylov_Get_MaxEigenvalue(kdata);
   ML_Krylov_Destroy(&kdata);
*/

   dataptr->sm_nodal->my_level->Amat = tmpmat;
   dataptr->sm_nodal->my_level->comm = tmpmat->comm;
   dataptr->TtATmat = tmpmat;

   ML_Smoother_HiptmairSubsmoother_Create(&(dataptr->ml_nodal),tmpmat,
					  nodal_smoother, nodal_args, 
					  dataptr->omega);

   /* Allocate some work vectors that are needed in the smoother. */

   dataptr->res_edge = (double *) ML_allocate( Amat->invec_leng*sizeof(double));
   dataptr->rhs_nodal = (double *)
                        ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->x_nodal = (double *)
                      ML_allocate(Tmat->invec_leng * sizeof(double));
   dataptr->edge_update = (double * )
                          ML_allocate(Amat->invec_leng * sizeof(double));
   /*
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
#endif       
   */
   return 0;
}
#endif /* ifdef MatrixProductHiptmair */

/*******************************************************************************
* The old ML_Smoother_Gen_Hiptmair_Data function.                              *
*******************************************************************************/

#ifdef MatrixFreeHiptmair
#ifdef debugSmoother
extern double *xxx;
#endif

int ML_Smoother_Gen_Hiptmair_Data(ML_Sm_Hiptmair_Data **data, ML_Operator *Amat,
                                 ML_Operator *Tmat, ML_Operator *Tmat_trans,
                                 double print)
{
   ML_Sm_Hiptmair_Data *dataptr;
   ML_Operator *tmpmat, *eyemat;
   int *cols = NULL, length, allocated_space;
   double *vals = NULL;
   double *tmpdiag;
#ifdef debugSmoother
   double yyy[5000];
#endif
   int i,j;
   ML_CommInfoOP *getrow_comm;
   int totalsize;
   struct ML_CSR_MSRdata *Amat_data, *eye_csr_data;
   ML_1Level *mylevel;

   dataptr = *data;
   dataptr->Tmat_trans = Tmat_trans;
   dataptr->Tmat = Tmat;
   dataptr->ATmat_trans = ML_Operator_Create(Tmat_trans->comm);

   /* We access AT by column.  So calculate the transpose of AT
      and access it by row. This assumes that A is symmetric. */


   if (Tmat_trans->invec_leng != Amat->outvec_leng)
   {
      printf("In ML_Smoother_Gen_Hiptmair_Data: Tmat_trans and Amat\n"
	         "\tdimensions do not agree:\n"
			 "\tTmat_trans->invec_leng = %d, Amat->outvec_leng = %d\n",
			 Tmat_trans->invec_leng, Amat->outvec_leng);
      exit(1);
   }
   ML_2matmult(Tmat_trans,Amat,dataptr->ATmat_trans, ML_CSR_MATRIX);


#ifdef debugSmoother
   ML_Operator_Apply(dataptr->ATmat_trans,
                     dataptr->ATmat_trans->invec_leng, xxx,
					 dataptr->ATmat_trans->outvec_leng,
					 yyy);
   printf("norm of yyy = %15.10e\n",
              sqrt(ML_gdot(dataptr->ATmat_trans->outvec_leng,
              yyy,yyy,dataptr->ATmat_trans->comm)));
#endif
   /*
   tmpmat = ML_Operator_Create(Amat->comm);
   ML_2matmult(Amat,Tmat,tmpmat, ML_CSR_MATRIX);

   ML_Operator_Transpose_byrow(tmpmat, dataptr->ATmat_trans);

   printf("%d: (Tmat_trans) outvec_leng = %d, invec_leng = %d\n",
          Tmat_trans->comm->ML_mypid,
		  dataptr->ATmat_trans->outvec_leng,
		  dataptr->ATmat_trans->invec_leng);

   getrow_comm= dataptr->ATmat_trans->getrow->pre_comm;
   if (getrow_comm != NULL) 
   {
      printf("%d: (ATmat_trans) num ghost rows = %d\n",
	         Tmat_trans->comm->ML_mypid,
	         getrow_comm->total_rcv_length);
      getrow_comm= dataptr->Tmat_trans->getrow->pre_comm;
      printf("%d: (Tmat_trans) num ghost rows = %d\n",
	         Tmat_trans->comm->ML_mypid,
	         getrow_comm->total_rcv_length);
      getrow_comm= Amat->getrow->pre_comm;
      printf("%d: (Amat) num ghost rows = %d\n",
	         Amat->comm->ML_mypid,
	         getrow_comm->total_rcv_length);
      fflush(stdout);
   }
   */

   /* Generate matrix triple product T^{*}AT. Use the already calculated
     product T^{*} A. */

   tmpmat = ML_Operator_Create(Tmat_trans->comm);
   if (dataptr->ATmat_trans->invec_leng != Tmat->outvec_leng)
   {
      printf("In ML_Smoother_Gen_Hiptmair_Data: ATmat_trans and Tmat\n"
	         "\tdimensions do not agree:\n"
			 "\tATmat_trans->invec_leng = %d, Tmat->outvec_leng = %d\n",
			 dataptr->ATmat_trans->invec_leng, Tmat->outvec_leng);
      exit(1);
   }
   ML_2matmult(dataptr->ATmat_trans, Tmat, tmpmat, ML_CSR_MATRIX);

   /* and pick off the diagonal entries. */
   dataptr->TtAT_diag = (double *)
                        ML_allocate( tmpmat->outvec_leng * sizeof(double) );
   ML_Operator_Get_Diag(tmpmat, tmpmat->outvec_leng, &tmpdiag );
   for (i=0; i< tmpmat->outvec_leng; i++)
      if (tmpdiag[i] == 0) dataptr->TtAT_diag[i] = 1.;
	  else dataptr->TtAT_diag[i] = tmpdiag[i];

   ML_Operator_Destroy(&tmpmat);

   /* Force the diagonal to not be zero */

   vals = NULL;
   if (Amat->getrow->external == MSR_getrows){
      Amat_data   = (struct ML_CSR_MSRdata *) Amat->data;
      vals  = Amat_data->values;
      cols = Amat_data->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &cols, &vals);
#endif
   if (vals != NULL) {
      for (i = 0; i < Amat->invec_leng; i++) 
         if (vals[i] == 0.0) vals[i] = 1.;
   }

   return 0;
}
#endif /* ifdef MatrixFreeHiptmair */

/* ************************************************************************* */
/* function to generate the factorizations of the diagonal blocks of A.      */
/* Factorizations are computed using lapack                                  */
/* ************************************************************************* */

int ML_Smoother_Gen_BGSFacts(ML_Sm_BGS_Data **data, ML_Operator *Amat, 
                             int blocksize) 
{
   int            i, j, *cols, allocated_space, length, Nrows, Nblocks, kkk;
   int            row_in_block, col_in_block, **perms, info, col;
   double         *vals, **blockfacts;
   ML_Sm_BGS_Data *dataptr;

   Nrows = Amat->getrow->Nrows;
   if (Nrows % blocksize != 0)
   {
      printf("Error: BGS requires an integer no. of blocks on each proc\n");
      printf("       Nrows, blocksize = %d %d \n", Nrows, blocksize);
      exit(1);
   }
   dataptr = (*data);

   Nblocks = Nrows / blocksize;
   dataptr->Nblocks = Nblocks;
   allocated_space = Amat->max_nz_per_row+2;
   dataptr->blocksize = blocksize;

   dataptr->blockfacts = (double **)ML_allocate(Nblocks*sizeof(double *));
   dataptr->perms = (int **)ML_allocate(Nblocks*sizeof(int *));
   blockfacts = dataptr->blockfacts;
   perms = dataptr->perms;
   for (j=0; j<Nblocks; j++) 
   {
      blockfacts[j]=(double *) ML_allocate(blocksize*blocksize*sizeof(double));
      for (kkk = 0; kkk < blocksize*blocksize; kkk++) blockfacts[j][kkk] = 0.;
      perms[j]=(int *)ML_allocate(blocksize*blocksize*sizeof(int));
   }
   cols = (int    *) ML_allocate(allocated_space*sizeof(int    ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double ));

   if (vals == NULL) 
      pr_error("Error in ML_Gen_BGSFacts(): Not enough space\n");

   for (i = 0; i < Nrows; i++) 
   {
      row_in_block=i%blocksize;
      ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,&length,0);
      for (j = 0; j < length; j++) 
      {
         col = cols[j];
         col_in_block=col%blocksize;
         if ((col < i+blocksize-row_in_block) && (col >= i-row_in_block))
            blockfacts[i/blocksize][col_in_block*blocksize+row_in_block]=vals[j];
      }
   }
   for (i = 0; i < Nblocks; i++) 
   {
      MLFORTRAN(dgetrf)(&blocksize, &blocksize, blockfacts[i], &blocksize,
	         perms[i], &info);
      if (info != 0)
         pr_error("Error in ML_Gen_BGSFacts:dgetrf returned a non-zero value\n");
   }
   ML_free(cols);
   ML_free(vals);
	
   return 0;
}

/* ************************************************************************* */
/* generate the block GS factorization (variable size block).                */
/* Blocking information taken from ML_Aggregate.                             */
/* ************************************************************************* */

int ML_Smoother_Gen_VBGSFacts(ML_Sm_BGS_Data **data, ML_Operator *Amat,
                              int Nblocks, int *blockIndices)
{
   int            i, j, *cols, allocated_space, length, Nrows;
   int            row_in_block, col_in_block, **perms, info, col, index;
   int            *block_sizes, *block_offset, block_num;
   double         *vals, **blockfacts;
   ML_Sm_BGS_Data *dataptr;
   int Nnz;

   Nrows   = Amat->getrow->Nrows;
   dataptr = (*data);
   allocated_space = Amat->max_nz_per_row+2;

   /* ----------------------------------------------------------- */
   /* error checking                                              */
   /* ----------------------------------------------------------- */

   dataptr->Nblocks = Nblocks;
   if ( Nblocks < 0 || Nblocks > Nrows )
   {
      printf("ML_Gen_VBGSFacts : invalid blocking information.\n");
      printf("ML_Gen_VBGSFacts : Nblocks = %d.\n", Nblocks);
      exit(1);
   }

   if ( blockIndices == NULL )
   { 
      printf("ML_Gen_VBGSFacts : blocking information not available.\n");
      exit(1);
   }
   dataptr->blockmap = (int *) ML_allocate( (Nrows+1) * sizeof(int));
   if (dataptr->blockmap == NULL) 
      pr_error("ML_Smoother_Gen_VBGSFacts: out of space\n");
   for (i = 0; i < Nrows; i++) dataptr->blockmap[i] = blockIndices[i];

   dataptr->blocklengths = (int*) ML_allocate( (Nblocks+1) * sizeof(int));
   if (dataptr->blocklengths == NULL) 
      pr_error("ML_Smoother_Gen_VBGSFacts: out of space\n");

   block_sizes = dataptr->blocklengths;

   /* ----------------------------------------------------------- */
   /* search for sizes of each block                              */
   /* ----------------------------------------------------------- */

   for ( i = 0; i < Nblocks; i++ ) block_sizes[i] = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      if ( blockIndices[i] < 0 || blockIndices[i] >= Nblocks )
      {
         if ( blockIndices[i] != -1 )
         {
            printf("ML_Gen_VBGSFacts : block index not valid %d. %d\n",
                                       blockIndices[i],i);
            exit(1);
         }
      } else
         block_sizes[blockIndices[i]]++;
   }
   
   /* ----------------------------------------------------------- */
   /* allocate memory for each block                              */
   /* ----------------------------------------------------------- */

   dataptr->blockfacts = (double **) ML_allocate(Nblocks*sizeof(double *));

   dataptr->perms = (int **)ML_allocate(Nblocks*sizeof(int *));
   blockfacts = dataptr->blockfacts;
   perms = dataptr->perms;
   for ( i = 0; i < Nblocks; i++) 
   {
      length = block_sizes[i] * block_sizes[i];
      blockfacts[i] = (double *)ML_allocate( length* sizeof(double) );
      for ( j = 0; j < length; j++) blockfacts[i][j] = 0.0; 
      perms[i] = (int *)ML_allocate( block_sizes[i] * sizeof(int));
   }

   /* ----------------------------------------------------------- */
   /* load the block matrices                                     */
   /* ----------------------------------------------------------- */

   block_offset = (int *) ML_allocate(Nrows*sizeof(int ));
   cols = (int    *) ML_allocate(allocated_space*sizeof(int    ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double ));
   if (vals == NULL) 
      pr_error("Error in ML_Smoother_Gen_VBGSFacts: Not enough space\n");

   for ( i = 0; i < Nblocks; i++) block_sizes[i] = 0; 
   for (i = 0; i < Nrows; i++) 
   {
      block_num       = blockIndices[i];
      if ( blockIndices[i] >= 0 && blockIndices[i] < Nblocks )
         block_offset[i] = block_sizes[block_num]++;
   }
   for (i = 0; i < Nrows; i++) 
   {
      block_num    = blockIndices[i];
      if ( blockIndices[i] >= 0 && blockIndices[i] < Nblocks )
      {
         row_in_block = block_offset[i];
         ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,&length,0);
	 Nnz = 0;
         for (j = 0; j < length; j++) {
            col = cols[j];
            if ( col < Nrows )
            {
               if ( blockIndices[col] == block_num )
               {
		 if (vals[j] != 0.) Nnz++;
		 col_in_block = block_offset[col];
		 index = col_in_block * block_sizes[block_num] + row_in_block;
		 blockfacts[block_num][index] = vals[j];
               }
            }
         }
	 /* Handle the case of a zero row. */
	 /* By just putting a 1 on the diagonal. */
	 if (Nnz == 0) {
	   index = row_in_block * block_sizes[block_num] + row_in_block;
	   blockfacts[block_num][index] = 1.;
	 }
      }
   }

   /* ----------------------------------------------------------- */
   /* perform factorization on each block                         */
   /* ----------------------------------------------------------- */

   for (i = 0; i < Nblocks; i++) 
   {
      length = block_sizes[i];
      MLFORTRAN(dgetrf)(&length, &length, blockfacts[i], &length, 
                        perms[i], &info);
      if (info != 0)
      {
         printf("Error in ML_Smoother_Gen_VBGSFacts: dgetrf returned %d (!=0)\n",info);
         printf("This was caused by block %d of size %d\n",i,length);
         exit(1);
      }
   }

   /* ----------------------------------------------------------- */
   /* clean up                                                    */
   /* ----------------------------------------------------------- */

   ML_free(cols);
   ML_free(vals);
   ML_free(block_offset);
	
   return 0;
}

/*****************************************************************************/
/* needed for overlapped smoothers                                           */
/*****************************************************************************/

int ML_Smoother_ComposeOverlappedMatrix(ML_Operator *Amat, ML_Comm *comm,
                 int *total_recv_leng, int **recv_lengths, int **int_buf,
                 double **dble_buf, int **sindex_array, int **sindex_array2,
                 int *offset)
{
   int           i, nprocs, mypid, Nrows, *proc_array, *proc_array2;
   int           extNrows, NrowsOffset, *index_array, *index_array2;
   double        *dble_array;
   ML_CommInfoOP *getrow_comm;

   nprocs = comm->ML_nprocs;
   mypid  = comm->ML_mypid;
   Nrows   = Amat->getrow->Nrows;

   /* ----------------------------------------------------------- */
   /* see if communicator is present                              */
   /* ----------------------------------------------------------- */

   if (Amat->getrow->ML_id == ML_EMPTY)
      pr_error("Error(ComposeOverlappedMatrix): Need getrow()\n");

   if (Amat->getrow->post_comm != NULL)
      pr_error("ComposeOverlappedmatrix Post communication not implemented\n");

   getrow_comm = Amat->getrow->pre_comm;
   if (getrow_comm != NULL) 
      extNrows = Nrows + getrow_comm->total_rcv_length;
   else extNrows = Nrows;

   /* ----------------------------------------------------------- */
   /* compose NrowsOffset and processor offsets proc_array        */
   /* ----------------------------------------------------------- */
 
   proc_array  = (int *) ML_allocate(nprocs * sizeof(int) );
   proc_array2 = (int *) ML_allocate(nprocs * sizeof(int) );
   for ( i = 0; i < nprocs; i++ ) proc_array[i] = 0;
   proc_array[mypid] = Nrows;
   ML_gsum_vec_int(&proc_array, &proc_array2, nprocs, comm);
   NrowsOffset = 0;
   for (i = 0; i < mypid; i++) NrowsOffset += proc_array[i];
   for (i = 1; i < nprocs; i++) proc_array[i] += proc_array[i-1];
   ML_free(proc_array2);

   /* ----------------------------------------------------------- */
   /* compose the column index map (index_array,index_array2)     */
   /* ----------------------------------------------------------- */

   dble_array  = (double *) ML_allocate(extNrows *sizeof(double));
   for (i = Nrows; i < extNrows; i++) dble_array[i] = 0.0;
   for (i = 0; i < Nrows; i++) dble_array[i] = 1.0 * ( i + NrowsOffset );
   if (getrow_comm != NULL)
      ML_exchange_bdry(dble_array,getrow_comm, Nrows,comm,ML_OVERWRITE,NULL);
   index_array = ( int *) ML_allocate((extNrows-Nrows) * sizeof(int));
   for (i = Nrows; i < extNrows; i++) index_array[i-Nrows] = dble_array[i];
   index_array2  = (int *) ML_allocate((extNrows-Nrows) *sizeof(int));
   for (i = 0; i < extNrows-Nrows; i++) index_array2[i] = i;
   ML_free( dble_array );

   /* ----------------------------------------------------------- */
   /* send the lengths of each row to remote processor            */
   /* at the end, additional row information should be given      */
   /* in total_recv_leng, recv_lengths, int_buf, dble_buf         */
   /* ----------------------------------------------------------- */

   ML_Smoother_GetRowLengths(getrow_comm, comm, Amat, total_recv_leng, 
                             recv_lengths); 
   ML_Smoother_GetOffProcRows(getrow_comm, comm, Amat, *total_recv_leng, 
                              *recv_lengths, NrowsOffset, index_array,
                               int_buf, dble_buf); 

   ML_free(proc_array);
   ML_az_sort(index_array, extNrows-Nrows, index_array2, NULL);
   (*sindex_array) = index_array;
   (*sindex_array2) = index_array2;
   (*offset) = NrowsOffset;
   return 0;
}
      
/*****************************************************************************/
/* needed for overlapped smoothers                                           */
/*****************************************************************************/

int ML_Smoother_GetRowLengths(ML_CommInfoOP *comm_info, ML_Comm *comm, 
                              ML_Operator *Amat, int *leng, int **recv_leng)
{
   int     N_neighbors, *neighbors, total_recv, mtype, msgtype, proc_id;
   int     i, j, index, *temp_list, nbytes, length, offset, allocated_space;
   int     *cols, m, nnz;
   double  *vals;
   USR_REQ *request;

   /* ----------------------------------------------------------- */
   /* fetch communication information                             */
   /* ----------------------------------------------------------- */

   N_neighbors = ML_CommInfoOP_Get_Nneighbors(comm_info);
   if ( N_neighbors <= 0 ) { (*leng) = 0; (*recv_leng) = NULL; return 0;}

   neighbors  = ML_CommInfoOP_Get_neighbors(comm_info);
   total_recv = 0;
   for ( i = 0; i < N_neighbors; i++ )
      total_recv += ML_CommInfoOP_Get_Nrcvlist(comm_info,neighbors[i]);

   (*leng) = total_recv;

   /* ----------------------------------------------------------- */
   /* Set up send messages                                        */
   /* ----------------------------------------------------------- */

   request      = (USR_REQ  *) ML_allocate(N_neighbors*sizeof(USR_REQ ));
   (*recv_leng) = (int  *)     ML_allocate(total_recv * sizeof(int));

   mtype = 2001;

   /* ----------------------------------------------------------- */
   /* post receives for all messages                              */
   /* ----------------------------------------------------------- */

   offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nbytes  = sizeof(int) * length;
      comm->USR_irecvbytes((void *) &((*recv_leng)[offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
   }

   /* ----------------------------------------------------------- */
   /* write out all messages                                      */
   /* ----------------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row + 2;
   cols = (int *) ML_allocate(allocated_space * sizeof(int));
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id   = neighbors[i];
      length    = ML_CommInfoOP_Get_Nsendlist(comm_info, proc_id);
      temp_list = ML_CommInfoOP_Get_sendlist(comm_info, proc_id);
      nnz = 0;
      for (j = 0; j < length; j++) 
      {
         index = temp_list[j];
         ML_get_matrix_row(Amat,1,&index,&allocated_space,&cols,&vals,&m,0);
         temp_list[j] = m;
         nnz += m;
      }
      msgtype = mtype;   
      nbytes  = sizeof(int) * length;
      comm->USR_sendbytes((void*) temp_list, nbytes, proc_id, msgtype, 
                          comm->USR_comm);
      ML_free( temp_list );
   }
   ML_free(cols);
   ML_free(vals);

   /* ----------------------------------------------------------- */
   /* wait for all messages                                       */
   /* ----------------------------------------------------------- */

   offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nbytes  = sizeof(int) * length;
      comm->USR_cheapwaitbytes((void *) &((*recv_leng)[offset]), nbytes, &proc_id, 
                          &msgtype, comm->USR_comm, request+i);
      offset += length;
   }
   ML_free(request);
   ML_free(neighbors);
   return 0;
}

/*****************************************************************************/
/* needed for overlapped smoothers                                           */
/*****************************************************************************/

int ML_Smoother_GetOffProcRows(ML_CommInfoOP *comm_info, ML_Comm *comm, 
                           ML_Operator *Amat, int leng, int *recv_leng,
                           int Noffset, int *map, int **int_buf, 
                           double **dble_buf)
{
   int     N_neighbors, *neighbors, total_recv, mtype, msgtype, proc_id;
   int     i, j, index, *temp_list, nbytes, length, offset, allocated_space;
   int     *cols, *isend_buf=NULL, Nrows, m, k, nnz, nnz_offset;
   double  *vals, *send_buf=NULL;
   USR_REQ *request;

   /* ----------------------------------------------------------- */
   /* fetch communication information                             */
   /* ----------------------------------------------------------- */

   N_neighbors = ML_CommInfoOP_Get_Nneighbors(comm_info);
   if ( N_neighbors <= 0 ) { (*int_buf) = NULL; (*dble_buf) = NULL; return 0;}

   neighbors  = ML_CommInfoOP_Get_neighbors(comm_info);
   total_recv = 0;
   for ( i = 0; i < leng; i++ ) total_recv += recv_leng[i];

   Nrows   = Amat->getrow->Nrows;

   /* ----------------------------------------------------------- */
   /* Set up communication                                        */
   /* ----------------------------------------------------------- */

   request     = (USR_REQ  *) ML_allocate(N_neighbors*sizeof(USR_REQ ));
   (*int_buf)  = (int  *)     ML_allocate(total_recv * sizeof(int));
   (*dble_buf) = (double  *)  ML_allocate(total_recv * sizeof(double));

   mtype = 2002;

   /* ----------------------------------------------------------- */
   /* post receives for all messages                              */
   /* ----------------------------------------------------------- */

   offset = 0;
   nnz_offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nnz = 0;
      for (j = 0; j < length; j++)  nnz += recv_leng[offset+j];

      nbytes  = sizeof(double) * nnz;
      comm->USR_irecvbytes((void *) &((*dble_buf)[nnz_offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
      nnz_offset += nnz;
   }

   /* ----------------------------------------------------------- */
   /* write out all messages                                      */
   /* ----------------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row + 2;
   cols = (int *) ML_allocate(allocated_space * sizeof(int));
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id   = neighbors[i];
      length    = ML_CommInfoOP_Get_Nsendlist(comm_info, proc_id);
      temp_list = ML_CommInfoOP_Get_sendlist(comm_info, proc_id);
      nnz       = 0;
      for (j = 0; j < length; j++) 
      {
         index = temp_list[j];
         ML_get_matrix_row(Amat,1,&index,&allocated_space,&cols,&vals,&m,0);
         nnz += m;
      }
      if ( nnz > 0 ) send_buf = (double *) ML_allocate( nnz * sizeof(double));
      offset = 0;
      for (j = 0; j < length; j++) 
      {
         index = temp_list[j];
         ML_get_matrix_row(Amat,1,&index,&allocated_space,&cols,&vals,&m,0);
         for (k = 0; k < m; k++) send_buf[offset+k] = vals[k]; 
         offset += m;
      }
      msgtype = mtype;   
      nbytes  = sizeof(double) * nnz;
      comm->USR_sendbytes((void*) send_buf, nbytes, proc_id, msgtype, 
                          comm->USR_comm);
      ML_free( temp_list );
      ML_free( send_buf );
   }
   ML_free(cols);
   ML_free(vals);

   /* ----------------------------------------------------------- */
   /* wait for all messages                                       */
   /* ----------------------------------------------------------- */

   offset = 0;
   nnz_offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nnz = 0;
      for (j = 0; j < length; j++)  nnz += recv_leng[offset+j];
      nbytes  = sizeof(double) * nnz;
      comm->USR_cheapwaitbytes((void *) &((*dble_buf)[nnz_offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
      nnz_offset += nnz;
   }

   mtype = 2003;

   /* ----------------------------------------------------------- */
   /* post receives for all messages                              */
   /* ----------------------------------------------------------- */

   offset = 0;
   nnz_offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nnz = 0;
      for (j = 0; j < length; j++)  nnz += recv_leng[offset+j];
      nbytes  = sizeof(int) * nnz;
      comm->USR_irecvbytes((void *) &((*int_buf)[nnz_offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
      nnz_offset += nnz;
   }

   /* ----------------------------------------------------------- */
   /* write out all messages                                      */
   /* ----------------------------------------------------------- */

   allocated_space = Amat->max_nz_per_row + 2;
   cols = (int *) ML_allocate(allocated_space * sizeof(int));
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id   = neighbors[i];
      length    = ML_CommInfoOP_Get_Nsendlist(comm_info, proc_id);
      temp_list = ML_CommInfoOP_Get_sendlist(comm_info, proc_id);
      nnz       = 0;
      for (j = 0; j < length; j++) 
      {
         index = temp_list[j];
         ML_get_matrix_row(Amat,1,&index,&allocated_space,&cols,&vals,&m,0);
         nnz += m;
      }
      if ( nnz > 0 ) isend_buf = (int *) ML_allocate( nnz * sizeof(int));
      offset = 0;
      for (j = 0; j < length; j++) 
      {
         index = temp_list[j];
         ML_get_matrix_row(Amat,1,&index,&allocated_space,&cols,&vals,&m,0);
         for (k = 0; k < m; k++) 
         {
            if ( cols[k] < Nrows ) isend_buf[offset+k] = cols[k] + Noffset; 
            else                   isend_buf[offset+k] = map[cols[k]-Nrows];
         }
         offset += m;
      }
      msgtype = mtype;   
      nbytes  = sizeof(int) * nnz;
      comm->USR_sendbytes((void*) isend_buf, nbytes, proc_id, msgtype, 
                          comm->USR_comm);
      ML_free( temp_list );
      ML_free( isend_buf );
   }
   ML_free(cols);
   ML_free(vals);

   /* ----------------------------------------------------------- */
   /* wait for all messages                                       */
   /* ----------------------------------------------------------- */

   offset = 0;
   nnz_offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      proc_id = neighbors[i];
      msgtype = mtype;   
      length  = ML_CommInfoOP_Get_Nrcvlist(comm_info, proc_id);
      nnz = 0;
      for (j = 0; j < length; j++)  nnz += recv_leng[offset+j];
      nbytes  = sizeof(int) * nnz;
      comm->USR_cheapwaitbytes((void *) &((*int_buf)[nnz_offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
      nnz_offset += nnz;
   }

   ML_free(request);
   ML_free(neighbors);
   return 0;
}

/*****************************************************************************/
/* function for doing ILUT decomposition                                     */
/*****************************************************************************/

int ML_Smoother_ILUTDecomposition(ML_Sm_ILUT_Data *data, ML_Operator *Amat, 
             ML_Comm *comm, int total_recv_leng, int *recv_lengths, int *ext_ja, 
             double *ext_aa, int *map, int *map2, int Noffset)
{
   int             fillin, *mat_ia, *mat_ja, i, m, allocated_space, *cols;
   int             index, first, Lcount, Ucount, j, k, total_nnz;
   int             sortcnt, colIndex, offset, nnz_count, Nrows, extNrows;
   int             track_leng, *track_array, *sortcols;
   double          *vals, ddata, tau, *mat_aa, *diagonal, *rowNorms;
   double          *dble_buf, *sortvals, absval, rel_tau;
   ML_Sm_ILUT_Data *ilut_ptr;
#ifdef ML_DEBUG_SMOOTHER
   int             mypid;
#endif

   /* ---------------------------------------------------------- */
   /* fetch ILUT parameters                                      */
   /* ---------------------------------------------------------- */

#ifdef ML_DEBUG_SMOOTHER
   mypid       = comm->ML_mypid;
#else
   if (Noffset == -47) ML_avoid_unused_param((void *) comm);
#endif
   ilut_ptr    = (ML_Sm_ILUT_Data *) data;
   fillin      = ilut_ptr->fillin;
   tau         = ilut_ptr->threshold;
   Nrows       = Amat->outvec_leng;
   extNrows    = Nrows + total_recv_leng;
   ilut_ptr->Nrows = extNrows;

   /* ---------------------------------------------------------- */
   /* allocate temporary storage space                           */
   /* ---------------------------------------------------------- */

   allocated_space = extNrows;
   cols = (int *) ML_allocate(allocated_space * sizeof(int));
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   sortcols = (int *)    ML_allocate(extNrows * sizeof(int));
   sortvals = (double *) ML_allocate(extNrows * sizeof(double));
   dble_buf = (double *) ML_allocate(extNrows * sizeof(double));
   diagonal = (double *) ML_allocate(extNrows * sizeof(double));
   rowNorms = (double *) ML_allocate(extNrows * sizeof(double));

   /* ---------------------------------------------------------- */
   /* compute the storage requirement for the ILU matrix         */
   /* ---------------------------------------------------------- */

   total_nnz   = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      rowNorms[i] = 0.0;
      ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,&m,0);
      total_nnz += m;
      for ( j = 0; j < m; j++ ) 
         if ( cols[j] < extNrows ) rowNorms[i] += ML_dabs(vals[j]);
      rowNorms[i] /= extNrows;
   }
   for ( i = 0; i < total_recv_leng; i++ ) total_nnz += recv_lengths[i];
   total_nnz *= (fillin + 1);
   ilut_ptr->mat_ia = (int *) ML_allocate( (extNrows + 1 ) * sizeof(int));
   ilut_ptr->mat_ja = (int *) ML_allocate( total_nnz * sizeof(int));
   ilut_ptr->mat_aa = (double *) ML_allocate( total_nnz * sizeof(double));
   mat_ia = ilut_ptr->mat_ia;
   mat_ja = ilut_ptr->mat_ja;
   mat_aa = ilut_ptr->mat_aa;

   offset = 0;
   for ( i = 0; i < total_recv_leng; i++ )
   {
      rowNorms[i+Nrows] = 0.0;
      for ( j = offset; j < offset+recv_lengths[i]; j++ ) 
      {
         index = ext_ja[j];
         if ( index >= Noffset && index < Noffset+Nrows )
            ext_ja[j] = index - Noffset; 
         else
         {
            m = ML_sorted_search(index, extNrows-Nrows, map); 
            if ( m >= 0 ) ext_ja[j] = map2[m] + Nrows;
            else          ext_ja[j] = -1;
         }
         if ( ext_ja[j] != -1 ) rowNorms[i+Nrows] += ML_dabs(ext_aa[j]);
      }
      rowNorms[i+Nrows] /= extNrows;
      offset += recv_lengths[i];
   }

   /* ---------------------------------------------------------- */
   /* process the first Nrows                                    */
   /* ---------------------------------------------------------- */

   nnz_count = 0;
   mat_ia[0] = 0;
   track_array = (int *) ML_allocate( extNrows * sizeof(int) );
   for ( i = 0; i < extNrows; i++ ) dble_buf[i] = 0.0;

   for ( i = 0; i < Nrows; i++ )
   {
#ifdef ML_DEBUG_SMOOTHER
      if ( i % 1000 == 0 )
         printf("%4d : ILUT Processing row %6d (%6d)\n", mypid, i, extNrows);
#endif
      track_leng = 0;
      ML_get_matrix_row(Amat,1,&i,&allocated_space,&cols,&vals,&m,0);
      for ( j = 0; j < m; j++ ) 
      {
         if ( cols[j] < extNrows ) 
         {
            dble_buf[cols[j]] = vals[j];
            track_array[track_leng++] = cols[j];
         }
      }
      Lcount = Ucount = first = 0;
      first  = extNrows;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( dble_buf[index] != 0 )
         {
            if ( index < i ) Lcount++;
            else if ( index > i ) Ucount++;
            else if ( index == i ) diagonal[i] = dble_buf[index];
            if ( index < first ) first = index;
         }
      }
      Lcount = Lcount * fillin;
      Ucount = Ucount * fillin;
      rel_tau = tau * rowNorms[i];
      for ( j = first; j < i; j++ ) 
      {
         if ( ML_dabs(dble_buf[j]) > rel_tau )
         {
            ddata = dble_buf[j] / diagonal[j];
            for ( k = mat_ia[j]; k < mat_ia[j+1]; k++ ) 
            {
               colIndex = mat_ja[k];
               if ( colIndex > j ) 
               {
                  if ( dble_buf[colIndex] != 0.0 )
                     dble_buf[colIndex] -= (ddata * mat_aa[k]);
                  else
                  {
                     dble_buf[colIndex] = - (ddata * mat_aa[k]);
                     if ( dble_buf[colIndex] != 0.0 )
                        track_array[track_leng++] = colIndex;
                  }
               }
            }
            dble_buf[j] = ddata;
         }
         else dble_buf[j] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] < extNrows )
         {
            vals[j] = dble_buf[cols[j]];
            if ( cols[j] != i ) dble_buf[cols[j]] = 0.0;
         }
      }
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index < i )
         {
            absval = ML_dabs( dble_buf[index] );
            if ( absval > rel_tau )
            {
               sortcols[sortcnt] = index;
               sortvals[sortcnt++] = absval * rowNorms[index];
            }
            else dble_buf[index] = 0.0;
         }
      }
      if ( sortcnt > Lcount ) 
      {
         ML_split_dsort(sortvals, sortcnt, sortcols, Lcount); 
         for ( j = Lcount; j < sortcnt; j++ ) dble_buf[sortcols[j]] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] < i && vals[j] != 0.0 )
         {
            mat_aa[nnz_count] = vals[j];
            mat_ja[nnz_count++] = cols[j];
         }
      }
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index < i && dble_buf[index] != 0.0 )
         {
            mat_aa[nnz_count] = dble_buf[index];
            mat_ja[nnz_count++] = index;
            dble_buf[index] = 0.0;
         }
      }
      diagonal[i] = dble_buf[i];
      if ( ML_dabs(diagonal[i]) < 1.0e-16 ) diagonal[i] = 1.0E-6;
      mat_aa[nnz_count] = diagonal[i]; 
      mat_ja[nnz_count++] = i;
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i )
         {
            absval = ML_dabs(dble_buf[index]);
            if ( absval > rel_tau )
            {
               sortcols[sortcnt] = index;
               sortvals[sortcnt++] = absval * rowNorms[index];
            }
            else dble_buf[index] = 0.0;
         }
      }
      if ( sortcnt > Ucount ) 
      {
         ML_split_dsort(sortvals, sortcnt, sortcols, Ucount); 
         for ( j = Ucount; j < sortcnt; j++ ) dble_buf[sortcols[j]] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] > i && vals[j] != 0.0 )
         {
            mat_aa[nnz_count] = vals[j];
            mat_ja[nnz_count++] = cols[j];
         }
      }
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i && dble_buf[index] != 0.0 )
         {
            mat_aa[nnz_count] = dble_buf[index];
            mat_ja[nnz_count++] = index;
            dble_buf[index] = 0.0;
         }
      }
      dble_buf[i] = 0.0;
      mat_ia[i+1] = nnz_count;
   } 

   /* ---------------------------------------------------------- */
   /* process the off-processor rows                             */
   /* ---------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < total_recv_leng; i++ )
   {
#ifdef ML_DEBUG_SMOOTHER
      if ( (i+Nrows) % 1000 == 0 )
         printf("%4d : ILUT Processing row %6d (%6d)\n",mypid,i+Nrows,extNrows);
#endif
      track_leng = m = 0;

      for ( j = offset; j < offset+recv_lengths[i]; j++ ) 
      {
         if ( ext_ja[j] != -1 ) 
         {
            dble_buf[ext_ja[j]] = ext_aa[j];
            track_array[track_leng++] = ext_ja[j];
            cols[m] = ext_ja[j];
            vals[m++] = ext_aa[j];
         }
      }
      Lcount = Ucount = 0;
      first  = extNrows;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( dble_buf[index] != 0.0 )
         {
            if ( index < i+Nrows ) Lcount++;
            else if ( index > i+Nrows ) Ucount++;
            else if ( i+Nrows == index ) diagonal[i+Nrows] = dble_buf[index];
            if ( index < first ) first = index;
         }
      }
      Lcount = Lcount * fillin;
      Ucount = Ucount * fillin;
      rel_tau = tau * rowNorms[i+Nrows];
      for ( j = first; j < i+Nrows; j++ )
      {
         if ( ML_dabs(dble_buf[j]) > rel_tau )
         {
            ddata = dble_buf[j] / diagonal[j];
            for ( k = mat_ia[j]; k < mat_ia[j+1]; k++ ) 
            {
               colIndex = mat_ja[k];
               if ( colIndex > j ) 
               {
                  if ( dble_buf[colIndex] != 0.0 )
                     dble_buf[colIndex] -= (ddata * mat_aa[k]);
                  else
                  {
                     dble_buf[colIndex] = - (ddata * mat_aa[k]);
                     if ( dble_buf[colIndex] != 0.0 )
                        track_array[track_leng++] = colIndex;
                  }
               }
            }
            dble_buf[j] = ddata;
         }
         else dble_buf[j] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] < extNrows )
         {
            vals[j] = dble_buf[cols[j]];
            if ( cols[j] != i+Nrows ) dble_buf[cols[j]] = 0.0;
         }
      }
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index < i+Nrows )
         {
            absval = ML_dabs( dble_buf[index] );
            if ( absval > rel_tau )
            {
               sortcols[sortcnt] = index;
               sortvals[sortcnt++] = absval * rowNorms[index];
            }
            else dble_buf[index] = 0.0;
         }
      }
      if ( sortcnt > Lcount ) 
      {
         ML_split_dsort(sortvals, sortcnt, sortcols, Lcount); 
         for ( j = Lcount; j < sortcnt; j++ ) dble_buf[sortcols[j]] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] < i+Nrows && vals[j] != 0.0 )
         {
            mat_aa[nnz_count] = vals[j];
            mat_ja[nnz_count++] = cols[j];
         }
      }
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index < i+Nrows && dble_buf[index] != 0.0 )
         {
            mat_aa[nnz_count] = dble_buf[index];
            mat_ja[nnz_count++] = index;
            dble_buf[index] = 0.0;
         }
      }
      diagonal[i+Nrows] = dble_buf[i+Nrows];
      if ( ML_dabs(diagonal[i+Nrows]) < 1.0e-16 ) diagonal[i+Nrows] = 1.0E-6;
      mat_aa[nnz_count] = diagonal[i+Nrows]; 
      mat_ja[nnz_count++] = i+Nrows;
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i+Nrows )
         {
            absval = ML_dabs( dble_buf[index] );
            if ( absval > rel_tau )
            {
               sortcols[sortcnt] = index;
               sortvals[sortcnt++] = absval * rowNorms[index];
            }
            else dble_buf[index] = 0.0;
         }
      }
      if ( sortcnt > Ucount ) 
      {
         ML_split_dsort(sortvals, sortcnt, sortcols, Ucount); 
         for ( j = Ucount; j < sortcnt; j++ ) dble_buf[sortcols[j]] = 0.0;
      }
      for ( j = 0; j < m; j++ )
      {
         if ( cols[j] > i+Nrows && cols[j] < extNrows && vals[j] != 0.0 )
         {
            mat_aa[nnz_count] = vals[j];
            mat_ja[nnz_count++] = cols[j];
         }
      }
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i+Nrows && dble_buf[index] != 0.0 )
         {
            mat_aa[nnz_count] = dble_buf[index];
            mat_ja[nnz_count++] = index;
            dble_buf[index] = 0.0;
         }
      }
      dble_buf[i+Nrows] = 0.0;
      mat_ia[i+Nrows+1] = nnz_count;
      offset += recv_lengths[i];
   } 
   if ( nnz_count > total_nnz )
   {
      printf("ERROR in ILUTDecomp : memory out of bound (consult ML developers)\n");
      exit(1);
   }
#ifdef ML_DEBUG_SMOOTHER
   printf("%4d : ILUT Smoother - nnz(ILUT) = \n", mypid, nnz_count);
#endif

   /* ---------------------------------------------------------- */
   /* deallocate temporary storage space                         */
   /* ---------------------------------------------------------- */

   ML_free(cols);
   ML_free(vals);
   ML_free(dble_buf);
   ML_free(diagonal);
   ML_free(rowNorms);
   ML_free(sortcols);
   ML_free(sortvals);
   ML_free(track_array);
   return 0;
}

#ifdef out
/*****************************************************************************/
/* function for setting up variable block overlapped Schwarz                 */
/*****************************************************************************/

int ML_Smoother_VBlockSchwarzDecomposition(ML_Sm_Schwarz_Data *data, 
             ML_Operator *Amat, ML_Comm *comm, int total_recv_leng, 
             int *recv_lengths, int *ext_ja, double *ext_aa, int *map, 
             int *map2, int Noffset)
{
#ifdef SUPERLU
   int                i, j, k, **bmat_ia, **bmat_ja, allocated_space;
   int                *blk_size, index, **blk_indices, **aux_bmat_ia;
   int                offset, nnz, Nrows, extNrows, **aux_bmat_ja;
   int                mypid, *tmp_blk_leng, *cols, blknum, ncnt, *blkinfo;
   int                rownum, rowleng, nblocks, col_ind, init_size, aux_nnz;
   int                *tmp_indices, cur_off_row, max_blk_size, *mat_ia, *mat_ja;
   double             **bmat_aa, *vals, *mat_aa, **aux_bmat_aa;
   ML_Sm_Schwarz_Data *schwarz_ptr;
   int                info, *perm_r, *perm_c, permc_spec, *etree, panel_size;
   int                lwork, nrows;
   double             *R, *C, *ferr, *berr, rpg, rcond, *trhs, *tsol;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* ---------------------------------------------------------- */
   /* fetch Schwarz parameters                                   */
   /* ---------------------------------------------------------- */

   mypid       = comm->ML_mypid;
   schwarz_ptr = (ML_Sm_Schwarz_Data *) data;
   Nrows       = Amat->outvec_leng;
   extNrows    = Nrows + total_recv_leng;
   schwarz_ptr->Nrows = extNrows;
   blkinfo     = schwarz_ptr->blk_info;
   nblocks     = schwarz_ptr->nblocks;

   /* ---------------------------------------------------------- */
   /* adjust the off-processor row data                          */
   /* ---------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < total_recv_leng; i++ )
   {
      for ( j = offset; j < offset+recv_lengths[i]; j++ ) 
      {
         index = ext_ja[j];
         if ( index >= Noffset && index < Noffset+Nrows )
            ext_ja[j] = index - Noffset; 
         else
         {
            col_ind = ML_sorted_search(index, extNrows-Nrows, map); 
            if ( col_ind >= 0 ) ext_ja[j] = map2[col_ind] + Nrows;
            else                ext_ja[j] = -1;
         }
      }
      offset += recv_lengths[i];
   }

   /* ---------------------------------------------------------- */
   /* compose the initial blk_size information                   */
   /* ---------------------------------------------------------- */

   schwarz_ptr->blk_size = (int *) ML_allocate(nblocks * sizeof(int) );
   schwarz_ptr->blk_indices = (int **) ML_allocate(nblocks * sizeof(int*) );
   blk_indices  = schwarz_ptr->blk_indices;
   blk_size     = schwarz_ptr->blk_size;
   tmp_blk_leng = (int *) ML_allocate(nblocks * sizeof(int) );
   for ( i = 0; i < nblocks; i++ ) blk_size[i] = 0; 
   for ( i = 0; i < Nrows; i++ )   blk_size[blkinfo[i]]++;
   for ( i = 0; i < nblocks; i++ ) 
   {
      if ( blk_size[i] == 0 )
      {
         printf("%4d : SchwarzDecomposition - block %d is empty\n",mypid,i);
         exit(1);
      }
   }
   for ( i = 0; i < nblocks; i++ ) 
   {
      tmp_blk_leng[i] = blk_size[i] * blk_size[i] + 5;
      blk_indices[i] = (int *) ML_allocate(tmp_blk_leng[i] * sizeof(int));
   }
   for ( i = 0; i < nblocks; i++ ) blk_size[i] = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      blknum = blkinfo[i];
      index = blk_size[blknum]++;
      blk_indices[blknum][index] = i;
   }

   /* ---------------------------------------------------------- */
   /* now extend the each block for the overlap                  */
   /* (at the end blk_indices and bli_size contains the info)    */
   /* ---------------------------------------------------------- */

   allocated_space = extNrows;
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   cols = (int *)    ML_allocate(allocated_space * sizeof(int));
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
   {
      init_size = blk_size[i];
      for ( j = 0; j < init_size; j++ ) 
      {
         rownum = blk_indices[i][j];
         ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,&vals,&rowleng,0);
         if ( blk_size[i] + rowleng > tmp_blk_leng[i] )
         {
            tmp_indices = blk_indices[i];
            tmp_blk_leng[i] = 2 * ( blk_size[i] + rowleng ) + 2;
            blk_indices[i] = (int *) ML_allocate(tmp_blk_leng[i] * sizeof(int));
            for (k = 0; k < blk_size[i]; k++) blk_indices[i][k] = tmp_indices[k]; 
            ML_free( tmp_indices );
         }   
         for ( k = 0; k < rowleng; k++ ) 
         {
            col_ind = cols[k];
            blk_indices[i][blk_size[i]++] = col_ind;
         }
      } 
      ML_az_sort(blk_indices[i], blk_size[i], NULL, NULL);
      ncnt = 0;
      for ( j = 1; j < blk_size[i]; j++ ) 
         if ( blk_indices[i][j] != blk_indices[i][ncnt] )
           blk_indices[i][++ncnt] = blk_indices[i][j];
      blk_size[i] = ncnt + 1;
      if ( blk_size[i] > max_blk_size ) max_blk_size = blk_size[i];
   }

   /* ---------------------------------------------------------- */
   /* compute the memory requirements for each block             */
   /* ---------------------------------------------------------- */

   schwarz_ptr->bmat_ia = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->bmat_ja = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->bmat_aa = (double **) ML_allocate(nblocks * sizeof(double*) );
   bmat_ia = schwarz_ptr->bmat_ia;
   bmat_ja = schwarz_ptr->bmat_ja;
   bmat_aa = schwarz_ptr->bmat_aa;
   schwarz_ptr->aux_bmat_ia = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->aux_bmat_ja = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->aux_bmat_aa = (double **) ML_allocate(nblocks * sizeof(double*) );
   aux_bmat_ia = schwarz_ptr->aux_bmat_ia;
   aux_bmat_ja = schwarz_ptr->aux_bmat_ja;
   aux_bmat_aa = schwarz_ptr->aux_bmat_aa;

   for ( i = 0; i < nblocks; i++ ) 
   {
      nnz = aux_nnz = offset = cur_off_row = 0;
      for ( j = 0; j < blk_size[i]; j++ ) 
      {
         rownum = blk_indices[i][j];
         if ( rownum < Nrows )
            ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,
                              &vals,&rowleng,0);
         else 
         {
            for ( k = cur_off_row; k < rownum-Nrows; k++ ) 
               offset += recv_lengths[k]; 
            cur_off_row = rownum - Nrows;
            rowleng = 0;
            for ( k = offset; k < offset+recv_lengths[cur_off_row]; k++ ) 
               if ( ext_ja[k] != -1 ) cols[rowleng++] = ext_ja[k];
         }
         for ( k = 0; k < rowleng; k++ )
         {
            index = ML_find_index( cols[k], blk_indices[i], blk_size[i]);
            if ( index >= 0 ) nnz++;
            else              aux_nnz++;
         }
      }
      bmat_ia[i] = (int *)    ML_allocate( (blk_size[i] + 1) * sizeof(int));
      bmat_ja[i] = (int *)    ML_allocate( nnz * sizeof(int));
      bmat_aa[i] = (double *) ML_allocate( nnz * sizeof(double));
      aux_bmat_ia[i] = (int *)    ML_allocate( (blk_size[i] + 1) * sizeof(int));
      aux_bmat_ja[i] = (int *)    ML_allocate( aux_nnz * sizeof(int));
      aux_bmat_aa[i] = (double *) ML_allocate( aux_nnz * sizeof(double));
   }

   /* ---------------------------------------------------------- */
   /* load the submatrices                                       */
   /* ---------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ ) 
   {
      nnz = aux_nnz     = offset = cur_off_row = 0;
      bmat_ia[i][0]     = 0;
      aux_bmat_ia[i][0] = 0;

      for ( j = 0; j < blk_size[i]; j++ ) 
      {
         rownum = blk_indices[i][j];
         if ( rownum < Nrows )
            ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,
                              &vals,&rowleng,0);
         else 
         {
            for ( k = cur_off_row; k < rownum-Nrows; k++ ) 
            {
               offset += recv_lengths[k]; 
            }
            cur_off_row = rownum - Nrows;
            rowleng = 0;
            for ( k = offset; k < offset+recv_lengths[cur_off_row]; k++ ) 
            {
               if ( ext_ja[k] != -1 ) 
               {
                  cols[rowleng] = ext_ja[k];
                  vals[rowleng++] = ext_aa[k];
               }
            }
         }
         for ( k = 0; k < rowleng; k++ )
         {
            index = ML_find_index( cols[k], blk_indices[i], blk_size[i]);
            if ( index >= 0 )
            {
               bmat_ja[i][nnz] = index;
               bmat_aa[i][nnz++] = vals[k];
            }
            else
            {
               aux_bmat_ja[i][aux_nnz] = cols[k];
               aux_bmat_aa[i][aux_nnz++] = vals[k];
            }
         }
         bmat_ia[i][j+1] = nnz;
         aux_bmat_ia[i][j+1] = aux_nnz;
      } 
   } 

   /* ---------------------------------------------------------- */
   /* call SuperLU to perform decomposition                      */
   /* ---------------------------------------------------------- */

   schwarz_ptr->slu_Amat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->slu_Lmat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->slu_Umat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->perm_r   = (int **) ML_allocate(nblocks*sizeof(int*));
   schwarz_ptr->perm_c   = (int **) ML_allocate(nblocks*sizeof(int*));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   tsol  = (double *) ML_allocate( max_blk_size * sizeof(double) );
   trhs  = (double *) ML_allocate( max_blk_size * sizeof(double) );
   for ( i = 0; i < max_blk_size; i ++ ) trhs[i] = 0.0;

   for ( i = 0; i < nblocks; i ++ )
   {
      schwarz_ptr->slu_Amat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      schwarz_ptr->slu_Lmat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      schwarz_ptr->slu_Umat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      A = schwarz_ptr->slu_Amat[i];
      L = schwarz_ptr->slu_Lmat[i];
      U = schwarz_ptr->slu_Umat[i];
      nrows  = blk_size[i];
      mat_ia = schwarz_ptr->bmat_ia[i];
      mat_ja = schwarz_ptr->bmat_ja[i];
      mat_aa = schwarz_ptr->bmat_aa[i];
      nnz    = mat_ia[nrows];
      dCreate_CompRow_Matrix(A,nrows,nrows,nnz,mat_aa,mat_ja,mat_ia,NR,_D,GE);
      schwarz_ptr->perm_r[i] = (int *) ML_allocate(nrows*sizeof(int));
      schwarz_ptr->perm_c[i] = (int *) ML_allocate(2*nrows*sizeof(int));
      perm_r = schwarz_ptr->perm_r[i];
      perm_c = schwarz_ptr->perm_c[i];
      permc_spec = 0;
      get_perm_c(permc_spec, A, perm_c);
      panel_size               = sp_ienv(1);
      iparam.panel_size        = panel_size;
      iparam.relax             = sp_ienv(2);
      iparam.diag_pivot_thresh = 1.0;
      iparam.drop_tol          = -1;
      lwork                    = 0;
      *fact                    = 'N';
      *equed                   = 'N';
      *trans                   = 'N';
      *refact                  = 'N';
      dCreate_Dense_Matrix(&B, nrows, 1, trhs, nrows, DN, _D, GE);
      dCreate_Dense_Matrix(&X, nrows, 1, tsol, nrows, DN, _D, GE);

      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
#ifdef ML_DEBUG_SMOOTHER
printf("Block %6d(%6d) : cond no. = %e\n", i, nblocks, 1.0 / rcond);
if ( info != 0 && info != (nrows+1) )
{
   for ( j = 0; j < nrows; j++ )
   {
      for ( k = mat_ia[j]; k < mat_ia[j+1]; k++ )
         printf("Block %4d : in data %6d(%6d) %6d = %e\n",i,j,
                nrows,blk_indices[i][mat_ja[k]],mat_aa[k]);
      for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
         printf("Block %4d : offdata %6d(%6d) %6d = %e\n",i,j,nrows,
                aux_bmat_ja[i][k],aux_bmat_aa[i][k]);
   }
}
#endif
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   /* ---------------------------------------------------------- */
   /* clean up                                                   */
   /* ---------------------------------------------------------- */

   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   ML_free(etree);
   ML_free(trhs);
   ML_free(tsol);
   ML_free(vals);
   ML_free(cols);
   ML_free(tmp_blk_leng);
   return 0;
#else
   printf("ML_Smoother_VBlockSchwarzDecomposition : not available.\n");
   exit(1);
   return 1;
#endif
}
#endif
/*****************************************************************************/
/*****************************************************************************/
/* AZTEC related smoothers                                                   */
/*****************************************************************************/

/*****************************************************************************/
/* Symmetric MSR Gauss-Seidel smoother with damping.                         */
/* ------------------------------------------------------------------------- */

int ML_Smoother_MSR_SGSdamping(void *sm,int inlen,double x[],int outlen,
				 double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum, *ptr_val;
   int             bindx_row, *bindx;
   double          *ptr_b, *val = NULL/*, omega2*/;
   struct ML_CSR_MSRdata *ptr = NULL;
   double omega, one_minus_omega;

   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
   if (val == NULL) {
     ML_Smoother_SGS(sm, inlen, x, outlen, rhs);
     return 0;
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO)
      {
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
      }
   }
   else x2 = x;
   omega = smooth_ptr->omega;
   one_minus_omega = 1.0 - omega;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      bindx_row = bindx[0];
      bindx_ptr = &bindx[bindx_row];
      ptr_val   = &val[bindx_row];
      ptr_b     = rhs;

      for (i = 0; i < Nrows; i++)
      {
         sum    = *ptr_b++;
	  
	     for (j = bindx[i]; j < bindx[i+1]; j++)
            sum -= *ptr_val++ * x2[*bindx_ptr++];
         if (val[i] != 0.0)
            x2[i] = omega * sum/val[i] + one_minus_omega * x2[i];
      }

      bindx_ptr--;
      ptr_val--;
      ptr_b--;

      for (i = Nrows - 1; i >= 0; i--)
      {
         sum    = *ptr_b--;

         for (j = bindx[i]; j < bindx[i+1]; j++)
            sum -= *ptr_val-- * x2[*bindx_ptr--];
         if (val[i] != 0.0)
            x2[i] = omega * sum/val[i] + one_minus_omega * x2[i];
      }

   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}

/*****************************************************************************/
/* Symmetric MSR Gauss-Seidel smoother with no damping.                      */
/* ------------------------------------------------------------------------- */

int ML_Smoother_MSR_SGSnodamping(void *sm,int inlen,double x[],int outlen,
				 double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum;
   int             *bindx;
   double          *val = NULL;
   struct ML_CSR_MSRdata *ptr = NULL;

   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
   if (val == NULL) {
     ML_Smoother_SGS(sm, inlen, x, outlen, rhs);
     return 0;
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO)
      {
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
      }
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      j = bindx[0];
      bindx_ptr = &bindx[j];
      for (i = 0; i < Nrows; i++) {
	sum =  rhs[i];
	while (j+10 < bindx[i+1]) {
	  sum -= (val[j+9]*x2[bindx_ptr[9]] +
		  val[j+8]*x2[bindx_ptr[8]] +
		  val[j+7]*x2[bindx_ptr[7]] +
		  val[j+6]*x2[bindx_ptr[6]] +
		  val[j+5]*x2[bindx_ptr[5]] +
		  val[j+4]*x2[bindx_ptr[4]] +
		  val[j+3]*x2[bindx_ptr[3]] +
		  val[j+2]*x2[bindx_ptr[2]] +
		  val[j+1]*x2[bindx_ptr[1]] +
		  val[j]*x2[*bindx_ptr]);
	  bindx_ptr += 10;
	  j += 10;
	}
	while (j < bindx[i+1]) {
	  sum -= val[j++] * x2[*bindx_ptr++];
	}
	if (val[i] != 0.0) x2[i] = sum/val[i];
/*	if (fabs(val[i]) > 1.0e-9)  x2[i] = sum/val[i]; These is needed on 
        nodal problems with very small mass matrix..... */
      }

      j--;
      bindx_ptr--;

      for (i = Nrows - 1; i >= 0; i--) {
	sum    = rhs[i];

	while (j-9 >= bindx[i]) {
	  bindx_ptr -= 10;
	  j -= 10;
	  sum -= (val[j+10]*x2[bindx_ptr[10]] +
		  val[j+9]*x2[bindx_ptr[9]] +
		  val[j+8]*x2[bindx_ptr[8]] +
		  val[j+7]*x2[bindx_ptr[7]] +
		  val[j+6]*x2[bindx_ptr[6]] +
		  val[j+5]*x2[bindx_ptr[5]] +
		  val[j+4]*x2[bindx_ptr[4]] +
		  val[j+3]*x2[bindx_ptr[3]] +
		  val[j+2]*x2[bindx_ptr[2]] +
		  val[j+1]*x2[bindx_ptr[1]]);
	}
	while (j >= bindx[i]) {
	  sum -= val[j--] * x2[*bindx_ptr--];
	}
	if (val[i] != 0.0) x2[i] = sum/val[i];
/*	if (fabs(val[i]) > 1.0e-9)  x2[i] = sum/val[i]; These is needed on 
        nodal problems with very small mass matrix..... */
      }
   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}

/*****************************************************************************/
/* Symmetric MSR Gauss-Seidel smoother                                       */
/* ------------------------------------------------------------------------- */

int ML_Smoother_MSR_SGS(void *sm,int inlen,double x[],int outlen,double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum, *ptr_val;
   int             bindx_row, *bindx;
   double          *ptr_b, *val /*, omega2*/;
   struct ML_CSR_MSRdata *ptr = NULL;
   double *omega_val, **data, *one_minus_omega;

   smooth_ptr = (ML_Smoother *) sm;
   data = (double **) smooth_ptr->smoother->data;
   omega_val= data[0];
   one_minus_omega = data[1];
/*
   omega = smooth_ptr->omega;
   omega2 = 1. - omega;
*/
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
      if (inlen == -47) ML_avoid_unused_param((void *) &outlen);
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO)
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++) {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      bindx_row = bindx[0];
      bindx_ptr = &bindx[bindx_row];
      ptr_val   = &val[bindx_row];
      ptr_b     = rhs;

      for (i = 0; i < Nrows; i++) {
         sum    = *ptr_b++;

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= *ptr_val++ * x2[*bindx_ptr++];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
/*
         x2[i] = omega2*x2[i] + sum * omega_val[i];
*/
      }

      bindx_ptr--;
      ptr_val--;
      ptr_b--;

      for (i = Nrows - 1; i >= 0; i--) {
         sum    = *ptr_b--;

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= *ptr_val-- * x2[*bindx_ptr--];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
/*
         x2[i] = omega2*x2[i] + sum * omega_val[i];
*/
      }

   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}

/*****************************************************************************/
/* Symmetric MSR Gauss-Seidel smoother (different version)                   */
/* ------------------------------------------------------------------------- */

int ML_MSR_SGSextra(void *sm, int inlen, double x[], int outlen, double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum, *ptr_val;
   int             bindx_row, *bindx;
   double          *ptr_b, *val /*, omega2 */;
   struct ML_CSR_MSRdata *ptr = NULL;
   double *omega_val, **data, *one_minus_omega;
int Nextra, *extra, ii;

   smooth_ptr = (ML_Smoother *) sm;
   data = (double **) smooth_ptr->smoother->data;
   omega_val= data[0];
   one_minus_omega = data[1];
Nextra = (int) data[2][0];
extra  = (int *) data[3];
/*
   omega = smooth_ptr->omega;
   omega2 = 1. - omega;
*/
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
      if (inlen == -47) ML_avoid_unused_param((void *) &outlen);
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO) {
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
      }
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++) {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      bindx_row = bindx[0];
      bindx_ptr = &bindx[bindx_row];
      ptr_val   = &val[bindx_row];
      ptr_b     = rhs;

      for (i = 0; i < Nrows; i++) {
         sum    = *ptr_b++;

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= *ptr_val++ * x2[*bindx_ptr++];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
/*
         x2[i] = omega2*x2[i] + sum * omega_val[i];
*/
      }

      ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);
      for (ii = 0; ii < Nextra; ii++) {
         i    = extra[ii];
         sum  = rhs[i];

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= val[j] * x2[bindx[j]];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
      }
      for (ii = Nextra-1; ii >= 0; ii--) {
         i    = extra[ii];
         sum  = rhs[i];

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= val[j] * x2[bindx[j]];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
      }
      ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      bindx_ptr--;
      ptr_val--;
      ptr_b--;

      for (i = Nrows - 1; i >= 0; i--) {
         sum    = *ptr_b--;

         for (j = bindx[i]; j < bindx[i+1]; j++) {
            sum -= *ptr_val-- * x2[*bindx_ptr--];
         }
         x2[i] = one_minus_omega[i]*x2[i] + sum * omega_val[i];
/*
         x2[i] = omega2*x2[i] + sum * omega_val[i];
*/
      }

   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}

/*****************************************************************************/
/* destructor                                                                */
/* ------------------------------------------------------------------------- */

void ML_Smoother_Clean_MSR_GS(void *data)
{
   double **ptr;

   ptr = (double **) data;
   ML_free(ptr[0]);
   ML_free(ptr[1]);
   ML_free(ptr);
}
void ML_MSR_GSextra_Clean(void *data)
{
   double **ptr;

   ptr = (double **) data;
   ML_free(ptr[0]);
   ML_free(ptr[1]);
   ML_free(ptr[2]);
   ML_free(ptr);
}

/*****************************************************************************/
/* ask Ray to find out what this function does                               */
/* ------------------------------------------------------------------------- */

int ML_Smoother_BackGS(void *sm,int inlen,double x[],int outlen,double rhs[])
{
   int iter, i, j, length, allocated_space, *cols, col;
   double dtemp, diag_term, *vals, omega;
   ML_Operator *Amat;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum, *ptr_val;
   int             bindx_row, j_last, *bindx = NULL;
   double          *val, omega2;
   struct ML_CSR_MSRdata *ptr = NULL;

   smooth_ptr = (ML_Smoother *) sm;

   omega = smooth_ptr->omega;
   omega2 = 1. - omega;
   Amat = smooth_ptr->my_level->Amat;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_SGS): Need getrow() for SGS smoother\n");

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
      if (inlen == -47) ML_avoid_unused_param((void *) &outlen);
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_SGS(): Not enough space\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
				   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;


   if (bindx == NULL) {
      for (iter = 0; iter < smooth_ptr->ntimes; iter++) {
         for (i = Nrows-1; i >= 0; i--) {
            dtemp = 0.0;
            diag_term = 0.0;
            ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
            for (j = 0; j < length; j++) {
               col = cols[j];
               dtemp += vals[j]*x2[col];
               if (col == i) diag_term = vals[j];
            }
            if (diag_term != 0.0)
               x2[i] += omega*(rhs[i] - dtemp)/diag_term;
         }
     }
   }
   else {
     ptr_val = val;
     for (i = 0; i < Nrows; i++) {
        (*ptr_val) = omega/(*ptr_val);
        ptr_val++;
     }
     for (iter = 0; iter < smooth_ptr->ntimes; iter++) {
        bindx_row = bindx[Nrows];
        bindx_ptr = &bindx[bindx_row-1];
        ptr_val   = &val[bindx_row-1];

        for (i = Nrows - 1; i >= 0; i--) {
           sum = rhs[i];
           j_last  = bindx[i+1] - bindx[i];

           for (j = 0; j < j_last; j++) {
              sum -= *ptr_val-- * x2[*bindx_ptr--];
           }
           x2[i] = omega2*x2[i] + sum * val[i];
        }
     }
     for (i = 0; i < Nrows; i++) val[i] = omega / val[i];
   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+2) 
      Amat->max_nz_per_row = allocated_space;

   ML_free(vals); ML_free(cols);

   return 0;
}

/*****************************************************************************/
/* ask Ray to find out what this function does                               */
/* ------------------------------------------------------------------------- */

int ML_Smoother_OrderedSGS(void *sm,int inlen,double x[],int outlen,
                           double rhs[])
{
   int iter, i, ii, j, length, allocated_space, *cols, col;
   double dtemp, diag_term, *vals, omega;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum, *ptr_val;
   int             bindx_row, j_last, *bindx = NULL, *ordering;
   double          *val, omega2;
   struct ML_CSR_MSRdata *ptr = NULL;

   smooth_ptr = (ML_Smoother *) sm;

   omega = smooth_ptr->omega;
   omega2 = 1. - omega;
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   ordering = (int *)smooth_ptr->smoother->data;

   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_SGS): Need getrow() for SGS smoother\n");

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
      if (inlen == -47) ML_avoid_unused_param((void *) &outlen);
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_SGS(): Not enough space\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
				   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
   }
   else x2 = x;

   if (bindx == NULL) {
      for (iter = 0; iter < smooth_ptr->ntimes; iter++) {
         if (getrow_comm != NULL)
            ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);
         for (ii = 0; ii < Nrows; ii++) {
            i = ordering[ii];
            dtemp = 0.0;
            diag_term = 0.0;
            ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
            for (j = 0; j < length; j++) {
               col = cols[j];
               dtemp += vals[j]*x2[col];
               if (col == i) diag_term = vals[j];
            }
/*
            if (diag_term == 0.0) {
               pr_error("Error: SGS() can not be used with a zero diagonal\n");
            }
*/
            if (diag_term != 0.0)
               x2[i] += omega*(rhs[i] - dtemp)/diag_term;
         }

         for (ii = Nrows-1; ii >= 0; ii--) {
            i = ordering[ii];
            dtemp = 0.0;
            diag_term = 0.0;
            ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
            for (j = 0; j < length; j++) {
               col = cols[j];
               dtemp += vals[j]*x2[col];
               if (col == i) diag_term = vals[j];
            }
/*
            if (diag_term == 0.0) {
               pr_error("Error: GS() can not be used with a zero diagonal\n");
            }
*/
            if (diag_term != 0.0)
               x2[i] += omega*(rhs[i] - dtemp)/diag_term;
         }
     }
   }
   else {
     ptr_val = val;
     for (i = 0; i < Nrows; i++) {
        (*ptr_val) = omega/(*ptr_val);
        ptr_val++;
     }
     for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
     {
        if (getrow_comm != NULL) 
           ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);

        for (ii = 0; ii < Nrows; ii++) {
           i = ordering[ii];
           sum    = rhs[i];
           bindx_row = bindx[i];
           j_last = bindx[i+1] - bindx_row;

           ptr_val  = &(val[bindx_row]);
           bindx_ptr= &(bindx[bindx_row]);
           for (j = 0; j < j_last; j++) {
              sum -= *ptr_val++ * x2[*bindx_ptr++];
           }
           x2[i] = omega2*x2[i] + sum * val[i];
        }

        for (ii = Nrows - 1; ii >= 0; ii--) {
           i = ordering[ii];
           sum = rhs[i];
           bindx_row = bindx[i];
           j_last    = bindx[i+1] - bindx_row;
           ptr_val   = &(val[bindx_row]);
           bindx_ptr = &(bindx[bindx_row]);
           for (j = 0; j < j_last; j++) {
              sum -= *ptr_val++ * x2[*bindx_ptr++];
           }
           x2[i] = omega2*x2[i] + sum * val[i];
        }
     }
     for (i = 0; i < Nrows; i++) val[i] = omega / val[i];
   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+2) 
      Amat->max_nz_per_row = allocated_space;

   ML_free(vals); ML_free(cols);

   return 0;
}

/*****************************************************************************/
/* ask Ray to find out what this function does                               */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Gen_Ordering(ML_Operator *Amat, int **data_ptr)
{
   int Nrows;
   int i,j, count = 0;
   char *not_done, *colorme;
   int  *ordering, *cols;
   int  allocated_space, length;
   double *vals;

   Nrows   = Amat->getrow->Nrows;
   allocated_space = Amat->max_nz_per_row+28;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) 
      pr_error("Error in Smoother_Gen_Ordering: Not enough space\n");

   not_done= (char *) ML_allocate(Nrows*sizeof(char));
   colorme = (char *) ML_allocate(Nrows*sizeof(char));
   ordering=  (int *) ML_allocate(Nrows*sizeof(int));
   if (ordering == NULL) pr_error("Out of spacing in Smoother_Gen_Order\n");

   for (i = 0; i < Nrows; i++) colorme[i] = 'y';
   for (i = 0; i < Nrows; i++) not_done[i] = 'y';

   while (count != Nrows) {
      for (i = 0; i < Nrows; i++) {
         if ( colorme[i] == 'y') {
            ordering[count++] = i;
            colorme[i] = 'n';
            not_done[i] = 'n';
            ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                              &length, 0);
            for ( j = 0; j < length; j++)
               if ( cols[j] < Nrows) colorme[cols[j]] = 'n';
         }
      }
      for (i = 0; i < Nrows; i++)
         colorme[i] = not_done[i];
   }
   ML_free(colorme);
   ML_free(not_done);
   ML_free(vals);
   ML_free(cols);
   *data_ptr = (int *) ordering;
   return(1);
}

/*****************************************************************************/
/* ask Ray to find out what this function does                               */
/* ------------------------------------------------------------------------- */

void ML_Smoother_Clean_OrderedSGS(void *data)
{
   ML_free(data);
}

#ifdef	MB_MODIF

/*****************************************************************************
 *
 * The MLS smoother. Implementation by Marian Brezina
 * The following section defines the functions relevant to the 
 * MLS-type smoother based on the Chebyshev polynomial of degree deg.
 * In the standard multigrid  context, value deg=1 is most relevant.
 *
 *****************************************************************************/ 


int ML_MLS_SandwPres(void *sm, int inlen, double x[], int outlen, double y[])
{
   /***************************************************************************	
    *
    * Apply the Chebyshev sandwich to x, return result in y.
    * Vectors x and y may not be aliased !
    * Vector  x is destroyed in the process, its output value is undefined.
    *
    ***************************************************************************/

    ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
    ML_Operator     *Amat = smooth_ptr->my_level->Amat;
    struct MLSthing *widget;
    int              i, deg, dg, n = outlen;
    double           *omV, om;      

    widget = (struct MLSthing *) smooth_ptr->smoother->data;

    deg    = widget->mlsDeg;
    omV    = widget->mlsOm;

    if (inlen != outlen) { 
	    pr_error("ML_MLS_Sandw: mtx. must be square\n");
    }

    for (i=0; i<n; i++) y[i] = x[i]; 

    for (dg=deg-1; dg>-1; dg--) { 

        ML_Operator_Apply(Amat, n, y, n, x);
        om = omV[dg]; 
        for (i=0; i<n; i++) y[i] -= om * x[i]; 

    }

    return 0; 
}

int ML_MLS_SandwPost(void *sm, int inlen, double x[], int outlen, double y[])
{
   /***************************************************************************	
    *
    * Apply the Chebyshev sandwich to x, return result in y.
    * Vectors x and y may not be aliased !
    * The contents of vector x are destroyed in the process, its output 
    * value is undefined.
    *
    ***************************************************************************/

    ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
    ML_Operator     *Amat = smooth_ptr->my_level->Amat;
    struct MLSthing *widget;
    int              i, deg, dg, n = outlen;
    double           *omV, om;      

    widget = (struct MLSthing *) smooth_ptr->smoother->data;

    deg    = widget->mlsDeg;
    omV    = widget->mlsOm;

    if (inlen != outlen) { 
	pr_error("ML_MLS_SandwPost: mtx. must be square\n");
    }

    for (i =0; i<n; i++) y[i] = x[i]; 

    for (dg=0; dg<deg; dg++) { 

        ML_Operator_Apply(Amat, n, y, n, x);
        om = omV[dg]; 
        for (i=0; i<n; i++) y[i] -= om * x[i]; 

    }

    return 0; 
}

int ML_MLS_SPrime_Apply(void *sm,int inlen,double x[],int outlen, double rhs[])
{
   /***************************************************************************	
    *
    * Approximate the solution of Ax=b using one itaration of the 
    * S_prime-based relaxation:   x <- S_prime(x, b)
    *
    ***************************************************************************/

    ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
    ML_Operator     *Amat = smooth_ptr->my_level->Amat;
    struct MLSthing *widget;
    int              i, n = outlen;
    double           cf, om2, over;
    double          *pAux, *y;      
#ifdef RST_MODIF
    double          *diagonal;
#endif

    widget = (struct MLSthing *) smooth_ptr->smoother->data;

    om2    = widget->mlsOm2;
    over   = widget->mlsOver;
    cf     = over * om2; 

    if (inlen != outlen) { 
	pr_error("ML_MLS_SPrime_Apply: mtx. must be square\n");
    }

    pAux   = (double *) ML_allocate(n*sizeof(double));
    y      = (double *) ML_allocate(n*sizeof(double));

    if (pAux == NULL) pr_error("ML_MLS_SPrime_Apply: allocation failed\n");
    if (y    == NULL) pr_error("ML_MLS_SPrime_Apply: allocation failed\n");

    ML_Operator_Apply(Amat, n, x, n, pAux);
#ifdef RST_MODIF
    ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);
    for (i=0; i<n; i++) pAux[i] -= rhs[i]/diagonal[i];
#else
    for (i=0; i<n; i++) pAux[i] -= rhs[i];
#endif
    ML_MLS_SandwPost(sm, n, pAux, n, y);
    ML_MLS_SandwPres(sm, n, y, n, pAux);
    for (i=0; i<n; i++) x[i] = x[i] - cf * pAux[i]; 

    ML_free(y);
    ML_free(pAux);

    return 0; 
}
int DinvA(void *data,  int in, double p[], int out, double ap[])
{
  void *olddata;
  struct DinvA_widget *DinvA_widget;
  ML_Operator *Amat;
  double *diagonal;
  int i;

  DinvA_widget = (struct DinvA_widget *) data;
  Amat = DinvA_widget->Amat;
  olddata = Amat->data;

  Amat->matvec->ML_id    = DinvA_widget->ML_id;
  Amat->matvec->internal = DinvA_widget->internal;
  Amat->matvec->external = DinvA_widget->external;
  Amat->data             = DinvA_widget->data;
  
  ML_Operator_Apply(Amat, in, p, out, ap);
  ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);
  for (i = 0; i < Amat->outvec_leng; i++) ap[i] = ap[i]/diagonal[i];

  Amat->matvec->ML_id    = ML_EXTERNAL;
  Amat->matvec->external = DinvA;
  Amat->data             = olddata;
  return 0;
}

int ML_Smoother_MLS_Apply(void *sm,int inlen,double x[],int outlen,
                          double rhs[])
{ /****************************************************************************	
   *
   *   Apply MLS smoother with the preselected degree deg:   x <- MLS(x, rhs)
   *
   ****************************************************************************/

   ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator     *Amat = smooth_ptr->my_level->Amat;
   struct MLSthing *widget;
   int              deg, dg, n, nn, i;
   double          *pAux, *res, *y, cf, over, *mlsCf;

#ifdef RST_MODIF
   int             *cols, allocated_space;
   double          *diagonal, *vals, *tdiag;
   struct DinvA_widget DinvA_widget;
   int j;
#endif
   n = outlen;
   widget = (struct MLSthing *) smooth_ptr->smoother->data;
   if (outlen == -47) ML_avoid_unused_param((void *) &inlen);

   deg    = widget->mlsDeg;
   mlsCf  = widget->mlsCf;
   over   = widget->mlsOver;

#define MB_FORNOW
#ifdef 	MB_FORNOW
   pAux  = (double *) ML_allocate((n+1)*sizeof(double));
   res   = (double *) ML_allocate((n+1)*sizeof(double));
   y     = (double *) ML_allocate((n+1)*sizeof(double));

   if (pAux == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
   if (res  == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
   if (y    == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
#endif

#ifdef RST_MODIF

   /* ----------------------------------------------------------------- */
   /* extract diagonal using getrow function if not found               */
   /* ----------------------------------------------------------------- */

   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(MLS_Apply): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&nn) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols); 
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < nn; j++) 
               if (cols[j] == i) tdiag[i] = vals[j];
	    if (tdiag[i] == 0.) tdiag[i] = 1.;
         }
         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);

   DinvA_widget.ML_id = Amat->matvec->ML_id;
   DinvA_widget.internal = Amat->matvec->internal;
   DinvA_widget.external = Amat->matvec->external;
   DinvA_widget.data     = Amat->data;
   DinvA_widget.Amat     = Amat;
   Amat->matvec->ML_id    = ML_EXTERNAL;
   Amat->matvec->external = DinvA;
   Amat->data             = &DinvA_widget;
#endif

   if (deg != 1) {
   if (smooth_ptr->init_guess == ML_NONZERO) {
     ML_Operator_Apply(Amat, n, x, n, pAux);
     for (i = 0; i < n; i++) pAux[i] = (rhs[i]
#ifdef RST_MODIF
			                /diagonal[i]
#endif
                                        - pAux[i]); 
   }
   else { 
     for (i = 0; i < n; i++) pAux[i] = rhs[i]
#ifdef RST_MODIF
			                /diagonal[i]
#endif
                                                    ;
   }
   }


   if (deg == 1) { 

       cf = over * mlsCf[0]; 

       /* for (i=0; i<n; i++) x[i] += cf * pAux[i];  */
       if (smooth_ptr->init_guess == ML_NONZERO) {
	 ML_Operator_Apply(Amat, n, x, n, pAux);
	 for (i = 0; i < n; i++) x[i] += cf * (rhs[i]
#ifdef RST_MODIF
			                /diagonal[i]
#endif
                                        - pAux[i]); 
       }
       else { 
	 for (i = 0; i < n; i++) x[i] = cf*rhs[i]
#ifdef RST_MODIF
			                /diagonal[i]
#endif
                                                    ;
       }
#ifdef 	MB_FORNOW
       if (y)    { ML_free(   y);    y = NULL; } 
       if (res)  { ML_free( res);  res = NULL; } 
       if (pAux) { ML_free(pAux); pAux = NULL; }
#endif
      /*
       * Apply the S_prime operator 
       */
       ML_MLS_SPrime_Apply(sm, n, x, n, rhs);
#ifdef RST_MODIF
       Amat->matvec->ML_id    = DinvA_widget.ML_id;
       Amat->matvec->internal = DinvA_widget.internal;
       Amat->matvec->external = DinvA_widget.external;
       Amat->data             = DinvA_widget.data;
#endif
       return 0;

   } else { 

       for (i=0;   i < n; i++) y[i]  = mlsCf[0] * pAux[i];
       for (dg=1; dg < deg; dg++) { 
            ML_Operator_Apply(Amat, n, pAux, n, res);
	    for (i=0; i < n; i++) pAux[i] = res[i];  
	    for (i=0; i < n; i++) y[i] += mlsCf[dg] * res[i];
       }

   }

   for (i=0; i < n; i++) x[i] += over * y[i];

#ifdef 	MB_FORNOW
   /* @@@ clean up later in destructor, right now clean in here .... */
   /* @@@ must deallocate here if we decide to stick with local allocations*/
   if (y)    { ML_free(   y);    y = NULL; } 
   if (res)  { ML_free( res);  res = NULL; } 
   if (pAux) { ML_free(pAux); pAux = NULL; }
#endif
  /*
   * Apply the S_prime operator 
   */
   ML_MLS_SPrime_Apply(sm, n, x, n, rhs);
#ifdef RST_MODIF
   Amat->matvec->ML_id    = DinvA_widget.ML_id;
   Amat->matvec->internal = DinvA_widget.internal;
   Amat->matvec->external = DinvA_widget.external;
   Amat->data             = DinvA_widget.data;
#endif
   return 0;
}
struct MLSthing *ML_Smoother_Create_MLS(void)
{
  struct MLSthing *widget;

  widget = (struct MLSthing *) ML_allocate(sizeof(struct MLSthing));
  widget->block_scaling   = NULL;
  widget->unscaled_matrix = NULL;
  widget->scaled_matrix   = NULL;
  widget->mlsBoost = 1.0;
  widget->mlsOver  = 1.1e0 ;
  widget->pAux     = NULL;   /* currently reserved */
  widget->res      = NULL;   /* currently reserved */
  widget->y        = NULL;   /* currently reserved */
  return(widget);

}

void ML_Smoother_Destroy_MLS(void *data)
{
   struct MLSthing *widget;

   widget = (struct  MLSthing *) data;

   if (widget->y)    ML_free(widget->y   );
   if (widget->res)  ML_free(widget->res );
   if (widget->pAux) ML_free(widget->pAux);

   if (widget->block_scaling) ML_Smoother_Destroy_BGS_Data((void *)
                                                       widget->block_scaling);

   if(widget->scaled_matrix) ML_Operator_Destroy(&widget->scaled_matrix);
   /* widget->unscaled_matrix was not create by   */
   /* MLS it is someone else's job to destroy it  */

   ML_free(widget);
}
/* ***************************************************************************
 *
 * end of the MLS-specific code.
 *
 * ***************************************************************************/ 

#endif/*MB_MODIF*/


/* ************************************************************************* */
/* Sparse approximate inverse smoother                                       */
/* ------------------------------------------------------------------------- */

#ifdef PARASAILS
#include "Matrix.h"
#include "ParaSails.h"

int ML_Smoother_ParaSails(void *sm,int inlen,double x[],int outlen,
                        double rhs[])
{
   ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator    *Amat = smooth_ptr->my_level->Amat;
   int            n = outlen, i;
   double         *res, *temp;
   struct widget { int parasails_factorized; ParaSails *ps;} *tptr;
   int            parasails_factorized;

   ParaSails      *ps;
	 
   tptr = (struct widget *) smooth_ptr->smoother->data;
   parasails_factorized = tptr->parasails_factorized;
   ps = tptr->ps;

   temp = (double *) ML_allocate(n*sizeof(double));
   res  = (double *) ML_allocate(n*sizeof(double));
   if (res == NULL) pr_error("ML_Smoother_ParaSails: out of space\n");

   ML_Operator_Apply(Amat, n, x, n, res);
   for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];

   if (!parasails_factorized)
   {
      MatrixMatvec(ps->M, res, temp);
      for (i = 0; i < n; i++) x[i] += temp[i];
   }
   else
   {
      MatrixMatvec(ps->M, res, temp);
      MatrixMatvecTrans(ps->M, temp, temp);
      for (i = 0; i < n; i++) x[i] += temp[i];
   }

   ML_free(temp);
   ML_free(res);

   (void) inlen;
   (void) outlen;
   return 0;
}

int ML_Smoother_ParaSailsTrans(void *sm,int inlen,double x[],int outlen,
                        double rhs[])
{
   ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator    *Amat = smooth_ptr->my_level->Amat;
   int            n = outlen, i;
   double         *res, *temp;
   struct widget { int parasails_factorized; ParaSails *ps;} *tptr;
   int            parasails_factorized;

   ParaSails      *ps;
	 
   tptr = (struct widget *) smooth_ptr->smoother->data;
   parasails_factorized = tptr->parasails_factorized;
   ps = tptr->ps;

   temp = (double *) ML_allocate(n*sizeof(double));
   res  = (double *) ML_allocate(n*sizeof(double));
   if (res == NULL) pr_error("ML_Smoother_ParaSails: out of space\n");

   ML_Operator_Apply(Amat, n, x, n, res);
   for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];

   if (!parasails_factorized)
   {
      MatrixMatvecTrans(ps->M, res, temp);
      for (i = 0; i < n; i++) x[i] += temp[i];
   }
   else
   {
      MatrixMatvec(ps->M, res, temp);
      MatrixMatvecTrans(ps->M, temp, temp);
      for (i = 0; i < n; i++) x[i] += temp[i];
   }

   ML_free(temp);
   ML_free(res);

   (void) inlen;
   (void) outlen;
   return 0;
}

int ML_Smoother_ParaSailsSym(void *sm,int inlen,double x[],int outlen,
                        double rhs[])
{
   ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator    *Amat = smooth_ptr->my_level->Amat;
   int            n = outlen, i;
   double         *res, *temp, *temp2;
   struct widget { int parasails_factorized; ParaSails *ps;} *tptr;
   int            parasails_factorized;

   ParaSails      *ps;
	 
   tptr = (struct widget *) smooth_ptr->smoother->data;
   parasails_factorized = tptr->parasails_factorized;
   ps = tptr->ps;

   temp = (double *) ML_allocate(n*sizeof(double));
   temp2= (double *) ML_allocate(n*sizeof(double));
   res  = (double *) ML_allocate(n*sizeof(double));
   if (res == NULL) pr_error("ML_Smoother_ParaSails: out of space\n");

   ML_Operator_Apply(Amat, n, x, n, res);
   for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];

   if (!parasails_factorized)
   {
      MatrixMatvec(ps->M, res, temp);
      MatrixMatvecTrans(ps->M, res, temp2);
      for (i = 0; i < n; i++) x[i] += 0.5 * (temp[i] + temp2[i]);
   }
   else
   {
      MatrixMatvec(ps->M, res, temp);
      MatrixMatvecTrans(ps->M, temp, temp);
      for (i = 0; i < n; i++) x[i] += temp[i];
   }

   ML_free(temp2);
   ML_free(temp);
   ML_free(res);

   (void) inlen;
   (void) outlen;
   return 0;
}

void ML_Smoother_Clean_ParaSails(void *data)
{
   struct widget { int parasails_factorized; ParaSails *ps;} *tptr;

   tptr = (struct widget *) data;
   ParaSailsDestroy(tptr->ps);
   ML_free(tptr);
}
#endif

/*******************************************************************************
 Reinitialize the smoother.  This is useful when just the operators have
 changed in the multigrid hierarchy (the interpolation and projection are
 the same) and the smoother needs to be regenerated on each level.

*******************************************************************************/

int ML_Smoother_Reinit(ML *ml )
{
    int             i;
    char            str[80];

    for (i = 0; i < ml->ML_num_levels; i++) {
      ML_Smoother_Clean( &(ml->pre_smoother[i]));
      ML_Smoother_Clean( &(ml->post_smoother[i]));
      ML_CSolve_Clean(&(ml->csolve[i]));

      ML_Smoother_Init( &(ml->pre_smoother[i]), &(ml->SingleLevel[i]) );
      ML_Smoother_Init( &(ml->post_smoother[i]), &(ml->SingleLevel[i]) );
      ML_CSolve_Init( &(ml->csolve[i]) );
      ML_CSolve_Set_1Level( &(ml->csolve[i]), &(ml->SingleLevel[i]) );

      sprintf(str,"PreS_%d",i);
      ML_Smoother_Set_Label( &(ml->pre_smoother[i]),str);
      sprintf(str,"PostS_%d",i);
      ML_Smoother_Set_Label( &(ml->post_smoother[i]),str);
      sprintf(str,"Solve_%d",i);
      ML_CSolve_Set_Label(&(ml->csolve[i]),str);
    }
    return 0;
}

ML *ML_Smoother_Get_Hiptmair_nodal(ML *ml, int level, int pre_or_post)
{
  ML_Sm_Hiptmair_Data *dataptr;

  if (level == ML_ALL_LEVELS) {
    printf("ML_Smoother_Get_Hiptmair_nodal: ML_ALL_LEVELS not supported.\n");
    ML_avoid_unused_param((void *) &pre_or_post);
    exit(1);
  }
  if (ml == NULL) {
    printf("ML_Smoother_Get_Hiptmair_nodal: ml parameter null?\n");
    exit(1);
  }
  if (level >= ml->ML_num_levels) {
    printf("ML_Smoother_Get_Hiptmair_nodal:Cannot set smoother on level %d\n",
	   level);
    printf("                               Only %d levels in the hierarchy.\n",
	   ml->ML_num_levels);
    exit(1);
  }
  if (ml->pre_smoother[level].smoother->internal != ML_Smoother_Hiptmair) {
    printf("ML_Smoother_Get_Hiptmair_nodal:Hiptmair not set on level %d\n",
	   level);
    exit(1);
  }
  
   /* pointer to private smoother data */

   dataptr = (ML_Sm_Hiptmair_Data *) ml->pre_smoother[level].smoother->data;
   if (dataptr->ml_nodal == NULL) {
     printf("ML_Smoother_Get_Hiptmair_nodal: Something wrong on level %d\n",
	   level);
     exit(1);
   }
   return(dataptr->ml_nodal);

}
void **ML_Smoother_Arglist_Create(int nargs)
{
  void **arglist;
  int *itmp, i;

  arglist = (void **) ML_allocate( (4+nargs)*sizeof(void *));
  for (i = 0; i < 4 + nargs; i++) arglist[i] = NULL;
  itmp = (int *) arglist;
  itmp[0] = ML_Set;
  itmp[1] = nargs;

  return arglist;
}
int ML_Smoother_Arglist_Delete(void ***arglist)
{
  int *itmp;
  if (*arglist == NULL) {
    /*
    printf("ML_Smoother_Arglist_Delete: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
    */
    return 0;
  }
  itmp = (int *) *arglist;
  if (itmp[0] != ML_Set) {
    printf("ML_Smoother_Arglist_Delete: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  ML_free(*arglist);
  *arglist = NULL;

  return 0;
}
int ML_Smoother_Arglist_Nargs(void **arglist)
{
  int *itmp;
  if (arglist == NULL) {
    printf("ML_Smoother_Arglist_Nargs: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  itmp = (int *) arglist;
  if (itmp[0] != ML_Set) {
    printf("ML_Smoother_Arglist_Nargs: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  return itmp[1];

}
void *ML_Smoother_Arglist_Get(void **arglist, int which_arg)
{
  int *itmp;
  if (arglist == NULL) {
    printf("ML_Smoother_Arglist_Get: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  itmp = (int *) arglist;
  if (itmp[0] != ML_Set) {
    printf("ML_Smoother_Arglist_Get: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  if (which_arg >= itmp[1]) {
    printf("ML_Smoother_Arglist_Get: argument #%d exceeds arglist length of",which_arg);
    printf(" %d allocated \n                         via ML_Smoother_Arglist_Create.",itmp[1]);
    printf(" Please remember that\n                         arguments are numbered starting from 0.\n");
    exit(1);
  }
  if (arglist[which_arg+4] == NULL) {
    printf("ML_Smoother_Arglist_Get: argument #%d has not been set.",which_arg);
    printf(" Please remember that\n                         arguments are numbered starting from 0.\n");
    exit(1);
  }

  return arglist[which_arg+4];
}
int ML_Smoother_Arglist_Set(void **arglist, int which_arg, void *ptr)
{
  int *itmp;

  if (arglist == NULL) {
    printf("ML_Smoother_Arglist_Set: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  itmp = (int *) arglist;
  if (itmp[0] != ML_Set) {
    printf("ML_Smoother_Arglist_Set: arglist not allocated via ML_Smoother_Arglist_Create\n");
    exit(1);
  }
  if (which_arg >= itmp[1]) {
    printf("ML_Smoother_Arglist_Set: argument #%d exceeds the arglist length of",which_arg);
    printf(" %d allocated \n                         via ML_Smoother_Arglist_Create.",itmp[1]);
    printf(" Please remember that\n                         arguments are numbered starting from 0.\n");
    exit(1);
  }
  arglist[which_arg+4] = ptr;

  return 0;
}

int ML_Smoother_HiptmairSubsmoother_Create(ML **ml_subproblem,
					   ML_Operator *Amat, void *smoother, 
					   void **args, double default_omega)
{
  double omega, *dbl_arg1;
  int *int_arg1, *int_arg2, *int_arg3;
#ifdef BLOCKMLS
  int nblocks, *block_list;
#endif

   ML_Create(ml_subproblem,1);
   ML_Operator_halfClone_Init( &((*ml_subproblem)->Amat[0]),
				   Amat);

   if (smoother == (void *) ML_Gen_Smoother_MLS) {


     if (ML_Smoother_Arglist_Nargs(args) != 2) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need 2 arguments for ML_Gen_Smoother_MLS() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1); /* eig ratio     */
     int_arg2 = (int *) ML_Smoother_Arglist_Get(args, 0);    /* poly degree   */
     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
       printf("Generating subsmoother MLS   %d\n",*int_arg2);
#ifdef BLOCKMLS
     printf("size is %d\n",(*ml_subproblem)->Amat[0].invec_leng);
     if (((*ml_subproblem)->Amat[0].invec_leng != 13322) &&
	 ((*ml_subproblem)->Amat[0].invec_leng != 655) &&
	 ((*ml_subproblem)->Amat[0].invec_leng != 4) ) {
       printf("CHANGING TO BLOCK CHEBY\n");
       nblocks = (int) (((double) (*ml_subproblem)->Amat[0].invec_leng)/50.);
       nblocks++;
       ML_Gen_Blocks_Metis(*ml_subproblem, 0, &nblocks, &block_list);
       ML_Gen_Smoother_BlockDiagScaledCheby(*ml_subproblem, 0, ML_PRESMOOTHER, 
					  *dbl_arg1,*int_arg2,
					  nblocks, block_list); 
     }
     else
       ML_Gen_Smoother_MLS(*ml_subproblem, 0, ML_PRESMOOTHER,*dbl_arg1,
			   *int_arg2*10);
#else
       ML_Gen_Smoother_MLS(*ml_subproblem, 0, ML_PRESMOOTHER,*dbl_arg1,
			   *int_arg2);
#endif

   }
   else if (smoother == (void *) ML_Gen_Smoother_ERF_1StepKrylov) {
     if (ML_Smoother_Arglist_Nargs(args) != 0) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need 0 arguments for ML_Gen_Smoother_ERG1StepKrylov() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
       printf("Generating subsmoother ERF_1StepKrylov\n");
     ML_Gen_Smoother_ERF_1StepKrylov(*ml_subproblem, 0, ML_PRESMOOTHER);
   }


   else if (smoother == (void *) ML_Gen_Smoother_Jacobi) {
     if (ML_Smoother_Arglist_Nargs(args) != 2) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need two arguments for ML_Gen_Smoother_Jacobi() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     int_arg1 = (int *) ML_Smoother_Arglist_Get(args, 0);
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1);
     omega = dbl_arg1[0];
     if (default_omega == 1.0) default_omega = .5;
     if ( ((int) omega) == ML_DEFAULT) omega= default_omega;

     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
        printf("Generating subsmoother Jacobi\n");
     ML_Gen_Smoother_Jacobi(*ml_subproblem, 0, ML_PRESMOOTHER,*int_arg1,
			    omega);
   }
   else if (smoother == (void *) ML_Gen_Smoother_GaussSeidel) {
     printf("Entering ML_Smoother_Gen_Hiptmair_Data (GS)\n");
     if (ML_Smoother_Arglist_Nargs(args) != 2) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need two arguments for ML_Gen_Smoother_GaussSeidel() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     int_arg1 = (int *) ML_Smoother_Arglist_Get(args, 0);
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1);
     omega = dbl_arg1[0];
     if ( ((int) omega) == ML_DEFAULT) omega= default_omega;

     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
        printf("Generating subsmoother Gauss Seidel\n");
     ML_Gen_Smoother_GaussSeidel(*ml_subproblem, 0, ML_PRESMOOTHER,*int_arg1,
			    omega);
   }
   else if (smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
     if (ML_Smoother_Arglist_Nargs(args) != 2) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need two arguments for ML_Gen_Smoother_SymGaussSeidel() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     int_arg1 = (int *) ML_Smoother_Arglist_Get(args, 0);
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1);
     omega = dbl_arg1[0];
     if ( ((int) omega) == ML_DEFAULT) omega= default_omega;

     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
       printf("Generating subsmoother symmetric Gauss Seidel\n");
     ML_Gen_Smoother_SymGaussSeidel(*ml_subproblem, 0, ML_PRESMOOTHER,*int_arg1,
			    omega);
   }
   else if (smoother == (void *) ML_Gen_Smoother_VBlockJacobi) {
     if (ML_Smoother_Arglist_Nargs(args) != 4) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need 4 arguments for ML_Gen_Smoother_VBlockJacobi() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     int_arg1 = (int *) ML_Smoother_Arglist_Get(args, 0);
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1);
     omega = dbl_arg1[0];
     if (default_omega == 1.0) default_omega = .5;
     if ( ((int) omega) == ML_DEFAULT) omega= default_omega;
     int_arg2 = (int *) ML_Smoother_Arglist_Get(args, 2);
     int_arg3 = (int *) ML_Smoother_Arglist_Get(args, 3);
     
     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
        printf("Generating subsmoother variable block Jacobi\n");
     ML_Gen_Smoother_VBlockJacobi(*ml_subproblem, 0, ML_PRESMOOTHER,
				  *int_arg1, omega, *int_arg2,int_arg3);
   }
   else if (smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel) {
     if (ML_Smoother_Arglist_Nargs(args) != 4) {
       printf("ML_Smoother_Gen_Hiptmair_Data: Need 4 arguments for ML_Gen_Smoother_VBlockSymGaussSeidel() got %d arguments\n", ML_Smoother_Arglist_Nargs(args));
       exit(1);
     }
     int_arg1 = (int *) ML_Smoother_Arglist_Get(args, 0);
     dbl_arg1 = (double *) ML_Smoother_Arglist_Get(args, 1);
     omega = dbl_arg1[0];
     if ( ((int) omega) == ML_DEFAULT) omega= default_omega;
     int_arg2 = (int *) ML_Smoother_Arglist_Get(args, 2);
     int_arg3 = (int *) ML_Smoother_Arglist_Get(args, 3);
     
     if (Amat->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() )
       printf("Generating subsmoother variable block symmetric Gauss Seidel\n");
     ML_Gen_Smoother_VBlockSymGaussSeidel(*ml_subproblem, 0, ML_PRESMOOTHER,
				  *int_arg1, omega, *int_arg2,int_arg3);
   }
   else {
   printf("ML_Smoother_Gen_Hiptmair_Data: Unknown smoother for Hiptmair subproblem\n");
     exit(1);
   }
   return 0;
}

 
 
int ML_Cheby(void *sm, int inlen, double x[], int outlen, double rhs[])
{

   ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator     *Amat = smooth_ptr->my_level->Amat;
   struct MLSthing *widget;
   int              deg, i, j, k, n, nn;
   double          *pAux, *dk;
   double beta, alpha, theta, delta, s1, rhok, rhokp1;
   int             *cols, allocated_space;
   double          *diagonal, *vals, *tdiag, dtemp1, dtemp2;
#ifdef GREG
   /* double tmp; */
#endif
   double lambda_min, lambda_max;

   n = outlen;
   widget = (struct MLSthing *) smooth_ptr->smoother->data;
   deg    = widget->mlsDeg;
   if (deg == 0) return 0;
   if (widget->block_scaling != NULL) {
     lambda_min = widget->scaled_matrix->lambda_min;
     lambda_max = widget->scaled_matrix->lambda_max;
   }
   else {
     lambda_min = Amat->lambda_min;
     lambda_max = Amat->lambda_max;
   }

   pAux  = (double *) ML_allocate((n+1)*sizeof(double));
   dk     = (double *) ML_allocate((n+1)*sizeof(double));

   if (pAux == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
   if (dk    == NULL) {
     pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
     ML_avoid_unused_param((void *) &inlen);
   }

   beta = 1.1*lambda_max;   /* try and bracket high */
   alpha = lambda_max/(widget->eig_ratio);
#ifdef GREG
#undef GREG
#endif
   delta = (beta - alpha)/2.;
   theta = (beta + alpha)/2.;
   s1 = theta/delta;
   rhok = 1./s1;
 

   /* ----------------------------------------------------------------- */
   /* extract diagonal using getrow function if not found               */
   /* ----------------------------------------------------------------- */

   if ((Amat->diagonal == NULL) && (widget->block_scaling == NULL))

   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(MLS_Apply): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&nn) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols); 
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
	    tdiag[i] = 0.;
            for (j = 0; j < nn; j++) 
               if (cols[j] == i) tdiag[i] = vals[j];
	    if (tdiag[i] == 0.) tdiag[i] = 1.;
         }
         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   if (widget->block_scaling == NULL)
     ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);

   /* This is meant for the case when the matrix is the identity.*/

   if ((lambda_min == 1.0) && (lambda_min == lambda_max)) {
     for (i = 0; i < n; i++) x[i] = rhs[i]/diagonal[i];
     return 0;
   }



   if (widget->block_scaling == NULL) { /* normal point scaling */

     if (smooth_ptr->init_guess == ML_NONZERO) {
       ML_Operator_Apply(Amat, n, x, n, pAux);
       for (i = 0; i < n; i++) {
	 dk[i] = (rhs[i] - pAux[i])/(theta*diagonal[i]);
	 x[i] += dk[i];
       }
     }
     else { 
       for (i = 0; i < n; i++) {
	 x[i] = dk[i] = rhs[i]/(theta*diagonal[i]);
       }
     }

     for (k = 0; k < deg-1; k++) {
       ML_Operator_Apply(Amat, n, x, n, pAux);
       rhokp1 = 1./(2.*s1 - rhok);
       dtemp1 = rhokp1*rhok;
       dtemp2 = 2.*rhokp1/delta;
       rhok = rhokp1;
       for (i = 0; i < n; i++) {
	 dk[i] = dtemp1 * dk[i] + dtemp2*(rhs[i]-pAux[i])/diagonal[i];
	 x[i] += dk[i];
       }
     }

   }
   else { /* block scaling */

     if (smooth_ptr->init_guess == ML_NONZERO) {
       ML_Operator_Apply(Amat, n, x, n, pAux);
       for (i = 0; i < n; i++) 	 dk[i] = (rhs[i] - pAux[i])/theta;
       ML_BlockDinv(widget->block_scaling, n, dk);
       for (i = 0; i < n; i++) 	 x[i] += dk[i];


     }
     else { 
       for (i = 0; i < n; i++) 	dk[i] = rhs[i]/theta;
       ML_BlockDinv(widget->block_scaling, n, dk);
       for (i = 0; i < n; i++) 	 x[i] = dk[i];
     }

     for (k = 0; k < deg-1; k++) {
       ML_Operator_Apply(Amat, n, x, n, pAux);
       rhokp1 = 1./(2.*s1 - rhok);
       dtemp1 = rhokp1*rhok;
       dtemp2 = 2.*rhokp1/delta;
       rhok = rhokp1;
       for (i = 0; i < n; i++) pAux[i] = rhs[i]-pAux[i];
       ML_BlockDinv(widget->block_scaling, n, pAux);
       for (i = 0; i < n; i++) {
	 dk[i] = dtemp1 * dk[i] + dtemp2*pAux[i];
	 x[i] += dk[i];
       }
     }
   }
   ML_free(dk);
   ML_free(pAux);

   return 0;	
}


/*****************************************************************************/
/* MSR Gauss-Seidel smoother with no damping.                      */
/* ------------------------------------------------------------------------- */

int ML_Smoother_MSR_GSforwardnodamping(void *sm,int inlen,double x[],
				       int outlen, double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum;
   int             *bindx;
   double          *val = NULL;
   struct ML_CSR_MSRdata *ptr = NULL;

   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
   if (val == NULL) {
     ML_Smoother_SGS(sm, inlen, x, outlen, rhs);
     return 0;
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO)
      {
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
      }
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      j = bindx[0];
      bindx_ptr = &bindx[j];
#ifndef rblack
      for (i = 0; i < Nrows; i++) {
#else
#ifdef fastrb
    if (iter == 0) {
      for (i = 0; i < Nrows; i += 2)
    if (val[i] != 0.0) x2[i] = rhs[i]/val[i];
    }
    else {
#endif
      for (i = 0; i < Nrows; i += 2) {
    j = bindx[i];             /* added for redblack */
    bindx_ptr = &bindx[j];    /* added for redblack */
#endif
	sum =  rhs[i];
	while (j+10 < bindx[i+1]) {
	  sum -= (val[j+9]*x2[bindx_ptr[9]] +
		  val[j+8]*x2[bindx_ptr[8]] +
		  val[j+7]*x2[bindx_ptr[7]] +
		  val[j+6]*x2[bindx_ptr[6]] +
		  val[j+5]*x2[bindx_ptr[5]] +
		  val[j+4]*x2[bindx_ptr[4]] +
		  val[j+3]*x2[bindx_ptr[3]] +
		  val[j+2]*x2[bindx_ptr[2]] +
		  val[j+1]*x2[bindx_ptr[1]] +
		  val[j]*x2[*bindx_ptr]);
	  bindx_ptr += 10;
	  j += 10;
	}
	while (j < bindx[i+1]) {
	  sum -= val[j++] * x2[*bindx_ptr++];
	}
	if (val[i] != 0.0) x2[i] = sum/val[i];
      }
#ifdef rblack
#ifdef fastrb
    }
#endif
      for (i = 1; i < Nrows; i += 2) {
    j = bindx[i];
    bindx_ptr = &bindx[j];
    sum =  rhs[i];
    while (j+10 < bindx[i+1]) {
      sum -= (val[j+9]*x2[bindx_ptr[9]] +
          val[j+8]*x2[bindx_ptr[8]] +
          val[j+7]*x2[bindx_ptr[7]] +
          val[j+6]*x2[bindx_ptr[6]] +
          val[j+5]*x2[bindx_ptr[5]] +
          val[j+4]*x2[bindx_ptr[4]] +
          val[j+3]*x2[bindx_ptr[3]] +
          val[j+2]*x2[bindx_ptr[2]] +
          val[j+1]*x2[bindx_ptr[1]] +
          val[j]*x2[*bindx_ptr]);
      bindx_ptr += 10;
      j += 10;
    }
    while (j < bindx[i+1]) {
      sum -= val[j++] * x2[*bindx_ptr++];
    }
    if (val[i] != 0.0) x2[i] = sum/val[i];
      }
#endif
   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}

int ML_Smoother_MSR_GSbackwardnodamping(void *sm,int inlen,double x[],
				       int outlen, double rhs[])
{
   int iter, i, j;
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2;
   ML_Smoother  *smooth_ptr;
   register int    *bindx_ptr;
   register double sum;
   int             *bindx;
   double          *val = NULL;
   struct ML_CSR_MSRdata *ptr = NULL;

   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
   if (val == NULL) {
     ML_Smoother_SGS(sm, inlen, x, outlen, rhs);
     return 0;
   }

   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                   *sizeof(double));
      if (x2 == NULL) {
         printf("Not enough space in Gauss-Seidel\n"); exit(1);
      }
      for (i = 0; i < inlen; i++) x2[i] = x[i];
      if (smooth_ptr->init_guess != ML_NONZERO)
      {
         for (i = inlen; i < inlen+getrow_comm->total_rcv_length+1; i++)
            x2[i] = 0.;
      }
   }
   else x2 = x;

   for (iter = 0; iter < smooth_ptr->ntimes; iter++)
   {

      if (((getrow_comm != NULL) && (smooth_ptr->init_guess == ML_NONZERO))
          || (iter != 0) )
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);


      j = bindx[Nrows];
      bindx_ptr = &(bindx[j]);
      j--;
      bindx_ptr--;
      for (i = Nrows - 1; i >= 0; i--) {
	sum    = rhs[i];

	while (j-9 >= bindx[i]) {
	  bindx_ptr -= 10;
	  j -= 10;
	  sum -= (val[j+10]*x2[bindx_ptr[10]] +
		  val[j+9]*x2[bindx_ptr[9]] +
		  val[j+8]*x2[bindx_ptr[8]] +
		  val[j+7]*x2[bindx_ptr[7]] +
		  val[j+6]*x2[bindx_ptr[6]] +
		  val[j+5]*x2[bindx_ptr[5]] +
		  val[j+4]*x2[bindx_ptr[4]] +
		  val[j+3]*x2[bindx_ptr[3]] +
		  val[j+2]*x2[bindx_ptr[2]] +
		  val[j+1]*x2[bindx_ptr[1]]);
	}
	while (j >= bindx[i]) {
	  sum -= val[j--] * x2[*bindx_ptr--];
	}
	if (val[i] != 0.0) x2[i] = sum/val[i];
/*	if (fabs(val[i]) > 1.0e-9)  x2[i] = sum/val[i]; These is needed on 
        nodal problems with very small mass matrix..... */
      }
   }

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }

   return 0;
}


int ML_complex_Cheby(void *sm, int inlen, double x[], int outlen, double rhs[])
{

  ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
  ML_Operator     *Amat = smooth_ptr->my_level->Amat;
  struct MLSthing *widget;
  int              deg, i, j, k, n, nn;
  double          *pAux, *dk;
  double beta_real, beta_img, alpha_real, alpha_img, theta_real, theta_img, delta_real, delta_img, s1_real, s1_img, rhok_real, rhok_img, rhokp1_real, rhokp1_img;
  int             *cols, allocated_space;

  double          *vals, *tdiag, dtemp1_real, dtemp1_img, dtemp2_real, dtemp2_img;
  struct ML_Operator_blockmat_data *blockmat;
  double t1, t2, t3, denom, denom_real, denom_img, numer_real, numer_img, ttti, tttipn, *d1, *d2;
  int status;
  double thenorm;
#ifdef GREG
  double tmp;
#endif

  n = outlen/2;
  widget = (struct MLSthing *) smooth_ptr->smoother->data;

  deg    = widget->mlsDeg;
  if (deg == 0)  return 0;

  blockmat = (struct ML_Operator_blockmat_data *) Amat->data;

  pAux  = (double *) ML_allocate(2*(n+1)*sizeof(double));
  dk     = (double *) ML_allocate(2*(n+1)*sizeof(double));

  if (pAux == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
  if (dk    == NULL) {
    pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
    ML_avoid_unused_param((void *) &inlen);
  }

  beta_real = 1.1*widget->beta_real;   /* try and bracket high */
  beta_img  = 1.1*widget->beta_img;    /* frequency errors.    */
  alpha_real= beta_real/(widget->eig_ratio);
  alpha_img = beta_img;

#ifdef GREG
if (Amat->comm->ML_mypid == 0) {
  printf("deg %d beta real %20.13e img %20.13e ratio %20.13e\n",deg,beta_real,beta_img,widget->eig_ratio);
   fflush(stdout);
}
#endif

  delta_real = (beta_real - alpha_real)/2.;
  delta_img = (beta_img - alpha_img)/2.;

  theta_real = (beta_real + alpha_real)/2.;
  theta_img = (beta_img + alpha_img)/2.;
  t1 = 1./(delta_real*delta_real + delta_img*delta_img);
  s1_real = (theta_real*delta_real+theta_img*delta_img)*t1;
  s1_img  = (theta_img*delta_real-theta_real*delta_img)*t1;
  t1 = 1./(s1_real*s1_real + s1_img*s1_img);
  rhok_real = s1_real*t1;
  rhok_img  = -s1_img*t1;

  /* ----------------------------------------------------------------- */
  /* extract diagonal using getrow function if not found               */
  /* ----------------------------------------------------------------- */

  if (blockmat->Ke_diag == NULL) {
    allocated_space = blockmat->N_Ke + blockmat->Nghost + 1;
    cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
    vals = (double *) ML_allocate(allocated_space*sizeof(double));
    tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
    for (i = 0; i < n; i++) {
      status = blockmat->Ke_getrow(blockmat->Ke_getrow_data,1,&i,
				   allocated_space, cols, vals, &nn);
      if (status == 0) {
	printf("ML_complex_Cheby: not enough space for getrow\n");
	exit(1);
      }
      tdiag[i] = 0.;
      for (j = 0; j < nn; j++) 
	if (cols[j] == i) tdiag[i] = vals[j];
      if (tdiag[i] == 0.) tdiag[i] = 1.;
    }
    ML_free(cols); ML_free(vals);
    blockmat->Ke_diag = tdiag;  tdiag = NULL;
  }
  d1 = blockmat->Ke_diag;
  if (blockmat->M_diag == NULL) {
    allocated_space = blockmat->N_Ke + blockmat->Nghost + 1;
    cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
    vals = (double *) ML_allocate(allocated_space*sizeof(double));
    tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
    for (i = 0; i < n; i++) {
      status = blockmat->M_getrow(blockmat->M_getrow_data,1,&i,
				   allocated_space, cols, vals, &nn);
      if (status == 0) {
	printf("ML_complex_Cheby: not enough space for getrow\n");
	exit(1);
      }
      tdiag[i] = 0.;
      for (j = 0; j < nn; j++) 
	if (cols[j] == i) tdiag[i] = vals[j];
    }
    ML_free(cols); ML_free(vals);
    blockmat->M_diag = tdiag;  tdiag = NULL;
  }
  d2 = blockmat->M_diag;
#ifdef GREG
tmp = sqrt(ML_gdot(2*n, x, x, Amat->comm));
if (Amat->comm->ML_mypid == 0) 
printf("x is %20.13e\n", tmp);
#endif

 if (smooth_ptr->init_guess == ML_NONZERO) {
#ifdef GREG
tmp = sqrt(ML_gdot(2*n, rhs, rhs, Amat->comm));
if (Amat->comm->ML_mypid == 0) 
printf("Res is %20.13e\n", tmp);
#endif
   ML_Operator_Apply(Amat, 2*n, x, 2*n, pAux);
   for (i = 0; i < 2*n; i++) {
     dk[i] = rhs[i] - pAux[i];
   }
   for (k = 1; k < n+1; k++) {
      denom_real = theta_real*d1[k-1] - theta_img*d2[k-1];
      denom_img  = theta_real*d2[k-1] + theta_img*d1[k-1];
      t1  = dk[k-1] * denom_real + dk[k+n-1]* denom_img;
      t2  = -dk[k-1] * denom_img + dk[k+n-1]* denom_real;
      t3 = 1./(denom_real*denom_real + denom_img*denom_img);
      dk[k  -1] = t1*t3;
      dk[k+n-1] = t2*t3;
      x[k  -1]  = x[k-1] + dk[k-1];
      x[k+n-1]  = x[k+n-1] + dk[k+n-1];
   }
}
  else { 
#ifdef GREG
tmp = sqrt(ML_gdot(2*n, rhs, rhs, Amat->comm));
if (Amat->comm->ML_mypid == 0) 
printf("res is %20.13e\n", tmp);
#endif
   for (k = 1; k < n+1; k++) {
       denom_real = theta_real*d1[k-1] - theta_img*d2[k-1];
       denom_img  = theta_real*d2[k-1] + theta_img*d1[k-1];
       t1  = rhs[k-1] * denom_real + rhs[k+n-1]* denom_img;
       t2  = -rhs[k-1] * denom_img + rhs[k+n-1]* denom_real;
       t3 = 1./(denom_real*denom_real + denom_img*denom_img);
       dk[k-1] = t1*t3;
       dk[k+n-1] = t2*t3;
       x[k-1] = dk[k-1];
       x[k+n-1] = dk[k+n-1];
     }
#ifdef GREG
tmp = sqrt(ML_gdot(2*n, x, x, Amat->comm));
if (Amat->comm->ML_mypid == 0) 
printf("x is %20.13e\n", tmp);
#endif
}
for (k = 0; k < deg-1; k++) {
 thenorm = 0.;
  ML_Operator_Apply(Amat, 2*n, x, 2*n, pAux);
   denom_real = 2*s1_real - rhok_real;
   denom_img  = 2*s1_img - rhok_img;
   t1 = 1./(denom_real*denom_real + denom_img*denom_img);
   rhokp1_real = denom_real*t1;
   rhokp1_img  = -denom_img*t1;
   dtemp1_real = rhokp1_real*rhok_real - rhokp1_img*rhok_img;
   dtemp1_img  = rhokp1_img*rhok_real + rhokp1_real*rhok_img;
   t1 = 1./(delta_real*delta_real + delta_img*delta_img);
   dtemp2_real = 2*(rhokp1_real*delta_real+rhokp1_img*delta_img)*t1;
   dtemp2_img  = 2*(rhokp1_img*delta_real-rhokp1_real*delta_img)*t1;
   rhok_real = rhokp1_real;
   rhok_img  = rhokp1_img;

   for (i = 1; i <= n; i++) {
    ttti = dtemp1_real*dk[i-1] - dtemp1_img*dk[i+n-1];
    dk[i+n-1] = dtemp1_img*dk[i-1] + dtemp1_real*dk[i+n-1];
    dk[i-1] = ttti;
    ttti = rhs[i-1] - pAux[i-1];
thenorm += (ttti*ttti);
    tttipn = rhs[i+n-1] - pAux[i+n-1];
thenorm += (tttipn*tttipn);
    numer_real = dtemp2_real*d1[i-1] + dtemp2_img*d2[i-1];
    numer_img  = dtemp2_img*d1[i-1] - dtemp2_real*d2[i-1];
    denom      = 1./(d1[i-1]*d1[i-1] + d2[i-1]*d2[i-1]);
    t1         = ttti*numer_real - tttipn*numer_img;
    t2         = ttti*numer_img + tttipn*numer_real;
    dk[i-1] = dk[i-1] + t1*denom;
    dk[i+n-1]= dk[i+n-1] + t2*denom;
    x[i  -1] = x[i  -1] + dk[i  -1];
    x[i+n-1] = x[i+n-1] + dk[i+n-1];
   }
#ifdef GREG
tmp = sqrt(ML_gsum_double(thenorm, Amat->comm));
if (Amat->comm->ML_mypid == 0) 
 printf("%d: the norm is %20.13e\n",k,tmp);
#endif
  }
#ifdef GREG
if (Amat->comm->ML_mypid == 0)  fflush(stdout);
#endif
  ML_free(dk);
  ML_free(pAux);
  return 0;
}


int ML_DiagScaled_1stepKrylov(void *sm, int inlen, double x[], int outlen, 
                              double rhs[])
{

  ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
  ML_Operator     *Amat = smooth_ptr->my_level->Amat;

  int    n, allocated_space, *cols, i, j, k, status, nn, flag_x = 0;
  double t1, t2, t3, aa = 0., bb = 0., ee = 0., alpha_real;
  double alpha_img, dks_real = 0., dks_img = 0.;
  double *vals, *tdiag, *Ke_diag, *M_diag, *pAux, *dk;
  double denom_real, denom_img, *dks, *pq, *pqAux;
  struct ML_Operator_blockmat_data *blockmat;


  n = outlen/2;

  blockmat = (struct ML_Operator_blockmat_data *) Amat->data;


  /* ----------------------------------------------------------------- */
  /* extract diagonal using getrow function if not found               */
  /* ----------------------------------------------------------------- */

  if (blockmat->Ke_diag == NULL) {
    allocated_space = blockmat->N_Ke + blockmat->Nghost + 1;
    cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
    vals = (double *) ML_allocate(allocated_space*sizeof(double));
    tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
    for (i = 0; i < n; i++) {
      status = blockmat->Ke_getrow(blockmat->Ke_getrow_data,1,&i,
				   allocated_space, cols, vals, &nn);
      if (status == 0) {
	printf("ML_complex_Cheby: not enough space for getrow\n");
	exit(1);
      }
      tdiag[i] = 0.;
      for (j = 0; j < nn; j++) 
	if (cols[j] == i) tdiag[i] = vals[j];
      if (tdiag[i] == 0.) tdiag[i] = 1.;
    }
    ML_free(cols); ML_free(vals);
    blockmat->Ke_diag = tdiag;  tdiag = NULL;
  }
  Ke_diag = blockmat->Ke_diag; /* stiffness diag */
  if (blockmat->M_diag == NULL) {
    allocated_space = blockmat->N_Ke + blockmat->Nghost + 1;
    cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
    vals = (double *) ML_allocate(allocated_space*sizeof(double));
    tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
    for (i = 0; i < n; i++) {
      status = blockmat->M_getrow(blockmat->M_getrow_data,1,&i,
				  allocated_space, cols, vals, &nn);
      if (status == 0) {
	printf("ML_complex_Cheby: not enough space for getrow\n");
	exit(1);
      }
      tdiag[i] = 0.;
      for (j = 0; j < nn; j++) 
	if (cols[j] == i) tdiag[i] = vals[j];
    }
    ML_free(cols); ML_free(vals);
    blockmat->M_diag = tdiag;  tdiag = NULL;
  }
  M_diag = blockmat->M_diag; /* mass diag */


  pAux  = (double *) ML_allocate(2*(n+1)*sizeof(double));
  dk     = (double *) ML_allocate(2*(n+1)*sizeof(double));

  if (pAux == NULL) pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
  if (dk    == NULL) {
    pr_error("ML_Smoother_MLS_Apply: allocation failed\n");
    ML_avoid_unused_param((void *) &inlen);
  }

  /* define previous UNSCALED residual 
   */
  flag_x = 0;
  if (smooth_ptr->init_guess == ML_NONZERO) {
    flag_x = 1;
    ML_Operator_Apply(Amat, 2*n, x, 2*n, pAux);
    for (i = 0; i < 2*n; i++) {
      dk[i] = rhs[i] - pAux[i];
    }
  }
    else {
      for (i = 0; i < 2*n; i++) {
      dk[i] = rhs[i];
      }
    }

    /* instantiate, initialize quantities 
     */
    dks = (double *) ML_allocate(2*(n+1)*sizeof(double));
    pqAux = (double *) ML_allocate(2*(n+1)*sizeof(double));
    pq = (double *) ML_allocate(2*(n+1)*sizeof(double));

    for (k = 1; k < n+1; k++) {
      denom_real = Ke_diag[k-1];
      denom_img  = M_diag[k-1];
      t3 = 1./(denom_real*denom_real + denom_img*denom_img);
      t1 = dk[k-1];
      t2 = dk[k+n-1];
      /* calculate dks = Dinv*dk without calculating Dinv, the equivalent real
         form of the inverse diagonal of the original complex matrix A (whose 
         equivalent real form is Amat)
      */
      dks[k-1] = (denom_real*t1 + denom_img*t2)*t3;
      dks[k+n-1] = (-denom_img*t1 + denom_real*t2)*t3;
    }    

    ML_Operator_Apply(Amat, 2*n, dks, 2*n, pqAux);
    for (k = 1; k < n+1; k++) {  
      /*      p     Dinv_r   Dinv_i   Ke_mat | -Me_mat
         pq = - =  --------|------- * -------|-------- * dks
              q    -Dinv_i   Dinv_r   Me_mat |  Ke_mat
      */
      denom_real = Ke_diag[k-1];
      denom_img  = M_diag[k-1];
      t3 = 1./(denom_real*denom_real + denom_img*denom_img);
      t1 = pqAux[k-1];
      t2 = pqAux[k+n-1];
      pq[k-1] = (denom_real*t1 +denom_img*t2)*t3;
      pq[k+n-1] = (-denom_img*t1 +denom_real*t2)*t3;
    }

    /* define step coefficient for scaled system based on 1-step GMRES 
     */
    for (k = 1; k < n+1; k++) {
      aa = aa + dks[k-1]*pq[k-1] + dks[k+n-1]*pq[k+n-1];
      bb = bb + dks[k+n-1]*pq[k-1] - dks[k-1]*pq[k+n-1];
      ee = ee + pq[k-1]*pq[k-1] + pq[k+n-1]*pq[k+n-1];
    }
    aa = ML_gsum_double(aa, Amat->comm);
    bb = ML_gsum_double(bb, Amat->comm);
    ee = ML_gsum_double(ee, Amat->comm);

    alpha_real = aa/ee; /*real part of step coefficient */
    alpha_img = bb/ee; /* imag part of step coefficient */

    for (k = 1; k < n+1; k++) {
      dks_real = alpha_real*dks[k-1] - alpha_img*dks[k+n-1];
      dks_img = alpha_img*dks[k-1] + alpha_real*dks[k+n-1];
      if (flag_x > 0) {
	x[k-1] = x[k-1] + dks_real;
	x[k+n-1] = x[k+n-1] + dks_img;
      }
      else {
	x[k-1] = dks_real;
	x[k+n-1] = dks_img;
      }
    }

  ML_free(dks);
  ML_free(pqAux);
  ML_free(pq);
  ML_free(dk);
  ML_free(pAux);
  return 0;

}




#ifdef WKC
/*  WKC
    BELOW IS OVERLOADED FUNCTIONS FOR EPETRA STUFF */

/* ************************************************************************* */
/* smoothing operation                                                       */
/* ************************************************************************* */



int ML_Smoother_Apply(ML_Smoother *pre, int inlen, Epetra_MultiVector &ep_sol, 
                      int outlen, Epetra_MultiVector &ep_rhs, int init_guess)
{
   double ** pp_sol;
   double ** pp_rhs;
   ep_sol.ExtractView ( &pp_sol );
   ep_rhs.ExtractView ( &pp_rhs );


   int         i, n;
   double      temp, *res, tol;
   ML_Operator *Amat;

#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   double      t0;
   t0 = GetClock();
#endif

   if (pre->smoother->ML_id == ML_EMPTY) return 1;
   pre->init_guess = init_guess;

   if (pre->smoother->ML_id == ML_EXTERNAL) {
      for ( int KK = 0 ; KK != ep_sol.NumVectors() ; KK++ ) {
         double *sol = pp_sol[KK];
         double *rhs = pp_rhs[KK];
        pre->smoother->external(pre->smoother->data,inlen,sol,outlen,rhs);
      }
   } else 
   {
      if (pre->ntimes == ML_CONVERGE) 
      {
         for ( int KK = 0 ; KK != ep_sol.NumVectors() ; KK++ ) {
            double *sol = pp_sol[KK];
            double *rhs = pp_rhs[KK];

         Amat = pre->my_level->Amat;
         n    = Amat->outvec_leng; 
         res  = (double *) ML_allocate( (n+1)*sizeof(double) );
         temp = sqrt(ML_gdot(n, rhs, rhs, pre->my_level->comm));
         tol  = temp*pre->tol;
         pre->ntimes = 100;
         while ( temp > tol ) 
         {
            pre->smoother->internal(pre,n,sol,n, rhs);
            ML_Operator_Apply(Amat, n, sol, n, res);
            for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];
            temp = sqrt(ML_gdot(n, res, res, pre->my_level->comm));
         }
         pre->ntimes = ML_CONVERGE;
         ML_free(res);
      }
      }
      else {

         if ( (void *)pre->smoother->internal == (void *)ML_Cheby )
            ML_Cheby_WKC ( pre , inlen , (double *)&ep_sol , outlen , 
                          (double *) &ep_rhs );

         else {
            for ( int KK = 0 ; KK != ep_sol.NumVectors() ; KK++ ) {
               double *sol = pp_sol[KK];
               double *rhs = pp_rhs[KK];

               pre->smoother->internal(pre,inlen,sol,outlen,rhs);
            }
         }
      }
   }
#if defined(ML_TIMING) || defined(ML_TIMING_DETAILED)
   pre->apply_time += (GetClock() - t0);
#endif
   return 1;
}



int ML_Cheby_WKC(void *sm, int inlen, double *pep_x, int outlen, double *pep_rhs)
{
   Epetra_MultiVector &ep_x ( *(Epetra_MultiVector *)pep_x );
   Epetra_MultiVector &ep_rhs ( *(Epetra_MultiVector *)pep_rhs );

   ML_Smoother     *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator     *Amat = smooth_ptr->my_level->Amat;
   struct MLSthing *widget;
   int              deg, i, j, k, n, nn;



   double beta, alpha, theta, delta, s1, rhok, rhokp1;
   int             *cols, allocated_space;
   double          *diagonal, *vals, *tdiag, dtemp1, dtemp2;

   n = outlen;
   widget = (struct MLSthing *) smooth_ptr->smoother->data;

   deg    = widget->mlsDeg;
   if (deg == 0) return 0;

   /* This is meant for the case when the matrix is the identity.*/
   if ((Amat->lambda_min == 1.0) && (Amat->lambda_min == Amat->lambda_max)) {
     for (i = 0; i < Amat->outvec_leng; i++) 
        for ( int KK = 0 ; KK != ep_x.NumVectors() ; KK++ )
           ep_x[KK][i] = ep_rhs[KK][i];
     return 0;
   }

   beta = 1.1*Amat->lambda_max;   /* try and bracket high */
   alpha = Amat->lambda_max/(widget->eig_ratio);

   delta = (beta - alpha)/2.;
   theta = (beta + alpha)/2.;
   s1 = theta/delta;
   rhok = 1./s1;
 

   /* ----------------------------------------------------------------- */
   /* extract diagonal using getrow function if not found               */
   /* ----------------------------------------------------------------- */

   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(MLS_Apply): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         tdiag = (double *) ML_allocate(Amat->outvec_leng*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&nn) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols); 
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
	    tdiag[i] = 0.;
            for (j = 0; j < nn; j++) 
               if (cols[j] == i) tdiag[i] = vals[j];
	    if (tdiag[i] == 0.) tdiag[i] = 1.;
         }
         ML_free(cols); ML_free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         ML_free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);

   Epetra_MultiVector ep_Aux ( ep_x );
   Epetra_MultiVector ep_dk ( ep_x );



   if (smooth_ptr->init_guess == ML_NONZERO) {
     ML_Operator_Apply(Amat, n, ep_x, n, ep_Aux);
     for (i = 0; i < n; i++) 
     for (int KK = 0 ; KK != ep_Aux.NumVectors() ; KK++ )
     {
       ep_dk[KK][i] = (ep_rhs[KK][i] - ep_Aux[KK][i])/(theta*diagonal[i]);
       ep_x[KK][i] += ep_dk[KK][i];
     }
   }
   else { 
     for (i = 0; i < n; i++) 
     for (int KK = 0 ; KK != ep_Aux.NumVectors() ; KK++ )
     {
       ep_x[KK][i] = ep_dk[KK][i] = ep_rhs[KK][i]/(theta*diagonal[i]);
     }
   }

   for (k = 0; k < deg-1; k++) {
     ML_Operator_Apply(Amat, n, ep_x, n, ep_Aux);
     rhokp1 = 1./(2.*s1 - rhok);
     dtemp1 = rhokp1*rhok;
     dtemp2 = 2.*rhokp1/delta;
     rhok = rhokp1;
     for (i = 0; i < n; i++) 
     for (int KK = 0 ; KK != ep_Aux.NumVectors() ; KK++ )
     {
       ep_dk[KK][i] = dtemp1 * ep_dk[KK][i] + dtemp2*
                             (ep_rhs[KK][i]-ep_Aux[KK][i])/diagonal[i];
       ep_x[KK][i] += ep_dk[KK][i];
     }
   }

   return 0;	
}


#endif
