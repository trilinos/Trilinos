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
/* ************************************************************************* */
/*  Methods available                                                        */
/*   - weighted Jacobi                                                       */
/*   - Gauss Seidel                                                          */
/*   - symmetric Gauss Seidel                                                */
/*   - block Gauss Seidel                                                    */
/*   - variable block Jacobi                                                 */
/*   - variable block Symmetric Gauss Seidel                                 */
/*   - variable block Symmetric Gauss Seidel (sequential)                    */
/*   - variable block Jacobi with Krylov                                     */
/*   - overlapped ILUT                                                       */
/* ************************************************************************* */

#include "ml_smoother.h"
#include "ml_lapack.h"

#define abs(x) (((x) > 0 ) ? x : -(x))

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
#ifdef ML_TIMING
   double          t1;
#endif
   ML_Sm_BGS_Data  *ml_data;
   ML_Sm_ILUT_Data *ilut_data;

#ifdef ML_TIMING
   if ( (ml_sm->label != NULL) && ( ml_sm->build_time != 0.0))
      printf(" Build time for %s\t= %e\n",ml_sm->label,ml_sm->build_time);
#endif

#ifdef ML_TIMING
   if (ml_sm->label != NULL) 
   {
       t1 = ML_gsum_double(ml_sm->apply_time, global_comm);
       t1 = t1/((double) global_comm->ML_nprocs);
       if ( (global_comm->ML_mypid == 0) && (t1 != 0.0))
          printf(" Apply time for %s\t= %e\n",ml_sm->label,t1);
   }
#endif

   ml_sm->ML_id = -1; 
   ml_sm->my_level = NULL;
   ml_sm->ntimes = 0;
   ml_sm->omega = 0;
   ml_sm->init_guess = ML_NONZERO;
   ml_sm->tol = 0;
   if (ml_sm->data_destroy != NULL) 
      ml_sm->data_destroy( ml_sm->smoother->data );
   ml_sm->data_destroy = NULL;
   if ( ml_sm->smoother->internal == ML_Smoother_VBlockSGS ||
        /* ml_sm->smoother->internal == ML_Smoother_BGS || */
        ml_sm->smoother->internal == ML_Smoother_VBlockJacobi )
   {
      ml_data = ml_sm->smoother->data;
      ML_Smoother_Destroy_BGS_Data(&(ml_data));
      ml_sm->smoother->data = NULL;
   }
   if ( ml_sm->smoother->internal == ML_Smoother_OverlappedILUT &&
        ml_sm->smoother->data != NULL )
   {
      ilut_data = ml_sm->smoother->data;
      ML_Smoother_Destroy_ILUT_Data(&(ilut_data));
      ml_sm->smoother->data = NULL;
   }
   ML_memory_free((void**)&(ml_sm->smoother));
   if (ml_sm->label != NULL) { free(ml_sm->label); ml_sm->label = NULL; }

   return 0;
}

/* ************************************************************************* */
/* set a label                                                               */
/* ************************************************************************* */

int ML_Smoother_Set_Label( ML_Smoother *smoo, char *label)
{
  int size;

   size = strlen(label) + 1;
   smoo->label = (char *) malloc(size*sizeof(char));
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
#ifdef ML_TIMING
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
         res  = (double *) malloc( (n+1)*sizeof(double) );
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
         free(res);
      }
      else pre->smoother->internal(pre,inlen,sol,outlen,rhs);
   }
#ifdef ML_TIMING
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
                    int ntimes, double omega)
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
   double *res,omega, *diagonal, *vals, *tdiag;
   ML_Smoother  *smooth_ptr;
#ifdef ML_SMOOTHER_DEBUG
   ML_Comm      *comm;
   double *res2, res_norm;
#endif

   smooth_ptr = (ML_Smoother *) sm;

#ifdef ML_SMOOTHER_DEBUG
   comm = smooth_ptr->my_level->comm;
#endif
   omega = smooth_ptr->omega;
   Amat = smooth_ptr->my_level->Amat;
   if (Amat->matvec->ML_id == ML_EMPTY) 
         pr_error("Error(ML_Jacobi): Need matvec\n");
   if (Amat->diagonal == NULL) 
   {
      if (Amat->getrow->ML_id == ML_EMPTY) 
         pr_error("Error(ML_Jacobi): Need diagonal\n");
      else 
      {
         allocated_space = 30;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         tdiag = (double *) malloc(Amat->outvec_leng*sizeof(double));
         for (i = 0; i < Amat->outvec_leng; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&n) == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols); 
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
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
         free(cols); free(vals);
         ML_Operator_Set_Diag(Amat, Amat->matvec->Nrows, tdiag);
         free(tdiag);
      } 
   }
   ML_DVector_GetDataPtr( Amat->diagonal, &diagonal);

   n     = Amat->outvec_leng;
   res   = (double *) malloc(n*sizeof(double));
#ifdef ML_SMOOTHER_DEBUG
   res2  = (double *) malloc(n*sizeof(double));
   printf("    ML_Jacobi, omega = %e\n", omega);
#endif

   for (j = 0; j < smooth_ptr->ntimes; j++) 
   {
      ML_Operator_Apply(Amat, n, x, n, res);
      for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];
      for (i = 0; i < n; i++) res[i] /= diagonal[i];
      for (i = 0; i < n; i++) x[i] += omega*res[i];
#ifdef ML_SMOOTHER_DEBUG
      ML_Operator_Apply(Amat, n, x, n, res2);
      for ( i = 0; i < n; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(n, res2, res2, comm));
      printf("      Jacobi : iter = %2d, rnorm = %e\n", j, res_norm);
#endif
   }
#ifdef ML_SMOOTHER_DEBUG
   free(res2);
#endif
   free(res);
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
#ifdef ML_SMOOTHER_DEBUG
   double *res2, res_norm;
#endif
   ML_Operator *Amat;
   ML_Comm *comm;
   ML_CommInfoOP *getrow_comm;
   int Nrows;
   double *x2, omega;
   ML_Smoother  *smooth_ptr;
   smooth_ptr = (ML_Smoother *) sm;

   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;
   omega = smooth_ptr->omega;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_GaussSeidel): Need getrow() for GS smoother\n");

   allocated_space = Amat->max_nz_per_row+1;
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
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

#ifdef ML_SMOOTHER_DEBUG
   res2 = (double*) malloc(Nrows * sizeof(double));
#endif
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);

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
#ifdef ML_SMOOTHER_DEBUG
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
   if (allocated_space != Amat->max_nz_per_row+1) 
      Amat->max_nz_per_row = allocated_space;

#ifdef ML_SMOOTHER_DEBUG
   free(res2);
#endif
   free(vals); free(cols);

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
#ifdef ML_SMOOTHER_DEBUG
   double *res2, res_norm, init_norm;
#endif
   ML_Smoother  *smooth_ptr;
   smooth_ptr = (ML_Smoother *) sm;

   omega = smooth_ptr->omega;
   Amat = smooth_ptr->my_level->Amat;
   comm = smooth_ptr->my_level->comm;
   Nrows = Amat->getrow->Nrows;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_SGS): Need getrow() for SGS smoother\n");

   allocated_space = Amat->max_nz_per_row+1;
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_SymGaussSeidel: Not enough space\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for SGS smoother\n");

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

#ifdef ML_SMOOTHER_DEBUG
   res2 = (double*) malloc(Nrows * sizeof(double));
#endif
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);

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
         /* ### Changes : C. Tong
         if (diag_term == 0.0)
            pr_error("Error: SGS() can not be used with a zero diagonal\n");
         */

         if (diag_term != 0.0)
            x2[i] += omega*(rhs[i] - dtemp)/diag_term;
      }
#ifdef ML_SMOOTHER_DEBUG
      ML_Operator_Apply(Amat, Nrows, x2, Nrows, res2);
      for ( i = 0; i < Nrows; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(Nrows, res2, res2, comm));
      printf("      SGS (for ) : iter = %2d, rnorm = %e\n", iter, res_norm);
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
         /* ### Changes : C. Tong
         if (diag_term == 0.0)
            pr_error("Error: GS() can not be used with a zero diagonal\n");
         if (diag_term != 0.0)
         */
            x2[i] += omega*(rhs[i] - dtemp)/diag_term;
      }
#ifdef ML_SMOOTHER_DEBUG
      ML_Operator_Apply(Amat, Nrows, x2, Nrows, res2);
      for ( i = 0; i < Nrows; i++ ) res2[i] = rhs[i] - res2[i];
      res_norm = sqrt(ML_gdot(Nrows, res2, res2, comm));
      printf("      SGS (back) : iter = %2d, rnorm = %e\n", iter, res_norm);
#endif
   }
#ifdef ML_SMOOTHER_DEBUG
   free(res2);
#endif

   if (getrow_comm != NULL) {
      for (i = 0; i < inlen; i++) x[i] = x2[i];
      ML_free(x2);
   }
   if (allocated_space != Amat->max_nz_per_row+1) {
      Amat->max_nz_per_row = allocated_space;
   }

   free(vals); free(cols);

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
   unsigned int   tmp;
	 
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

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_blockGaussSeidel): Need getrow() for smoother\n");

   allocated_space = Amat->max_nz_per_row+1;
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
   /* this is very space inefficient, but since lapack doesn't do sparsity,
      I'm not sure how to do better without rewriting their solver routine */
   Atimesx = (double *) malloc(blocksize*sizeof(double));
   correc = (double *) malloc(blocksize*sizeof(double));
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
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);

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
			   correc, &blocksize, &info, tmp);
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
                           correc, &blocksize, &info, tmp);
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
   if (allocated_space != Amat->max_nz_per_row+1) {
      Amat->max_nz_per_row = allocated_space;
   }
	 
   free(vals); free(cols); free(Atimesx); free(correc);
	 
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
   int            *do_update;
   ML_Comm        *comm;
   ML_CommInfoOP  *getrow_comm;
   ML_Operator    *Amat;
   double         *x_ext, **blockdata, *Ax, **res, *vals, omega;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_VBlockJacobi): Need getrow() for smoother\n");
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
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) malloc( Nblocks * sizeof(int) );
   aggr_group  = (int *) malloc( Nrows   * sizeof(int) );
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
      x_ext = (double *) malloc((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      Ax      = (double *) malloc( maxBlocksize * sizeof(double) );
   }
   if ( Nblocks > 0 )
   {
      res     = (double **) malloc( Nblocks * sizeof(double*) );
      for ( i = 0; i < Nblocks; i++ )
         res[i] = (double *) malloc( maxBlocksize * sizeof(double) );
      do_update = (int *) malloc( Nblocks * sizeof(int) );
      if ( do_update == NULL ) 
      {
         printf("ERROR : memory allocation.\n");
         exit(1);
      }
   }

   /* ----------------------------------------------------- */
   /* iterate                                               */
   /* ----------------------------------------------------- */

   strcpy(N,"N");
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE);

      for (i = 0; i < Nblocks; i++) 
      {
         do_update[i] = 0;
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
               if ( col == row ) do_update[i]++;
               Ax[k] += (vals[j]*x_ext[col]);
            }
            res[i][k] = rhs[row] - Ax[k];
         }
      }
 
      for (i = 0; i < Nblocks; i++) 
      {
         blocksize = blocklengths[i];
         if ( do_update[i] == blocksize && blocksize != 0 )
         {
            MLFORTRAN(dgetrs)(N,&blocksize,&one,blockdata[i],&blocksize,
                              perms[i], res[i], &blocksize, &info, itmp);
            if ( info != 0 ) 
            {
               printf("dgetrs returns with %d at block %d\n",info,i); 
               exit(1);
            }
            for (k = 0; k < blocksize; k++)
            {
               x_ext[aggr_group[aggr_offset[i]+k]] += (omega * res[i][k]);
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
      free(x_ext);
   }

   if ( maxBlocksize > 0 ) free( Ax );
   free( vals ); 
   free( cols ); 
   for (i = 0; i < Nblocks; i++) free( res[i] );
   if ( Nblocks > 0 ) free( res );
   free( aggr_offset );
   free( aggr_group );
   if ( Nblocks > 0 ) free( do_update );
	 
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
   double         *x_ext, **blockdata, *Ax, *res, *vals, omega;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_VBlockSymGS): Need getrow() for smoother\n");
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
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) malloc( Nblocks * sizeof(int) );
   aggr_group  = (int *) malloc( Nrows   * sizeof(int) );
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
      x_ext = (double *) malloc((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      Ax  = (double *) malloc( maxBlocksize * sizeof(double) );
      res = (double *) malloc( maxBlocksize * sizeof(double) );
   }

   /* ----------------------------------------------------- */
   /* iterate                                               */
   /* ----------------------------------------------------- */

   strcpy(N,"N");
   for (iter = 0; iter < smooth_ptr->ntimes; iter++) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE);

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
   }
   for (iter = smooth_ptr->ntimes; iter > 0; iter--) 
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE);

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
      free(x_ext);
   }
   if ( maxBlocksize > 0 )
   {
      free( Ax );
      free( res );
   }
   free( vals ); 
   free( cols ); 
   if ( Nblocks > 0 ) free( aggr_offset );
   if ( Nrows > 0 ) free( aggr_group );
	 
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
   double         *x_ext, **blockdata, *Ax, *res, *vals, omega;
   ML_Smoother    *smooth_ptr;
   ML_Sm_BGS_Data *dataptr;
   char           N[2];
   unsigned int   itmp;
	 
   /* ----------------------------------------------------- */
   /* fetch parameters                                      */
   /* ----------------------------------------------------- */

   smooth_ptr    = (ML_Smoother *) sm;
   comm          = smooth_ptr->my_level->comm;
   Amat          = smooth_ptr->my_level->Amat;
   omega         = smooth_ptr->omega;
   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_VBlockSymGSSeq): Need getrow() for smoother\n");
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
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
   maxBlocksize = 0;
   for ( i = 0; i < Nblocks; i++ )
      if ( blocklengths[i] > maxBlocksize ) 
         maxBlocksize = blocklengths[i];
   aggr_offset = (int *) malloc( Nblocks * sizeof(int) );
   aggr_group  = (int *) malloc( Nrows   * sizeof(int) );
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
      x_ext = (double *) malloc((inlen+getrow_comm->total_rcv_length+1)*
                                 sizeof(double));
      for (i = 0; i < inlen; i++) x_ext[i] = x[i];
   }
   else x_ext = x;
   if ( maxBlocksize > 0 )
   {
      Ax  = (double *) malloc( maxBlocksize * sizeof(double) );
      res = (double *) malloc( maxBlocksize * sizeof(double) );
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
            ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE);

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
            ML_exchange_bdry(x_ext,getrow_comm, inlen,comm,ML_OVERWRITE);

         if ( token == mypid )
         {
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
      free(x_ext);
   }

   if ( maxBlocksize > 0 )
   {
      free( Ax );
      free( res );
   }
   free( vals ); 
   free( cols ); 
   if ( Nblocks > 0 ) free( aggr_offset );
   if ( Nrows > 0 ) free( aggr_group );
	 
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

   dbuffer = (double *) malloc(extNrows * sizeof(double));
   idiag   = (int *)    malloc(extNrows * sizeof(int));
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
         ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE);
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
   free(dbuffer);
   free(idiag);
      
   return 0;
}

/* ******************************************************************** */
/* ******************************************************************** */
/* setup routines for various smoothers                                 */
/* ******************************************************************** */
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
   ml_data->blocksize = 1;
   ml_data->blocklengths = NULL;
   ml_data->blockmap = NULL;
   return 0;
}

/* ******************************************************************** */
/* Destructor for Sm_BGS_Data                                          */
/* ******************************************************************** */

int ML_Smoother_Destroy_BGS_Data(ML_Sm_BGS_Data **data)
{
   int i;
   ML_Sm_BGS_Data *ml_data;

   ml_data = (*data);
   if ( ml_data->blockfacts != NULL )
   {
      for ( i = 0; i < ml_data->Nblocks; i++ )
         if ( ml_data->blockfacts[i] != NULL ) free(ml_data->blockfacts[i]);
      free( ml_data->blockfacts );
   }
   if ( ml_data->perms != NULL )
   {
      for ( i = 0; i < ml_data->Nblocks; i++ )
         if ( ml_data->perms[i] != NULL ) free(ml_data->perms[i]);
      free( ml_data->perms );
   }
   if ( ml_data->blocklengths != NULL )
      free( ml_data->blocklengths );
   ml_data->blockmap = NULL;
   ML_memory_free((void**) data);
   return 0;
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
      free(perms[i]);
      free(blockfacts[i]);
   }

   free(perms);
   free(blockfacts);
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

int ML_Smoother_Destroy_ILUT_Data(ML_Sm_ILUT_Data **data)
{
   ML_Sm_ILUT_Data *ml_data;

   ml_data = (ML_Sm_ILUT_Data *) (*data);
   if ( ml_data->mat_ia != NULL ) free(ml_data->mat_ia);
   if ( ml_data->mat_ia != NULL ) free(ml_data->mat_ja);
   if ( ml_data->mat_ia != NULL ) free(ml_data->mat_aa);
   ML_memory_free( (void **) data );
   (*data) = NULL;
   return 0;
}

/* ************************************************************************* */
/* function to generate the factorizations of the diagonal blocks of A.      */
/* Factorizations are computed using lapack                                  */
/* ************************************************************************* */

int ML_Smoother_Gen_BGSFacts(ML_Sm_BGS_Data **data, ML_Operator *Amat, 
                             int blocksize) 
{
   int            i, j, *cols, allocated_space, length, Nrows, Nblocks;
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
   allocated_space = Amat->max_nz_per_row+1;
   dataptr->blocksize = blocksize;

   dataptr->blockfacts = (double **)malloc(Nblocks*sizeof(double *));
   dataptr->perms = (int **)malloc(Nblocks*sizeof(int *));
   blockfacts = dataptr->blockfacts;
   perms = dataptr->perms;
   for (j=0; j<Nblocks; j++) 
   {
      blockfacts[j]=(double *)calloc(blocksize*blocksize,sizeof(double));
      perms[j]=(int *)malloc(blocksize*blocksize*sizeof(int));
   }
   cols = (int    *) malloc(allocated_space*sizeof(int    ));
   vals = (double *) malloc(allocated_space*sizeof(double ));

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
   free(cols);
   free(vals);
	
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

   Nrows   = Amat->getrow->Nrows;
   dataptr = (*data);
   allocated_space = Amat->max_nz_per_row+1;

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
   dataptr->blockmap = blockIndices;
   if ( blockIndices == NULL )
   { 
      printf("ML_Gen_VBGSFacts : blocking information not available.\n");
      exit(1);
   }
   dataptr->blocklengths = (int*) malloc( Nblocks * sizeof(int));
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
            printf("ML_Gen_VBGSFacts : block index not valid %d. \n",
                                       blockIndices[i]);
            exit(1);
         }
      } else
         block_sizes[blockIndices[i]]++;
   }
   
   block_offset = (int *) malloc(Nrows*sizeof(int ));
   cols = (int    *) malloc(allocated_space*sizeof(int    ));
   vals = (double *) malloc(allocated_space*sizeof(double ));
   
   /* ----------------------------------------------------------- */
   /* allocate memory for each block                              */
   /* ----------------------------------------------------------- */

   dataptr->blockfacts = (double **) malloc(Nblocks*sizeof(double *));
   dataptr->perms = (int **)malloc(Nblocks*sizeof(int *));
   blockfacts = dataptr->blockfacts;
   perms = dataptr->perms;
   for ( i = 0; i < Nblocks; i++) 
   {
      length = block_sizes[i] * block_sizes[i];
      blockfacts[i] = (double *)malloc( length* sizeof(double) );
      for ( j = 0; j < length; j++) blockfacts[i][j] = 0.0; 
      perms[i] = (int *)malloc( block_sizes[i] * sizeof(int));
   }

   /* ----------------------------------------------------------- */
   /* load the block matrices                                     */
   /* ----------------------------------------------------------- */

   block_offset = (int *) malloc(Nrows*sizeof(int ));
   cols = (int    *) malloc(allocated_space*sizeof(int    ));
   vals = (double *) malloc(allocated_space*sizeof(double ));
   if (vals == NULL) 
      pr_error("Error in ML_Gen_VBJacobiFacts: Not enough space\n");

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
         for (j = 0; j < length; j++) {
            col = cols[j];
            if ( col < Nrows )
            {
               if ( blockIndices[col] == block_num )
               {
                  col_in_block = block_offset[col];
                  index = col_in_block * block_sizes[block_num] + row_in_block;
                  blockfacts[block_num][index] = vals[j];
               }
            }
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
         printf("Error in ML_Gen_VBJacobi: dgetrf returned %d (!=0)\n",info);
         exit(1);
      }
   }

   /* ----------------------------------------------------------- */
   /* clean up                                                    */
   /* ----------------------------------------------------------- */

   free(cols);
   free(vals);
   free(block_offset);
	
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
 
   proc_array  = (int *) malloc(nprocs * sizeof(int) );
   proc_array2 = (int *) malloc(nprocs * sizeof(int) );
   for ( i = 0; i < nprocs; i++ ) proc_array[i] = 0;
   proc_array[mypid] = Nrows;
   ML_gsum_vec_int(proc_array, proc_array2, nprocs, comm);
   NrowsOffset = 0;
   for (i = 0; i < mypid; i++) NrowsOffset += proc_array[i];
   for (i = 1; i < nprocs; i++) proc_array[i] += proc_array[i-1];
   free(proc_array2);

   /* ----------------------------------------------------------- */
   /* compose the column index map (index_array,index_array2)     */
   /* ----------------------------------------------------------- */

   dble_array  = (double *) malloc(extNrows *sizeof(double));
   for (i = Nrows; i < extNrows; i++) dble_array[i] = 0.0;
   for (i = 0; i < Nrows; i++) dble_array[i] = 1.0 * ( i + NrowsOffset );
   if (getrow_comm != NULL)
      ML_exchange_bdry(dble_array,getrow_comm, Nrows,comm,ML_OVERWRITE);
   index_array = ( int *) malloc((extNrows-Nrows) * sizeof(int));
   for (i = Nrows; i < extNrows; i++) index_array[i-Nrows] = dble_array[i];
   index_array2  = (int *) malloc((extNrows-Nrows) *sizeof(int));
   for (i = 0; i < extNrows-Nrows; i++) index_array2[i] = i;
   free( dble_array );

   /* ----------------------------------------------------------- */
   /* send the lengths of each row to remote processor            */
   /* at the end, additional row information should be given      */
   /* in total_recv_leng, recv_lengths, int_buf, dble_buf         */
   /* ----------------------------------------------------------- */

   ML_Smoother_GetRowLengths(getrow_comm, comm, Amat, total_recv_leng, 
                             recv_lengths); 
   ML_Smoother_GetOffProcRows(getrow_comm, comm, Amat, *total_recv_leng, 
                              *recv_lengths, NrowsOffset, index_array,
                              index_array2, int_buf, dble_buf); 

   free(proc_array);
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

   request      = (USR_REQ  *) malloc(N_neighbors*sizeof(USR_REQ ));
   (*recv_leng) = (int  *)     malloc(total_recv * sizeof(int));

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

   allocated_space = Amat->max_nz_per_row + 1;
   cols = (int *) malloc(allocated_space * sizeof(int));
   vals = (double *) malloc(allocated_space * sizeof(double));
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
      free( temp_list );
   }
   free(cols);
   free(vals);

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
      comm->USR_waitbytes((void *) &((*recv_leng)[offset]), nbytes, &proc_id, 
                          &msgtype, comm->USR_comm, request+i);
      offset += length;
   }
   free(request);
   free(neighbors);
   return 0;
}

/*****************************************************************************/
/* needed for overlapped smoothers                                           */
/*****************************************************************************/

int ML_Smoother_GetOffProcRows(ML_CommInfoOP *comm_info, ML_Comm *comm, 
                           ML_Operator *Amat, int leng, int *recv_leng,
                           int Noffset, int *map, int *map2, int **int_buf, 
                           double **dble_buf)
{
   int     N_neighbors, *neighbors, total_recv, mtype, msgtype, proc_id;
   int     i, j, index, *temp_list, nbytes, length, offset, allocated_space;
   int     *cols, *isend_buf, Nrows, m, k, nnz, nnz_offset;
   double  *vals, *send_buf;
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

   request     = (USR_REQ  *) malloc(N_neighbors*sizeof(USR_REQ ));
   (*int_buf)  = (int  *)     malloc(total_recv * sizeof(int));
   (*dble_buf) = (double  *)  malloc(total_recv * sizeof(double));

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

   allocated_space = Amat->max_nz_per_row + 1;
   cols = (int *) malloc(allocated_space * sizeof(int));
   vals = (double *) malloc(allocated_space * sizeof(double));
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
      if ( nnz > 0 ) send_buf = (double *) malloc( nnz * sizeof(double));
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
      free( temp_list );
      free( send_buf );
   }
   free(cols);
   free(vals);

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
      comm->USR_waitbytes((void *) &((*dble_buf)[nnz_offset]), nbytes, 
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

   allocated_space = Amat->max_nz_per_row + 1;
   cols = (int *) malloc(allocated_space * sizeof(int));
   vals = (double *) malloc(allocated_space * sizeof(double));
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
      if ( nnz > 0 ) isend_buf = (int *) malloc( nnz * sizeof(int));
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
      free( temp_list );
      free( isend_buf );
   }
   free(cols);
   free(vals);

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
      comm->USR_waitbytes((void *) &((*int_buf)[nnz_offset]), nbytes, 
 		   &proc_id, &msgtype, comm->USR_comm, request+i);
      offset += length;
      nnz_offset += nnz;
   }

   free(request);
   free(neighbors);
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
   int             mypid, track_leng, *track_array, *sortcols;
   double          *vals, ddata, tau, *mat_aa, *diagonal, *rowNorms;
   double          *dble_buf, *sortvals, absval, rel_tau;
   ML_Sm_ILUT_Data *ilut_ptr;

   /* ---------------------------------------------------------- */
   /* fetch ILUT parameters                                      */
   /* ---------------------------------------------------------- */

   mypid       = comm->ML_mypid;
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
   cols = (int *) malloc(allocated_space * sizeof(int));
   vals = (double *) malloc(allocated_space * sizeof(double));
   sortcols = (int *)    malloc(extNrows * sizeof(int));
   sortvals = (double *) malloc(extNrows * sizeof(double));
   dble_buf = (double *) malloc(extNrows * sizeof(double));
   diagonal = (double *) malloc(extNrows * sizeof(double));
   rowNorms = (double *) malloc(extNrows * sizeof(double));

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
         if ( cols[j] < extNrows ) rowNorms[i] += abs(vals[j]);
      rowNorms[i] /= extNrows;
   }
   for ( i = 0; i < total_recv_leng; i++ ) total_nnz += recv_lengths[i];
   total_nnz *= (fillin + 1);
   ilut_ptr->mat_ia = (int *) malloc( (extNrows + 1 ) * sizeof(int));
   ilut_ptr->mat_ja = (int *) malloc( total_nnz * sizeof(int));
   ilut_ptr->mat_aa = (double *) malloc( total_nnz * sizeof(double));
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
         if ( ext_ja[j] != -1 ) rowNorms[i+Nrows] += abs(ext_aa[j]);
      }
      rowNorms[i+Nrows] /= extNrows;
      offset += recv_lengths[i];
   }

   /* ---------------------------------------------------------- */
   /* process the first Nrows                                    */
   /* ---------------------------------------------------------- */

   nnz_count = 0;
   mat_ia[0] = 0;
   track_array = (int *) malloc( extNrows * sizeof(int) );
   for ( i = 0; i < extNrows; i++ ) dble_buf[i] = 0.0;

   for ( i = 0; i < Nrows; i++ )
   {
#ifdef ML_SMOOTHER_DEBUG
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
         if ( abs(dble_buf[j]) > rel_tau )
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
            absval = abs( dble_buf[index] );
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
      if ( abs(diagonal[i]) < 1.0e-16 ) diagonal[i] = 1.0E-6;
      mat_aa[nnz_count] = diagonal[i]; 
      mat_ja[nnz_count++] = i;
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i )
         {
            absval = dble_buf[index];
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
#ifdef ML_SMOOTHER_DEBUG
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
         if ( abs(dble_buf[j]) > rel_tau )
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
            absval = abs( dble_buf[index] );
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
      if ( abs(diagonal[i+Nrows]) < 1.0e-16 ) diagonal[i+Nrows] = 1.0E-6;
      mat_aa[nnz_count] = diagonal[i+Nrows]; 
      mat_ja[nnz_count++] = i+Nrows;
      sortcnt = 0;
      for ( j = 0; j < track_leng; j++ )
      {
         index = track_array[j];
         if ( index > i+Nrows )
         {
            absval = abs( dble_buf[index] );
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
#ifdef ML_SMOOTHER_DEBUG
   printf("%4d : ILUT Smoother - nnz(ILUT) = \n", mypid, nnz_count);
#endif

   /* ---------------------------------------------------------- */
   /* deallocate temporary storage space                         */
   /* ---------------------------------------------------------- */

   free(cols);
   free(vals);
   free(dble_buf);
   free(diagonal);
   free(rowNorms);
   free(sortcols);
   free(sortvals);
   free(track_array);
   return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/* AZTEC related smoothers                                                   */
/*****************************************************************************/

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
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);


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
         ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);


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

      ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);
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
      ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);


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
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   allocated_space = Amat->max_nz_per_row+1;
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
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
   if (allocated_space != Amat->max_nz_per_row+1) 
      Amat->max_nz_per_row = allocated_space;

   free(vals); free(cols);

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
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif

   allocated_space = Amat->max_nz_per_row+1;
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
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
            ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);
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
           ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE);

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
   if (allocated_space != Amat->max_nz_per_row+1) 
      Amat->max_nz_per_row = allocated_space;

   free(vals); free(cols);

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
   cols = (int    *) malloc(allocated_space*sizeof(int   ));
   vals = (double *) malloc(allocated_space*sizeof(double));
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
   free(data);
}

