/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the CG solver                                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_cg.h"


/* ******************************************************************** */
/* ML_CG_Solve                                                          */
/* ******************************************************************** */

int ML_CG_Solve(ML_Krylov *data, int length, double *rhs, double *sol)
{
   int         i, its, maxiter, converged=0, print_freq;
   double      alpha, beta, rho, rhom1, init_norm, sigma, eps1, tol;
   double      *r, *p, *ap, *z, res_norm;
   void        *precon;
   int         (*precfcn)(void*,int,double*,int,double*);
   ML_Operator *matrix;
   ML_Comm     *comm;

   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object                           */
   /* ----------------------------------------------------------------*/

   maxiter    = ML_Krylov_Get_MaxIterations(data);
   tol        = ML_Krylov_Get_Tolerance(data);
   precon     = ML_Krylov_Get_Precon(data);
   precfcn    = data->ML_precfcn;
   matrix     = ML_Krylov_Get_Amatrix(data);
   comm       = ML_Krylov_Get_Comm(data);
   print_freq = ML_Krylov_Get_PrintFreq(data);

   /* ----------------------------------------------------------------*/
   /* allocate temporary memory */
   /* ----------------------------------------------------------------*/

   r  = (double *) ML_allocate(length * sizeof(double));
   p  = (double *) ML_allocate(length * sizeof(double));
   ap  = (double *) ML_allocate(length * sizeof(double));
   if ( precfcn == NULL ) z = r;
   else {
      z = (double *) ML_allocate(length * sizeof(double));
      for ( i = 0; i < length; i++ ) z[i] = 0.0;
   }

   /* ----------------------------------------------------------------*/
   /* compute initial residual vector and norm */
   /* ----------------------------------------------------------------*/

   ML_Operator_Apply(matrix, length, sol, length, r);
   for ( i = 0; i < length; i++ ) {r[i] = rhs[i] - r[i]; p[i] = 0.0;}
   res_norm = sqrt(ML_gdot(length, r, r, comm));
   init_norm = res_norm;
   if ( comm->ML_mypid == 0 && print_freq > 0 && print_freq < 1000 )
       printf("ML_CG initial residual norm = %e \n", init_norm);
   if ( init_norm == 0.0 ) 
   {
      ML_free(r); ML_free(p); ML_free(ap);
      if ( precon != NULL ) ML_free(z);
      return 1;
   }

   /* ----------------------------------------------------------------*/
   /* initialization */
   /* ----------------------------------------------------------------*/

   its  = 0;
   eps1 = tol * init_norm;
   rho  = 0.0;

   /* ----------------------------------------------------------------*/
   /* loop until convergence is achieved or maxit has been exceeded*/
   /* ----------------------------------------------------------------*/

   while (converged == 0) 
   {
      while (res_norm > eps1 && its < maxiter) 
      {
         its++;
         if ( precfcn != NULL ) precfcn(precon,length,z,length,r);
         if ( its > 1 && rho == 0.0 ) 
         {
            printf("ML_CG breakdown (1).\n");
            exit(1);
         }
         rhom1 = rho;
         rho = ML_gdot(length, r, z, comm);
         if (its == 1) beta = 0.0;
         else          beta = rho / rhom1;
         for (i=0; i<length; i++) p[i] = beta * p[i] + z[i];
         ML_Operator_Apply(matrix, length, p, length, ap);
         sigma = ML_gdot(length, p, ap, comm);
         if ( sigma == 0.0 ) 
         {
            printf("ML_CG breakdown (2).\n");
            exit(1);
         }
         alpha  = rho / sigma;
         for (i=0; i<length; i++) 
         {
            sol[i] += alpha * p[i];
            r[i]   -= alpha * ap[i];
         }
         res_norm = sqrt(ML_gdot(length, r, r, comm));
         if ( print_freq != 0 && its % print_freq == 0 && comm->ML_mypid == 0)
            printf("ML_CG : iter %d - res. norm = %e \n",its,res_norm);
      }
      ML_Operator_Apply(matrix, length, sol, length, r);
      for ( i = 0; i < length; i++ ) r[i] = rhs[i] - r[i];
      res_norm = sqrt(ML_gdot(length, r, r, comm));
      if (comm->ML_mypid == 0 && print_freq > 0 && print_freq < 1000 )
         printf("ML_CG final residual norm = %e \n", res_norm);
      if (res_norm < eps1 || its >= maxiter) converged = 1;
   }
   if (comm->ML_mypid == 0 && print_freq > 0 && print_freq < 1000 )
      printf("ML_CG : total number of iterations = %d \n",its);

   /* ----------------------------------------------------------------*/
   /* de-allocate storage for temporary vectors */
   /* ----------------------------------------------------------------*/

   ML_free(r); ML_free(p); ML_free(ap);
   if ( precon != NULL) ML_free(z);
   return 1;
}

/* ******************************************************************** */
/* ML_CG_ComputeEigenvalues                                             */
/*                                                                      */
/* NOTE: if scale_by_diag == ML_TRUE compute lambda (D^{1/2} A D^{1/2}) */
/*       if scale_by_diag == ML_FALSE compute  lambda ( A )             */
/* ******************************************************************** */

int ML_CG_ComputeEigenvalues(ML_Krylov *data, int length, int scale_by_diag)
{
   int         i, j, k, its, maxiter, ncnt, *colInd, allocated, print_freq;
   int         *offset_array, myoffset, *itmp_array, nprocs, mypid, totallength;
   int         ext_leng, *index_array, total_length, original_maxiter, Nbc;
   double      alpha, beta, rho, rhom1, sigma, offdiag_norm;
   double      *r = NULL, *p = NULL, *ap = NULL, res_norm, *alpha_array, *colVal, *diag=NULL;
   double      *rhs=NULL, *rnorm_array, **Tmat, init_offdiag_norm;
   double      app, aqq, arr, ass, apq, sign, tau, t, c, s, *u;
   double      max_row_sum, min_row_sum, sum;
   ML_Operator *matrix;
   ML_Comm     *comm;
   ML_CommInfoOP *getrow_comm;
#ifdef PRINTMAT
   FILE        *fp;
#endif
   char        fname[100];

   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* ----------------------------------------------------------------*/

   matrix      = ML_Krylov_Get_Amatrix(data);
   comm        = ML_Krylov_Get_Comm(data);
   totallength = ML_Comm_GsumInt(comm, length);
   print_freq  = ML_Krylov_Get_PrintFreq(data);
   nprocs      = comm->ML_nprocs;
   mypid       = comm->ML_mypid;
   maxiter     = 10;
   if ( totallength < maxiter ) maxiter = totallength;

   /* There may be a bug in this function if A is 1x1. */
/*
   maxiter = - maxiter;
   maxiter = ML_gmax_int(maxiter, comm);
   maxiter = - maxiter;
*/
/* Is this really what we want ... if the length on one processor is small */
/* we could get 0 for maxiter? */

   /* ----------------------------------------------------------------*/
   /* set up to write matrix to a file */
   /* ----------------------------------------------------------------*/

   offset_array = (int *) ML_allocate(nprocs * sizeof(int));
   itmp_array   = (int *) ML_allocate(nprocs * sizeof(int));
   for ( i = 0; i < nprocs; i++ ) offset_array[i] = 0;
   offset_array[mypid] = length;
   ML_gsum_vec_int(&offset_array, &itmp_array, nprocs, comm);
   ML_free(itmp_array);
   myoffset = 0;
   for ( i = 0; i < mypid; i++ ) myoffset += offset_array[i];
   total_length = 0;
   for ( i = 0; i < nprocs; i++ ) total_length += offset_array[i];
   ML_free(offset_array);

   getrow_comm = matrix->getrow->pre_comm;
   if (getrow_comm != NULL) {
      ext_leng = length + getrow_comm->total_rcv_length;
   } else ext_leng = length;
   u = (double *) ML_allocate((ext_leng+1)*sizeof(double));
   for (i = 0; i < length; i++) u[i] = myoffset + i;
   for (i = length; i <= ext_leng; i++) u[i] = 0.0;
   if (getrow_comm != NULL) {
      ML_exchange_bdry(u,getrow_comm, length,comm,ML_OVERWRITE,NULL);
   }
   index_array = (int *) ML_allocate((ext_leng+1)*sizeof(int));
   for (i = 0; i <= ext_leng; i++) index_array[i] = (int) u[i];
   ML_free(u);

   /* ----------------------------------------------------------------*/
   /* allocate temporary memory  */
   /* ----------------------------------------------------------------*/

   if ( length > 0 )
   {
      u    = (double *) ML_allocate(length * sizeof(double));
      r    = (double *) ML_allocate(length * sizeof(double));
      p    = (double *) ML_allocate(length * sizeof(double));
      ap   = (double *) ML_allocate(length * sizeof(double));
      rhs  = (double *) ML_allocate(length * sizeof(double));
      diag = (double *) ML_allocate(length * sizeof(double));
      if ( diag == NULL )
      {
         printf("ML : ERROR in allocating memory.\n");
         exit(1);
      }
   }
   ML_random_vec(rhs, length, comm); 

   res_norm = sqrt(ML_gdot(length, rhs, rhs, comm));
   for (i = 0; i < length; i++) rhs[i] /= res_norm;

   alpha_array = (double *) ML_allocate((maxiter+1) * sizeof(double));
   rnorm_array = (double *) ML_allocate((maxiter+1) * sizeof(double));
   Tmat = (double **) ML_allocate((maxiter+1) * sizeof(double*));
   original_maxiter = maxiter;
   for ( i = 0; i <= maxiter; i++ ) 
   { 
      Tmat[i] = (double *) ML_allocate((maxiter+1) * sizeof(double));
      for ( j = 0; j <= maxiter; j++ ) Tmat[i][j] = 0.0;
      Tmat[i][i] = 1.0;
   }

   /* ----------------------------------------------------------------*/
   /* retrieve the diagonals */
   /* ----------------------------------------------------------------*/

   allocated = 100;
   colInd = (int    *) ML_allocate( allocated * sizeof(int) );
   colVal = (double *) ML_allocate( allocated * sizeof(double) );
   sprintf(fname, "mat_%d", mypid);
#ifdef PRINTMAT
   fp = fopen(fname, "w");
#endif
   Nbc = 0;  /* rst to handle nonsymmetric (due to BCs) matrices */
   if (scale_by_diag) {
   for ( i = 0; i < length; i++ )
   {
      while(ML_Operator_Getrow(matrix,1,&i,allocated,colInd,colVal,&ncnt)==0)
      {
         allocated *= 2;
         ML_free(colInd); ML_free(colVal);
         colInd = (int    *) ML_allocate( allocated * sizeof(int) );
         colVal = (double *) ML_allocate( allocated * sizeof(double) );
      }

      sum = 0.0;
      for ( j = 0; j < ncnt; j++ )
      {
         if ( colInd[j] == i ) diag[i] = colVal[j];
         else sum += ML_dabs(colVal[j]);
#ifdef PRINTMAT
         fprintf(fp,"A(%d,%d) = %e;\n",myoffset+i+1,index_array[colInd[j]]+1,
            colVal[j]);
#endif
      }
      /* kludging this in to handle Dirichlet BC's */
      if ( sum == 0.0) { rhs[i] = 0.; Nbc++; diag[i] = 1.;}

      else {
         if ( diag[i] == 0.0 ) 
         {
            printf("%d : diagonal[%d] == 0.0.\n", comm->ML_mypid, i);
	    /*            exit(1); */
	    diag[i] = 1.;
         }
         else if ( diag[i] < 0.0 )
         {
            printf("%d : diagonal[%d] = %e < 0.0.\n", comm->ML_mypid, i, diag[i]);
         }
	 else {
	   diag[i] = 1.0 / sqrt(ML_dabs(diag[i]));
	 }
      }
   }
   }
   else {
     for (i = 0; i < length; i++) diag[i] = 1.;
   }
#ifdef PRINTMAT
   fclose(fp);
#endif


   max_row_sum = -1.0E10;
   min_row_sum =  1.0E10;
   for (i=0; i<length; i++) u[i] = 1.0;
   ML_Operator_Apply(matrix, length, u, length, p);

   for (i=0; i<length; i++) p[i] *= diag[i];
   for (i=0; i<length; i++) 
   {
      if ( p[i] > max_row_sum ) max_row_sum = p[i];
      if ( p[i] < min_row_sum ) min_row_sum = p[i];
   }
   Nbc = ML_Comm_GsumInt(comm, Nbc);
   max_row_sum = ML_gmax_double(max_row_sum, comm);
   min_row_sum = - min_row_sum;
   min_row_sum = ML_gmax_double(min_row_sum, comm);
   min_row_sum = - min_row_sum;
   if ( comm->ML_mypid == 0 && print_freq > 0 )
   {
      printf("Max Row Sum = %e\n", max_row_sum);
      printf("Min Row Sum = %e\n", min_row_sum);
      if (Nbc != 0) printf("no. of BC's = %d\n", Nbc);
   }
   ML_free(colInd); ML_free(colVal);
   if ( totallength-Nbc < maxiter ) maxiter = totallength-Nbc;
   /* if maxiter is zero we will assume that we have a diagonal matrix? */
   if (maxiter == 0) {
     if (totallength == 0) {
       data->ML_eigen_max = 0.;
       data->ML_eigen_min = 0.;
     }
     else {
       data->ML_eigen_max = 1.;
       data->ML_eigen_min = 1.;
     }
     if (r  != NULL) ML_free(r);
     if (p  != NULL) ML_free(p);
     if (ap != NULL) ML_free(ap);
     return 1;
   }

   /* ----------------------------------------------------------------*/
   /* compute initial residual vector and norm */
   /* ----------------------------------------------------------------*/

   for ( i = 0; i < length; i++ ) {r[i] = rhs[i]; p[i] = 0.0;}
   res_norm = sqrt(ML_gdot(length, r, r, comm));
   rnorm_array[0] = res_norm;
   rho = 0.0;
   if ( res_norm == 0.0 ) 
   {
     if (Nbc != matrix->invec_leng) {
       printf("ML_CG_ComputeEigenvalues : cannot compute (res = 0).\n");
       data->ML_eigen_max = 1.0;
       data->ML_eigen_min = 1.0;
     }
     else {
       data->ML_eigen_max = max_row_sum;
       data->ML_eigen_min = min_row_sum;
     }
     if (r  != NULL) ML_free(r);
     if (p  != NULL) ML_free(p);
     if (ap != NULL) ML_free(ap);

     return 1;
   }

   /* ----------------------------------------------------------------*/
   /* loop until convergence is achieved or maxit has been exceeded   */
   /* ----------------------------------------------------------------*/

   for ( its = 0; its < maxiter; its++ )
   {
      rhom1 = rho;
      rho = ML_gdot(length, r, r, comm);
      if (its == 0) beta = 0.0;
      else 
      {
         beta = rho / rhom1;
         Tmat[its-1][its] = -beta;
      }
      for (i=0; i<length; i++) p[i] = beta * p[i] + r[i];
      for (i=0; i<length; i++) u[i] = p[i] * diag[i];
      ML_Operator_Apply(matrix, length, u, length, ap);
      for (i=0; i<length; i++) ap[i] = ap[i] * diag[i];
      sigma = ML_gdot(length, p, ap, comm);
      if ( fabs(sigma) < 1.0E-12 )
      {
	alpha_array[its] = sigma;
	maxiter = its + 1;
	break;
      }
      alpha  = rho / sigma;
      alpha_array[its] = sigma;
      for (i=0; i<length; i++) r[i] -= alpha * ap[i];
      res_norm = sqrt(ML_gdot(length, r, r, comm));
      rnorm_array[its+1] = res_norm;
      if ( rnorm_array[its+1] < 1.0E-8 * rnorm_array[0] )
      {
         maxiter = its + 1;
         break;
      }
   }
if (maxiter == 0) {
  for (i=0; i<length; i++) u[i] = p[i] * diag[i];
  ML_Operator_Apply(matrix, length, u, length, ap);
  for (i=0; i<length; i++) ap[i] = ap[i] * diag[i];
  sigma = ML_gdot(length, p, ap, comm);
  alpha_array[0] = sigma;
}

   /* ----------------------------------------------------------------*/
   /* construct T */
   /* ----------------------------------------------------------------*/

   Tmat[0][0] = alpha_array[0];
   for ( i = 1; i < maxiter; i++ )
   {
      Tmat[i][i]=alpha_array[i]+alpha_array[i-1]*Tmat[i-1][i]*Tmat[i-1][i];
   }
   for ( i = 0; i < maxiter; i++ )
   {
      Tmat[i][i+1] *= alpha_array[i];
      Tmat[i+1][i] = Tmat[i][i+1];
      rnorm_array[i] = 1.0 / rnorm_array[i];
   }
   for ( i = 0; i < maxiter; i++ )
     for ( j = 0; j < maxiter; j++ ) {
         Tmat[i][j] = Tmat[i][j] * rnorm_array[i] * rnorm_array[j];
     }
   
   /* ----------------------------------------------------------------*/
   /* diagonalize T using Jacobi iteration */
   /* ----------------------------------------------------------------*/

   offdiag_norm = 0.0;
   for ( i = 0; i < maxiter; i++ )
      for ( j = 0; j < i; j++ ) offdiag_norm += (Tmat[i][j] * Tmat[i][j]);
   offdiag_norm *= 2.0;
   init_offdiag_norm = offdiag_norm;

   while ( offdiag_norm > init_offdiag_norm * 1.0E-8 )
   {
      for ( i = 1; i < maxiter; i++ )
      {
         for ( j = 0; j < i; j++ )
         {
            apq = Tmat[i][j];
            if ( apq != 0.0 )
            {
               app = Tmat[j][j];
               aqq = Tmat[i][i];
               tau = ( aqq - app ) / (2.0 * apq);
               sign = (tau >= 0.0) ? 1.0 : -1.0;
               t  = sign / (tau * sign + sqrt(1.0 + tau * tau));
               c  = 1.0 / sqrt( 1.0 + t * t );
               s  = t * c;
               for ( k = 0; k < maxiter; k++ )
               {
                  arr = Tmat[j][k];
                  ass = Tmat[i][k];
                  Tmat[j][k] = c * arr - s * ass;
                  Tmat[i][k] = s * arr + c * ass;
               }
               for ( k = 0; k < maxiter; k++ )
               {
                  arr = Tmat[k][j];
                  ass = Tmat[k][i];
                  Tmat[k][j] = c * arr - s * ass;
                  Tmat[k][i] = s * arr + c * ass;
               }
            }
         }
      }
      offdiag_norm = 0.0;
      for ( i = 0; i < maxiter; i++ )
         for ( j = 0; j < i; j++ ) offdiag_norm += (Tmat[i][j] * Tmat[i][j]);
      offdiag_norm *= 2.0;
      /*for ( i = 0; i < maxiter; i++ )
           printf("%13.6e %13.6e %13.6e %13.6e %13.6e\n", Tmat[i][0],
                  Tmat[i][1], Tmat[i][2], Tmat[i][3], Tmat[i][4]);
      */
   }
   
   /* ----------------------------------------------------------------*/
   /* search for max eigenvalue*/
   /* ----------------------------------------------------------------*/
 
   t = Tmat[0][0];
   for ( i = 1; i < maxiter; i++ )
   {
      t = (Tmat[i][i] > t) ? Tmat[i][i] : t;
   }
   if ( comm->ML_mypid == 0 && print_freq > 0 )
      printf("max eigenvalue = %e\n", t);
   data->ML_eigen_max = t;
   t = Tmat[0][0];
   for ( i = 1; i < maxiter; i++ )
   {
      t = (Tmat[i][i] < t) ? Tmat[i][i] : t;
   }
   if ( comm->ML_mypid == 0 && print_freq > 0 )
      printf("min eigenvalue = %e\n", t);

   data->ML_eigen_min = t;   
   /* ----------------------------------------------------------------*/
   /* de-allocate storage for temporary vectors */
   /* ----------------------------------------------------------------*/

   if ( length > 0 )
   {
      ML_free(u); ML_free(r); ML_free(p); ML_free(ap);
      ML_free(rhs);
      ML_free(diag);
   }
   ML_free(alpha_array);
   ML_free(rnorm_array);
   ML_free(index_array);
   for ( i = 0; i <= original_maxiter; i++ ) ML_free(Tmat[i]);
   ML_free(Tmat);
   return 1;
}


/* ******************************************************************** */
/* ML_Power_ComputeEigenvalues                                             */
/*                                                                      */
/* NOTE: if scale_by_diag == ML_TRUE compute lambda (D^{1/2} A D^{1/2}) */
/*       if scale_by_diag == ML_FALSE compute  lambda ( A )             */
/* ******************************************************************** */

int ML_Power_ComputeEigenvalues(ML_Krylov *data, int length, int scale_by_diag)
{
   int         totallength, print_freq, nprocs, mypid, maxiter;
   int         i, j, ncnt, Nbc, *colInd = NULL, allocated;
   double      *p = NULL, *ap = NULL, *colVal = NULL, norm, *diag = NULL, sum;
   ML_Operator *matrix;
   ML_Comm     *comm;
   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* ----------------------------------------------------------------*/

   matrix      = ML_Krylov_Get_Amatrix(data);
   comm        = ML_Krylov_Get_Comm(data);
   totallength = ML_Comm_GsumInt(comm, length);
   print_freq  = ML_Krylov_Get_PrintFreq(data);
   nprocs      = comm->ML_nprocs;
   mypid       = comm->ML_mypid;
   maxiter     = 10;
   if ( totallength < maxiter ) maxiter = totallength;


   /* ----------------------------------------------------------------*/
   /* allocate temporary memory  */
   /* ----------------------------------------------------------------*/

   if ( length > 0 )
   {
     ap   = (double *) ML_allocate(length * sizeof(double));
     p    = (double *) ML_allocate(length * sizeof(double));
     diag = (double *) ML_allocate(length * sizeof(double));
      if ( diag == NULL )
      {
         printf("ML : ERROR in allocating memory.\n");
         exit(1);
      }
   }
   ML_random_vec(p, length, comm); 

   /* ----------------------------------------------------------------*/
   /* retrieve the diagonals */
   /* ----------------------------------------------------------------*/

   allocated = 100;
   colInd = (int    *) ML_allocate( allocated * sizeof(int) );
   colVal = (double *) ML_allocate( allocated * sizeof(double) );
   Nbc = 0;  /* rst to handle nonsymmetric (due to BCs) matrices */
   if (scale_by_diag) {
   for ( i = 0; i < length; i++ )
   {
      while(ML_Operator_Getrow(matrix,1,&i,allocated,colInd,colVal,&ncnt)==0)
      {
         allocated *= 2;
         ML_free(colInd); ML_free(colVal);
         colInd = (int    *) ML_allocate( allocated * sizeof(int) );
         colVal = (double *) ML_allocate( allocated * sizeof(double) );
      }

      sum = 0.0;
      for ( j = 0; j < ncnt; j++ )
      {
         if ( colInd[j] == i ) diag[i] = colVal[j];
         else sum += ML_dabs(colVal[j]);
      }
      /* kludging this in to handle Dirichlet BC's */
      if ( sum == 0.0) { p[i] = 0.; Nbc++; diag[i] = 1.;}

      else {
         if ( diag[i] == 0.0 ) 
         {
            printf("%d : diagonal[%d] == 0.0.\n", comm->ML_mypid, i);
	    /*            exit(1); */
	    diag[i] = 1.;
         }
         else if ( diag[i] < 0.0 )
         {
            printf("%d : diagonal[%d] = %e < 0.0.\n", comm->ML_mypid, i, diag[i]);
         }
	 else {
	   diag[i] = 1.0 / (ML_dabs(diag[i]));
	 }
      }
   }
   }
   else {
     for (i = 0; i < length; i++) diag[i] = 1.;
   }

   norm = sqrt(ML_gdot(length, p, p, comm));
   if (norm == 0.0) {
     data->ML_eigen_max = 1.;
     data->ML_eigen_min = 1.;
     if ( length > 0 ) {
      ML_free(ap); 
      ML_free(p);
      ML_free(diag);
     }
     return 1;
   }
   else norm = 1./norm;

   for (j = 0; j < length; j++) p[j] *= norm;
   for (i = 0; i < maxiter; i++) {
     ML_Operator_Apply(matrix, length, p, length, ap);
     for (j = 0; j < length; j++) ap[j] = ap[j]*diag[j];
     norm = 1./sqrt(ML_gdot(length, ap, ap, comm));
     for (j = 0; j < length; j++) p[j] = norm*ap[j];
   }
   data->ML_eigen_max = 1.0/norm;
   data->ML_eigen_min = 0.; /* no estimate for smallest eigenvalue */
   if ( length > 0 )    {
      ML_free(ap); 
      ML_free(p);
      ML_free(diag);
   }

   return 1;
}


