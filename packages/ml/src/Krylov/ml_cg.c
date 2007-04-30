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
#include "ml_lapack.h"


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
   double      *u = NULL, *p = NULL, *ap = NULL, *r = NULL, *scale = NULL;
   int         totallength, print_freq, maxiter, info, one = 1, allocated;
   int         NegDefinite = 0, Nbc =0, NegDiag = 0, Nignored = 0;
   double      *Tdiag = NULL, *Tsubdiag = NULL, *colVal = NULL;
   double      *alpha_array = NULL, *rnorm_array = NULL;
   double      rho, rhom1, alpha, beta, sigma, sum;
   int         i, j, *colInd = NULL, ncnt;
   ML_Operator *matrix = NULL;
   ML_Comm     *comm = NULL;
   char        EigsOnly[2];

   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* ----------------------------------------------------------------*/
   
   matrix      = ML_Krylov_Get_Amatrix(data);
   comm        = ML_Krylov_Get_Comm(data);
   totallength = ML_Comm_GsumInt(comm, length);
   print_freq  = ML_Krylov_Get_PrintFreq(data);
   maxiter     = ML_Krylov_Get_MaxIterations(data);

   if (totallength == 0) {
       data->ML_eigen_max = 0.;
       data->ML_eigen_min = 0.;
       return 1;
   }
   if ( (totallength == 1) && (scale_by_diag)) {
       data->ML_eigen_max = 1.;
       data->ML_eigen_min = 1.;
       return 1;
   }

   if ( length > 0 ) {
      r     = (double *) ML_allocate(length * sizeof(double));
      scale = (double *) ML_allocate(length * sizeof(double));
      if ( scale == NULL ) {
         printf("ML : ERROR in allocating memory.\n");
         exit(1);
      }
   }
   ML_random_vec(r, length, comm); 

   /* ----------------------------------------------------------------*/
   /* retrieve diagonals                                              */
   /* ----------------------------------------------------------------*/

   if (scale_by_diag) {
     allocated = 100;
     colInd = (int    *) ML_allocate( allocated * sizeof(int) );
     colVal = (double *) ML_allocate( allocated * sizeof(double) );
     for ( i = 0; i < length; i++ ) {
        while(ML_Operator_Getrow(matrix,1,&i,allocated,colInd,colVal,&ncnt)==0){
           allocated *= 2;
           ML_free(colInd); ML_free(colVal);
           colInd = (int    *) ML_allocate( allocated * sizeof(int) );
           colVal = (double *) ML_allocate( allocated * sizeof(double) );
        }
  
        sum = 0.0;
        for ( j = 0; j < ncnt; j++ ) {
           if ( colInd[j] == i ) scale[i]  = colVal[j];
           else sum += ML_dabs(colVal[j]);
        }

        if ( sum == 0.0) { /* kludging this to handle Dirichlet BC's */
           scale[i] = 0.;  r[i] = 0; Nbc++;
        } 
        else {
           if ( scale[i] > 0.0 ) scale[i]  = 1.0 / sqrt(scale[i]);
           else {  /* Something is strange. The diagonal should not be 0 or*/
                   /* negative for SPD matrices. Set scale = -sqrt(-scale).*/
                   /*                                                      */
                   /* Later check                                          */
                   /*   if all scales negative (negative definite)         */
                   /*        Take D to be negative diagonal, compute       */
                   /*        eigenvalue, and then flip sign.               */
                   /*   else                                               */
                   /*      Ignore this row and column by setting scaling   */
                   /*      to zero.                                        */
               if (scale[i] != 0.) scale[i]  = -1.0 / sqrt(-scale[i]);
           }
        } /* if ( sum == 0.0) */
     }
     ML_free(colInd); ML_free(colVal);

     /* Now see how many diagonals are negative and zero */

     NegDiag = 0; Nignored = 0;
     for (i = 0; i < length; i++) {
        if (scale[i]  < 0.) NegDiag++;
        if (scale[i] == 0.) Nignored++;
     }
     NegDiag  = ML_Comm_GsumInt(comm, NegDiag);
     Nignored = ML_Comm_GsumInt(comm, Nignored);
     if (NegDiag+Nignored== totallength) {
         if ( NegDiag != 0) { /* matrix is negative definite */
             NegDefinite = 1;
             for (i = 0; i < length; i++) scale[i] = -scale[i];
         }
         else { /* matrix has no rows and columns that are not ignored */
           //           data->ML_eigen_max = 0.;//cms
           //           data->ML_eigen_min = 0.;//cms
           data->ML_eigen_max = 1.;
           data->ML_eigen_min = 1.;
           if (scale != NULL) ML_free(scale);
            if (r != NULL) ML_free(r);
            return 1;
         }
     }
     else {  /* ignore rows & columns with negative diagonals */
        for (i = 0; i < length; i++) if (scale[i] < 0.) {scale[i] = 0.;r[i]=0.;}
        Nignored += NegDiag;
        if ((NegDiag > 0) && (comm->ML_mypid == 0) && (print_freq > 0)) {
           printf("%d diagonals are negative while computing\n",NegDiag);
           printf("eigenvalues. This should not occur for SPD\n");
           printf("systems. Is the original matrix being solved SPD?\n");
           printf("If so, send this example to the ml-users list so we\n");
           printf("figure out what is wrong.\n");
        }
     }
   }
   else {
       /* Try and detect Dirichlet BC's by doing a matrix vector product */

      ML_Operator_Apply(matrix, length, r, length, scale);
      Nbc = 0;
      for (i = 0; i < length; i++) {
         if (scale[i] == r[i] ) {
            Nbc++;  Nignored++; scale[i] = 0.; r[i] = 0.;
         }
         else if (scale[i] == 0) { Nignored++; r[i] = 0.;}
         else scale[i] = 1.;
      }
      Nignored = ML_Comm_GsumInt(comm, Nignored);
   }

   Nbc     = ML_Comm_GsumInt(comm, Nbc);
   if ( (comm->ML_mypid == 0) && (print_freq > 0) && (Nbc != 0) )
      printf("no. of BC's = %d\n", Nbc);

   if ( maxiter > totallength - Nignored)  maxiter = totallength - Nignored;


   /* All of the equations are boundary conditions and we have scaled */
   /* by the matrix diagonal => all the eigenvalues are 1.            */

   if (totallength == Nbc) {
     data->ML_eigen_max = 1.;
     data->ML_eigen_min = 1.;
     if (scale != NULL) ML_free(scale);
     if (r != NULL) ML_free(r);
     return 1;
   }

   rho = 0.0;

   alpha_array = (double *) ML_allocate((maxiter+1) * sizeof(double));
   rnorm_array = (double *) ML_allocate((maxiter+1) * sizeof(double));
   if ( length > 0 ) {
      u    = (double *) ML_allocate(length * sizeof(double));
      p    = (double *) ML_allocate(length * sizeof(double));
      ap   = (double *) ML_allocate(length * sizeof(double));
      if ( ap == NULL ) {
         printf("ML : ERROR in allocating memory.\n");
         exit(1);
      }
   }
   rnorm_array[0] = sqrt(ML_gdot(length, r, r, comm));
   for ( i = 0; i < length; i++ ) {p[i] = 0.0;}

   /* ----------------------------------------------------------------*/
   /* loop until convergence is achieved or maxit has been exceeded   */
   /* ----------------------------------------------------------------*/
   Tdiag       = (double *) ML_allocate((maxiter+1) * sizeof(double));
   Tsubdiag    = (double *) ML_allocate((maxiter+1) * sizeof(double));
   for ( i = 0; i <= maxiter; i++ ) { 
      Tdiag[i] = 1.; Tsubdiag[i] = 0.;
   }

   for ( i = 0; i < maxiter; i++ )
   {
      rhom1 = rho;
      rho = ML_gdot(length, r, r, comm);
      if (i == 0) beta = 0.0;
      else 
      {
         beta = rho / rhom1;
         Tsubdiag[i-1] = -beta;
      }
      for (j=0; j<length; j++) p[j] = beta * p[j] + r[j];
      for (j=0; j<length; j++) u[j] = p[j] * scale[j];
      ML_Operator_Apply(matrix, length, u, length, ap);
      for (j=0; j<length; j++) ap[j] = ap[j] * scale[j];
      sigma = ML_gdot(length, p, ap, comm);
      if ( fabs(sigma) < 1.0E-12 )
      {
    alpha_array[i] = sigma;
    maxiter = i + 1;
    break;
      }
      alpha  = rho / sigma;
      alpha_array[i] = sigma;
      for (j=0; j<length; j++) r[j] -= alpha * ap[j];
      rnorm_array[i+1] = sqrt(ML_gdot(length, r, r, comm));
      if ( rnorm_array[i+1] < 1.0E-8 * rnorm_array[0] )
      {
         maxiter = i + 1;
         break;
      }
   }

   /* ----------------------------------------------------------------*/
   /* construct T */
   /* ----------------------------------------------------------------*/



   Tdiag[0] = alpha_array[0];
   for ( i = 1; i < maxiter; i++ )
   {
      Tdiag[i] = alpha_array[i]+alpha_array[i-1]*Tsubdiag[i-1]*Tsubdiag[i-1];
   }
   for ( i = 0; i < maxiter; i++ )
   {
      Tsubdiag[i] *= alpha_array[i];
      rnorm_array[i] = 1.0 / rnorm_array[i];
   }
   for ( i = 0; i < maxiter; i++ ) {
     Tdiag[i] = Tdiag[i] * rnorm_array[i] * rnorm_array[i];
     if (i != maxiter-1) 
        Tsubdiag[i] = Tsubdiag[i] * rnorm_array[i] * rnorm_array[i+1];
   }
   
   strcpy(EigsOnly,"N");
   DSTEQR_F77(EigsOnly,&maxiter,Tdiag,Tsubdiag, NULL, &one, NULL, &info);
   if (NegDefinite) for (i = 0; i < maxiter; i++) Tdiag[i] = -Tdiag[i];

   if (maxiter > 0) {
      data->ML_eigen_max = Tdiag[0];
      data->ML_eigen_min = Tdiag[0];
   }
   else {
      data->ML_eigen_max = 0.;
      data->ML_eigen_min = 0.;
   }
   for (i = 1; i < maxiter; i++) {
      if (Tdiag[i] > data->ML_eigen_max) data->ML_eigen_max = Tdiag[i];
      if (Tdiag[i] < data->ML_eigen_min) data->ML_eigen_min = Tdiag[i];
   }
   if ( comm->ML_mypid == 0 && print_freq > 0 ) {
      printf("max eigenvalue = %e\n", data->ML_eigen_max);
      printf("min eigenvalue = %e\n", data->ML_eigen_min);
   }
   /* ----------------------------------------------------------------*/
   /* de-allocate storage for temporary vectors */
   /* ----------------------------------------------------------------*/

   if (u           != NULL) ML_free(u);
   if (p           != NULL) ML_free(p);
   if (ap          != NULL) ML_free(ap);
   if (r           != NULL) ML_free(r);
   if (scale       != NULL) ML_free(scale);
   if (alpha_array != NULL) ML_free(alpha_array);
   if (rnorm_array != NULL) ML_free(rnorm_array);
   if (Tdiag       != NULL) ML_free(Tdiag);
   if (Tsubdiag    != NULL) ML_free(Tsubdiag);

   return 1;
}


/* ******************************************************************** */
/* ML_SubspaceIteration_ComputeEigenvalues                              */
/*      This is essentially the power method to compute eigenvalues     */
/*      with subspace dimension = 2. The initial guess for the first    */
/*      vector is real and the initial guess for the 2nd vector is      */
/*      pure imaginary. I don't know how good this routine really is    */
/*      but it appears a lot better than just doing the straight        */
/*      power method on a nonsymmetric system.                          */
/*                                                                      */
/* NOTE: if scale_by_diag == ML_TRUE compute lambda (D^{1/2} A D^{1/2}) */
/*       if scale_by_diag == ML_FALSE compute  lambda ( A )             */
/* ******************************************************************** */

int ML_SubspaceIteration_ComputeEigenvalues(ML_Krylov *data, int length, int scale_by_diag)
{
   int         totallength, maxiter;
   int         i, j, ncnt, Nbc, *colInd = NULL, allocated, level;
   double      *colVal = NULL, *diag = NULL, sum;
   ML_Operator *matrix;
   ML_Comm     *comm;
   double      *v1real=NULL,*v2imag=NULL,*y1real=NULL,*y2imag=NULL;
   double      norm1, norm2, alpha, Toneone, Ttwotwo, b, c;

   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* ----------------------------------------------------------------*/

   matrix      = ML_Krylov_Get_Amatrix(data);
   level = -1;
   if (matrix->to != NULL) level = matrix->to->levelnum;

   comm        = ML_Krylov_Get_Comm(data);
   totallength = ML_Comm_GsumInt(comm, length);
   maxiter     = ML_Krylov_Get_MaxIterations(data);
   if ( totallength < maxiter ) maxiter = totallength;

   if (totallength == 0) {
       data->ML_eigen_max = 0.;
       data->ML_eigen_min = 0.;
       return 1;
   }
   if ( (totallength == 1) && (scale_by_diag)) {
       data->ML_eigen_max = 1.;
       data->ML_eigen_min = 1.;
       return 1;
   }


   /* ----------------------------------------------------------------*/
   /* allocate temporary memory  */
   /* ----------------------------------------------------------------*/

   if ( length > 0 )
   {
     v1real   = (double *) ML_allocate(length * sizeof(double));
     v2imag   = (double *) ML_allocate(length * sizeof(double));
     y1real   = (double *) ML_allocate(length * sizeof(double));
     y2imag   = (double *) ML_allocate(length * sizeof(double));
     diag = (double *) ML_allocate(length * sizeof(double));
      if ( diag == NULL )
      {
         printf("ML : ERROR in allocating memory.\n");
         exit(1);
      }
      for (i = 0; i < length; i++) diag[i] = 1.;
   }
   ML_random_vec(v1real, length, comm); 
   ML_random_vec(v2imag, length, comm); 

   /* ----------------------------------------------------------------*/
   /* retrieve the diagonals */
   /* ----------------------------------------------------------------*/

   allocated = 100;
   colInd = (int    *) ML_allocate( allocated * sizeof(int) );
   colVal = (double *) ML_allocate( allocated * sizeof(double) );
   Nbc = 0;  /* rst to handle nonsymmetric (due to BCs) matrices */
   if (scale_by_diag)
   {
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
        for ( j = 0; j < ncnt; j++ ) {
           if ( colInd[j] == i ) diag[i] = colVal[j];
           else sum += ML_dabs(colVal[j]);
        }
        /* kludging this in to handle Dirichlet BC's */
        if ( sum == 0.0) {
          v1real[i] = 0.; v2imag[i] = 0.; Nbc++; diag[i] = 1.;
        } else {
           if ( diag[i] == 0.0 ) {
             if (ML_Get_PrintLevel() > 0) {
               if (level != -1)
                 printf("%d : diagonal[%d] == 0.0 for matrix stored on level %d within MG hierarchy\n", comm->ML_mypid, i, level);
               else
                 printf("%d : diagonal[%d] == 0.0\n", comm->ML_mypid, i);
             }
             diag[i] = 1.;
           }
           /* MS * added on 01-Mar-06 */
           else
             diag[i] = 1.0 / diag[i];
        } /*if ( sum == 0.0) */
     } /*for ( i = 0; i < length; i++ )*/
   }
   else {
     for (i = 0; i < length; i++) diag[i] = 1.;
   }

   ML_free(colInd);
   ML_free(colVal);

   norm1 = sqrt(ML_gdot(length, v1real, v1real, comm));
   norm2 = sqrt(ML_gdot(length, v2imag, v2imag, comm));
   if ( (norm1 == 0.0) || (norm2 == 0.0)) {
     data->ML_eigen_max = 1.;
     data->ML_eigen_min = 1.;
     if ( diag   != NULL) ML_free(diag);
     if ( v1real != NULL ) ML_free(v1real);
     if ( y1real != NULL ) ML_free(y1real);
     if ( v2imag != NULL ) ML_free(v2imag);
     if ( y2imag != NULL ) ML_free(y2imag);
     return 1;
   }
   norm1 = 1./norm1;
   norm2 = 1./norm2;
   for (j = 0; j < length; j++) v1real[j] *= norm1;
   for (j = 0; j < length; j++) v2imag[j] *= norm2;

   for (i = 0; i < maxiter; i++) {
     ML_Operator_Apply(matrix, length, v1real, length, y1real);
     for (j = 0; j < length; j++) y1real[j] *= diag[j];
     ML_Operator_Apply(matrix, length, v2imag, length, y2imag);
     for (j = 0; j < length; j++) y2imag[j] *= diag[j];
     alpha = sqrt(ML_gdot(length, y1real, y1real, comm));
     for (j = 0; j < length; j++) v1real[j] = y1real[j]/alpha;
     alpha = ML_gdot(length, v1real, y2imag, comm);
     for (j = 0; j < length; j++) v2imag[j] = y2imag[j] - alpha*v1real[j];
     alpha = sqrt(ML_gdot(length, v2imag, v2imag, comm));
     for (j = 0; j < length; j++) v2imag[j] = v2imag[j]/alpha;
   }
   ML_Operator_Apply(matrix, length, v1real, length, y1real);
   for (j = 0; j < length; j++) y1real[j] *= diag[j];
   ML_Operator_Apply(matrix, length, v2imag, length, y2imag);
   for (j = 0; j < length; j++) y2imag[j] *= diag[j];
   Toneone = ML_gdot(length, v1real, y1real, comm);
   Ttwotwo = ML_gdot(length, v2imag, y2imag, comm);
   b = -( Toneone + Ttwotwo);
   c = Toneone*Ttwotwo - ML_gdot(length, v1real, y2imag, comm)*
       ML_gdot(length, v2imag, y1real, comm);
   if (  b*b-4.*c > 0) {
      data->ML_eigen_max = (ML_abs(b) + sqrt(b*b - 4.*c))/2;
   }
   else {
      data->ML_eigen_max = sqrt(b*b + ML_abs(b*b - 4.*c))/2.;
   }

   data->ML_eigen_min = 0.; /* no estimate for smallest eigenvalue */
   if ( diag   != NULL) ML_free(diag);
   if ( v1real != NULL ) ML_free(v1real);
   if ( y1real != NULL ) ML_free(y1real);
   if ( v2imag != NULL ) ML_free(v2imag);
   if ( y2imag != NULL ) ML_free(y2imag);

   return 1;
}

/* ******************************************************************** */
/* ML_Power_ComputeEigenvalues                                          */
/*   NOTE: It is BETTER to use ML_SubspaceIteration_ComputeEigenvalues  */
/*         for nonsymmetric systems. This keeps a subspace dimension of */
/*         2 and better handles complex conjugates.                     */
/*                                                                      */
/* NOTE: if scale_by_diag == ML_TRUE compute lambda (D^{1/2} A D^{1/2}) */
/*       if scale_by_diag == ML_FALSE compute  lambda ( A )             */
/* ******************************************************************** */

int ML_Power_ComputeEigenvalues(ML_Krylov *data, int length, int scale_by_diag)
{
   int         totallength, maxiter;
   int         i, j, ncnt, Nbc, *colInd = NULL, allocated, level;
   double      *p = NULL, *ap = NULL, *colVal = NULL, norm, *diag = NULL, sum;
   ML_Operator *matrix;
   ML_Comm     *comm;
   /* ----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* ----------------------------------------------------------------*/

   matrix      = ML_Krylov_Get_Amatrix(data);
   level = -1;
   if (matrix->to != NULL) level = matrix->to->levelnum;

   comm        = ML_Krylov_Get_Comm(data);
   totallength = ML_Comm_GsumInt(comm, length);
   maxiter     = ML_Krylov_Get_MaxIterations(data);
   if ( totallength < maxiter ) maxiter = totallength;

   if (totallength == 0) {
       data->ML_eigen_max = 0.;
       data->ML_eigen_min = 0.;
       return 1;
   }
   if ( (totallength == 1) && (scale_by_diag)) {
       data->ML_eigen_max = 1.;
       data->ML_eigen_min = 1.;
       return 1;
   }

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
   if (scale_by_diag)
   {
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
        for ( j = 0; j < ncnt; j++ ) {
           if ( colInd[j] == i ) diag[i] = colVal[j];
           else sum += ML_dabs(colVal[j]);
        }
        /* kludging this in to handle Dirichlet BC's */
        if ( sum == 0.0) {
          p[i] = 0.; Nbc++; diag[i] = 1.;
        } else {
           if ( diag[i] == 0.0 ) {
             if (ML_Get_PrintLevel() > 0) {
               if (level != -1)
                 printf("%d : diagonal[%d] == 0.0 for matrix stored on level %d within MG hierarchy\n", comm->ML_mypid, i, level);
               else
                 printf("%d : diagonal[%d] == 0.0\n", comm->ML_mypid, i);
             }
             diag[i] = 1.;
           }
           /* MS * added on 01-Mar-06 */
           else
             diag[i] = 1.0 / diag[i];
#if 0
           else if ( diag[i] < 0.0 ) {
             if (ML_Get_PrintLevel() > 0) {
               if (level != -1) 
                 printf("%d : diagonal[%d] = %e < 0 for matrix stored on level %d within MG hierarchy\n", comm->ML_mypid, i, diag[i], level);
               else
                 printf("%d : diagonal[%d] = %e < 0.0.\n", comm->ML_mypid, i, diag[i]);
               }
           }
           else
             diag[i] = 1.0 / (ML_dabs(diag[i]));
#endif
        } /*if ( sum == 0.0) */
     } /*for ( i = 0; i < length; i++ )*/
   }
   else {
     for (i = 0; i < length; i++) diag[i] = 1.;
   }

   ML_free(colInd);
   ML_free(colVal);

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


