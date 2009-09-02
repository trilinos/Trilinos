/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the ML_Operator structure                              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#include <math.h>
#include "ml_bicgstabl.h"

/* ******************************************************************** */
/* ML_BICGSTABL_Solve                                                   */
/* ******************************************************************** */

int ML_BICGSTABL_Solve(ML_Krylov *data,int length,double *rhs,double *sol)
{
   int         offset, i, j, k, its, maxiter, converged=0, blen, print_freq;
   double      dtmp, omega, beta, rho, rho1, init_norm, eps1, tol;
   double      *r, *rh, *xh, *ut, *rt, *t, *tt, res_norm;
   double      alpha, *sigma, *gammap, *gammanp, darray[2];
   double      *gammapp, **mat, **tau, gamma;
   void        *precon;
   int         (*precfcn)(void*,int,double*,int,double*);
   ML_Operator *matrix;
   ML_Comm     *comm;

   /* -----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* -----------------------------------------------------------------*/

   maxiter    = ML_Krylov_Get_MaxIterations(data);
   tol        = ML_Krylov_Get_Tolerance(data);
   precon     = ML_Krylov_Get_Precon(data);
   precfcn    = data->ML_precfcn;
   matrix     = ML_Krylov_Get_Amatrix(data);
   comm       = ML_Krylov_Get_Comm(data);
   blen       = ML_Krylov_Get_BICGSTABLSize(data);
   print_freq = ML_Krylov_Get_PrintFreq(data);

   /* -----------------------------------------------------------------*/
   /* allocate temporary memory*/
   /* -----------------------------------------------------------------*/

   r  = (double *) ML_allocate( ((blen+2)*2+6) * length * sizeof(double));
   rh = r + length;
   xh = rh + length;
   t  = xh + length;
   tt = t + length;
   ut = tt + length;
   rt = ut + (blen+2)*length;

   /* -----------------------------------------------------------------*/
   /* compute initial residual vector and norm*/
   /* -----------------------------------------------------------------*/

   ML_Operator_Apply(matrix, length, sol, length, r);
   for ( i = 0; i < length; i++ ) r[i] = rhs[i] - r[i]; 
   res_norm = sqrt(ML_gdot(length, r, r, comm));
   init_norm = res_norm;
   if (comm->ML_mypid == 0 && print_freq < 1000 )
       printf("ML_BICGSTABL initial residual norm = %e \n", init_norm);
   if ( init_norm == 0.0 )
   {
      ML_free(r); 
      return 1;
   }

   /* -----------------------------------------------------------------*/
   /* initialization*/
   /* -----------------------------------------------------------------*/

   eps1 = tol * init_norm;
   sigma   = (double *) ML_allocate( sizeof(double) * (blen+1));
   gammap  = (double *) ML_allocate( sizeof(double) * (blen+1));
   gammanp = (double *) ML_allocate( sizeof(double) * (blen+1));
   gammapp = (double *) ML_allocate( sizeof(double) * (blen+1));
   mat     = (double **) ML_allocate( sizeof(double*) * (blen+1));
   tau     = (double **) ML_allocate( sizeof(double*) * (blen+1));
   for (i=1; i<=blen; i++) 
   {
      mat[i] = (double *) ML_allocate(sizeof(double)*(blen+1));
      tau[i] = (double *) ML_allocate(sizeof(double)*(blen+1));
   }

   /* -----------------------------------------------------------------*/
   /* loop until convergence is achieved or maxiter is exceeded */
   /* -----------------------------------------------------------------*/

   its = 0;
   while (converged == 0) 
   {
      for ( i = 0; i < length; i++ ) 
      {
         rh[i] = r[i]; 
         ut[i] = 0.0;
         xh[i] = sol[i];
         rt[i] = r[i]; 
      }
      omega = rho = 1.0; alpha = 0.0; 
      while (res_norm > eps1 && its < maxiter) {
         its = its + blen;
         for (i=0; i<length; i++)
         {
            ut[length+i] = ut[i];
            rt[length+i] = rt[i];
         }
         rho = -omega * rho;
         for (j=0; j<=blen-1; j++) {
            rho1 = ML_gdot(length, rh, rt+length*(j+1), comm);
            beta = alpha * rho1 / rho;
            rho = rho1;
            dtmp = -beta;
            for (k=0; k<=j; k++) 
            {
               offset = (k+1) * length;
               for (i=0; i<length; i++)
               {
                  ut[offset+i] = dtmp * ut[offset+i] + rt[offset+i];
               }
            }
            if ( precon != NULL ) 
               precfcn(precon,length,t, length, ut+(j+1)*length);
            else
               for (i=0; i<length; i++) t[i] = ut[(j+1)*length+i]; 
            ML_Operator_Apply(matrix, length, t, length, ut+(j+2)*length);
            gamma = ML_gdot(length, rh, ut+length*(j+2), comm);
            alpha = rho / gamma; dtmp = -alpha;
            for (k=0; k<=j; k++) 
            {
               offset = (k+1)*length;
               for (i=0; i<length; i++)
               {
                  rt[offset+i] += (dtmp * ut[offset+length+i]);
               }
            }
            if (precon != NULL) precfcn(precon,length,t,length,rt+(j+1)*length);
            else
               for (i=0; i<length; i++) t[i] = rt[(j+1)*length+i]; 
            ML_Operator_Apply(matrix, length, t, length, rt+(j+2)*length);
            for (i=0; i<length; i++) xh[i] += (alpha * ut[length+i]);
         }
         for (j=1; j<=blen; j++) 
            for (k=1; k<=blen; k++) mat[k][j] = 0.0; 
         for (j=1; j<=blen; j++) 
         { 
            for (k=1; k<=j-1; k++) {
               dtmp = ML_gdot(length, rt+(k+1)*length, rt+length*(j+1),comm);
               tau[k][j] = dtmp / sigma[k];
               mat[k][j] = tau[k][j] * sigma[k];
               dtmp = -tau[k][j];
               for (i=0; i<length; i++) 
                  rt[(j+1)*length+i] += (dtmp * rt[(k+1)*length+i]);
            }
            darray[0] = ML_gdot(length, rt+(j+1)*length,rt+(j+1)*length,comm);
            darray[1] = ML_gdot(length, rt+length,rt+(j+1)*length,comm);
            sigma[j] = darray[0];
            mat[j][j] = sigma[j];
            gammap[j] = darray[1] / sigma[j];
         }
         gammanp[blen] = gammap[blen];
         omega = gammanp[blen];
         for (j=blen-1; j>=1; j--) { 
            gammanp[j] = gammap[j];
            for (k=j+1; k<=blen; k++) 
               gammanp[j] = gammanp[j] - tau[j][k] * gammanp[k];
         }
         for (j=1; j<=blen-1; j++) { 
            gammapp[j] = gammanp[j+1];
            for (k=j+1; k<=blen-1; k++) 
               gammapp[j] = gammapp[j] + tau[j][k] * gammanp[k+1];
         }
         dtmp = gammanp[1];
         for (i=0; i<length; i++) xh[i] += (dtmp * rt[length+i]);
         dtmp = - gammap[blen];
         for (i=0; i<length; i++) 
            rt[length+i] += (dtmp * rt[(1+blen)*length+i]);
         dtmp = - gammanp[blen];
         for (i=0; i<length; i++) 
            ut[length+i] += (dtmp * ut[(1+blen)*length+i]);
         for (j=1; j<=blen-1; j++) { 
            dtmp = - gammanp[j]; 
            for (i=0; i<length; i++) 
               ut[length+i] += (dtmp * ut[(1+j)*length+i]);
            dtmp = gammapp[j]; 
            for (i=0; i<length; i++) xh[i] += (dtmp * rt[(1+j)*length+i]);
            dtmp = - gammap[j];  
            for (i=0; i<length; i++) 
               rt[length+i] += (dtmp * rt[(1+j)*length+i]);
         }
         for (i=0; i<length; i++)
         {
            ut[i]  = ut[length+i];
            rt[i]  = rt[length+i];
            sol[i] = xh[i];
         }
         res_norm = sqrt(ML_gdot(length, rt+length, rt+length, comm));
         if (comm->ML_mypid == 0 && its % print_freq == 0 )
            printf(" BiCGSTAB(L) : iter %4d - res. norm = %e \n", its, res_norm);
      }
      if ( precon != NULL ) 
      {
         precfcn(precon,length,t,length,sol);
         for (i=0; i<length; i++) sol[i] = t[i];
      }

      /* ---------------------------------------------------------------*/
      /* compute final residual vector and norm*/
      /* ---------------------------------------------------------------*/

      ML_Operator_Apply(matrix, length, sol, length, r);
      for ( i = 0; i < length; i++ ) r[i] = rhs[i] - r[i]; 
      res_norm = sqrt(ML_gdot(length, r, r, comm));

      if (comm->ML_mypid == 0 && print_freq < 1000 )
         printf("ML_BiCGSTAB(L) final residual norm = %e \n",res_norm);
      if (res_norm < eps1 || its >= maxiter) converged = 1;
   }
   if (comm->ML_mypid == 0 && print_freq < 1000 )
      printf("ML_BiCGSTAB(L) : total number of iterations = %d \n", its);

  /* -----------------------------------------------------------------*/
  /* de-allocate storage for temporary vectors*/
  /* -----------------------------------------------------------------*/

   ML_free(r);
   for (i=1; i<=blen; i++) {
      ML_free( mat[i] );
      ML_free( tau[i] );
   }
   ML_free( gammap );
   ML_free( gammanp );
   ML_free( gammapp );
   ML_free( mat );
   ML_free( tau );
   ML_free( sigma );
   return 1;
}

