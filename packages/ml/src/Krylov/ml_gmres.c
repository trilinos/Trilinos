/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for GMRES solver                                           */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#include <math.h>
#include "ml_gmres.h"
#include "ml_utils.h"

/* ******************************************************************** */
/* ML_CG_Solve                                                          */
/* ******************************************************************** */

int ML_GMRES_Solve(ML_Krylov *data, int length, double *rhs, double *sol)
{
   int         i, j, k, k1, icnt, icnt2, its, maxiter, mlen, *indlist;
   int         print_freq;
   double      init_norm, eps1, tol, **ws, res_norm, rnorm2;
   double      **HH, *RS, *S, *C, ror, *darray, gam, epsmac=1.0e-10, t;
   void        *precon;
   int         (*precfcn)(void*,int,double*,int,double*);
   ML_Operator *matrix;
   ML_Comm     *comm;

   /* -----------------------------------------------------------------*/
   /* get all parameters from parent object*/
   /* -----------------------------------------------------------------*/

   maxiter    = ML_Krylov_Get_MaxIterations(data);
   tol        = ML_Krylov_Get_Tolerance(data);
   mlen       = ML_Krylov_Get_GMRESSize(data);
   precon     = ML_Krylov_Get_Precon(data);
   precfcn    = data->ML_precfcn;
   matrix     = ML_Krylov_Get_Amatrix(data);
   comm       = ML_Krylov_Get_Comm(data);
   print_freq = ML_Krylov_Get_PrintFreq(data);
   if ( mlen > maxiter ) mlen = maxiter;

   /* -----------------------------------------------------------------*/
   /* allocate temporary memory */
   /* -----------------------------------------------------------------*/

   darray = (double*) ML_allocate((mlen+1)*sizeof(double));
   HH = (double**) ML_allocate((mlen+2)*sizeof(double*));
   for (i=0; i<=mlen+1; i++) 
      HH[i] = (double*) ML_allocate((mlen+2)*sizeof(double));
   RS = (double*) ML_allocate((mlen+2)*sizeof(double));
   S  = (double*) ML_allocate((mlen+2)*sizeof(double));
   C  = (double*) ML_allocate((mlen+2)*sizeof(double));
   indlist  = (int*) ML_allocate((2*mlen+2)*sizeof(int));
   ws = (double**) ML_allocate((mlen+3)*sizeof(double*));
   for (i=0; i<=mlen+2; i++) 
      ws[i] = (double*) ML_allocate(length*sizeof(double));

   /* -----------------------------------------------------------------*/
   /* compute initial residual vector and norm */
   /* -----------------------------------------------------------------*/

   ML_Operator_Apply(matrix, length, sol, length, ws[0]);
   for ( i = 0; i < length; i++ ) ws[0][i] = rhs[i] - ws[0][i];
   res_norm = sqrt(ML_gdot(length, ws[0], ws[0], comm));
   init_norm = res_norm;
   if (comm->ML_mypid == 0 && print_freq < 1000 )
       printf("ML_GMRES initial residual norm = %e \n", init_norm);
   if ( init_norm == 0.0 ) 
   {
      for (i=0; i<=mlen+2; i++) ML_free(ws[i]);
      ML_free(ws);
      ML_free(darray);
      for (i=1; i<=mlen+1; i++) ML_free( HH[i] );
      ML_free(HH);
      ML_free(indlist);
      ML_free(RS);
      ML_free(S);
      ML_free(C);
      return 1;
   }

   /* -----------------------------------------------------------------*/
   /* initialization */
   /* -----------------------------------------------------------------*/

   its  = 0;
   eps1 = tol * init_norm;

   /* -----------------------------------------------------------------*/
   /* loop until convergence is achieved or maxit has been exceeded*/
   /* -----------------------------------------------------------------*/

   while (res_norm > eps1 && its < maxiter) 
   {
      ror = 1.0 / res_norm;
      for (i=0; i<length; i++) ws[0][i] *= ror;
      RS[1] = init_norm;
      icnt = 0; 
      rnorm2 = res_norm;
      while (icnt < mlen && (rnorm2/init_norm) > eps1) 
      {
         icnt++;
         its++;
         icnt2 = icnt + 1;
         if (precon != NULL) precfcn(precon,length,ws[icnt+1],length,ws[icnt-1]);
         else 
            for (i=0; i<length; i++) ws[icnt+1][i] = ws[icnt-1][i];
         ML_Operator_Apply(matrix, length, ws[icnt+1], length, ws[icnt]);
         for (j=1; j<=icnt; j++) 
         {
            darray[j-1] = ML_gdot(length, ws[j-1], ws[icnt2-1], comm);
            t = darray[j-1];
            HH[j][icnt] = t;  t = - t;
            for (i=0; i<length; i++) ws[icnt2-1][i] += (t*ws[j-1][i]);
         }
         t = sqrt(ML_gdot(length, ws[icnt2-1], ws[icnt2-1], comm));
         HH[icnt2][icnt] = t;
         if (t != 0.0) {
            t = 1.0 / t;
            for (i=0; i<length; i++) ws[icnt2-1][i] *= t;
         }       
         if (icnt != 1) {
            for (k=2; k<=icnt; k++) {
               k1 = k - 1;
               t = HH[k1][icnt];
               HH[k1][icnt] =  C[k1] * t + S[k1] * HH[k][icnt];
               HH[k][icnt]  = -S[k1] * t + C[k1] * HH[k][icnt];
            }
         }
         gam=sqrt(HH[icnt][icnt]*HH[icnt][icnt]+
                  HH[icnt2][icnt]*HH[icnt2][icnt]);
         if (gam == 0.0) gam = epsmac;
         C[icnt] = HH[icnt][icnt] / gam;
         S[icnt] = HH[icnt2][icnt] / gam;
         RS[icnt2] = -S[icnt] * RS[icnt];
         RS[icnt]  = C[icnt] * RS[icnt];
         HH[icnt][icnt] = C[icnt] * HH[icnt][icnt] + 
                          S[icnt] * HH[icnt2][icnt];
         rnorm2 = ML_dabs(RS[icnt2]);
         if (its % print_freq == 0 && comm->ML_mypid == 0)
            printf("ML_GMRES : iter %4d - res. norm = %e (%e)\n",its,
                    rnorm2, eps1);
      }  
      res_norm = rnorm2;
      RS[icnt] = RS[icnt] / HH[icnt][icnt];
      for (i=2; i<=icnt; i++) {
         k = icnt - i + 1;
         k1 = k + 1;
         t = RS[k];
         for (j=k1; j<=icnt; j++) t = t - HH[k][j] * RS[j];
         RS[k] = t / HH[k][k];
      }
      t = RS[1];
      for (i=0; i<length; i++) ws[0][i] *= t;
      for (j=2; j<=icnt; j++) 
      {
         t = RS[j];
         for (i=0; i<length; i++) ws[0][i] += (t * ws[j-1][i]);
      }
      if (precon != NULL) {
         precfcn(precon,length,ws[1],length,ws[0]);
         for (i=0; i<length; i++) ws[0][i] = ws[1][i];
      }
      for (i=0; i<length; i++) sol[i] += ws[0][i];
      ML_Operator_Apply(matrix, length, sol, length, ws[0]);
      for ( i = 0; i < length; i++ ) ws[0][i] = rhs[i] - ws[0][i];
      res_norm = sqrt(ML_gdot(length, ws[0], ws[0], comm));
      if (res_norm > eps1 && print_freq < 1000 && comm->ML_mypid == 0)
         printf("ML_GMRES : iter %4d - true res. norm = %e (%e)\n",its,
                    res_norm, eps1);
   }
   if (comm->ML_mypid == 0 && print_freq < 1000 )
      printf("ML_GMRES : total number of iterations = %d %e\n",its,res_norm);

  /* -----------------------------------------------------------------*/
  /* de-allocate storage for temporary vectors*/
  /* -----------------------------------------------------------------*/

   for (i=0; i<=mlen+2; i++) ML_free(ws[i]);
   ML_free(ws);
   ML_free(darray);
   for (i=1; i<=mlen+1; i++) ML_free( HH[i] );
   ML_free(HH);
   ML_free(indlist);
   ML_free(RS);
   ML_free(S);
   ML_free(C);
   return 1;
}

