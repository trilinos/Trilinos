/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Some tools for two grid analysis.                                    */
/* Comments use matlab notation.                                        */
/* -------------------------------------------------------------------- */
/* Author       : Jonathan Hu (SNL)                                     */
/* Date         : September, 2001                                       */
/* ML_gdot_H0                                                           */
/* ML_gdot_H1                                                           */
/* ML_gdot_H2                                                           */
/* ML_GetCoarseGridConst                                                */
/* ML_GetSmoothingConst                                                 */
/* ML_GetTwoLevelConvergenceFactor                                      */
/* ******************************************************************** */

#include <stdio.h>
#include "ml_twogrid_analysis.h"

/*******************************************************************************

 Calculate the H0 inner product of vectors u and v defined by

        <u,v>_H0 = <D \ u,v>_2.

*******************************************************************************/

double ML_gdot_H0(ML_Operator *Amat, double *vec1, double *vec2)
{
   double dtemp, *diagonal,*tempvec;
   int i;

   tempvec = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
   ML_allocate_check(tempvec);
   ML_Operator_Get_Diag(Amat, Amat->outvec_leng, &diagonal);
   for (i = 0; i < Amat->outvec_leng; i++)
      tempvec[i] = diagonal[i] * vec1[i];
   dtemp = ML_gdot(Amat->outvec_leng, tempvec, vec2, Amat->comm);
   ML_free(tempvec);

   return dtemp;
}

/*******************************************************************************

 Calculate the H1 (energy) inner product of vectors u and v defined by

        <u,v>_H1 = <Au,v>_2.

*******************************************************************************/

double ML_gdot_H1(ML_Operator *Amat, double *vec1, double *vec2)
{
   double dtemp, *tempvec;

   tempvec = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
   ML_allocate_check(tempvec);
   ML_Operator_Apply(Amat, Amat->invec_leng, vec1, Amat->outvec_leng, tempvec);
   dtemp = ML_gdot(Amat->outvec_leng, tempvec, vec2, Amat->comm);
   ML_free(tempvec);

   return dtemp;
}

/*******************************************************************************

 Calculate the H2 inner product of vectors u and v defined by

        <u,v>_H2 = <D \ Au,Av>_2.

*******************************************************************************/

double ML_gdot_H2(ML_Operator *Amat, double *vec1, double *vec2)
{
   double dtemp, *diagonal, *tempvec1, *tempvec2;
   int i;

   tempvec1 = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
   ML_allocate_check(tempvec1);
   ML_Operator_Apply(Amat, Amat->invec_leng, vec1, Amat->outvec_leng, tempvec1);
   ML_Operator_Get_Diag(Amat, Amat->outvec_leng, &diagonal);
   for (i = 0; i < Amat->outvec_leng; i++)
      tempvec1[i] = tempvec1[i] / diagonal[i];

   tempvec2 = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
   ML_allocate_check(tempvec2);
   ML_Operator_Apply(Amat, Amat->invec_leng, vec2, Amat->outvec_leng, tempvec2);

   dtemp = ML_gdot(Amat->outvec_leng, tempvec1, tempvec2, Amat->comm);

   ML_free(tempvec1);
   ML_free(tempvec2);

   return dtemp;
}

/*******************************************************************************

 Find the constant alpha associated with the smoothing property

        |Se|^2_H1  <=  |e|^2_H1  - alpha *|e|^2_H2.

 (See Introduction to Algebraic Multigrid, Christian Wagner, U. of Heidelberg.)

*******************************************************************************/

double ML_GetSmoothingConst(ML_Operator *Amat, double *err_h, ML_Smoother *sm)
{
   double alpha, eH1norm, eH2norm, smerrH1norm, *sm_err;
   int ntimes;

   /* Calculate ||e||^2_H1 and ||e||^2_H2. */
   eH1norm = fabs(ML_gdot_H1(Amat,err_h,err_h));
   eH2norm = fabs(ML_gdot_H2(Amat,err_h,err_h));

   /* Calculate ||Se||^2_H1 where S is the smoother. */

   sm_err = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
   ML_allocate_check(sm_err);

   /* Only do one iteration of the smoother. */
   ntimes = sm->ntimes;
   sm->ntimes = 1;
   ML_Smoother_Apply(sm, Amat->invec_leng, sm_err,
                     Amat->outvec_leng, err_h, ML_ZERO);
   sm->ntimes = ntimes;;
   smerrH1norm = fabs(ML_gdot_H1(Amat, sm_err, sm_err));
   ML_free(sm_err);

   alpha = (eH1norm - smerrH1norm) / eH2norm;

   return alpha;
}

/*******************************************************************************

 Find the constant beta associated with the coarse grid correction

        min(e_H) ||e - Pe_H ||^2_H0 <= beta * ||e||^2_H1.

 We use the fact that e_H = (P'P) \ P'e.

*******************************************************************************/

double ML_GetCoarseGridConst(ML_Operator *Amat, ML_Operator *Rmat,
                             ML_Operator *Pmat, double *err_h)
{
   ML_Operator  *RPmat;
   ML_Krylov    *kdata;
   double       *rhs, *err_H, *tempvec, dtemp, dtemp2, beta;
   int          i;

   /* Solve for e_H. */
   RPmat = ML_Operator_Create(Amat->comm);
   ML_2matmult(Rmat,Pmat,RPmat, ML_CSR_MATRIX);
   rhs = (double *) ML_allocate( Rmat->outvec_leng * sizeof(double) );
   ML_allocate_check(rhs);
   ML_Operator_Apply(Rmat, Rmat->invec_leng, err_h, Rmat->outvec_leng,rhs);
   err_H = (double *) ML_allocate( RPmat->invec_leng * sizeof(double) );
   ML_allocate_check(err_H);

   kdata = ML_Krylov_Create( RPmat->comm );
   ML_Krylov_Set_PrintFreq( kdata, 0 );
   ML_Krylov_Set_Amatrix(kdata, RPmat);
   ML_Krylov_Solve(kdata, RPmat->outvec_leng, rhs, err_H);
   ML_Krylov_Destroy(&kdata);

   /* Calculate ||e - Pe_H||^2_H0. */
   tempvec = (double *) ML_allocate( Pmat->outvec_leng * sizeof(double) );
   ML_allocate_check(tempvec);
   ML_Operator_Apply(Pmat, Pmat->invec_leng, err_H, Pmat->outvec_leng, tempvec);
   for (i = 0; i < Pmat->outvec_leng; i++)
      tempvec[i] = err_h[i] - tempvec[i];

   dtemp = fabs(ML_gdot_H0(Amat, tempvec, tempvec));

   /* Calculate ||e||^2_H1 and beta. */
   dtemp2 = fabs(ML_gdot_H1(Amat, err_h, err_h));
   beta = dtemp / dtemp2;

   ML_free(rhs);
   ML_free(tempvec);
   ML_free(err_H);

   return beta;
}

/*******************************************************************************

 Calculate two level convergence factor associated with post smoothing only.

*******************************************************************************/

double ML_GetTwoLevelConvergenceFactor(ML *ml,
                  /*                  ML_Operator *Amat, ML_Operator *Rmat,
                                    ML_Operator *Pmat, ML_Smoother *sm, */
                                    double *approx_soln, double *exact_soln)
{
   double *err_h, alpha, beta, conv_factor;
   int i;
   ML_Operator *Amat, *Rmat, *Pmat;
   ML_Smoother *sm;

   Amat = ml->Amat;
   Rmat = ml->Rmat;
   Pmat = ml->Pmat;
   sm   = ml->post_smoother;


   if (exact_soln == NULL)
      err_h = approx_soln;
   else
   {
      err_h = (double *) ML_allocate( Amat->outvec_leng * sizeof(double) );
      ML_allocate_check(err_h);
      for (i=0; i< Amat->outvec_leng; i++)
         err_h[i] = exact_soln[i] - approx_soln[i];
   }
      
   alpha =  ML_GetSmoothingConst(Amat, err_h, sm);
   beta = ML_GetCoarseGridConst(Amat, Rmat, Pmat, err_h);

   if (beta != 0)
   {
      conv_factor = sqrt(1 - alpha/beta);
   }
   else
   {
      printf("In ML_GetTwoLevelConvergenceFactor: Could not calculate convergence factor\n");
      conv_factor = -1.0;
   }
   return conv_factor;
}
