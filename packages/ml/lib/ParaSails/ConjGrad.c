/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * ConjGrad - Preconditioned conjugate gradient algorithm using the
 * ParaSails preconditioner.
 *
 *****************************************************************************/

#include "math.h"
#include "Common.h"
#include "Matrix.h"
#include "ParaSails.h"

double ddot_(int *, double *, int *, double *, int *);
void dcopy_(int *, double *, int *, double *, int *);
void dscal_(int *, double *, double *, int *);
void daxpy_(int *, double *, double *, int *, double *, int *);

static double InnerProd(int n, double *x, double *y, MPI_Comm comm)
{
    double local_result, result;

    int one = 1;
    local_result = ddot_(&n, x, &one, y, &one);

    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

    return result;
}

static void CopyVector(int n, double *x, double *y)
{
    int one = 1;
    dcopy_(&n, x, &one, y, &one);
}

static void ScaleVector(int n, double alpha, double *x)
{
    int one = 1;
    dscal_(&n, &alpha, x, &one);
}

static void Axpy(int n, double alpha, double *x, double *y)
{
    int one = 1;
    daxpy_(&n, &alpha, x, &one, y, &one);
}


/*--------------------------------------------------------------------------
 * PCG_ParaSails - PCG solver using ParaSails.
 * Use NULL for ps if to get unpreconditioned solve.
 * Solver will stop at step 500 if rel. resid. norm reduction is not less 
 * than 0.1 at that point.
 *--------------------------------------------------------------------------*/

void PCG_ParaSails(Matrix *mat, ParaSails *ps, double *b, double *x,
   double tol, int max_iter)
{
   double *p, *s, *r;
   double alpha, beta;
   double gamma, gamma_old;
   double bi_prod, i_prod, eps;
   int i = 0;
   int mype;

   /* local problem size */
   int n = mat->end_row - mat->beg_row + 1;

   MPI_Comm comm = mat->comm;
   MPI_Comm_rank(comm, &mype);

   /* compute square of absolute stopping threshold  */
   /* bi_prod = <b,b> */
   bi_prod = InnerProd(n, b, b, comm);
   eps = (tol*tol)*bi_prod;

   /* Check to see if the rhs vector b is zero */
   if (bi_prod == 0.0)
   {
      /* Set x equal to zero and return */
      CopyVector(n, b, x);
      return;
   }

   p = (double *) malloc(n * sizeof(double));
   s = (double *) malloc(n * sizeof(double));
   r = (double *) malloc(n * sizeof(double));

   /* r = b - Ax */
   MatrixMatvec(mat, x, r);  /* r = Ax */
   ScaleVector(n, -1.0, r);  /* r = -r */
   Axpy(n, 1.0, b, r);       /* r = r + b */
 
   /* p = C*r */
   if (ps != NULL)
      ParaSailsApply(ps, r, p);
   else
      CopyVector(n, r, p);

   /* gamma = <r,p> */
   gamma = InnerProd(n, r, p, comm);

   while ((i+1) <= max_iter)
   {
      i++;

      /* s = A*p */
      MatrixMatvec(mat, p, s);

      /* alpha = gamma / <s,p> */
      alpha = gamma / InnerProd(n, s, p, comm);

      gamma_old = gamma;

      /* x = x + alpha*p */
      Axpy(n, alpha, p, x);

      /* r = r - alpha*s */
      Axpy(n, -alpha, s, r);
         
      /* s = C*r */
      if (ps != NULL)
         ParaSailsApply(ps, r, s);
      else
         CopyVector(n, r, s);

      /* gamma = <r,s> */
      gamma = InnerProd(n, r, s, comm);

      /* set i_prod for convergence test */
      i_prod = InnerProd(n, r, r, comm);

#ifdef PARASAILS_CG_PRINT
      if (mype == 0 && i % 100 == 0)
         printf("Iter (%d): rel. resid. norm: %e\n", i, sqrt(i_prod/bi_prod));
#endif

      /* check for convergence */
      if (i_prod < eps)
         break;

      /* non-convergence test */
      if (i >= 1000 && i_prod/bi_prod > 0.01)
      {
         if (mype == 0)
            printf("Aborting solve due to slow or no convergence.\n");
         break;
      }
 
      /* beta = gamma / gamma_old */
      beta = gamma / gamma_old;

      /* p = s + beta p */
      ScaleVector(n, beta, p);   
      Axpy(n, 1.0, s, p);
   }

   free(p);
   free(s);

   /* compute exact relative residual norm */
   MatrixMatvec(mat, x, r);  /* r = Ax */
   ScaleVector(n, -1.0, r);  /* r = -r */
   Axpy(n, 1.0, b, r);       /* r = r + b */
   i_prod = InnerProd(n, r, r, comm);

   free(r);

   if (mype == 0)
      printf("Iter (%4d): computed rrn    : %e\n", i, sqrt(i_prod/bi_prod));
}
