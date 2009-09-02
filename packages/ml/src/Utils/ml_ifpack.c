/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/*       User Interface Functions                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS)
#include "ml_struct.h"
#include "ml_lapack.h"

/* ------------------------------------------------------------------------- */
/* generate the Ifpack smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
		       double rhs[])
{
  ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
  void *Ifpack_Handle = smooth_ptr->smoother->data;
  double* x2 = NULL,* rhs2 = NULL;
  /*int i;*/
  int n, kk;
  int one_int = 1;
  double minus_one_double = -1.0;

  if (sm->init_guess == ML_NONZERO) 
  {
    n = sm->my_level->Amat->invec_leng;
    assert (n == sm->my_level->Amat->outvec_leng);

    rhs2 = (double*) ML_allocate(sizeof(double) * (n + 1));
    x2   = (double*) ML_allocate(sizeof(double) * (n + 1));

    ML_Operator_Apply(sm->my_level->Amat, n, x, n, rhs2);
    DCOPY_F77(&n, x, &one_int, x2, &one_int);
    DAXPY_F77(&n, &minus_one_double, rhs, &one_int, rhs2, &one_int);
    ML_Ifpack_Solve(Ifpack_Handle, x2, rhs2);
    DAXPY_F77(&n, &minus_one_double, x2, &one_int, x, &one_int);

    ML_free(rhs2);
    ML_free(x2);
  }
  else
    ML_Ifpack_Solve(Ifpack_Handle, x, rhs);

  for (kk = 1; kk < sm->ntimes; kk++) {
    n = sm->my_level->Amat->invec_leng;
    assert (n == sm->my_level->Amat->outvec_leng);

    rhs2 = (double*) malloc(sizeof(double) * (n + 1));
    x2 = (double*) malloc(sizeof(double) * (n + 1));

    ML_Operator_Apply(sm->my_level->Amat, n, x, n, rhs2);
    DCOPY_F77(&n, x, &one_int, x2, &one_int);
    DAXPY_F77(&n, &minus_one_double, rhs, &one_int, rhs2, &one_int);
    ML_Ifpack_Solve(Ifpack_Handle, x2, rhs2);
    DAXPY_F77(&n, &minus_one_double, x2, &one_int, x, &one_int);

    ML_free(rhs2);
    ML_free(x2);
  }
  return 0;
} /* ML_Smoother_Ifpack */

/* ------------------------------------------------------------------------- */
/* clean the Ifpack smoother                                                 */
/* ------------------------------------------------------------------------- */

void ML_Smoother_Clean_Ifpack(void *Ifpack_Handle)
{

  ML_Ifpack_Destroy(Ifpack_Handle);
  return;
  
} /* ML_Smoother_Clean_Ifpack */

#endif
