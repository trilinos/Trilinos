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

/* ------------------------------------------------------------------------- */
/* generate the Ifpack smoother                                              */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
		       double rhs[])
{

  ML_Smoother    *smooth_ptr = (ML_Smoother *) sm;
  void *Ifpack_Handle = smooth_ptr->smoother->data;

  ML_Ifpack_Solve(Ifpack_Handle, x, rhs);
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
