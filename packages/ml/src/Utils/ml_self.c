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
#include "ml_self.h"

/* ------------------------------------------------------------------------- */
/* generate the ML's self smoother                                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_Self(ML_Smoother *sm,int inlen,double x[],int outlen,
                     double rhs[])
{
  int iter;
  void *Self_Handle = sm->smoother->data;
  ML_Comm *comm = sm->my_level->comm;
  ML_CommInfoOP *getrow_comm = sm->my_level->Amat->getrow->pre_comm;

  /* No data exchange necessary, as this takes place in Ifpack_Additive_Schwarz. */

  for (iter = 0; iter < sm->ntimes; iter++) {
    ML_Self_Solve(Self_Handle, x, rhs);
  }

  return 0;

} /* ML_Smoother_Self */

/* ------------------------------------------------------------------------- */
/* clean the ML's self smoother                                              */
/* ------------------------------------------------------------------------- */

void ML_Smoother_Clean_Self(void *Self_Handle)
{

  ML_Self_Destroy(Self_Handle);
  return;
  
} /* ML_Smoother_Clean_Self */

#endif
