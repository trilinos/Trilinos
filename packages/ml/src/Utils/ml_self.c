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
  int iter,i;
  void *Self_Handle = sm->smoother->data;
  double *x2;
  ML_Comm *comm = sm->my_level->comm;
  ML_CommInfoOP *getrow_comm = sm->my_level->Amat->getrow->pre_comm;

  /*Exchange data.*/
  if (getrow_comm != NULL) {
     x2 = (double *) ML_allocate((inlen+getrow_comm->total_rcv_length+1)
                                 *sizeof(double));
     if (x2 == NULL) {
       if (comm->ML_mypid == 0)
         fprintf(stderr,"Not enough space in ML_Smoother_Self()\n");
#      ifdef HAVE_MPI
       MPI_Finalize();
#      endif
       exit(EXIT_FAILURE);
     }
     for (i = 0; i < inlen; i++) x2[i] = x[i];
  }
  else x2 = x;

  for (iter = 0; iter < sm->ntimes; iter++) {
    if (getrow_comm != NULL)
      ML_exchange_bdry(x2,getrow_comm, inlen,comm,ML_OVERWRITE,NULL);
    ML_Self_Solve(Self_Handle, x2, rhs);
  }

  if (getrow_comm != NULL) {
    for (i = 0; i < inlen; i++) x[i] = x2[i];
    ML_free(x2);
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
