/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRA)
#include "ml_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_VbrMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"
// converter from ML_Operator to Epetra_RowMatrix (only wraps)
#include "ml_RowMatrix.h"
// IFPACK factory class
#include "Ifpack.h"

using namespace ML_Epetra;

static map<void*, bool> MemoryManager;

int ML_Ifpack_Gen(ML *ml, const char* Type, int Overlap, int curr_level, 
                  Teuchos::ParameterList& List, 
                  const Epetra_Comm& Comm, 
                  void ** Ifpack_Handle);

// ====================================================================== 
// MS // This does not work yet with ML_ALL_LEVELS
// MS // I also ask for the Epetra_Comm, as IFPACK only works
// MS // with the Epetra interface.
int ML_Gen_Smoother_Ifpack(ML *ml, const char* Type, int Overlap,
                           int nl, int pre_or_post,
                           Teuchos::ParameterList& List,
                           const Epetra_Comm& Comm)
{

   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int status = 1;
   char str[80];
   void *Ifpack_Handle ;

   fun = ML_Smoother_Ifpack;

   /* Creates IFPACK objects */

   status = ML_Ifpack_Gen(ml, Type, Overlap, nl, List, Comm, &Ifpack_Handle) ; 
   assert (status == 0); 

   /* Sets function pointers */

   if (pre_or_post == ML_PRESMOOTHER) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]), (void*)Ifpack_Handle,
			      fun, 1, 0.0, str);
     ml->pre_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void*)Ifpack_Handle, fun, 1, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else if (pre_or_post == ML_BOTH) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]),
			      (void*)Ifpack_Handle,
			      fun, 1,  0.0, str);
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]),
			      (void*)Ifpack_Handle, fun, 1, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
   }
   else 
     return(pr_error("ML_Gen_Smoother_Ifpack: unknown pre_or_post choice\n"));

   return(status);

}
// ================================================ ====== ==== ==== == =

int ML_Ifpack_Gen(ML *ml, const char* Type, int Overlap, int curr_level, 
                  Teuchos::ParameterList& List, 
                  const Epetra_Comm& Comm, 
                  void ** Ifpack_Handle)
{
  ML_Operator *Ke = &(ml->Amat[curr_level]);

  Epetra_RowMatrix* Ifpack_Matrix;

  if (Ke->type == ML_TYPE_ROW_MATRIX)
  {
    Ifpack_Matrix = (Epetra_RowMatrix*) Ke->data;
    // I have to remember not to delete this guy
    MemoryManager[(void*)Ifpack_Matrix] = false;
  }
  else if(Ke->type == ML_TYPE_CRS_MATRIX)
  {
    Ifpack_Matrix = (Epetra_CrsMatrix*) Ke->data;
    // I have to remember not to delete this guy
    MemoryManager[(void*)Ifpack_Matrix] = false;
  }
  else if(Ke->type == ML_TYPE_VBR_MATRIX)
  {
    Ifpack_Matrix = (Epetra_VbrMatrix*) Ke->data;
    // I have to remember not to delete this guy
    MemoryManager[(void*)Ifpack_Matrix] = false;
  }
  else
  {
    // creates the wrapper from ML_Operator to Epetra_RowMatrix
    // (ML_Epetra::RowMatrix). This is a cheap conversion
    Ifpack_Matrix = new RowMatrix(Ke, &Comm);
    assert (Ifpack_Matrix != 0);
    // this guy has to be deleted
    MemoryManager[(void*)Ifpack_Matrix] = true;
  }

  // we enter the IFPACK world through the factory only
  Ifpack Factory;
  Ifpack_Preconditioner* Prec;

  // create the preconditioner 
  Prec = Factory.Create(Type, Ifpack_Matrix, Overlap);
  assert (Prec != 0);

  Prec->SetParameters(List);
  ML_CHK_ERR(Prec->Initialize());
  ML_CHK_ERR(Prec->Compute());

  *Ifpack_Handle = (void *)Prec;

  return 0;
  
} /* ML_Ifpack_Gen */

// ================================================ ====== ==== ==== == =

int ML_Ifpack_Solve(void * Ifpack_Handle, double * x, double * rhs )
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle;

  Epetra_Vector Erhs(View, Prec->OperatorRangeMap(), rhs);
  Epetra_Vector Ex(View, Prec->OperatorDomainMap(), x);
  Prec->ApplyInverse(Erhs,Ex); 

  return 0;

} /* ML_Ifpack_Solve */

// ================================================ ====== ==== ==== == =

void ML_Ifpack_Destroy(void * Ifpack_Handle)
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle;
  if (ML_Get_PrintLevel() > 8)
    cout << *Prec;

  if (MemoryManager[(void*)(&(Prec->Matrix()))])
  {
    delete &(Prec->Matrix());
    MemoryManager[(void*)(&(Prec->Matrix()))] = false;
  }
  delete Prec;

} /* ML_Ifpack_Destroy */

#endif
