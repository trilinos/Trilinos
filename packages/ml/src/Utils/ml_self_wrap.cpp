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
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_self.h"
#include "ml_self_wrap.h"
#include "ml_RowMatrix.h"
#include "ml_Ifpack_ML.h"
#include "Ifpack_AdditiveSchwarz.h"

using namespace ML_Epetra;

int ML_Self_Gen(ML *ml, int Overlap, int curr_level, 
                  Teuchos::ParameterList& List, 
                  const Epetra_Comm& Comm, 
                  void ** Self_Handle);

// ====================================================================== 
// This file is a conceptual copy of ml_ifpack_wrap.cpp.
// Copied on 08-Mar-05 (women's day).
// ====================================================================== 

int ML_Gen_Smoother_Self(ML *ml, int Overlap, int nl, int pre_or_post,
                         int niters, Teuchos::ParameterList& List,
                         const Epetra_Comm& Comm)
{

   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int status = 1;
   char str[80];
   void *Self_Handle ;

   fun = ML_Smoother_Self;

   /* Creates IFPACK objects */

   status = ML_Self_Gen(ml, Overlap, nl, List, Comm, &Self_Handle) ; 
   assert (status == 0); 

   /* Sets function pointers */

   if (pre_or_post == ML_PRESMOOTHER) {
     sprintf(str,"self_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]), (void*)Self_Handle,
			      fun, niters, 0.0, str);
     ml->pre_smoother[nl].data_destroy = ML_Smoother_Clean_Self;
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
     sprintf(str,"self_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void*)Self_Handle, fun, niters, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Self;
   }
   else if (pre_or_post == ML_BOTH) {
     sprintf(str,"self_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]),
			      (void*)Self_Handle,
			      fun, niters,  0.0, str);
     sprintf(str,"self_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]),
			      (void*)Self_Handle, fun, niters, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Self;
   }
   else 
     pr_error("ML_Gen_Smoother_Self: unknown pre_or_post choice\n");

   return(status);

}
// ================================================ ====== ==== ==== == =

#ifdef IFPACK_NODE_AWARE_CODE
int ML_NODE_ID = -1;  //FIXME delete this
#endif

int ML_Self_Gen(ML *ml, int Overlap, int curr_level, 
                Teuchos::ParameterList& List, const Epetra_Comm& Comm, 
                void ** Self_Handle)
{

  ML_Operator *Ke = &(ml->Amat[curr_level]);

  // creates the wrapper from ML_Operator to Epetra_RowMatrix
  // (ML_Epetra::RowMatrix). This is a cheap conversion
  RowMatrix* Self_Matrix = new RowMatrix(Ke, &Comm);
  assert (Self_Matrix != 0);

  Ifpack_Preconditioner* Prec;

  Prec = new Ifpack_AdditiveSchwarz<Ifpack_ML>(Self_Matrix, Overlap);
  assert (Prec != 0);

  List.set("zero starting solution", true);
  List.set("schwarz: compute condest", false);
#ifdef IFPACK_NODE_AWARE_CODE
  ML_NODE_ID = List.get("ML node id",-1);
#endif
  Prec->SetParameters(List);
  ML_CHK_ERR(Prec->Initialize());
  ML_CHK_ERR(Prec->Compute());

  *Self_Handle = (void *)Prec;

  return 0;
  
} /* ML_Self_Gen */

// ================================================ ====== ==== ==== == =

int ML_Self_Solve(void * Self_Handle, double * x, double * rhs )
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Self_Handle;

  Epetra_Vector Erhs(View, Prec->OperatorRangeMap(), rhs);
  Epetra_Vector Ex(View, Prec->OperatorDomainMap(), x);
  Prec->ApplyInverse(Erhs,Ex); 

  return 0;

} /* ML_Self_Solve */

// ================================================ ====== ==== ==== == =

void ML_Self_Destroy(void * Self_Handle)
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Self_Handle;
  if (ML_Get_PrintLevel() > 8)
    cout << *Prec;

  delete &(Prec->Matrix());
  delete Prec;

} /* ML_Self_Destroy */

#endif
