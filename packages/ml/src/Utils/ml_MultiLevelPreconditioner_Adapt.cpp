/*!
 * 
 *  \file ml_MultiLevelPreconditioner_Adapt.cpp
 *
 *  \brief Methods to define adaptive smoothed aggregation.
 *
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_include.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"

// ================================================ ====== ==== ==== == =
// the tentative null space is in input because the user
// has to remember to allocate and fill it, and then to delete
// it after calling this method.
int ML_Epetra::MultiLevelPreconditioner::
ComputeAdaptivePreconditioner(int TentativeNullSpaceSize,
                              double* TentativeNullSpace)
{

  // 0.- get parameters from the input list
  int MaxSweeps = List_.get("adaptive: max sweeps", 10);
  int NumAdaptiveVectors = List_.get("adaptive: num vectors", 1);

  // 1.- destroy what we may already have
  if (IsComputePreconditionerOK_ == true) {
    DestroyPreconditioner();
  }

  // 2.- compute the preconditioner, set null space from user
  //     (who will have to delete vector TentativeNullSpace)
  List_.set(Prefix_ + "null space: type", "pre-computed");
  List_.set(Prefix_ + "null space: dimension", TentativeNullSpaceSize);
  List_.set(Prefix_ + "null space: vectors", TentativeNullSpace);
  ComputePreconditioner();

  // 3.- look for "bad" modes by simple relaxation
  Epetra_MultiVector RHS(RowMatrixRowMap(),NumAdaptiveVectors);
  Epetra_MultiVector LHS(RowMatrixRowMap(),NumAdaptiveVectors);
  LHS.Random();
  vector<double> Norm(NumAdaptiveVectors);

  for (int i = 0 ; i < MaxSweeps ; ++i) {
    // RHS = (I - ML^{-1} A) LHS
    ML_CHK_ERR(RowMatrix_->Multiply(false,LHS,RHS));
    // FIXME: can do something slightly better here
    ML_CHK_ERR(ApplyInverse(RHS,RHS));
    ML_CHK_ERR(LHS.Update(-1.0,RHS,1.0));
    LHS.Norm2(&Norm[0]);
    cout << Norm[0] << endl;
  }

  // 4.- copy old and new null space into NewNullSpace, which
  //     now becomes the standard null space.
  int NewNullSpaceSize = NumAdaptiveVectors + TentativeNullSpaceSize;
  double* NewNullSpace = new double[NumMyRows() * NewNullSpaceSize];
  assert (NewNullSpace != 0);
  int itmp = TentativeNullSpaceSize * NumMyRows();
  for (int i = 0 ; i < itmp ; ++i) {
      NewNullSpace[i] = TentativeNullSpace[i];
  }
  for (int i = 0 ; i < NumAdaptiveVectors ; ++i) {
    for (int j = 0 ; j < NumMyRows() ; ++j) {
      NewNullSpace[itmp + i * NumMyRows() + j] = LHS[i][j];
    }
  }
   
  // 5.- Destroy the old preconditioner
  DestroyPreconditioner();
  
  // 6.- build the new preconditioner with the new null space
  List_.set(Prefix_ + "null space: type", "pre-computed");
  List_.set(Prefix_ + "null space: dimension", NewNullSpaceSize);
  List_.set(Prefix_ + "null space: vectors", NewNullSpace);
  ComputePreconditioner();

  // 7.- keep trace of this pointer, it will be delete'd later
  NullSpaceToFree_ = NewNullSpace;

  return(0);

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
