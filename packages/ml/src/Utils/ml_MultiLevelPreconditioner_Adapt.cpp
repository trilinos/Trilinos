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
#include "ml_viz_xyz.h"
#include "ml_agg_info.h"

// ================================================ ====== ==== ==== == =
// the tentative null space is in input because the user
// has to remember to allocate and fill it, and then to delete
// it after calling this method.
int ML_Epetra::MultiLevelPreconditioner::
ComputeAdaptivePreconditioner(int TentativeNullSpaceSize,
                              double* TentativeNullSpace)
{

  if ((TentativeNullSpaceSize == 0) || (TentativeNullSpace == 0))
    ML_CHK_ERR(-1);
   
  // ================================== //
  // get parameters from the input list //
  // ================================== //
  
  // maximum number of relaxation sweeps
  int MaxSweeps = List_.get("adaptive: max sweeps", 10);
  // number of vector to be added to the tentative null space
  int NumAdaptiveVectors = List_.get("adaptive: num vectors", 1);
  // set to zero all vectors whose \infty norm is below this value
  double threshold = List_.get("adaptive: threshold", 0.0);

  // ==================================================== //
  // compute the preconditioner, set null space from user //
  // (who will have to delete vector TentativeNullSpace)  //
  // ==================================================== //
  
  // destroy what we may already have
  if (IsComputePreconditionerOK_ == true) {
    DestroyPreconditioner();
  }

  List_.set(Prefix_ + "null space: type", "pre-computed");
  List_.set(Prefix_ + "null space: dimension", TentativeNullSpaceSize);
  List_.set(Prefix_ + "null space: vectors", TentativeNullSpace);
  ComputePreconditioner();

  // ========================================= //
  // look for "bad" modes by simple relaxation //
  // ========================================= //
  
  Epetra_MultiVector RHS(RowMatrixRowMap(),NumAdaptiveVectors);
  Epetra_MultiVector LHS(RowMatrixRowMap(),NumAdaptiveVectors);
  LHS.Random();
  vector<double> Norm(NumAdaptiveVectors);

  if (verbose_) 
    cout << PrintMsg_ << "*** Relaxing up to " << MaxSweeps << " times," << endl
         << PrintMsg_ << "to compute " << NumAdaptiveVectors 
         << " additional vectors" << endl;

  for (int i = 0 ; i < MaxSweeps ; ++i) {
    // RHS = (I - ML^{-1} A) LHS
    ML_CHK_ERR(RowMatrix_->Multiply(false,LHS,RHS));
    // FIXME: can do something slightly better here
    ML_CHK_ERR(ApplyInverse(RHS,RHS));
    ML_CHK_ERR(LHS.Update(-1.0,RHS,1.0));
    LHS.Norm2(&Norm[0]);
    if (verbose_) {
      cout << PrintMsg_ << "*** Iter " << i << ", ||x||_2 = ";
      for (int j = 0 ; j < NumAdaptiveVectors ; ++j)
        cout << Norm[j] << " ";
      cout << endl;
    }
  }

  // scaling vectors
  for (int i = 0 ; i < NumAdaptiveVectors ; ++i) {
    double NormInf;
    LHS(i)->NormInf(&NormInf);
    LHS(i)->Scale(1.0 / NormInf);
  }

  // ========================================================= //
  // copy tentative and computed null space into NewNullSpace, //
  // which now becomes the standard null space                 //
  // ========================================================= //
  
  int NewNullSpaceSize = NumAdaptiveVectors + TentativeNullSpaceSize;
  double* NewNullSpace = new double[NumMyRows() * NewNullSpaceSize];
  assert (NewNullSpace != 0);
  int itmp = TentativeNullSpaceSize * NumMyRows();
  for (int i = 0 ; i < itmp ; ++i) {
      NewNullSpace[i] = TentativeNullSpace[i];
  }

  for (int i = 0 ; i < NumAdaptiveVectors ; ++i) {
    double Norm1;
    LHS(i)->Norm1(&Norm1);
    if (Norm1 >= threshold) {
      for (int j = 0 ; j < NumMyRows() ; ++j) {
        NewNullSpace[itmp + i * NumMyRows() + j] = LHS[i][j];
      }
    }
  }
   
  // Destroy the old preconditioner
  DestroyPreconditioner();
  
  // ==================================================== //
  // build the new preconditioner with the new null space //
  // ==================================================== //
  
  List_.set(Prefix_ + "null space: type", "pre-computed");
  List_.set(Prefix_ + "null space: dimension", NewNullSpaceSize);
  List_.set(Prefix_ + "null space: vectors", NewNullSpace);
  ComputePreconditioner();

  // keep trace of this pointer, it will be delete'd later
  NullSpaceToFree_ = NewNullSpace;

  // =============== //
  // visualize modes //
  // =============== //

  if (List_.get("adaptive: visualize", false)) {

    double* x_coord = List_.get("viz: x-coordinates", (double*)0);
    double* y_coord = List_.get("viz: y-coordinates", (double*)0);
    double* z_coord = List_.get("viz: z-coordinates", (double*)0);
    assert (x_coord != 0);

    vector<double> plot_me(NumMyRows()/NumPDEEqns_);
    ML_Aggregate_Viz_Stats info;
    info.Amatrix = &(ml_->Amat[LevelID_[0]]);
    info.x = x_coord;
    info.y = y_coord;
    info.z = z_coord;
    info.Nlocal = NumMyRows() / NumPDEEqns_;
    info.Naggregates = 1;
    
    for (int ieqn = 0 ; ieqn < NumPDEEqns_ ; ++ieqn) {
      for (int i = 0 ; i < NumAdaptiveVectors ; ++i) {
        for (int j = 0 ; j < NumMyRows() ; j+=NumPDEEqns_) {
          plot_me[j / NumPDEEqns_] = LHS[i][j + ieqn];
          char FileName[80];
          sprintf(FileName,"nullspace-mode%d-eq%d.xyz", i, ieqn);
          ML_Aggregate_VisualizeXYZ(info,FileName,
                                    ml_->comm,&plot_me[0]);
        }
      }
    }
  }

  return(0);

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
