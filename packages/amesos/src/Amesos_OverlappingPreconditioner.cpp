#include "Amesos_ConfigDefs.h"
#include "Amesos_OverlappingPreconditioner.h"
#include "Amesos_InverseFactory.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_LocalRowMatrix.h"
#include "Amesos.h"
#include "Amesos_Utils.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Amesos_OverlappingPreconditioner::
Amesos_OverlappingPreconditioner(char* SolverType,
				 Epetra_CrsMatrix* Matrix,
				 Amesos_InverseFactory& Factory,
				 int OverlappingLevel,
				 bool IsSymmetric)
{

  Teuchos::ParameterList List;
  AMESOS_CHK_ERRV(Compute(SolverType,Matrix,OverlappingLevel,
			  IsSymmetric,Factory,List));

}

//==============================================================================
Amesos_OverlappingPreconditioner::
Amesos_OverlappingPreconditioner(char* SolverType,
				 Epetra_CrsMatrix* Matrix,
				 Amesos_InverseFactory& Factory,
				 int OverlappingLevel,
				 bool IsSymmetric,
				 Teuchos::ParameterList& List)
{
  AMESOS_CHK_ERRV(Compute(SolverType,Matrix,OverlappingLevel,
			  IsSymmetric,Factory,List));
}

int Amesos_OverlappingPreconditioner::
Compute(char* SolverType, Epetra_CrsMatrix* Matrix,
	int OverlappingLevel, bool IsSymmetric,
	Amesos_InverseFactory& Factory, Teuchos::ParameterList& List)
{

  Matrix_ = Matrix;
  OverlappingLevel_ = OverlappingLevel;
  IsSymmetric_ = IsSymmetric;
  Inverse_ = 0;

  // =========================== //
  // build the overlapping graph //
  // =========================== //
 
  Label_ = string(SolverType) + " preconditioner, overlap = " +
    Amesos_toString(OverlappingLevel);

  if ((OverlappingLevel_ > 0) && (Comm().NumProc() > 1)) {

    OverlappingMatrix_ = Amesos_CreateOverlappingCrsMatrix(Matrix, 
							   OverlappingLevel_);
    assert(OverlappingMatrix_ != 0);

  }
  else 
    OverlappingMatrix_ = Matrix_;

    Importer_ = new Epetra_Import(OverlappingMatrix_->RowMatrixRowMap(), 
				  Matrix->RowMatrixRowMap());
    Exporter_ = new Epetra_Export(Matrix->RowMatrixRowMap(),
				  OverlappingMatrix_->RowMatrixRowMap());

  OverlappingX_ = 
    new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);
  OverlappingY_ = 
    new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);

  // to import off-process data before the application of the
  // preconditioner
  // FIXME: can I save one?
  // wrap up the overlapping matrix into a local matrix
  LocalOverlappingMatrix_ = 
    new Amesos_LocalRowMatrix(OverlappingMatrix_);
  assert (LocalOverlappingMatrix_ != 0);

  // ================= //
  // build the inverse //
  // ================= //

  Inverse_ = Factory.Create(SolverType,LocalOverlappingMatrix_,
			    Matrix_,List);
 
  if (Inverse_ == 0)
    AMESOS_CHK_ERR(-1);
    
  return(0);
}

//==============================================================================
Amesos_OverlappingPreconditioner::~Amesos_OverlappingPreconditioner()
{

  if (OverlappingMatrix_ && OverlappingLevel_)
    delete OverlappingMatrix_;

  if (OverlappingX_)
    delete OverlappingX_;

  if (OverlappingY_)
    delete OverlappingY_;

  if (LocalOverlappingMatrix_)
    delete LocalOverlappingMatrix_;

  if (Importer_)
    delete Importer_;

  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

  if (Inverse_)
    delete Inverse_;
}

//==============================================================================
int Amesos_OverlappingPreconditioner::SetUseTranspose(bool UseTranspose)
{
  AMESOS_CHK_ERR(-99); // not implemented
  return(-99);
}

//==============================================================================
int Amesos_OverlappingPreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // apply the original (non-overlapping) matrix
  AMESOS_CHK_ERR(Matrix_->Apply(X,Y));
}

//==============================================================================
int Amesos_OverlappingPreconditioner::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (X.NumVectors() != Y.NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong

  if (X.NumVectors() != OverlappingX_->NumVectors()) {
    delete OverlappingX_;
    delete OverlappingY_;
    OverlappingX_ = 
      new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),
			     X.NumVectors());
    OverlappingY_ = 
      new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),
			     X.NumVectors());
    assert (OverlappingX_ != 0);
    assert (OverlappingY_ != 0);
  }

  AMESOS_CHK_ERR(OverlappingX_->Export(X,*Exporter_,Insert));

  Epetra_MultiVector LocalX(View,LocalOverlappingMatrix_->RowMatrixRowMap(),
			    OverlappingX_->Pointers(),X.NumVectors());
  Epetra_MultiVector LocalY(View,LocalOverlappingMatrix_->RowMatrixRowMap(),
			    OverlappingY_->Pointers(),Y.NumVectors());

  AMESOS_CHK_ERR(Inverse_->ApplyInverse(LocalX,LocalY));

  if (IsSymmetric()) {
    AMESOS_CHK_ERR(Y.Export(*OverlappingY_,*Importer_,Add));
  }
  else {
    AMESOS_CHK_ERR(Y.Export(*OverlappingY_,*Importer_,Zero));
  }

  return(0);
}

//==============================================================================
double Amesos_OverlappingPreconditioner::NormInf() const
{
  return(-1.0);
}

//==============================================================================
char * Amesos_OverlappingPreconditioner::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Amesos_OverlappingPreconditioner::UseTranspose() const
{
  return(false);
}

//==============================================================================
bool Amesos_OverlappingPreconditioner::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Amesos_OverlappingPreconditioner::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Amesos_OverlappingPreconditioner::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Amesos_OverlappingPreconditioner::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}
