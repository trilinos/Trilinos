#include "Amesos_ConfigDefs.h"
#include "Amesos_OverlappingPreconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Amesos_OverlappingPreconditioner::
Amesos_OverlappingPreconditioner(char* SolverType,
				 Epetra_CrsMatrix* Matrix,
				 int OverlappingLevel) 
{

  Teuchos::ParameterList List;
  Amesos_OverlappingPreconditioner(SolverType, Matrix, OverlappingLevel,
				   List);

}

//==============================================================================
Amesos_OverlappingPreconditioner::
Amesos_OverlappingPreconditioner(char* SolverType,
				 Epetra_CrsMatrix* Matrix,
				 int OverlappingLevel,
				 Teuchos::ParameterList& List) :
  Matrix_(Matrix),
  OverlappingLevel_(OverlappingLevel),
  Label_("Amesos_OverlappingPreconditioner")
{

  // FIXME: to add an option for symmetric/non-symmetric
  // FIXME:
  // not I don't consider no-overlap. Should I ???
  assert(OverlappingLevel > 0);

  OverlappingLevel_ = 1;
  if (OverlappingLevel_ > 0) {
    OverlappingMatrix_ = CreateOverlappingMatrix(Matrix, 
						 OverlappingMatrix_);
    assert(OverlappingMatrix_ != 0);


    OverlappingX_ = 
      new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);
    OverlappingY_ = 
      new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),1);

    // to import off-process data before the application of the
    // preconditioner
    Importer_ = new Epetra_Importer(Matrix->RowMatrixRowMap(),
				    OverlappingMatrix->RowMatrixRowMap());

    // wrap up the overlapping matrix into a local matrix
    LocalOverlappingMatrix_ = 
      new Amesos_LocalRowMatrix(OverlappingMatrix_);

    Problem_ = new Epetra_LinearProblem;
    Problem_->SetOperator(LocalOverlappingMatrix_);

    Amesos Factory;
    Solver_ = Factory.Create(SolverType,*Problem_);

    if (Solver_ == 0) {
      // try to create KLU, it is generally enabled
      Solver_ = Factory.Create("Amesos_Klu",*Problem_);

      if (Solver_ == 0) {
	AMESOS_CHK_ERRV(-1);
      }
    }

    Solver_->SetParameters(List);

    // factorize matrices
    AMESOS_CHK_ERRV(Solver_->SymbolicFactorization());
    AMESOS_CHK_ERRV(Solver_->NumericFactorization());

  } // fix the case without overlap
}

//==============================================================================
Amesos_OverlappingPreconditioner::~Amesos_OverlappingPreconditioner()
{
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

  // FIXME
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

  X.Import(Overlapping_X,*Importer_,Insert);

  Problem_->SetLHS(OverlappingY_);
  Problem_->SetRHS(OverlappingX_);

  AMESOS_CHK_ERR(Solver_->Solve());

  Y.Export(OverlappingY_,*Importer_,Insert);
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

