#include "Amesos_ConfigDefs.h"
#include "Amesos_Preconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_LocalRowMatrix.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Amesos_Preconditioner::Amesos_Preconditioner(char* SolverType,
					     Epetra_RowMatrix* Matrix,
					     const bool LocalizeMatrix)
{

  Teuchos::ParameterList List;
  AMESOS_CHK_ERRV(Compute(SolverType,Matrix,List,LocalizeMatrix));

}

//==============================================================================
Amesos_Preconditioner::Amesos_Preconditioner(char* SolverType,
					     Epetra_RowMatrix* Matrix,
					     Teuchos::ParameterList& List,
					     const bool LocalizeMatrix)
{
  AMESOS_CHK_ERRV(Compute(SolverType,Matrix,List,LocalizeMatrix));
}

//==============================================================================
int Amesos_Preconditioner::Compute(char* SolverType,
				   Epetra_RowMatrix* Matrix,
				   Teuchos::ParameterList& List,
				   const bool LocalizeMatrix)
{

  Matrix_ = Matrix;
  LocalizedMatrix_ = 0;
  IsLocalized_ = LocalizeMatrix;

  Label_ = string(SolverType) + " preconditioner, overlap = 0";

  Problem_ = new Epetra_LinearProblem;

  if (IsLocalized() == true) {

    LocalizedMatrix_ = new Amesos_LocalRowMatrix(Matrix);
    assert (LocalizedMatrix_ != 0);
    Problem_->SetOperator(LocalizedMatrix_);

  }
  else
    Problem_->SetOperator(Matrix_);

  Amesos Factory;
  Solver_ = Factory.Create(SolverType,*Problem_);
  
  if (Solver_ == 0) {
    // try to create KLU, it is generally enabled
    Solver_ = Factory.Create("Amesos_Klu",*Problem_);

    if (Solver_ == 0) {
      AMESOS_CHK_ERR(-1);
    }
  }

  Solver_->SetParameters(List);

  AMESOS_CHK_ERR(Solver_->SymbolicFactorization());
  AMESOS_CHK_ERR(Solver_->NumericFactorization());

  return(0);
}

//==============================================================================
Amesos_Preconditioner::~Amesos_Preconditioner()
{
  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

  if (LocalizedMatrix_)
    delete LocalizedMatrix_;

}

//==============================================================================
int Amesos_Preconditioner::SetUseTranspose(bool UseTranspose)
{
  AMESOS_CHK_ERR(-99); // not implemented
  return(-99);
}

//==============================================================================
int Amesos_Preconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsLocalized() == false) {
    AMESOS_CHK_ERR(Matrix_->Apply(X,Y));
  } 
  else {
    AMESOS_CHK_ERR(LocalizedMatrix_->Apply(X,Y));
  }
}

//==============================================================================
int Amesos_Preconditioner::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();

  if (NumVectors != Y.NumVectors())
    AMESOS_CHK_ERR(-1); // wrong input
  
  if (IsLocalized()) {
    
    Epetra_MultiVector LocalizedX(Copy,LocalizedMatrix_->RowMatrixRowMap(),
				  X.Pointers(), NumVectors);
    Epetra_MultiVector LocalizedY(View,LocalizedMatrix_->RowMatrixRowMap(),
				  Y.Pointers(), NumVectors);
    
    Problem_->SetLHS(&LocalizedY);
    Problem_->SetRHS(&LocalizedX);

    AMESOS_CHK_ERR(Solver_->Solve());

  }
  else {

    Epetra_MultiVector Xtmp(X);
    Problem_->SetLHS(&Y);
    Problem_->SetRHS((Epetra_MultiVector*)&Xtmp);
    AMESOS_CHK_ERR(Solver_->Solve());

  }
    

  return(0);
}

//==============================================================================
double Amesos_Preconditioner::NormInf() const
{
  return(-1.0);
}

//==============================================================================
char * Amesos_Preconditioner::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Amesos_Preconditioner::UseTranspose() const
{
  return(false);
}

//==============================================================================
bool Amesos_Preconditioner::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Amesos_Preconditioner::Comm() const
{
  if (IsLocalized() == false) {
    return(Matrix_->Comm());
  }
  else {
    return(LocalizedMatrix_->Comm());
  }
}

//==============================================================================
const Epetra_Map & Amesos_Preconditioner::OperatorDomainMap() const
{
  if (IsLocalized() == false) {
    return(Matrix_->OperatorDomainMap());
  }
  else {
    return(LocalizedMatrix_->OperatorDomainMap());
  }
}

//==============================================================================
const Epetra_Map & Amesos_Preconditioner::OperatorRangeMap() const
{
  if (IsLocalized() == false)
    return(Matrix_->OperatorRangeMap());
  else
    return(LocalizedMatrix_->OperatorRangeMap());
}
