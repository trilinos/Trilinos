#include "Amesos_ConfigDefs.h"
#include "Amesos_Preconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Amesos_Preconditioner::Amesos_Preconditioner(char* SolverType,
				 Epetra_RowMatrix* Matrix) 
{

  Teuchos::ParameterList List;
  Amesos_Preconditioner(SolverType, Matrix, List);

}

//==============================================================================
Amesos_Preconditioner::Amesos_Preconditioner(char* SolverType,
				 Epetra_RowMatrix* Matrix,
				 Teuchos::ParameterList& List) :
  Matrix_(Matrix),
  Label_("Amesos_Preconditioner")
{

  Problem_ = new Epetra_LinearProblem;
  Problem_->SetOperator(Matrix_);

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

  AMESOS_CHK_ERRV(Solver_->SymbolicFactorization());
  AMESOS_CHK_ERRV(Solver_->NumericFactorization());

}

//==============================================================================
Amesos_Preconditioner::~Amesos_Preconditioner()
{
  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

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
  AMESOS_CHK_ERR(Matrix_->Apply(X,Y));
}

//==============================================================================
int Amesos_Preconditioner::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // FIXME
  Epetra_MultiVector Xtmp(X);
  Problem_->SetLHS(&Y);
  Problem_->SetRHS((Epetra_MultiVector*)&Xtmp);

  AMESOS_CHK_ERR(Solver_->Solve());

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
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Amesos_Preconditioner::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Amesos_Preconditioner::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}
