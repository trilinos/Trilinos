#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Amesos.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Ifpack_Amesos::Ifpack_Amesos(Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  Solver_(0),
  Problem_(0),
  IsInitialized_(false),
  IsComputed_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  Label_("Amesos_Klu")
{

  Problem_ = new Epetra_LinearProblem;

  Problem_->SetOperator(Matrix_);

}

//==============================================================================
Ifpack_Amesos::~Ifpack_Amesos()
{
  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

}

//==============================================================================
int Ifpack_Amesos::SetParameters(Teuchos::ParameterList& List)
{

  Label_ = List.get("amesos: solver type", "Amesos_Klu");

  return(0);
}

//==============================================================================
int Ifpack_Amesos::Initialize()
{

  IsInitialized_ = false;
  IsComputed_ = false;

  if (Matrix_ == 0)
    IFPACK_CHK_ERR(-1);

  // reallocate the solver. This allows the user to call
  // or not SetParameters()
  if (Solver_)
    delete Solver_;

  Amesos Factory;
  Solver_ = Factory.Create((char*)Label_.c_str(),*Problem_);
  
  if (Solver_ == 0) {
    // try to create KLU, it is generally enabled
    Solver_ = Factory.Create("Amesos_Klu",*Problem_);
  }
  if (Solver_ == 0)
    IFPACK_CHK_ERR(-1);

  Solver_->SetParameters(List_);
  IFPACK_CHK_ERR(Solver_->SymbolicFactorization());

  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_Amesos::Compute()
{

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;

  if (Matrix_ == 0)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Solver_->NumericFactorization());

  IsComputed_ = true;
  return(0);
}

//==============================================================================
int Ifpack_Amesos::SetUseTranspose(bool UseTranspose)
{
  IFPACK_CHK_ERR(-99); // not implemented
}

//==============================================================================
int Ifpack_Amesos::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check for maps ???
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}

//==============================================================================
int Ifpack_Amesos::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);

  // check for maps??

  int NumVectors = X.NumVectors();

  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input
  
  Problem_->SetLHS(&Y);
  Problem_->SetRHS((Epetra_MultiVector*)&X);
  IFPACK_CHK_ERR(Solver_->Solve());

  return(0);
}

//==============================================================================
double Ifpack_Amesos::NormInf() const
{
  return(-1.0);
}

//==============================================================================
char * Ifpack_Amesos::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Ifpack_Amesos::UseTranspose() const
{
  return(false);
}

//==============================================================================
bool Ifpack_Amesos::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Ifpack_Amesos::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Ifpack_Amesos::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Ifpack_Amesos::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
std::ostream& Ifpack_Amesos::Print(std::ostream& os) const
{
  if (Matrix().Comm().MyPID())
    return(os);

  os << "*** Ifpack_Amesos" << endl;
  return(os);
}
#endif // HAVE_IFPACK_AMESOS && HAVE_IFPACK_TEUCHOS
