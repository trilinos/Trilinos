// how to keep the AztecOO's data up to the destructor??
#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_TEUCHOS) && defined(HAVE_IFPACK_AZTECOO)
#include "Ifpack_Preconditioner.h"
#include "Ifpack_AztecOO.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

#include "Ifpack_LocalRowMatrix.h"
// FIXME
#include "Ifpack_Jacobi.h"

//==============================================================================
Ifpack_AztecOO::Ifpack_AztecOO(Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  Solver_(0),
  Problem_(0),
  IsComputed_(false)
{
}

//==============================================================================
Ifpack_AztecOO::~Ifpack_AztecOO()
{
  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

}

//==============================================================================
int Ifpack_AztecOO::SetParameters(Teuchos::ParameterList& List)
{

  Problem_ = new Epetra_LinearProblem;
  Problem_->SetOperator(Matrix_);
  Solver_ = new AztecOO(*Problem_);

  // FIXME
  Ifpack_Preconditioner* Prec = new Ifpack_Jacobi(Matrix_);
  Prec->Compute();

  Solver_->SetPrecOperator(Prec);
  return(0);
}

//==============================================================================
int Ifpack_AztecOO::Compute()
{

  if (Solver_ == 0)
    IFPACK_CHK_ERR(-1);

  Solver_->SetAztecOption(AZ_solver,AZ_gmres);
  Solver_->SetAztecOption(AZ_output,1);

  Solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  Solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilu);

  IsComputed_ = true;

  return(0);
}

//==============================================================================
int Ifpack_AztecOO::SetUseTranspose(bool UseTranspose)
{
  IFPACK_CHK_ERR(-99); // not implemented
  return(-99);
}

//==============================================================================
int Ifpack_AztecOO::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check for maps ???
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
}

//==============================================================================
int Ifpack_AztecOO::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);

  // check for maps??

  int NumVectors = X.NumVectors();

  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input
  
  Epetra_MultiVector Xtmp(X);
  Epetra_MultiVector Ytmp(Y);
  Solver_->SetLHS(&Ytmp);
  Solver_->SetRHS((Epetra_MultiVector*)&Xtmp);
  IFPACK_CHK_ERR(Solver_->Iterate(1550, 1e-9));

  return(0);
}

//=============================================================================
double Ifpack_AztecOO::NormInf() const
{
  return(-1.0);
}

//==============================================================================
char * Ifpack_AztecOO::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Ifpack_AztecOO::UseTranspose() const
{
  return(false);
}

//==============================================================================
bool Ifpack_AztecOO::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Ifpack_AztecOO::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Ifpack_AztecOO::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Ifpack_AztecOO::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}
#endif // HAVE_IFPACK_AZTECOO && HAVE_IFPACK_TEUCHOS
