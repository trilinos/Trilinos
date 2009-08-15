// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

#include "Epetra_config.h"
#include "EpetraExt_TimedEpetraOperator.hpp"

EpetraExt::Epetra_Timed_Operator::Epetra_Timed_Operator(const Teuchos::RCP<Epetra_Operator>& A_) 
  : A(A_)
{
ApplyTimer = Teuchos::rcp(new Teuchos::Time("apply timer",false));
ApplyInverseTimer = Teuchos::rcp(new Teuchos::Time("apply inverse timer",false));
}

EpetraExt::Epetra_Timed_Operator::~Epetra_Timed_Operator()
{
}

int 
EpetraExt::Epetra_Timed_Operator::SetUseTranspose(bool useTranspose) 
{
  int success;
  success = A->SetUseTranspose(useTranspose);
  return success;
}

int 
EpetraExt::Epetra_Timed_Operator::Apply(const Epetra_MultiVector& Input, 
				   Epetra_MultiVector& Result) const
{
  int success;
  ApplyTimer->start();
  success = A->Apply(Input,Result);
  ApplyTimer->stop();
  return success;
}

int 
EpetraExt::Epetra_Timed_Operator::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  int success;
  ApplyInverseTimer->start();
  success = A->ApplyInverse(Input,Result);
  ApplyInverseTimer->stop();
  return success;
}

double 
EpetraExt::Epetra_Timed_Operator::NormInf() const
{
  return A->NormInf();
}


const char* 
EpetraExt::Epetra_Timed_Operator::Label () const
{
  return A->Label();
}
  
bool 
EpetraExt::Epetra_Timed_Operator::UseTranspose() const
{
  return A->UseTranspose();
}

bool 
EpetraExt::Epetra_Timed_Operator::HasNormInf() const
{
  return A->HasNormInf();
}

const Epetra_Comm & 
EpetraExt::Epetra_Timed_Operator::Comm() const
{
  return A->Comm();
}
const Epetra_Map& 
EpetraExt::Epetra_Timed_Operator::OperatorDomainMap() const
{
  return A->OperatorDomainMap();
}

const Epetra_Map& 
EpetraExt::Epetra_Timed_Operator::OperatorRangeMap() const
{
  return A->OperatorRangeMap();
}
