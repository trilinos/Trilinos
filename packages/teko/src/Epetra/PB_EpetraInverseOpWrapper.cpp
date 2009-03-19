#include "PB_EpetraInverseOpWrapper.hpp"

using namespace Teuchos;

namespace PB {
namespace Epetra {

const Epetra_Map& EpetraInverseOpWrapper::OperatorDomainMap() const
{return EpetraOperatorWrapper::OperatorRangeMap();}

const Epetra_Map& EpetraInverseOpWrapper::OperatorRangeMap() const
{return EpetraOperatorWrapper::OperatorDomainMap();}

int EpetraInverseOpWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{ EpetraOperatorWrapper::ApplyInverse(X,Y);}

int EpetraInverseOpWrapper::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{ EpetraOperatorWrapper::Apply(X,Y);}

} // end namespace Epetra
} // end namespace PB
