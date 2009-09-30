#include "PB_EpetraInverseOpWrapper.hpp"

using namespace Teuchos;

namespace PB {
namespace Epetra {

int EpetraInverseOpWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{ return EpetraOperatorWrapper::ApplyInverse(X,Y);}

int EpetraInverseOpWrapper::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{ return EpetraOperatorWrapper::Apply(X,Y);}

} // end namespace Epetra
} // end namespace PB
