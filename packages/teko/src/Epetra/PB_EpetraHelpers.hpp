#ifndef __PB_EpetraHelpers_hpp__
#define __PB_EpetraHelpers_hpp__

// stl includes
#include <string>

// Epetra includes
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_VectorBase.hpp"

namespace PB {

namespace Epetra {

// build an epetra operator from a block 2x2 matrix
Epetra_Operator * block2x2(const Epetra_Operator * sub00,const Epetra_Operator * sub01,
                           const Epetra_Operator * sub10,const Epetra_Operator * sub11,
                           const std::string & str="ANYM");

// swaps the Apply,ApplyInverse to ApplyInverse,Apply
Epetra_Operator * mechanicalInverse(const Epetra_Operator * inverse);


// here we assume that Epetra_Vector is blocked and each block is distributed across all processors
void thyraToEpetra(const Teuchos::RCP<const Thyra::VectorBase<double> > & t, Epetra_MultiVector & e);

// here we assume that Epetra_Vector is blocked and each block is distributed across all processors
void epetraToThyra(const Epetra_MultiVector & e,const Teuchos::Ptr<Thyra::VectorBase<double> > & t);

} // end namespace Epetra
} // end namespace PB

#endif
