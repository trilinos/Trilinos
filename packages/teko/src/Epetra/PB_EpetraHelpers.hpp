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

/** \brief Builds an epetra operator from a block 2x2 matrix.
  *
  * Builds an epetra operator from a block 2x2 matrix.
  * <em>*** The calling user is responsible for deleting the resulting Epetra_Operator! ***</em>
  */
Epetra_Operator * block2x2(const Epetra_Operator * sub00,const Epetra_Operator * sub01,
                           const Epetra_Operator * sub10,const Epetra_Operator * sub11,
                           const std::string & str="ANYM");

/** \brief Swaps the Apply/ApplyInverse to ApplyInverse/Apply.
  *
  * Swaps the Apply/ApplyInverse operations to ApplyInverse/Apply.
  * <em>*** The calling user is responsible for deleting the resulting Epetra_Operator! ***</em>
  */
Epetra_Operator * mechanicalInverse(const Epetra_Operator * inverse);

} // end namespace Epetra
} // end namespace PB

#endif
