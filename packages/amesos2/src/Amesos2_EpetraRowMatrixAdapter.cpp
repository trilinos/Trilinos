#ifndef AMESOS2_EPETRA_ROWMATRIX_ADAPTER_CPP
#define AMESOS2_EPETRA_ROWMATRIX_ADAPTER_CPP

#include "Amesos2_EpetraRowMatrixAdapter_decl.hpp"

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
#  include "Amesos2_EpetraRowMatrixAdapter_def.hpp"
#  include "Teuchos_ExplicitInstantiationHelpers.hpp"
namespace Amesos {
/* Need to figure out the explicit instantiation system for Amesos2, since we
 * are instantiating not only on the Scalar types, but Matrices and Vectors
 */

// TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES( )
}
#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#endif  // AMESOS2_EPETRA_ROWMATRIX_ADAPTER_CPP
