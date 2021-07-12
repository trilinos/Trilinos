#ifndef TPETRA_CREATEDEEPCOPY_CRSMATRIX_DECL_HPP
#define TPETRA_CREATEDEEPCOPY_CRSMATRIX_DECL_HPP

#include "Tpetra_createDeepCopy.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_RowMatrix_fwd.hpp"

#ifdef TPETRA_ENABLE_DEPRECATED_CODE

namespace Tpetra {

template<class SC, class LO, class GO, class NT>
CrsMatrix<SC, LO, GO, NT>
TPETRA_DEPRECATED
createDeepCopy (const RowMatrix<SC, LO, GO, NT>& in);

} // namespace Tpetra

#endif // TPETRA_ENABLE_DEPRECATED_CODE

#endif // TPETRA_CREATEDEEPCOPY_CRSMATRIX_DECL_HPP
