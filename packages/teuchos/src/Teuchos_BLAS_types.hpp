// /////////////////////////////////////////////
// Teuchos_BLAS_types.hpp

#ifndef TEUCHOS_BLAS_TYPES_HPP
#define TEUCHOS_BLAS_TYPES_HPP

namespace Teuchos {

enum ESide{ LEFT_SIDE, RIGHT_SIDE };
enum ETransp { NO_TRANS, TRANS, CONJ_TRANS };
enum EUplo { UPPER_TRI, LOWER_TRI };
enum EDiag { UNIT_DIAG, NON_UNIT_DIAG };

} // namespace Teuchos

#endif // TEUCHOS_BLAS_TYPES_HPP
