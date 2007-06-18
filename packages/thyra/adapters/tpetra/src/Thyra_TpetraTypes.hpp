
#ifndef THYRA_TPETRA_TYPES_HPP
#define THYRA_TPETRA_TYPES_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Tpetra {
template<typename Ordinal, typename Scalar> class VectorSpace;
template<typename Ordinal, typename Scalar> class Vector;
template<typename Ordinal, typename Scalar> class Operator;
} // namespace Tpetra

namespace Thyra {

/** \defgroup Tpetra_Thyra_Op_Vec_support_code_grp Tpetra to Thyra Operator/Vector Adapter Support Code

\ingroup Tpetra_Thyra_Op_Vec_adapters_grp

This is some basic support code that the Tpetra to %Thyra operator/vector adapter Code is built on.

*/

/** \brief Determine if adjoints are supported on Tpetra_Opeator or not.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_support_code_grp
 */
enum EAdjointTpetraOp {
  TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED    ///< Adjoint not supported
  ,TPETRA_OP_ADJOINT_SUPPORTED               ///< Adjoint (conjugate transpose) supported
  ,TPETRA_OP_TRANSPOSE_SUPPORTED             ///< Tranpose (non-conjugate transpose) supported
  ,TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED     ///< Adjoint not supported
};

/** \brief . 
 *
 * \ingroup Tpetra_Thyra_Op_Vec_support_code_grp
 */
inline
const char* toString(const EAdjointTpetraOp adjointTpetraOp)
{
  switch(adjointTpetraOp) {
    case TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED:
      return "TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED";
    case TPETRA_OP_ADJOINT_SUPPORTED:
      return "TPETRA_OP_ADJOINT_SUPPORTED";
    case TPETRA_OP_TRANSPOSE_SUPPORTED:
      return "TPETRA_OP_TRANSPOSE_SUPPORTED";
    case TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED:
      return "TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return NULL;
}

template<class Ordinal, class Scalar> class TpetraLinearOpBase;
template<class Ordinal, class Scalar> class TpetraLinearOp;

} // namespace Thyra

#endif // THYRA_TPETRA_TYPES_HPP
