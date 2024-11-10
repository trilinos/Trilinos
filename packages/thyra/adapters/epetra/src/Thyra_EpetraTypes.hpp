// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_TYPES_HPP
#define THYRA_EPETRA_TYPES_HPP

#include "Thyra_OperatorVectorTypes.hpp"

// Define this to see selected timers
//#define EPETRA_THYRA_TEUCHOS_TIMERS

class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Operator;


namespace Thyra {


/** \brief Determine if adjoints are supported on Epetra_Opeator or not.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
enum EAdjointEpetraOp {
  EPETRA_OP_ADJOINT_SUPPORTED      ///< Adjoint supported
  ,EPETRA_OP_ADJOINT_UNSUPPORTED   ///< Adjoint not supported
};


/** \brief . 
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
inline
const std::string toString(const EAdjointEpetraOp adjointEpetraOp)
{
  switch(adjointEpetraOp) {
    case EPETRA_OP_ADJOINT_SUPPORTED:
      return "EPETRA_OP_ADJOINT_SUPPORTED";
    case EPETRA_OP_ADJOINT_UNSUPPORTED:
      return "EPETRA_OP_ADJOINT_UNSUPPORTED";
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  TEUCHOS_UNREACHABLE_RETURN("");
}


/** \brief Determine how the apply an Epetra_Operator as a linear operator
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
enum EApplyEpetraOpAs {
  EPETRA_OP_APPLY_APPLY            ///< Apply using Epetra_Operator::Apply(...)
  ,EPETRA_OP_APPLY_APPLY_INVERSE   ///< Apply using Epetra_Operator::ApplyInverse(...)
};


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
inline
const std::string toString(const EApplyEpetraOpAs applyEpetraOpAs)
{
  switch(applyEpetraOpAs) {
    case EPETRA_OP_APPLY_APPLY:
      return "EPETRA_OP_APPLY_APPLY";
    case EPETRA_OP_APPLY_APPLY_INVERSE:
      return "EPETRA_OP_APPLY_APPLY_INVERSE";
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  TEUCHOS_UNREACHABLE_RETURN("");
}


/** \brief . */
class EpetraLinearOp;


} // namespace Thyra

#endif // THYRA_EPETRA_TYPES_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

