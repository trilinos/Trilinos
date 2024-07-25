// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_GET_EPETRA_OPERATOR_HPP
#define THYRA_GET_EPETRA_OPERATOR_HPP

#include "Thyra_EpetraTypes.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <stdexcept> // std::invalid_argument

namespace Thyra {

/** \brief Get smart pointer to non-<tt>const</tt>
 * <tt>Epetra_Operator</tt> object from reference to a
 * non-<tt>const</tt> <tt>EpetraLinearOp</tt> accessed through its
 * <tt>LinearOpBase</tt> interface.
 *
 * \param op [in] Reference to operator to extract <tt>Epetra_Operator</tt>
 * out of.
 *
 * Preconditions:<ul>
 * <li><tt>dynamic_cast<EpetraLinearOp*>(&op) != NULL</tt>
 * </ul>
 *
 * This function is designed to provide an easy way for non-C++ experts
 * to get at the <tt>Epetra_Operator</tt> object that was stuck into
 * an <tt>EpetraLinearOp</tt> object.
 *
 * If the dynamic cast fails then a <tt>std::bad_cast</tt> exception
 * is thrown containing a very detailed error message as to why the
 * cast failed.
 *
 * This function is simple enough and developers can see what needs to
 * be done to accomplish this type of access by looking at the source
 * code by clicking on:
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
template<class Scalar>
Teuchos::RCP<Epetra_Operator>
get_Epetra_Operator (LinearOpBase<Scalar>& /* op */)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::invalid_argument, "Thyra::get_Epetra_Operator: This function "
     "only works if Scalar=double, because for Epetra objects, the only Scalar"
     " type is double.  Instead, Scalar = " <<
     Teuchos::TypeNameTraits<Scalar>::name () << ".");
}

//! Full specialization for Scalar=double.
template<>
Teuchos::RCP<Epetra_Operator>
get_Epetra_Operator( LinearOpBase<double> &op );

/** \brief Get smart pointer to <tt>const</tt>
 * <tt>Epetra_Operator</tt> object from reference to a <tt>const</tt>
 * <tt>EpetraLinearOp</tt> accessed through its <tt>LinearOpBase</tt>
 * interface.
 *
 * \param op [in] Reference to operator to extract <tt>Epetra_Operator</tt>
 * out of.
 *
 * Preconditions:<ul>
 * <li><tt>dynamic_cast<const EpetraLinearOp*>(&op) != NULL</tt>
 * </ul>
 *
 * This function is designed to provide an easy way for non-C++ experts
 * to get at the <tt>Epetra_Operator</tt> object that was stuck into
 * an <tt>EpetraLinearOp</tt> object.
 *
 * If the dynamic cast fails then a <tt>std::bad_cast</tt> exception
 * is thrown containing a very detailed error message as to why the
 * cast failed.
 *
 * This function is simple enough and developers can see what needs to
 * be done to accomplish this type of access by looking at the source
 * code by clicking on:
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
template<class Scalar>
Teuchos::RCP<const Epetra_Operator>
get_Epetra_Operator( const LinearOpBase<Scalar> & /* op */ )
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::invalid_argument, "Thyra::get_Epetra_Operator: This function "
     "only works if Scalar=double, because for Epetra objects, the only Scalar"
     " type is double.  Instead, Scalar = " <<
     Teuchos::TypeNameTraits<Scalar>::name () << ".");
}

//! Full specialization for Scalar=double.
template<>
Teuchos::RCP<const Epetra_Operator>
get_Epetra_Operator( const LinearOpBase<double> &op );

} // namespace Thyra


#endif // THYRA_GET_EPETRA_OPERATOR_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

