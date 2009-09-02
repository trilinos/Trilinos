// @HEADER
// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_GET_EPETRA_OPERATOR_HPP
#define THYRA_GET_EPETRA_OPERATOR_HPP

#include "Thyra_EpetraTypes.hpp"


namespace Thyra {


/** \defgroup Epetra_Thyra_Op_Vec_get_Epetra_Operator_grp Epetra_Operator extraction utility functions

\ingroup Epetra_Thyra_Op_Vec_adapters_grp

These function allow the extraction of an <tt>Epetra_Operator</tt> from a <tt>Thyra::LinearOpBase</tt> object.

*/


/** \brief Get smart pointer to non-<tt>const</tt>
 * <tt>Epetra_Operator</tt> object from reference to a
 * non-<tt>const</tt> <tt>EpetraLinearOp</tt> accessed through its
 * <tt>LinearOpBase</tt> interface.
 *
 * @param op [in] Reference to operator to extract <tt>Epetra_Operator</tt> out of.
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
 * \ingroup Epetra_Thyra_Op_Vec_get_Epetra_Operator_grp
 */
Teuchos::RCP<Epetra_Operator>
get_Epetra_Operator( LinearOpBase<double> &op );


/** \brief Get smart pointer to <tt>const</tt>
 * <tt>Epetra_Operator</tt> object from reference to a <tt>const</tt>
 * <tt>EpetraLinearOp</tt> accessed through its <tt>LinearOpBase</tt>
 * interface.
 *
 * @param op [in] Reference to operator to extract <tt>Epetra_Operator</tt> out of.
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
 * \ingroup Epetra_Thyra_Op_Vec_get_Epetra_Operator_grp
 */
Teuchos::RCP<const Epetra_Operator>
get_Epetra_Operator( const LinearOpBase<double> &op );


} // namespace Thyra


#endif // THYRA_GET_EPETRA_OPERATOR_HPP
