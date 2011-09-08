// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_GET_EPETRA_OPERATOR_HPP
#define THYRA_GET_EPETRA_OPERATOR_HPP

#include "Thyra_EpetraTypes.hpp"


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
Teuchos::RCP<const Epetra_Operator>
get_Epetra_Operator( const LinearOpBase<double> &op );


} // namespace Thyra


#endif // THYRA_GET_EPETRA_OPERATOR_HPP
