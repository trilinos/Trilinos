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


/** \defgroup Epetra_Thyra_Op_Vec_support_code_grp Epetra to Thyra Operator/Vector Adapter Support Code

\ingroup Epetra_Thyra_Op_Vec_adapters_grp

This is some basic support code that the Epetra to %Thyra operator/vector adapter Code is built on.

*/


/** \brief Determine if adjoints are supported on Epetra_Opeator or not.
 *
 * \ingroup Epetra_Thyra_Op_Vec_support_code_grp
 */
enum EAdjointEpetraOp {
  EPETRA_OP_ADJOINT_SUPPORTED      ///< Adjoint supported
  ,EPETRA_OP_ADJOINT_UNSUPPORTED   ///< Adjoint not supported
};


/** \brief . 
 *
 * \ingroup Epetra_Thyra_Op_Vec_support_code_grp
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
      TEST_FOR_EXCEPT(true);
  }
  return "";
}


/** \brief Determine how the apply an Epetra_Operator as a linear operator
 *
 * \ingroup Epetra_Thyra_Op_Vec_support_code_grp
 */
enum EApplyEpetraOpAs {
  EPETRA_OP_APPLY_APPLY            ///< Apply using Epetra_Operator::Apply(...)
  ,EPETRA_OP_APPLY_APPLY_INVERSE   ///< Apply using Epetra_Operator::ApplyInverse(...)
};


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_support_code_grp
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
      TEST_FOR_EXCEPT(true);
  }
  return "";
}


/** \brief . */
class EpetraLinearOp;


} // namespace Thyra

#endif // THYRA_EPETRA_TYPES_HPP
