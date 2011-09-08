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
      TEST_FOR_EXCEPT(true);
  }
  return "";
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
      TEST_FOR_EXCEPT(true);
  }
  return "";
}


/** \brief . */
class EpetraLinearOp;


} // namespace Thyra

#endif // THYRA_EPETRA_TYPES_HPP
