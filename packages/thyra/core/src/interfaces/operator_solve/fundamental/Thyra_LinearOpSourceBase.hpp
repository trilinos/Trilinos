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

#ifndef THYRA_LINEAR_OP_SOURCE_BASE_HPP
#define THYRA_LINEAR_OP_SOURCE_BASE_HPP

#include "Thyra_SolveSupportTypes.hpp"
#include "Teuchos_Describable.hpp"


namespace Thyra {


/** \brief Base interface for objects that can return a linear operator.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
template<class Scalar>
class LinearOpSourceBase : virtual public Teuchos::Describable
{
public:

  /** @name Pure virtual public functions that must be overridden in subclasses */
  //@{

  /** \brief Return if the underlying linear operator is const-only or not.
   */
  virtual bool isOpConst() const = 0;

  /** \brief Return a non-const reference to the underlying linear operator.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isOpConst()==true</tt>] <tt>getOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstOp() = 0;

  /** \brief Return a const left preconditioner linear operator if one is
   * designed or targeted to be applied on the left.
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getOp()const = 0;
  
  //@}

};

} // namespace Thyra

#endif // THYRA_LINEAR_OP_SOURCE_BASE_HPP
