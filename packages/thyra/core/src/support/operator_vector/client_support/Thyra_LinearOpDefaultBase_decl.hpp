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

#ifndef THYRA_LINEAR_OP_DEFAULT_DECL_HPP
#define THYRA_LINEAR_OP_DEFAULT_DECL_HPP

#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

/** \brief Node subclass that provides a good default implementation for
 * the <tt>describe()</tt> function.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class LinearOpDefaultBase : virtual public LinearOpBase<Scalar> {
public:
  
#ifdef THYRA_INJECT_USING_DECLARATIONS
  using LinearOpBase<Scalar>::apply;
  using LinearOpBase<Scalar>::describe;
#endif

  /** @name Public functions overridden from Teuchos::Describable */
  //@{

  /** \brief Default description that gives the label, type, and dimenstion . */
  std::string description() const;

  /** \brief Generates a default outputting for all linear operators.
   *
   * Calls on the <tt>this->description()</tt> function for the name of the
   * class (and possibly its instance name) and then if <tt>verbLevel >=
   * VERB_EXTREME</tt>, then the linear operators elements themselves are
   * printed as well.  The format of the output is as follows:
   *
   \verbatim

   type = 'this->description()', rangeDim = m, domainDim = n
     1:1:a11 1:2:a12 ... 1:n:a1n
     2:1:a21 2:2:a22 ... 1:n:a2n
     .       .           .
     .       .           .
     .       .           .
     m:1:am1 m:2:am2 ... m:n:amn
   \endverbatim
   *
   * The above matrix coefficients are with respect to the natural basis as
   * defined by the scalar products.
   *
   * Before <tt>type = 'this->description()'</tt> is printed and after
   * each newline, <tt>leadingIndent</tt> is output.  The
   * <tt>index:value</tt> lines are offset an additional
   * <tt>indentSpacer</tt> amount.  A newline is printed after the
   * last <tt>m:n:amn</tt> entry.
   */
  void describe(
    Teuchos::FancyOStream                &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;

  //@}

};	// end class LinearOpDefaultBase

}	// end namespace Thyra

#endif	// THYRA_LINEAR_OP_DEFAULT_DECL_HPP
