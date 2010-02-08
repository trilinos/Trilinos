// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
