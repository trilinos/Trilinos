// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
