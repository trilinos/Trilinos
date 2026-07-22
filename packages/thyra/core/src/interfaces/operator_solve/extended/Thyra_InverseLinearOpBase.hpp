// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_INVERSE_LINEAR_OP_BASE_HPP
#define THYRA_INVERSE_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"


namespace Thyra {


/** \brief Base interface for <ttLinearOpBase</tt> objects that are implemented
 * in terms of the solve function on a <tt>LinearOpWithSolveBase</tt> object.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 */
template<class Scalar>
class InverseLinearOpBase : virtual public LinearOpBase<Scalar>
{
public:

  /** \brief Determine if the underlying <tt>LinearOpWithSolveBase</tt> is
   * const-only or not.
   */
  virtual bool isLowsConst() const = 0;

  /** \brief Extra a non-const view of the underlying
   * <tt>LinearOpWithSolveBase</tt> object.
   */
  virtual Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLows() = 0; 

  /** \brief Extra a const view of the underlying
   * <tt>LinearOpWithSolveBase</tt> object.
   */
  virtual Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
  getLows() const = 0; 

};


} // namespace Thyra


#endif // THYRA_INVERSE_LINEAR_OP_BASE_HPP
