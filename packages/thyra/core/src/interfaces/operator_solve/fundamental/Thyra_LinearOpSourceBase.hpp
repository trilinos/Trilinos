// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
