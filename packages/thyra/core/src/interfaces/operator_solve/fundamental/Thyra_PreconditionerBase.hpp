// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_PRECONDITIONER_BASE_HPP
#define THYRA_PRECONDITIONER_BASE_HPP

#include "Thyra_OperatorSolveTypes.hpp"
#include "Teuchos_Describable.hpp"


namespace Thyra {


/** \brief Simple interface class to access a precreated preconditioner as one
 * or more linear operators objects and information on how they are meant to
 * be applied.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 *
 * This class supports four different kinds of preconditioners:
 * <ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the left:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()!=NULL</tt>
 *   <li><tt>getRightPrecOp().get()==NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Single preconditioner linear operator designed or targeted to be applied on the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()==NULL</tt>
 *   <li><tt>getRightPrecOp().get()!=NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Split two-sided preconditioner with linear operators designed or
 *   targeted to be applied on the left and the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()!=NULL</tt>
 *   <li><tt>getRightPrecOp().get()!=NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()==NULL</tt>
 *   </ul>
 * <li>Single preconditioner linear operator not designed or targeted to be
 *   applied on the left or the right:
 *   <ul>
 *   <li><tt>getLeftPrecOp().get()==NULL</tt>
 *   <li><tt>getRightPrecOp().get()==NULL</tt>
 *   <li><tt>getUnspecifiedPrecOp().get()!=NULL</tt>
 *   </ul>
 * </ul>
 */
template<class Scalar>
class PreconditionerBase : virtual public Teuchos::Describable
{
public:

  /** @name Pure virtual public functions that must be overridden in subclasses */
  //@{

  /** \brief Return if the underlying left preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isLeftPrecOpConst() const = 0;

  /** \brief Return a non-const left preconditioner linear operator if one is
   * designed or targeted to be applied on the left.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isLeftPrecOpConst()==true</tt>] <tt>getLeftPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp() = 0;

  /** \brief Return a const left preconditioner linear operator if one is
   * designed or targeted to be applied on the left.
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp()const = 0;

  /** \brief Return if the underlying right preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isRightPrecOpConst() const = 0;

  /** \brief Return a non-const right preconditioner linear operator if one is
   * designed or targeted to be applied on the right.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isRightPrecOpConst()==true</tt>] <tt>getRightPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp() = 0;

  /** \brief Return a const right preconditioner linear operator if one is
   * designed or targeted to be applied on the right.
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const = 0;

  /** \brief Return if the underlying unspecified preconditioner operator is
   * const-only or allows non-const access.
   */
  virtual bool isUnspecifiedPrecOpConst() const = 0;

  /** \brief Return a non-const generic preconditioner linear operator that is
   * not designed or targeted to be applied on the left or on the right.
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp() = 0;

  /** \brief Return a const generic preconditioner linear operator that is not
   * designed or targeted to be applied on the left or on the right.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>isUnspecifiedPrecOpConst()==true</tt>] <tt>getUnspecifiedPrecOp().get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const = 0;
  
  //@}

private:
  
  // Not defined and not to be called
  PreconditionerBase<Scalar>&
  operator=(const PreconditionerBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_PRECONDITIONER_BASE_HPP
