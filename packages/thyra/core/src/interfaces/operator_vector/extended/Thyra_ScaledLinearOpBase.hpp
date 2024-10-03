// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SCALED_LINEAR_OP_BASE_HPP
#define THYRA_SCALED_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase_decl.hpp"


namespace Thyra {


/** \brief Applies left or right sclaing to the linear operator.
 *
 * This interface represents a linear operator <tt>M</tt> that explicitly
 * applies left or right scaling by a diagonal (vector) operator <tt>d</tt>.
 #
 * Left scaling:
 * 
 \verbatim
 M = dM
 \endverbatim
 *
 * or Right scaling:
 *
 \verbatim
 M = Md
 \endverbatim
 *
 * where:
 *
 * <ul>
 *
 * <li> <tt>M</tt> is the <tt>LinearOp</tt> object,
 *
 * <li> <tt>d</tt> is the <tt>VectorBase</tt> object representing the diagonal
 * scaling operator.
 *
 * </ul>
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class ScaledLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** @name Non-virtual public interface functions. */
  //@{

  /** \brief Determines if this objects supports left scaling.
   */
  bool supportsScaleLeft() const
    { return supportsScaleLeftImpl(); }

  /** \brief Determines if this objects supports right scaling.
   */
  bool supportsScaleRight() const
    { return supportsScaleRightImpl(); }

  /** \brief Left scales operator with diagonal scaling operator.
   *
   * \precondtion <tt>supportsScaleLeft()==true</tt>
   */
  void scaleLeft(const VectorBase<Scalar> &row_scaling)
    { scaleLeftImpl(row_scaling); }

  /** \brief Right scales operator with diagonal scaling operator.
   *
   * \precondtion <tt>supportsScaleRight()==true</tt>
   */
  void scaleRight(const VectorBase<Scalar> &col_scaling)
    { scaleRightImpl(col_scaling); }

  //@}

protected:

  /** \name Protected virtual functions to be overridden by subclasses. */
  //@{

  /** \brief . */
  virtual bool supportsScaleLeftImpl() const = 0;

  /** \brief . */
  virtual bool supportsScaleRightImpl() const = 0;

  /** \brief . */
  virtual void scaleLeftImpl(const VectorBase<Scalar> &row_scaling) = 0;

  /** \brief . */
  virtual void scaleRightImpl(const VectorBase<Scalar> &col_scaling) = 0;

  //@}

};


}	// end namespace Thyra


#endif	// THYRA_SCALED_LINEAR_OP_BASE_HPP
