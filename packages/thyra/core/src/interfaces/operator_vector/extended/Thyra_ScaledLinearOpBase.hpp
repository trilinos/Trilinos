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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
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
