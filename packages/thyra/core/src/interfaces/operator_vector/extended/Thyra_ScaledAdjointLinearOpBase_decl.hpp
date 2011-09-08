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

#ifndef THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_DECL_HPP
#define THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_DECL_HPP

#include "Thyra_LinearOpBase.hpp"


namespace Thyra {


/** \brief Base class for <tt>LinearOpBase</tt> decorator subclasses that wrap
a <tt>LinearOpBase</tt> object and adds on an extra scaling factor and/or a
new transpose enum.

This interface class represents a scaled, adjointed (transposed) linear
operator <tt>M</tt> of the form:

\verbatim
 
  M = scalar * op(Op)
\endverbatim

where <tt>Op</tt> is another <tt>LinearOpBase</tt> object, <tt>scalar</tt> is
a <tt>Scalar</tt>, and the operation <tt>op(Op)</tt> is specified by a
<tt>EOpTransp</tt> and is given as <tt>op(Op) = Op</tt> (<tt>NOTRANS</tt>), or
<tt>op(Op) = Op^T</tt> (<tt>TRANS</tt>), or <tt>op(Op) = Op^H</tt>
(<tt>CONJTRANS</tt>).

\ingroup Thyra_Op_Vec_extended_interfaces_code_grp

*/
template<class Scalar>
class ScaledAdjointLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

#ifdef THYRA_INJECT_USING_DECLARATIONS
  using LinearOpBase<Scalar>::apply;
#endif

  /** @name Pure virtual functions to be overridden in subclasses */
  //@{

  /** \brief Return the overall scale factor.
   */
  virtual Scalar overallScalar() const = 0;

  /** \brief Return the overall transpose (adjoint) enum.
   */
  virtual EOpTransp overallTransp() const = 0;

  /** \brief Return the non-const original linear operator <tt>origOp</tt>.
   *
   * Note that <tt>*this</tt> is only guaranteed to be fully modified once the
   * returned RCP goes away.
   */
  virtual RCP<LinearOpBase<Scalar> > getNonconstOrigOp() = 0;

  /** \brief Return the const original linear operator <tt>origOp</tt>.
   */
  virtual RCP<const LinearOpBase<Scalar> > getOrigOp() const = 0;

  //@}

};


/** \brief Extract the <tt>overallScalar</tt>, <tt>overallTransp</tt> and
 * <tt>const</tt> <tt>origOp</tt> from a <tt>const</tt>
 * <tt>LinearOpBase</tt> object.
 *
 * \param Op [in] The input, possibly scaled and/or adjoined, linear operator
 *
 * \param scalar [out] The overall scaling factor.
 *
 * \param transp [out] The overall adjoint (transposition) enum.
 *
 * \param origOp [out] The underlying, non-scaled and non-adjoined linear
 * operator.  This pointer returns a non-persisting relationship that is to be
 * used and then immediately forgotten.
 *
 * Preconditions:<ul>
 * <li><tt>scalar!==NULL</tt>
 * <li><tt>transp!==NULL</tt>
 * <li><tt>origOp!==NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>*origOp!==NULL</tt>
 * </ul>
 *
 * The purpose of this function is to strip off the
 * <tt>ScaledAdjointLinearOpBase</tt> wrapper and get at the underlying linear
 * operator for the purpose of further dynamic casting to some more derived
 * interface.
 *
 * The implementation of this function is not too complicated and is
 * appropriate for study.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
void unwrap(
  const LinearOpBase<Scalar> &Op,
  Scalar *scalar,
  EOpTransp *transp,
  const LinearOpBase<Scalar>* *origOp
  );


/** \brief Extract the <tt>overallScalar</tt>, <tt>overallTransp</tt> and
 * <tt>RCP</tt> wrapped <tt>const</tt> <tt>origOp</tt> from a
 * <tt>RCP</tt> wrapped <tt>const</tt> <tt>LinearOpBase</tt>
 * object.
 *
 * \param Op [in] The input, possibly scaled and/or adjoined, linear operator
 *
 * \param scalar [out] The overall scaling factor.
 *
 * \param transp [out] The overall adjoint (transposition) enum.
 *
 * \param origOp [out] The underlying, non-scaled and non-adjoined linear
 * operator.  This pointer returns a non-persisting relationship that is to be
 * used and then immediately forgotten.
 *
 * Preconditions:<ul>
 * <li><tt>scalar!==NULL</tt>
 * <li><tt>transp!==NULL</tt>
 * <li><tt>origOp!==NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>*origOp!==NULL</tt>
 * </ul>
 *
 * The purpose of this function is to strip off the
 * <tt>ScaledAdjointLinearOpBase</tt> wrapper and get at the underlying linear
 * operator for the purpose of further dynamic casting to some more derived
 * interface.
 *
 * The implementation of this function is not too complicated and is
 * appropriate for study.
 *
 * \ingroup Thyra_Op_Vec_ScaledAdjointedLinearOp_helpers_grp
 */
template<class Scalar>
void unwrap(
  const RCP<const LinearOpBase<Scalar> > &Op,
  Scalar *scalar,
  EOpTransp *transp,
  RCP<const LinearOpBase<Scalar> > *origOp
  );


} // namespace Thyra


#endif	// THYRA_SCALED_ADJOINT_LINEAR_OP_BASE_DECL_HPP
