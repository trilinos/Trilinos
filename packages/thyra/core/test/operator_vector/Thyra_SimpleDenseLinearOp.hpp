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

#ifndef THYRA_SIMPLE_DENSE_LINEAR_OP_HPP
#define THYRA_SIMPLE_DENSE_LINEAR_OP_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_MultiVectorBase.hpp"


namespace Thyra {


/** \brief Simple test LinearOpBase subclass implemented in terms of a MultiVectorBase object.
 */
template<class Scalar>
class SimpleDenseLinearOp
  : public Thyra::LinearOpDefaultBase<Scalar>
{
public:

  /** \brief Construct to uninitialized. */
  SimpleDenseLinearOp() {}
  
  /** Initialize given a fully formed MultiVectorBase object.
   */
  void initialize(const RCP<MultiVectorBase<Scalar> > &mv)
    {
      mv_ = mv;
    }

protected:

  // Overridden from LinearOpBase

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > range() const
    { return mv_->range(); }
 
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const
    { return mv_->domain(); }

  /** \brief . */
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const
    { return true; }

  /** \brief . */
  void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X_in,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta
    ) const
    {
      mv_->apply(M_trans, X_in, Y_inout, alpha, beta);
    }

private:

  RCP<MultiVectorBase<Scalar> > mv_;

};


/** Non-member constructor.
 *
 * \relates SimpleDenseLinearOp
 */
template<class Scalar>
RCP<SimpleDenseLinearOp<Scalar> >
createNonconstSimpleDenseLinearOp(
  const RCP<MultiVectorBase<Scalar> > &mv)
{
  const RCP<SimpleDenseLinearOp<Scalar> > sdlo =
    Teuchos::rcp(new SimpleDenseLinearOp<Scalar>);
  sdlo->initialize(mv);
  return sdlo;
}


} // namespace Thyra


#endif // THYRA_SIMPLE_DENSE_LINEAR_OP_HPP
