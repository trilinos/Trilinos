/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DEF_HPP
#define OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DEF_HPP


#include "OptiPack_Types.hpp"
#include "OptiPack_LineSearchPointEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "GlobiPack_MeritFunc1DBase.hpp"
#include "Teuchos_Assert.hpp"


namespace OptiPack {


// Constructor/Initializers/Accessors


template<typename Scalar>
UnconstrainedOptMeritFunc1D<Scalar>::UnconstrainedOptMeritFunc1D()
  : paramIndex_(-1),
    responseIndex_(-1)
{}


template<typename Scalar>
void UnconstrainedOptMeritFunc1D<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const int paramIndex,
  const int responseIndex
  )
{
#ifdef TEUCHOS_DEBUG
  model.assert_not_null();
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( paramIndex, 0, model->Np() );
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( responseIndex, 0, model->Ng() );
#endif
  model_ = model;
  paramIndex_ = paramIndex;
  responseIndex_ = responseIndex;
}


template<typename Scalar>
void UnconstrainedOptMeritFunc1D<Scalar>::setEvaluationQuantities(
  const RCP<const LineSearchPointEvaluatorBase<Scalar> > &pointEvaluator,
  const RCP<Thyra::VectorBase<Scalar> > &p,
  const RCP<Thyra::VectorBase<Scalar> > &g_vec,
  const RCP<Thyra::VectorBase<Scalar> > &g_grad_vec
  )
{
#ifdef TEUCHOS_DEBUG
  pointEvaluator.assert_not_null();
  p.assert_not_null();
  g_vec.assert_not_null();
  // ToDo: Check that pointEvaluator, p, g_vec, and g_grad_vec are compatible!
#endif
  pointEvaluator_ = pointEvaluator;
  p_ = p;
  g_vec_ = g_vec;
  g_grad_vec_ = g_grad_vec; // Can be null and that is okay
}


// Overridden from MeritFunc1DBase


template<typename Scalar>
bool UnconstrainedOptMeritFunc1D<Scalar>::supportsDerivEvals() const
{
  return !is_null(g_grad_vec_);
}


template<typename Scalar>
void UnconstrainedOptMeritFunc1D<Scalar>::eval(
  const ScalarMag &alpha, const Ptr<ScalarMag> &phi,
  const Ptr<ScalarMag> &Dphi ) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::get_ele;
  using Thyra::eval_g;
  using Thyra::eval_g_DgDp;
  pointEvaluator_->computePoint(alpha, p_.ptr());
  if (!is_null(g_grad_vec_)) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
      "Error, g_grad_vec has not been implemented yet!.");
  }
  else {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(is_null(Dphi));
#endif
    eval_g( *model_, paramIndex_, *p_, responseIndex_, g_vec_.ptr() );
    *phi = get_ele(*g_vec_, 0);
  }
}


} // namespace OptiPack


#endif // OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DEF_HPP
