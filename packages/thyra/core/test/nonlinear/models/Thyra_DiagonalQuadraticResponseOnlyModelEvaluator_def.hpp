/*
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
*/

#ifndef OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP
#define OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP


#include "Thyra_DiagonalQuadraticResponseOnlyModelEvaluator_decl.hpp"
#include "Thyra_DiagonalScalarProd.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


//
// DiagonalQuadraticResponseOnlyModelEvaluator
//


// Constructors, Initialization, Misc.


template<class Scalar>
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::DiagonalQuadraticResponseOnlyModelEvaluator(
  const int localDim,
  const RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm
  )
  :Np_(1), Ng_(1), comm_(comm), localDim_(localDim),
   nonlinearTermFactor_(0.0), g_offset_(0.0)
{

  typedef ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  TEUCHOS_ASSERT( localDim > 0 );

  // Get the comm
  if (is_null(comm_)) {
    comm_ = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();
  }

  // Locally replicated space for g
  g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, 1, 1);

  // Distributed space for p
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, localDim, -1);

  // Default solution
  const RCP<Thyra::VectorBase<Scalar> > ps = createMember<Scalar>(p_space_);
  V_S(ps.ptr(), ST::zero());
  ps_ = ps;

  // Default diagonal
  const RCP<Thyra::VectorBase<Scalar> > diag = createMember<Scalar>(p_space_);
  V_S(diag.ptr(), ST::one());
  diag_ = diag;
  diag_bar_ = diag;

  // Default inner product
  const RCP<Thyra::VectorBase<Scalar> > s_bar = createMember<Scalar>(p_space_);
  V_S(s_bar.ptr(), ST::one());
  s_bar_ = s_bar;

  // Default response offset
  g_offset_ = ST::zero();

}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setSolutionVector(
  const RCP<const Thyra::VectorBase<Scalar> > &ps)
{
  ps_ = ps.assert_not_null();
}


template<class Scalar>
const RCP<const Thyra::VectorBase<Scalar> >
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::getSolutionVector() const
{
  return ps_;
}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setDiagonalVector(
  const RCP<const Thyra::VectorBase<Scalar> > &diag)
{
  diag_ = diag;
  diag_bar_ = diag;
}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setDiagonalBarVector(
  const RCP<const Thyra::VectorBase<Scalar> > &diag_bar)
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::createMember;
  using Thyra::ele_wise_divide;
  using Thyra::V_S;

  diag_bar_ = diag_bar.assert_not_null();

  // Reset the scalar product for p_space!

  RCP<Thyra::VectorBase<Scalar> > s_bar = createMember<Scalar>(p_space_->clone());
  // NOTE: We have to clone the vector space in order to avoid creating a
  // circular reference between the space and the vector that defines the
  // scalar product for the vector space.

  // s_bar[i] = diag[i] / diag_bar[i] 
  V_S( s_bar.ptr(), ST::zero() );
  ele_wise_divide( ST::one(), *diag_, *diag_bar_, s_bar.ptr() );
  s_bar_ = s_bar;
  
  const RCP<Thyra::ScalarProdVectorSpaceBase<Scalar> > sp_p_space =
    rcp_dynamic_cast<Thyra::ScalarProdVectorSpaceBase<Scalar> >(p_space_, true);
  sp_p_space->setScalarProd(diagonalScalarProd<Scalar>(s_bar_));

}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setNonlinearTermFactor(
  const Scalar &nonlinearTermFactor)
{
  nonlinearTermFactor_ = nonlinearTermFactor;
}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setScalarOffset(
  const Scalar &g_offset)
{
  g_offset_ = g_offset;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
int DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::Np() const
{
  return Np_;
}


template<class Scalar>
int DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::Ng() const
{
  return Ng_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::get_p_space(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::get_g_space(int j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(Np_);
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeSupport DS;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np_,Ng_);
  outArgs.setSupports(MEB::OUT_ARG_DgDp, 0 ,0, MEB::DERIV_TRANS_MV_BY_ROW);
  return outArgs;
}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::get_mv;
  using Thyra::ConstDetachedSpmdVectorView;
  using Thyra::DetachedSpmdVectorView;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeMultiVector<Scalar> DMV;

  const ConstDetachedSpmdVectorView<Scalar> p(inArgs.get_p(0));
  const ConstDetachedSpmdVectorView<Scalar> ps(ps_);
  const ConstDetachedSpmdVectorView<Scalar> diag(diag_);
  const ConstDetachedSpmdVectorView<Scalar> s_bar(s_bar_);

  // g(p)
  if (!is_null(outArgs.get_g(0))) {
    Scalar g_val = ST::zero();
    for (Ordinal i = 0; i < p.subDim(); ++i) {
      const Scalar p_ps = p[i] - ps[i];
      g_val += diag[i] * p_ps*p_ps;
      if (nonlinearTermFactor_ != ST::zero()) {
        g_val += nonlinearTermFactor_ * p_ps * p_ps * p_ps;
      }
    }
    Scalar global_g_val;
    Teuchos::reduceAll<Ordinal, Scalar>(*comm_, Teuchos::REDUCE_SUM, 
      g_val, outArg(global_g_val) );
    DetachedSpmdVectorView<Scalar>(outArgs.get_g(0))[0] =
      as<Scalar>(0.5) * global_g_val + g_offset_;
  }

  // DgDp[i]
  if (!outArgs.get_DgDp(0,0).isEmpty()) {
    const RCP<Thyra::MultiVectorBase<Scalar> > DgDp_trans_mv =
      get_mv<Scalar>(outArgs.get_DgDp(0,0), "DgDp^T", MEB::DERIV_TRANS_MV_BY_ROW);
    const DetachedSpmdVectorView<Scalar> DgDp_grad(DgDp_trans_mv->col(0));
    for (Thyra::Ordinal i = 0; i < p.subDim(); ++i) {
      const Scalar p_ps = p[i] - ps[i];
      Scalar DgDp_grad_i = diag[i] * p_ps;
      if (nonlinearTermFactor_ != ST::zero()) {
        DgDp_grad_i += as<Scalar>(1.5) * nonlinearTermFactor_ * p_ps * p_ps;
      }
      DgDp_grad[i] = DgDp_grad_i / s_bar[i];

    }
  }
  
}


} // namespace Thyra


#endif // OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP
