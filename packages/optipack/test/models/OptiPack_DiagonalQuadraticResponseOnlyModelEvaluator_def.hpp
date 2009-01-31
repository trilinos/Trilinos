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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP
#define OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP


#include "OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator_decl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace OptiPack {


//
// Implementation
//


// Constructors, Initialization, Misc.


template<class Scalar>
DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::DiagonalQuadraticResponseOnlyModelEvaluator(
  const int localDim
  )
  :Np_(1), Ng_(1)
{

  typedef ScalarTraits<Scalar> ST;

  // Get the comm
  comm_ = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

  // Parallel space for p
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, localDim, -1);

  // Locally replicated space for g
  g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, 1, 1);

  // Default solution
  const RCP<Thyra::VectorBase<Scalar> > ps = createMember(p_space_);
  V_S(ps.ptr(), ST::zero());
  ps_ = ps;

  // Default offset
  g_offset_ = ST::zero();

}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setSolutionVector(
  const RCP<const Thyra::VectorBase<Scalar> > &ps)
{
  ps_ = ps.assert_not_null();
}


template<class Scalar>
void DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>::setScalarOffset(
  const Scalar &s)
{
  g_offset_ = s;
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

  if (!is_null(outArgs.get_g(0))) {
    Scalar g_val = ST::zero();
    for (Ordinal i = 0; i < p.subDim(); ++i) {
      const Scalar p_ps = p[i] - ps[i];
      g_val += p_ps*p_ps;
    }
    Scalar global_g_val;
    Teuchos::reduceAll<Ordinal, Scalar>(*comm_, Teuchos::REDUCE_SUM, g_val,
      outArg(global_g_val));
    DetachedSpmdVectorView<Scalar>(outArgs.get_g(0))[0] =
      as<Scalar>(0.5) * global_g_val + g_offset_;
  }

  if (!outArgs.get_DgDp(0,0).isEmpty()) {
    const RCP<Thyra::MultiVectorBase<Scalar> > DgDp_trans_mv =
      get_mv<Scalar>(outArgs.get_DgDp(0,0), "DgDp^T", MEB::DERIV_TRANS_MV_BY_ROW);
    const DetachedSpmdVectorView<Scalar> DgDp_grad(DgDp_trans_mv->col(0));
    for (Thyra::Ordinal i = 0; i < p.subDim(); ++i) {
      DgDp_grad[i] = p[i] - ps[i];
    }
  }
  
}


} // namespace OptiPack


#endif // OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DEF_HPP
