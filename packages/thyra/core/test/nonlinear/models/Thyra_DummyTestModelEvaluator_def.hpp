// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP
#define DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP


#include "Thyra_DummyTestModelEvaluator_decl.hpp"
#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {

// Nonmember constuctors


template<class Scalar>
Teuchos::RCP<DummyTestModelEvaluator<Scalar> >
dummyTestModelEvaluator(
  const Ordinal x_size,
  const ArrayView<const Ordinal> &p_sizes,
  const ArrayView<const Ordinal> &g_sizes,
  const bool supports_x_dot,
  const bool supports_x_dot_dot,
  const bool supports_extended_inargs,
  const bool supports_extended_outargs,
  const bool supports_derivatives
  )
{
  return Teuchos::rcp(new DummyTestModelEvaluator<Scalar>(x_size, p_sizes, g_sizes, supports_x_dot, supports_x_dot_dot,supports_extended_inargs,supports_extended_outargs,supports_derivatives));
}


// Initializers/Accessors


template<class Scalar>
DummyTestModelEvaluator<Scalar>::DummyTestModelEvaluator(
  const Ordinal x_size,
  const ArrayView<const Ordinal> &p_sizes,
  const ArrayView<const Ordinal> &g_sizes,
  const bool supports_x_dot,
  const bool supports_x_dot_dot,
  const bool supports_extended_inargs,
  const bool supports_extended_outargs,
  const bool supports_derivatives
  )
{
  
  typedef ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  x_space_ = defaultSpmdVectorSpace<Scalar>(x_size);

  p_space_.resize(p_sizes.size());
  for (Ordinal l = 0; l < p_sizes.size(); ++l) {
    p_space_[l] = defaultSpmdVectorSpace<Scalar>(p_sizes[l]);
  }
  
  f_space_ = x_space_;
  
  g_space_.resize(g_sizes.size());
  for (Ordinal j = 0; j < g_sizes.size(); ++j) {
    g_space_[j] = defaultSpmdVectorSpace<Scalar>(g_sizes[j]);
  }
  
  W_factory_ = defaultSerialDenseLinearOpWithSolveFactory<Scalar>();

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np_Ng(p_space_.size(),g_sizes.size());
  inArgs.setSupports(MEB::IN_ARG_x);
  if (supports_x_dot)
    inArgs.setSupports(MEB::IN_ARG_x_dot);
  if (supports_x_dot_dot)
    inArgs.setSupports(MEB::IN_ARG_x_dot_dot);
  inArgs.setSupports(MEB::IN_ARG_step_size);
  inArgs.setSupports(MEB::IN_ARG_stage_number);
  inArgs.template setSupports<Thyra::MockExtendedInArgs<Scalar> >(true);
  // test the removal of support
  if (!supports_extended_inargs)
    inArgs.template setSupports<Thyra::MockExtendedInArgs<Scalar> >(false);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(p_space_.size(), g_space_.size());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  outArgs.template setSupports<Thyra::MockExtendedOutArgs<Scalar> >(true);
  // test the removal of support
  if (!supports_extended_outargs)
    outArgs.template setSupports<Thyra::MockExtendedOutArgs<Scalar> >(false);
  if (supports_derivatives)
  {
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    outArgs.setHessianSupports(true);
#endif

    MEB::DerivativeSupport derivativeSupport(MEB::DERIV_MV_BY_COL);

    for (int j=0; j<g_space_.size(); ++j)
      outArgs.setSupports(MEB::OUT_ARG_DgDx,j,derivativeSupport);
  }
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  const RCP<VectorBase<Scalar> > x0 = createMember(x_space_);
  V_S(x0.ptr(), ST::zero());
  nominalValues_.set_x(x0);

}


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_p_space(int l) const
{
  return p_space_[l];
}


template<class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
DummyTestModelEvaluator<Scalar>::get_p_names(int /* l */) const
{
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_g_space(int j) const
{
  return g_space_[j];
}


template<class Scalar>
Teuchos::ArrayView<const std::string>
DummyTestModelEvaluator<Scalar>::get_g_names(int /* j */) const
{
  return g_names_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getLowerBounds() const
{
  return ModelEvaluatorBase::InArgs<Scalar>();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getUpperBounds() const
{
  return ModelEvaluatorBase::InArgs<Scalar>();
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DummyTestModelEvaluator<Scalar>::create_W_op() const
{
  return createNonconstSimpleDenseLinearOp<Scalar>(
    createMembers<Scalar>(f_space_, x_space_->dim())
    );
}


template<class Scalar>
Teuchos::RCP<PreconditionerBase<Scalar> >
DummyTestModelEvaluator<Scalar>::create_W_prec() const
{
  return nonconstUnspecifiedPrec<Scalar>(
    createNonconstSimpleDenseLinearOp<Scalar>(
      createMembers<Scalar>(f_space_, x_space_->dim())
      )
    );
}


template<class Scalar>
Teuchos::RCP<const LinearOpWithSolveFactoryBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &/* finalPoint */,
  const bool /* wasSolved */
  )
{
  // ToDo: Capture the final point and then provide in interface.
}

template<class Scalar>
void
DummyTestModelEvaluator<Scalar>::
change_p_size_incorrectly(const Ordinal new_size)
{
  using MEB = ModelEvaluatorBase;
  {
    MEB::InArgsSetup<Scalar> inArgs(prototypeInArgs_);
    inArgs.set_Np(new_size);
    prototypeInArgs_ = inArgs;
  }
  {
    MEB::OutArgsSetup<Scalar> outArgs(prototypeOutArgs_);
    outArgs.set_Np_Ng(new_size,g_space_.size());
    prototypeOutArgs_ = outArgs;
  }
  // forgot to call initializeDefaultBase() or resetDefaultBase()
}

template<class Scalar>
void
DummyTestModelEvaluator<Scalar>::
change_p_size_correctly(const Ordinal new_size)
{
  this->change_p_size_incorrectly(new_size);
  this->resetDefaultBase();
}

// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DummyTestModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &/* inArgs */,
  const ModelEvaluatorBase::OutArgs<Scalar> &/* outArgs */
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement to just copy inArgs and outArgs!
}


} // namespace Thyra


//
// Explicit instantiation macro
//
// Must be expanded from within the global namespace!
//

#define DUMMY_TEST_MODEL_EVALUATOR_INSTANT(SCALAR) \
  \
  template class DummyTestModelEvaluator<SCALAR >; \
  \
  template Teuchos::RCP<DummyTestModelEvaluator<SCALAR > > \
  dummyTestModelEvaluator( \
    const Ordinal x_size, \
    const ArrayView<const Ordinal> &p_sizes, \
    const ArrayView<const Ordinal> &g_sizes, \
    const bool supports_x_dot, \
    const bool supports_x_dot_dot,              \
    const bool supports_extended_inargs,        \
    const bool supports_extended_outargs,       \
    const bool supports_derivatives             \
    ); \


#endif // DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP
