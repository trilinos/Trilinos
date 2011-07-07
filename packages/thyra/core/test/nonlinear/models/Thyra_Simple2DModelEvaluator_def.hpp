
#ifndef THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP
#define THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP


#include "Thyra_Simple2DModelEvaluator_decl.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


// Nonmember constuctors


template<class Scalar>
Teuchos::RCP<Simple2DModelEvaluator<Scalar> >
simple2DModelEvaluator()
{
  return Teuchos::rcp(new Simple2DModelEvaluator<Scalar>);
}


// Initializers/Accessors


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_d(const Scalar &d)
{
  d_ = d;
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Simple2DModelEvaluator<Scalar>::create_W_op() const
{
  return Thyra::createMembers(f_space_, x_space_->dim());
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
Simple2DModelEvaluator<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Simple2DModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void Simple2DModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  TEST_FOR_EXCEPT(true);
}


// private


template<class Scalar>
Simple2DModelEvaluator<Scalar>::Simple2DModelEvaluator()
  : x_space_(Thyra::defaultSpmdVectorSpace<Scalar>(2)),
    f_space_(x_space_),
    W_factory_(Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>()),
    d_(0.0),
    p_(x_space_->dim(), Teuchos::ScalarTraits<Scalar>::zero()),
    showGetInvalidArg_(false)
{

  using Teuchos::RCP;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  x0_ = createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());
  nominalValues_.set_x(x0_);

  set_d(10.0);
  set_p(Teuchos::tuple<Scalar>(2.0, 0.0)());
  set_x0(Teuchos::tuple<Scalar>(1.0, 1.0)());

}


} // namespace Thyra


//
// Explicit instantiation macro
//
// Must be expanded from within the global namespace!
//

#define SIMPLE_2D_MODEL_EVALUATOR_INSTANT(SCALAR) \
  \
  template class Simple2DModelEvaluator<SCALAR >; \
  \
  template Teuchos::RCP<Simple2DModelEvaluator<SCALAR > > \
  simple2DModelEvaluator(); \


#endif // THYRA_SIMPLE_2D_MODEL_EVALUATOR_DEF_HPP
