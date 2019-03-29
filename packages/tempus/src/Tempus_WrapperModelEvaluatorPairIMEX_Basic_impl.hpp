// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
#define Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"


namespace Tempus {

template <typename Scalar>
void
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
initialize()
{
  using Teuchos::RCP;

  wrapperImplicitInArgs_  = this->createInArgs();
  wrapperImplicitOutArgs_ = this->createOutArgs();

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION( !(explicitModel_->get_x_space()->isCompatible(
                                       *(implicitModel_->get_x_space()))),
    std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
    "  Explicit and Implicit vector x spaces are incompatible!\n"
    "  Explicit vector x space = " << *(explicitModel_->get_x_space())<< "\n"
    "  Implicit vector x space = " << *(implicitModel_->get_x_space())<< "\n");

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION( !(explicitModel_->get_f_space()->isCompatible(
                                       *(implicitModel_->get_f_space()))),
    std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
    "  Explicit and Implicit vector f spaces are incompatible!\n"
    "  Explicit vector f space = " << *(explicitModel_->get_f_space())<< "\n"
    "  Implicit vector f space = " << *(implicitModel_->get_f_space())<< "\n");
}

template <typename Scalar>
void
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
setAppModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & /* me */)
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::setAppModel\n"
    "  should not be used.  One should instead use setExplicitModel,\n"
    "  setImplicitModel, or create a new WrapperModelEvaluatorPairIMEX.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
getAppModel() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::getAppModel\n"
    "  should not be used.  One should instead use getExplicitModel,\n"
    "  getImplicitModel, or directly use this WrapperModelEvaluatorPairIMEX\n"
    "  object.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_x_space() const
{
  return this->implicitModel_->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_g_space(int i) const
{
  return this->implicitModel_->get_g_space(i);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_p_space(int i) const
{
  return this->implicitModel_->get_p_space(i);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs = this->createInArgs();
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  //MEB::InArgsSetup<Scalar> inArgs(implicitModel_->createInArgs());
  MEB::InArgsSetup<Scalar> inArgs(implicitModel_->getNominalValues());
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs(implicitModel_->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}

template <typename Scalar>
void WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar>  & inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> >   x_dot = Thyra::createMember(get_x_space());
  timeDer_->compute(x, x_dot);

  MEB::InArgs<Scalar>  appImplicitInArgs (wrapperImplicitInArgs_);
  MEB::OutArgs<Scalar> appImplicitOutArgs(wrapperImplicitOutArgs_);
  appImplicitInArgs.set_x(x);
  appImplicitInArgs.set_x_dot(x_dot);
  for (int i=0; i<implicitModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      appImplicitInArgs.set_p(i, inArgs.get_p(i));
  }

  appImplicitOutArgs.set_f(outArgs.get_f());
  appImplicitOutArgs.set_W_op(outArgs.get_W_op());

  implicitModel_->evalModel(appImplicitInArgs,appImplicitOutArgs);
}

} // end namespace Tempus

#endif // Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
