//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_NONAUTOSRC_MODEL_IMPL_HPP
#define TEMPUS_TEST_NONAUTOSRC_MODEL_IMPL_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

#include <iostream>
#include <cmath>

namespace Tempus_Test {

template <class Scalar>
NonAutoSrcModel<Scalar>::NonAutoSrcModel(Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_     = false;
  dim_               = 1;
  Np_                = 0;
  np_                = 0;
  Ng_                = 1;
  ng_                = dim_;
  acceptModelParams_ = false;
  useDfDpAsTangent_  = false;
  haveIC_            = true;

  gamma_       = 1.0;
  xs_          = 1.0;
  avogadro_    = 6.02214076e23;

  flux_amp_    = 1.0;
  a_           = 1.0;
  b_           = 1.0;
  kappa_       = 1.0;
  flux_beta_   = 1.0;

  initial_amp_ = 0.0;
  t0_ic_       = 0.0;

  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);

  setParameterList(pList_);
}

template <class Scalar>
Scalar NonAutoSrcModel<Scalar>::sourceFlux_(const Scalar time) const
{
  return flux_amp_ / Scalar(2.0)
         * ((std::tanh(a_ * time - b_) + Scalar(1.0))
            - (std::tanh(kappa_ * time - flux_beta_) + Scalar(1.0)));
}

template <class Scalar>
Scalar NonAutoSrcModel<Scalar>::nonautoSource_(const Scalar time) const
{
  return sourceFlux_(time) * gamma_ * xs_ / avogadro_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
NonAutoSrcModel<Scalar>::getExactSolution(double t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  inArgs.set_t(t);

  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x =
      createMember(x_space_);

  {
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);

    const Scalar time = Scalar(t);

    const Scalar const_s = flux_amp_ * gamma_ * xs_ / avogadro_;

    const Scalar const_c =
        -const_s * std::log(std::cosh(b_)) / (Scalar(2.0) * a_)
        + const_s * std::log(std::cosh(-flux_beta_)) / (Scalar(2.0) * kappa_)
        + initial_amp_;

    exact_x_view[0] =
        const_s * std::log(std::cosh(b_ - a_ * time)) / (Scalar(2.0) * a_)
        - const_s * std::log(std::cosh(kappa_ * time - flux_beta_))
              / (Scalar(2.0) * kappa_)
        + const_c;
  }

  inArgs.set_x(exact_x);

  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x_dot =
      createMember(x_space_);

  {
    Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
    exact_x_dot_view[0] = nonautoSource_(Scalar(t));
  }

  inArgs.set_x_dot(exact_x_dot);

  return inArgs;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
NonAutoSrcModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
NonAutoSrcModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
NonAutoSrcModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
NonAutoSrcModel<Scalar>::create_W() const
{
  using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<Scalar> > multivec =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,
                                                                   true);
    {
      // Set matrix = I
      for (int j = 0; j < dim_; ++j) {
        RCP<Thyra::VectorBase<Scalar> > e = Thyra::createMember(x_space_);
        {
          Thyra::DetachedVectorView<Scalar> ev(*e);
          for (int i = 0; i < dim_; ++i) ev[i] = Scalar(0.0);
          ev[j] = Scalar(1.0);
        }
        V_V(multivec->col(j).ptr(), *e);
      }
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
      Thyra::linearOpWithSolve<Scalar>(*W_factory, matrix);
  return W;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> > NonAutoSrcModel<Scalar>::create_W_op()
    const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);

  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
NonAutoSrcModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> NonAutoSrcModel<Scalar>::createInArgs()
    const
{
  setupInOutArgs_();
  return inArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
NonAutoSrcModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void NonAutoSrcModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::VectorBase;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  const RCP<VectorBase<Scalar> > f_out          = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();

  const Scalar time = inArgs.get_t();
  const Scalar source = nonautoSource_(time);

  if (!inArgs.get_x_dot().is_null()) {
    RCP<const VectorBase<Scalar> > x_dot_in =
        inArgs.get_x_dot().assert_not_null();

    Scalar alpha = inArgs.get_alpha();

    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);

      f_out_view[0] = x_dot_in_view[0] - source;
    }

    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out, true);

      Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);

      matrix_view(0, 0) = alpha;
    }
  }
}

// private

template <class Scalar>
void NonAutoSrcModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) {
    return;
  }

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;
  {
    // Set up prototypical InArgs
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_t);
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs.setSupports(MEB::IN_ARG_beta);
    inArgs.setSupports(MEB::IN_ARG_x_dot);
    inArgs.setSupports(MEB::IN_ARG_alpha);
    if (acceptModelParams_) {
      inArgs.set_Np(Np_);
    }
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    outArgs_ = outArgs;
  }

  // Set up nominal values
  nominalValues_ = inArgs_;
  if (haveIC_) {
    nominalValues_.set_t(t0_ic_);
    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      x_ic_view[0] = initial_amp_;
    }
    nominalValues_.set_x(x_ic);
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      x_dot_ic_view[0] = nonautoSource_(t0_ic_);
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

template <class Scalar>
void NonAutoSrcModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("NonAutoSrcModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();
  bool acceptModelParams         = get<bool>(*pl, "Accept model parameters");
  bool haveIC                    = get<bool>(*pl, "Provide nominal values");
  bool useDfDpAsTangent          = get<bool>(*pl, "Use DfDp as Tangent");
  if ((acceptModelParams != acceptModelParams_) || (haveIC != haveIC_)) {
    isInitialized_ = false;
  }
  acceptModelParams_ = acceptModelParams;
  haveIC_            = haveIC;
  useDfDpAsTangent_  = useDfDpAsTangent;
  gamma_       = get<Scalar>(*pl, "Gamma");
  xs_          = get<Scalar>(*pl, "Cross Section");
  avogadro_    = get<Scalar>(*pl, "Avogadro Number");

  flux_amp_    = get<Scalar>(*pl, "Flux Amplitude");
  a_           = get<Scalar>(*pl, "Flux a");
  b_           = get<Scalar>(*pl, "Flux b");
  kappa_       = get<Scalar>(*pl, "Flux kappa");
  flux_beta_   = get<Scalar>(*pl, "Flux beta");

  initial_amp_ = get<Scalar>(*pl, "Initial Species");
  t0_ic_       = get<Scalar>(*pl, "IC t0");
  setupInOutArgs_();
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
NonAutoSrcModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters", false);
    pl->set("Provide nominal values", true);
    pl->set("Use DfDp as Tangent", false);
    pl->set<std::string>("Output File Name", "Tempus_EPI_NonAutoSrc");
    Teuchos::setDoubleParameter("Gamma", 1.0,
                            "Gamma coefficient", &*pl);

Teuchos::setDoubleParameter("Cross Section", 1.0,
                            "Area cross section", &*pl);

Teuchos::setDoubleParameter("Avogadro Number", 6.02214076e23,
                            "Avogadro number", &*pl);

Teuchos::setDoubleParameter("Flux Amplitude", 1.0,
                            "Source flux amplitude", &*pl);

Teuchos::setDoubleParameter("Flux a", 1.0,
                            "First tanh ramp coefficient a", &*pl);

Teuchos::setDoubleParameter("Flux b", 1.0,
                            "First tanh shift b", &*pl);

Teuchos::setDoubleParameter("Flux kappa", 1.0,
                            "Second tanh ramp coefficient kappa", &*pl);

Teuchos::setDoubleParameter("Flux beta", 1.0,
                            "Second tanh shift beta", &*pl);

Teuchos::setDoubleParameter("Initial Species", 0.0,
                            "Initial value of the species", &*pl);

Teuchos::setDoubleParameter("IC t0", 0.0,
                            "Initial time t0", &*pl);
    Teuchos::setIntParameter("Number of Time Step Sizes", 1,
                             "Number time step sizes for convergence study",
                             &*pl);
    validPL = pl;
  }
  return validPL;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_NONAUTOSRC_MODEL_IMPL_HPP
