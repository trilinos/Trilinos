//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_REACTION_MODEL_IMPL_HPP
#define TEMPUS_TEST_REACTION_MODEL_IMPL_HPP

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

namespace Tempus_Test {

template <class Scalar>
ReactionModel<Scalar>::ReactionModel(Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_     = false;
  dim_               = 3;
  Np_                = 3;     // Number of parameter vectors (p, dx/dp, dx_dot/dp)
  np_                = 3;     // Number of parameters in this vector (3)
  Ng_                = 1;     // Number of observation functions (1)
  ng_                = dim_;  // Number of elements in this observation function ( == x )
  acceptModelParams_ = false;
  useDfDpAsTangent_  = false;
  haveIC_            = true;
  lambda0_         = 1.0;
  lambda1_         = 0.5;
  lambda2_         = 0.25;
  x0_ic_             = 1.0;
  x1_ic_             = 0.0;
  x2_ic_             = 0.0;
  t0_ic_             = 0.0;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);

  setParameterList(pList_);
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ReactionModel<Scalar>::getExactSolution(
    double t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  double exact_t                                   = t;
  inArgs.set_t(exact_t);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    Scalar sum_k = Scalar(0.0);
    Scalar prod_l = Scalar(1.0);
    Scalar lmbd_prod = Scalar(1.0);

    for (int num = 1; num < dim_ + 1; ++num)
    {
        lmbd_prod = 1.0;
        for (int i = 0; i < num - 1; i++)
        {
            lmbd_prod *= -Lmat_[i][i];
        }
        sum_k = 0;
        for (int i = 0; i < num; i++)
        {
            prod_l = 1.0;
            for (int j = 0; j < num; j++)
            {
                if (j != i)
                    prod_l *= (-Lmat_[j][j]
                                + Lmat_[i][i]);
            }
            sum_k += exp(+Lmat_[i][i] * t) / prod_l;
        }
            exact_x_view[num - 1] = x0_ic_ * lmbd_prod * sum_k;
    }
  }
  inArgs.set_x(exact_x);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x_dot = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    for (int i = 0; i < dim_; ++i) 
    {
      Scalar sum = Scalar(0.0);
      for (int j = 0; j < dim_; ++j) 
      {
        sum += Lmat_[i][j] * exact_x_view[j];
      }
      exact_x_dot_view[i] = sum;
    }
  }
  inArgs.set_x_dot(exact_x_dot);
  return (inArgs);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ReactionModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ReactionModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ReactionModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
ReactionModel<Scalar>::create_W() const
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
Teuchos::RCP<Thyra::LinearOpBase<Scalar> > ReactionModel<Scalar>::create_W_op()
    const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);

  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ReactionModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ReactionModel<Scalar>::createInArgs()
    const
{
  setupInOutArgs_();
  return inArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ReactionModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void ReactionModel<Scalar>::evalModelImpl(
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

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view(*x_in);

  // double t = inArgs.get_t();
  // Scalar a = a_;
  // Scalar f = f_;
  // Scalar L = L_;
  // if (acceptModelParams_) {
  //   const RCP<const VectorBase<Scalar> > p_in =
  //       inArgs.get_p(0).assert_not_null();
  //   Thyra::ConstDetachedVectorView<Scalar> p_in_view(*p_in);
  //   a = p_in_view[0];
  //   f = p_in_view[1];
  //   L = p_in_view[2];
  // }

  Scalar beta = inArgs.get_beta();

  const RCP<VectorBase<Scalar> > f_out          = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();

  if (!inArgs.get_x_dot().is_null()) {
    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const VectorBase<Scalar> > x_dot_in;
    x_dot_in     = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      for (int i = 0; i < 3; ++i) {
        Scalar Lu_i = Scalar(0.0);
        for (int j = 0; j < 3; ++j)
          Lu_i += Lmat_[i][j] * x_in_view[j];
        f_out_view[i] = x_dot_in_view[i] - Lu_i;
      }
    }
    if (!is_null(W_out)) {
  RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out, true);
  Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          matrix_view(i,j) = -beta * Lmat_[i][j];
        matrix_view(i,i) += alpha;
      }
}
  }
}

// private

template <class Scalar>
void ReactionModel<Scalar>::setupInOutArgs_() const
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
      x_ic_view[0] = x0_ic_;
      x_ic_view[1] = x1_ic_;
      x_ic_view[2] = x2_ic_;
    }
    nominalValues_.set_x(x_ic);
    if (acceptModelParams_) {
      const RCP<Thyra::VectorBase<Scalar> > p_ic = createMember(p_space_);
      {
        Thyra::DetachedVectorView<Scalar> p_ic_view(*p_ic);
        p_ic_view[0] = lambda0_;
        p_ic_view[1] = lambda1_;
        p_ic_view[2] = lambda2_;
      }
      nominalValues_.set_p(0, p_ic);
    }
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      for (int i = 0; i < dim_; ++i) {
        Scalar sum = Scalar(0.0);
        for (int j = 0; j < dim_; ++j) {
          sum += Lmat_[i][j] * x_ic_view[j];
        }
        x_dot_ic_view[i] = sum;
      }
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

template <class Scalar>
void ReactionModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("ReactionModel"));
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
  lambda0_                 = get<Scalar>(*pl, "Coeff L1");
  lambda1_                 = get<Scalar>(*pl, "Coeff L2");
  lambda2_                 = get<Scalar>(*pl, "Coeff L3");
  x0_ic_             = get<Scalar>(*pl, "IC x0");
  x1_ic_             = get<Scalar>(*pl, "IC x1");
  x2_ic_             = get<Scalar>(*pl, "IC x2");
  t0_ic_             = get<Scalar>(*pl, "IC t0");
  buildDecayMatrix_();
  setupInOutArgs_();
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
ReactionModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters", false);
    pl->set("Provide nominal values", true);
    pl->set("Use DfDp as Tangent", false);
    pl->set<std::string>("Output File Name", "Tempus_BDF2_Reaction");
    Teuchos::setDoubleParameter("Coeff L1", 1.0, "Coefficient L1 in model", &*pl);
    Teuchos::setDoubleParameter("Coeff L2", 0.1, "Coefficient L2 in model", &*pl);
    Teuchos::setDoubleParameter("Coeff L3", 0.01, "Coefficient L3 in model", &*pl);
    Teuchos::setDoubleParameter("IC x0", 1.0, "Initial Condition for x0", &*pl);
    Teuchos::setDoubleParameter("IC x1", 0.0, "Initial Condition for x1", &*pl);
    Teuchos::setDoubleParameter("IC x2", 0.0, "Initial Condition for x2", &*pl);
    Teuchos::setDoubleParameter("IC t0", 0.0, "Initial time t0", &*pl);
    Teuchos::setIntParameter("Number of Time Step Sizes", 1,
                             "Number time step sizes for convergence study",
                             &*pl);
    validPL = pl;
  }
  return validPL;
}

template <class Scalar>
void ReactionModel<Scalar>::buildDecayMatrix_()
{
  Lmat_[0][0] = -lambda0_;  Lmat_[0][1] = 0.0;      Lmat_[0][2] = 0.0;
  Lmat_[1][0] =  lambda0_;  Lmat_[1][1] = -lambda1_; Lmat_[1][2] = 0.0;
  Lmat_[2][0] = 0.0;       Lmat_[2][1] =  lambda1_; Lmat_[2][2] = -lambda2_;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_REACTION_MODEL_IMPL_HPP
