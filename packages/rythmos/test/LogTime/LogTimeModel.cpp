//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "LogTimeModel.hpp"

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"

#ifdef LOGTIMEMODEL_DEBUG
#include <iostream>
#endif // LOGTIMEMODEL_DEBUG

namespace {
  const std::string Implicit_name = "Implicit model formulation";
  const bool Implicit_default = false;

  const std::string AcceptModelParams_name = "Accept model parameters";
  const bool AcceptModelParams_default = false;

  const std::string HaveIC_name = "Provide nominal values";
  const bool HaveIC_default = true;

  const std::string Coeff_a_name = "Coeff a";
  const double Coeff_a_default = 1.4;

  const std::string Coeff_b_name = "Coeff b";
  const double Coeff_b_default = 0.0001;

  const std::string Coeff_c_name = "Coeff c";
  const double Coeff_c_default = 0.1;

  const std::string Coeff_d_name = "Coeff d";
  const double Coeff_d_default = 1.0e-36;

  const std::string IC_x_name = "IC x_0";
  const double IC_x_default = 0.0;

  const std::string IC_t0_name = "IC t_0";
  const double IC_t0_default = 0.0;

} // namespace

namespace Rythmos {

// non-member Constructor
RCP<LogTimeModel> logTimeModel()
{
  RCP<LogTimeModel> model = rcp(new LogTimeModel);
  return(model);
}

// non-member Constructor
RCP<LogTimeModel> logTimeModel(bool implicit)
{
  RCP<LogTimeModel> model = logTimeModel();
  model->setImplicitFlag(implicit);
  return(model);
}

// Constructor
LogTimeModel::LogTimeModel()
{
  isInitialized_ = false;
  dim_ = 1;
  Np_ = 1; // Number of parameter vectors (1)
  np_ = 1; // Number of parameters in this vector (1)
  Ng_ = 1; // Number of observation functions (1)
  ng_ = 1; // Number of elements in this observation function (1)
  isImplicit_ = Implicit_default;
  acceptModelParams_ = AcceptModelParams_default;
  haveIC_ = HaveIC_default;
  a_ = Coeff_a_default;
  b_ = Coeff_b_default;
  c_ = Coeff_c_default;
  d_ = Coeff_d_default;
  x_ic_ = IC_x_default;
  t0_ic_ = IC_t0_default;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  // Create p_space and g_space
  p_space_ = Thyra::defaultSpmdVectorSpace<double>(np_);
  g_space_ = Thyra::defaultSpmdVectorSpace<double>(ng_);
}

void LogTimeModel::setImplicitFlag(bool implicit)
{
  if (isImplicit_ != implicit) {
    isInitialized_ = false;
  }
  isImplicit_ = implicit;
  setupInOutArgs_();
}

ModelEvaluatorBase::InArgs<double> LogTimeModel::getExactSolution(double t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_x_view(*exact_x);
    exact_x_view[0] = a_*(b_*pow(t,4) + c_*pow(t,4.5))
                      *pow(b_ + pow(t,0.5),-1)*pow(d_ + pow(t,4),-1);
  }
  inArgs.set_x(exact_x);
  if (isImplicit_) {
    RCP<VectorBase<double> > exact_x_dot = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_x_dot_view(*exact_x_dot);
      exact_x_dot_view[0] =
        ( a_*pow(t,3)*(8*c_*d_*t + 8*d_*pow(b_,2)
          + b_*pow(t,0.5)*((7 + 9*c_)*d_ + (-1 + c_)*pow(t,4)))
        *pow(b_ + pow(t,0.5),-2)*pow(d_ + pow(t,4),-2))/2.0;
    }
    inArgs.set_x_dot(exact_x_dot);
  }
  return(inArgs);
}


RCP<const Thyra::VectorSpaceBase<double> >
LogTimeModel::get_x_space() const
{
  return x_space_;
}


RCP<const Thyra::VectorSpaceBase<double> >
LogTimeModel::get_f_space() const
{
  return f_space_;
}


ModelEvaluatorBase::InArgs<double>
LogTimeModel::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  return nominalValues_;
}




RCP<Thyra::LinearOpWithSolveBase<double> >
LogTimeModel::create_W() const
{
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = this->get_W_factory();
  RCP<Thyra::LinearOpBase<double> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<double> > multivec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(matrix,true);
    {
      RCP<Thyra::VectorBase<double> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        vec_view[0] = 1.0;
      }
      V_V(multivec->col(0).ptr(),*vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<double> > W =
    Thyra::linearOpWithSolve<double>(
      *W_factory,
      matrix
      );
  return W;
}
//RCP<Thyra::LinearOpWithSolveBase<double> >
//LogTimeModel::create_W() const
//{
//  return Thyra::multiVectorLinearOpWithSolve<double>();
//}





RCP<Thyra::LinearOpBase<double> >
LogTimeModel::create_W_op() const
{
  RCP<Thyra::MultiVectorBase<double> > matrix = Thyra::createMembers(x_space_, dim_);
  return(matrix);
}




RCP<const Thyra::LinearOpWithSolveFactoryBase<double> >
LogTimeModel::get_W_factory() const
{
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > W_factory =
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  return W_factory;
}
//  RCP<const Thyra::LinearOpBase<double> > fwdOp = this->create_W_op();
//  RCP<Thyra::LinearOpWithSolveBase<double> > W =
//    Thyra::linearOpWithSolve<double>(
//      *W_factory,
//      fwdOp
//      );
//  W_factory->initializeOp(
//      Thyra::defaultLinearOpSource<double>(fwdOp),
//      &*W,
//      Thyra::SUPPORT_SOLVE_UNSPECIFIED
//      );



ModelEvaluatorBase::InArgs<double>
LogTimeModel::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


ModelEvaluatorBase::OutArgs<double>
LogTimeModel::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}


void LogTimeModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );

  const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<double> x_in_view( *x_in );

  //double t = inArgs.get_t();
  double a = a_;
  double b = b_;
  double c = c_;
  double d = d_;
  if (acceptModelParams_) {
    const RCP<const VectorBase<double> > p_in = inArgs.get_p(0).assert_not_null();
    Thyra::ConstDetachedVectorView<double> p_in_view( *p_in );
    a = p_in_view[0];
    b = p_in_view[1];
    c = p_in_view[2];
    d = p_in_view[3];
  }

  RCP<const VectorBase<double> > x_dot_in;
  //double beta = inArgs.get_beta();
  double alpha = -1.0;
  if (isImplicit_) {
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    alpha = inArgs.get_alpha();
  }

  const RCP<VectorBase<double> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<double> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<double> > DfDp_out;
  if (acceptModelParams_) {
    Derivative<double> DfDp = outArgs.get_DfDp(0);
    DfDp_out = DfDp.getMultiVector();
  }

  double t = inArgs.get_t();

  if (!isImplicit_) { // isImplicit_ == false
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out );
      f_out_view[0] =
       ( a*pow(t,3)*(8*c*d*t + 8*d*pow(b,2)
          + b*pow(t,0.5)*((7 + 9*c)*d + (-1 + c)*pow(t,4)))
        *pow(b + pow(t,0.5),-2)*pow(d + pow(t,4),-2))/2.0;
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = 0.0;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<double> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
    }
  } else { // isImplicit_ == true
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out );
      Thyra::ConstDetachedVectorView<double> x_dot_in_view( *x_dot_in );
      f_out_view[0] = x_dot_in_view[0] -
       ( a*pow(t,3)*(8*c*d*t + 8*d*pow(b,2)
          + b*pow(t,0.5)*((7 + 9*c)*d + (-1 + c)*pow(t,4)))
        *pow(b + pow(t,0.5),-2)*pow(d + pow(t,4),-2))/2.0;
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = alpha;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<double> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
    }
  }
}

RCP<const Thyra::VectorSpaceBase<double> > LogTimeModel::get_p_space(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_;
}

RCP<const Teuchos::Array<std::string> > LogTimeModel::get_p_names(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  RCP<Teuchos::Array<std::string> > p_strings =
    Teuchos::rcp(new Teuchos::Array<std::string>());
  p_strings->push_back("Model Coefficient:  a");
  p_strings->push_back("Model Coefficient:  b");
  p_strings->push_back("Model Coefficient:  c");
  p_strings->push_back("Model Coefficient:  d");
  return p_strings;
}

RCP<const Thyra::VectorSpaceBase<double> > LogTimeModel::get_g_space(int j) const
{
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_;
}

// private

void LogTimeModel::setupInOutArgs_() const
{
  if (isInitialized_) {
    return;
  }

  {
    // Set up prototypical InArgs
    ModelEvaluatorBase::InArgsSetup<double> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x );
    inArgs.setSupports( ModelEvaluatorBase::IN_ARG_beta );
    if (isImplicit_) {
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x_dot );
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_alpha );
    }
    if (acceptModelParams_) {
      inArgs.set_Np(Np_);
    }
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    ModelEvaluatorBase::OutArgsSetup<double> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_f );
    //if (isImplicit_) { // Thyra_ModelEvaluatorBase requires this
      outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_W_op );
    //}
    if (acceptModelParams_) {
      outArgs.set_Np_Ng(Np_,Ng_);
      outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_DfDp,0,DERIV_MV_BY_COL );
    }
    outArgs_ = outArgs;
  }

  // Set up nominal values
  nominalValues_ = inArgs_;
  if (haveIC_)
  {
    nominalValues_.set_t(t0_ic_);
    const RCP<VectorBase<double> > x_ic = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
      x_ic_view[0] = a_*(b_*pow(t0_ic_,4) + c_*pow(t0_ic_,4.5))
                      *pow(b_ + pow(t0_ic_,0.5),-1)*pow(d_ + pow(t0_ic_,4),-1);
    }
    nominalValues_.set_x(x_ic);
    if (acceptModelParams_) {
      const RCP<VectorBase<double> > p_ic = createMember(p_space_);
      {
        Thyra::DetachedVectorView<double> p_ic_view( *p_ic );
        p_ic_view[0] = a_;
        p_ic_view[1] = b_;
        p_ic_view[2] = c_;
        p_ic_view[3] = d_;
      }
      nominalValues_.set_p(0,p_ic);
    }
    if (isImplicit_) {
      const RCP<VectorBase<double> > x_dot_ic = createMember(x_space_);
      { // scope to delete DetachedVectorView
        Thyra::DetachedVectorView<double> x_dot_ic_view( *x_dot_ic );
        x_dot_ic_view[0] =
        ( a_*pow(t0_ic_,3)*(8*c_*d_*t0_ic_ + 8*d_*pow(b_,2)
          + b_*pow(t0_ic_,0.5)*((7 + 9*c_)*d_ + (-1 + c_)*pow(t0_ic_,4)))
          *pow(b_ + pow(t0_ic_,0.5),-2)*pow(d_ + pow(t0_ic_,4),-2))/2.0;

      }
      nominalValues_.set_x_dot(x_dot_ic);
    }
  }

  isInitialized_ = true;

}

void LogTimeModel::setParameterList(RCP<ParameterList> const& paramList)
{
  using Teuchos::get;
  TEUCHOS_TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  // 06/16/09 tscoffe:  TODO:  Only set the parameters that explicitely show up
  // in the new parameter list I.e.  Save all the previous options that have
  // been set in the past rather than resetting everything to defaults, even if
  // its not mentioned in the new parameter list.
  this->setMyParamList(paramList);
  RCP<ParameterList> pl = this->getMyNonconstParamList();
  bool isImplicit = get<bool>(*pl,Implicit_name);
  bool acceptModelParams = get<bool>(*pl,AcceptModelParams_name);
  bool haveIC = get<bool>(*pl,HaveIC_name);
  if ( (isImplicit != isImplicit_) ||
       (acceptModelParams != acceptModelParams_) ||
       (haveIC != haveIC_)
     ) {
    isInitialized_ = false;
  }
  isImplicit_ = isImplicit;
  acceptModelParams_ = acceptModelParams;
  haveIC_ = haveIC;
  a_ = get<double>(*pl,Coeff_a_name);
  b_ = get<double>(*pl,Coeff_b_name);
  c_ = get<double>(*pl,Coeff_c_name);
  d_ = get<double>(*pl,Coeff_d_name);
  x_ic_ = get<double>(*pl,IC_x_name);
  t0_ic_ = get<double>(*pl,IC_t0_name);
  setupInOutArgs_();
}

RCP<const ParameterList> LogTimeModel::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(Implicit_name,Implicit_default);
    pl->set(AcceptModelParams_name, AcceptModelParams_default);
    pl->set(HaveIC_name, HaveIC_default);
    Teuchos::setDoubleParameter(
        Coeff_a_name, Coeff_a_default, "Coefficient a in model",
        &*pl
        );
    Teuchos::setDoubleParameter(
        Coeff_b_name, Coeff_b_default, "Coefficient b in model",
        &*pl
        );
    Teuchos::setDoubleParameter(
        Coeff_c_name, Coeff_c_default, "Coefficient c in model",
        &*pl
        );
    Teuchos::setDoubleParameter(
        Coeff_d_name, Coeff_d_default, "Coefficient d in model",
        &*pl
        );
    Teuchos::setDoubleParameter(
        IC_x_name, IC_x_default, "Initial Condition for x",
        &*pl
        );
    Teuchos::setDoubleParameter(
        IC_t0_name, IC_t0_default, "Initial time t0",
        &*pl
        );
    validPL = pl;
  }
  return validPL;
}

} // namespace Rythmos

