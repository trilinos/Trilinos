#include "SinCosModel.hpp"

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"

#include <iostream>

namespace {
  const std::string Implicit_name = "Implicit model formulation";
  const bool Implicit_default = false;

  const std::string AcceptModelParams_name = "Accept model parameters";
  const bool AcceptModelParams_default = false;

  const std::string HaveIC_name = "Provide nominal values";
  const bool HaveIC_default = true;

  const std::string Coeff_a_name = "Coeff a";
  const double Coeff_a_default = 0.0;

  const std::string Coeff_f_name = "Coeff f";
  const double Coeff_f_default = 1.0;

  const std::string Coeff_L_name = "Coeff L";
  const double Coeff_L_default = 1.0;

  const std::string IC_x0_name = "IC x_0";
  const double IC_x0_default = 0.0;

  const std::string IC_x1_name = "IC x_1";
  const double IC_x1_default = 1.0;

  const std::string IC_t0_name = "IC t_0";
  const double IC_t0_default = 0.0;

} // namespace

namespace Tempus_Test {


// non-member Constructor
RCP<SinCosModel> sinCosModel(RCP<ParameterList> pList_)
{
  RCP<SinCosModel> model = rcp(new SinCosModel(pList_));
  return(model);
}


SinCosModel::SinCosModel(RCP<ParameterList> pList_)
{
  isInitialized = false;
  dim = 2;
  Np = 1; // Number of parameter vectors (1)
  np = 3; // Number of parameters in this vector (3)
  Ng = 1; // Number of observation functions (1)
  ng = 1; // Number of elements in this observation function (1)
  isImplicit = Implicit_default;
  acceptModelParams = AcceptModelParams_default;
  haveIC = HaveIC_default;
  a_ = Coeff_a_default;
  f_ = Coeff_f_default;
  L_ = Coeff_L_default;
  x0_ic = IC_x0_default;
  x1_ic = IC_x1_default;
  t0_ic = IC_t0_default;
  setParameterList(pList_);

  // Create x_space and f_space
  x_space = Thyra::defaultSpmdVectorSpace<double>(dim);
  f_space = Thyra::defaultSpmdVectorSpace<double>(dim);
  // Create p_space and g_space
  p_space = Thyra::defaultSpmdVectorSpace<double>(np);
  g_space = Thyra::defaultSpmdVectorSpace<double>(ng);
}

void SinCosModel::setImplicitFlag(bool implicit)
{
  if (isImplicit != implicit) {
    isInitialized = false;
  }
  isImplicit = implicit;
  setupInOutArgs();
}

ModelEvaluatorBase::InArgs<double> SinCosModel::getExactSolution(double t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs;
  double exact_t = t;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_x = createMember(x_space);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_x_view(*exact_x);
    exact_x_view[0] = a_+b_*sin((f_/L_)*t+phi_);
    exact_x_view[1] = b_*(f_/L_)*cos((f_/L_)*t+phi_);
  }
  inArgs.set_x(exact_x);
  if (isImplicit) {
    RCP<VectorBase<double> > exact_x_dot = createMember(x_space);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_x_dot_view(*exact_x_dot);
      exact_x_dot_view[0] = b_*(f_/L_)*cos((f_/L_)*t+phi_);
      exact_x_dot_view[1] = -b_*(f_/L_)*(f_/L_)*sin((f_/L_)*t+phi_);
    }
    inArgs.set_x_dot(exact_x_dot);
  }
  return(inArgs);
}

//
// 06/24/09 tscoffe:
// These are the exact sensitivities for the problem assuming the initial conditions are specified as:
//    s(0) = [1;0]    s(1) = [0;b/L]                 s(2) = [0;-b*f/(L*L)]
// sdot(0) = [0;0] sdot(1) = [0;-3*f*f*b/(L*L*L)] sdot(2) = [0;3*f*f*f*b/(L*L*L*L)]
//
ModelEvaluatorBase::InArgs<double> SinCosModel::getExactSensSolution(int j, double t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs;
  if (!acceptModelParams) {
    return inArgs;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, np );
  double exact_t = t;
  double b = b_;
  double f = f_;
  double L = L_;
  double phi = phi_;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_s = createMember(x_space);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_s_view(*exact_s);
    if (j == 0) { // dx/da
      exact_s_view[0] = 1.0;
      exact_s_view[1] = 0.0;
    } else if (j == 1) { // dx/df
      exact_s_view[0] = (b/L)*t*cos((f/L)*t+phi);
      exact_s_view[1] = (b/L)*cos((f/L)*t+phi)-(b*f*t/(L*L))*sin((f/L)*t+phi);
    } else { // dx/dL
      exact_s_view[0] = -(b*f*t/(L*L))*cos((f/L)*t+phi);
      exact_s_view[1] = -(b*f/(L*L))*cos((f/L)*t+phi)+(b*f*f*t/(L*L*L))*sin((f/L)*t+phi);
    }
  }
  inArgs.set_x(exact_s);
  if (isImplicit) {
    RCP<VectorBase<double> > exact_s_dot = createMember(x_space);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_s_dot_view(*exact_s_dot);
      if (j == 0) { // dxdot/da
        exact_s_dot_view[0] = 0.0;
        exact_s_dot_view[1] = 0.0;
      } else if (j == 1) { // dxdot/df
        exact_s_dot_view[0] = (b/L)*cos((f/L)*t+phi)-(b*f*t/(L*L))*sin((f/L)*t+phi);
        exact_s_dot_view[1] = -(2.0*b*f/(L*L))*sin((f/L)*t+phi)-(b*f*f*t/(L*L*L))*cos((f/L)*t+phi);
      } else { // dxdot/dL
        exact_s_dot_view[0] = -(b*f/(L*L))*cos((f/L)*t+phi)+(b*f*f*t/(L*L*L))*sin((f/L)*t+phi);
        exact_s_dot_view[1] = (2.0*b*f*f/(L*L*L))*sin((f/L)*t+phi)+(b*f*f*f*t/(L*L*L*L))*cos((f/L)*t+phi);
      }
    }
    inArgs.set_x_dot(exact_s_dot);
  }
  return(inArgs);
}

RCP<const Thyra::VectorSpaceBase<double> >
SinCosModel::get_x_space() const
{
  return x_space;
}


RCP<const Thyra::VectorSpaceBase<double> >
SinCosModel::get_f_space() const
{
  return f_space;
}


ModelEvaluatorBase::InArgs<double>
SinCosModel::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  return nominalValues;
}


RCP<Thyra::LinearOpWithSolveBase<double> >
SinCosModel::create_W() const
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
      RCP<Thyra::VectorBase<double> > vec = Thyra::createMember(x_space);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        vec_view[0] = 0.0;
        vec_view[1] = 1.0;
      }
      V_V(multivec->col(0).ptr(),*vec);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        vec_view[0] = 1.0;
        vec_view[1] = 0.0;
      }
      V_V(multivec->col(1).ptr(),*vec);
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
//SinCosModel::create_W() const
//{
//  return Thyra::multiVectorLinearOpWithSolve<double>();
//}


RCP<Thyra::LinearOpBase<double> >
SinCosModel::create_W_op() const
{
  RCP<Thyra::MultiVectorBase<double> > matrix = Thyra::createMembers(x_space, dim);
  return(matrix);
}


RCP<const Thyra::LinearOpWithSolveFactoryBase<double> >
SinCosModel::get_W_factory() const
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
SinCosModel::createInArgs() const
{
  setupInOutArgs();
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


ModelEvaluatorBase::OutArgs<double>
SinCosModel::createOutArgsImpl() const
{
  setupInOutArgs();
  return outArgs;
}


void SinCosModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );

  const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<double> x_in_view( *x_in );

  //double t = inArgs.get_t();
  double a = a_;
  double f = f_;
  double L = L_;
  if (acceptModelParams) {
    const RCP<const VectorBase<double> > p_in = inArgs.get_p(0).assert_not_null();
    Thyra::ConstDetachedVectorView<double> p_in_view( *p_in );
    a = p_in_view[0];
    f = p_in_view[1];
    L = p_in_view[2];
  }

  RCP<const VectorBase<double> > x_dot_in;
  double beta = inArgs.get_beta();
  double alpha = -1.0;
  if (isImplicit) {
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    alpha = inArgs.get_alpha();
  }

  const RCP<VectorBase<double> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<double> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<double> > DfDp_out;
  if (acceptModelParams) {
    Derivative<double> DfDp = outArgs.get_DfDp(0);
    DfDp_out = DfDp.getMultiVector();
  }

  if (!isImplicit) { // isImplicit == false
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out );
      f_out_view[0] = x_in_view[1];
      f_out_view[1] = (f/L)*(f/L)*(a-x_in_view[0]);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = 0.0;
      matrix_view(0,1) = +beta;
      matrix_view(1,0) = -beta*(f/L)*(f/L);
      matrix_view(1,1) = 0.0;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<double> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
      DfDp_out_view(0,1) = 0.0;
      DfDp_out_view(0,2) = 0.0;
      DfDp_out_view(1,0) = (f/L)*(f/L);
      DfDp_out_view(1,1) = (2.0*f/(L*L))*(a-x_in_view[0]);
      DfDp_out_view(1,2) = -(2.0*f*f/(L*L*L))*(a-x_in_view[0]);
    }
  } else { // isImplicit == true
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out );
      Thyra::ConstDetachedVectorView<double> x_dot_in_view( *x_dot_in );
      f_out_view[0] = x_dot_in_view[0] - x_in_view[1];
      f_out_view[1] = x_dot_in_view[1] - (f/L)*(f/L)*(a-x_in_view[0]);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = alpha;
      matrix_view(0,1) = -beta;
      matrix_view(1,0) = +beta*(f/L)*(f/L);
      matrix_view(1,1) = alpha;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<double> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
      DfDp_out_view(0,1) = 0.0;
      DfDp_out_view(0,2) = 0.0;
      DfDp_out_view(1,0) = -(f/L)*(f/L);
      DfDp_out_view(1,1) = -(2.0*f/(L*L))*(a-x_in_view[0]);
      DfDp_out_view(1,2) = +(2.0*f*f/(L*L*L))*(a-x_in_view[0]);
    }
  }
}

RCP<const Thyra::VectorSpaceBase<double> > SinCosModel::get_p_space(int l) const
{
  if (!acceptModelParams) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np );
  return p_space;
}

RCP<const Teuchos::Array<std::string> > SinCosModel::get_p_names(int l) const
{
  if (!acceptModelParams) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np );
  RCP<Teuchos::Array<std::string> > p_strings =
    Teuchos::rcp(new Teuchos::Array<std::string>());
  p_strings->push_back("Model Coefficient:  a");
  p_strings->push_back("Model Coefficient:  f");
  p_strings->push_back("Model Coefficient:  L");
  return p_strings;
}

RCP<const Thyra::VectorSpaceBase<double> > SinCosModel::get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng );
  return g_space;
}

// private

void SinCosModel::setupInOutArgs() const
{
  if (isInitialized) {
    return;
  }

  {
    // Set up prototypical InArgs
    ModelEvaluatorBase::InArgsSetup<double> inArgs_;
    inArgs_.setModelEvalDescription(this->description());
    inArgs_.setSupports( ModelEvaluatorBase::IN_ARG_t );
    inArgs_.setSupports( ModelEvaluatorBase::IN_ARG_x );
    inArgs_.setSupports( ModelEvaluatorBase::IN_ARG_beta );
    if (isImplicit) {
      inArgs_.setSupports( ModelEvaluatorBase::IN_ARG_x_dot );
      inArgs_.setSupports( ModelEvaluatorBase::IN_ARG_alpha );
    }
    if (acceptModelParams) {
      inArgs_.set_Np(Np);
    }
    inArgs = inArgs_;
  }

  {
    // Set up prototypical OutArgs
    ModelEvaluatorBase::OutArgsSetup<double> outArgs_;
    outArgs_.setModelEvalDescription(this->description());
    outArgs_.setSupports( ModelEvaluatorBase::OUT_ARG_f );
    //if (isImplicit) { // Thyra_ModelEvaluatorBase requires this
      outArgs_.setSupports( ModelEvaluatorBase::OUT_ARG_W_op );
    //}
    if (acceptModelParams) {
      outArgs_.set_Np_Ng(Np,Ng);
      outArgs_.setSupports( ModelEvaluatorBase::OUT_ARG_DfDp,0,DERIV_MV_BY_COL );
    }
    outArgs = outArgs_;
  }

  // Set up nominal values
  nominalValues = inArgs;
  if (haveIC)
  {
    nominalValues.set_t(t0_ic);
    const RCP<VectorBase<double> > x_ic = createMember(x_space);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
      x_ic_view[0] = a_+b_*sin((f_/L_)*t0_ic_+phi_);
      x_ic_view[1] = b_*(f_/L_)*cos((f_/L_)*t0_ic_+phi_);
    }
    nominalValues.set_x(x_ic);
    if (acceptModelParams) {
      const RCP<VectorBase<double> > p_ic = createMember(p_space);
      {
        Thyra::DetachedVectorView<double> p_ic_view( *p_ic );
        p_ic_view[0] = a_;
        p_ic_view[1] = f_;
        p_ic_view[2] = L_;
      }
      nominalValues_.set_p(0,p_ic);
    }
    if (isImplicit) {
      const RCP<VectorBase<double> > x_dot_ic = createMember(x_space);
      { // scope to delete DetachedVectorView
        Thyra::DetachedVectorView<double> x_dot_ic_view( *x_dot_ic );
        x_dot_ic_view[0] = b_*(f_/L_)*cos((f_/L_)*t0_ic_+phi_);
        x_dot_ic_view[1] = -b_*(f_/L_)*(f_/L_)*sin((f_/L_)*t0_ic_+phi_);
      }
      nominalValues_.set_x_dot(x_dot_ic);
    }
  }

  isInitialized = true;

}

void SinCosModel::setParameterList(RCP<ParameterList> const& paramList)
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
  bool isImplicit_ = get<bool>(*pl,Implicit_name);
  bool acceptModelParams = get<bool>(*pl,AcceptModelParams_name);
  bool haveIC_ = get<bool>(*pl,HaveIC_name);
  if ( (isImplicit != isImplicit_) ||
       (acceptModelParams != acceptModelParams) ||
       (haveIC_ != haveIC)
     ) {
    isInitialized_ = false;
  }
  isImplicit = isImplicit_;
  acceptModelParams_ = acceptModelParams;
  haveIC = haveIC_;
  a_ = get<double>(*pl,Coeff_a_name);
  f_ = get<double>(*pl,Coeff_f_name);
  L_ = get<double>(*pl,Coeff_L_name);
  x0_ic_ = get<double>(*pl,IC_x0_name);
  x1_ic_ = get<double>(*pl,IC_x1_name);
  t0_ic_ = get<double>(*pl,IC_t0_name);
  calculateCoeffFromIC();
  setupInOutArgs();
}

RCP<const ParameterList> SinCosModel::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(Implicit_name,Implicit_default);
    pl->set(AcceptModelParams_name, AcceptModelParams_default);
    pl->set(HaveIC_name, HaveIC_default);
    Teuchos::setDoubleParameter(
        Coeff_a_name, Coeff_a_default, "Coefficient a in model", &*pl);
    Teuchos::setDoubleParameter(
        Coeff_f_name, Coeff_f_default, "Coefficient f in model", &*pl);
    Teuchos::setDoubleParameter(
        Coeff_L_name, Coeff_L_default, "Coefficient L in model", &*pl);
    Teuchos::setDoubleParameter(
        IC_x0_name, IC_x0_default, "Initial Condition for x0", &*pl);
    Teuchos::setDoubleParameter(
        IC_x1_name, IC_x1_default, "Initial Condition for x1", &*pl);
    Teuchos::setDoubleParameter(
        IC_t0_name, IC_t0_default, "Initial time t0", &*pl);
    validPL = pl;
  }
  return validPL;
}

void SinCosModel::calculateCoeffFromIC()
{
  phi_ = atan(((f_/L_)/x1_ic_)*(x0_ic_-a_))-(f_/L_)*t0_ic_;
  b_ = x1_ic_/((f_/L_)*cos((f_/L_)*t0_ic_+phi_));
}

} // namespace Tempus_Test
