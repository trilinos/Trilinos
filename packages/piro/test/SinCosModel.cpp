// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SinCosModel.hpp"

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

using Teuchos::RCP;
using Teuchos::rcp;

SinCosModel::
SinCosModel(Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_ = false;
  dim_ = 2;
  Np_ = 3; // Number of parameter vectors (p, dx/dp, dx_dot/dp)
  np_ = 3; // Number of parameters in this vector (3)
  Ng_ = 1; // Number of observation functions (1)
  ng_ = dim_; // Number of elements in this observation function ( == x )
  acceptModelParams_ = false;
  useDfDpAsTangent_ = false;
  haveIC_ = true;
  a_ = 0.0;
  f_ = 1.0;
  L_ = 1.0;
  x0_ic_ = 0.0;
  x1_ic_ = 1.0;
  t0_ic_ = 0.0;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  // Create p_space and g_space
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(ng_);

  setParameterList(pList_);

  // Create DxDp product space
  DxDp_space_ = Thyra::multiVectorProductVectorSpace(x_space_, np_);
}

//
// 06/24/09 tscoffe:
// These are the exact sensitivities for the problem assuming the initial conditions are specified as:
//    s(0) = [1;0]    s(1) = [0;b/L]                 s(2) = [0;-b*f/(L*L)]
// sdot(0) = [0;0] sdot(1) = [0;-3*f*f*b/(L*L*L)] sdot(2) = [0;3*f*f*f*b/(L*L*L*L)]
//


Thyra::ModelEvaluatorBase::InArgs<Scalar> 
SinCosModel::
getExactSolution(Scalar t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setupInOutArgs_ must be called first!\n");
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  Scalar exact_t = t;
  inArgs.set_t(exact_t);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    exact_x_view[0] = a_+b_*sin((f_/L_)*t+phi_);
    exact_x_view[1] = b_*(f_/L_)*cos((f_/L_)*t+phi_);
  }
  inArgs.set_x(exact_x);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x_dot = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
    exact_x_dot_view[0] = b_*(f_/L_)*cos((f_/L_)*t+phi_);
    exact_x_dot_view[1] = -b_*(f_/L_)*(f_/L_)*sin((f_/L_)*t+phi_);
  }
  inArgs.set_x_dot(exact_x_dot);
  return(inArgs);
}

//
// 06/24/09 tscoffe:
// These are the exact sensitivities for the problem assuming the initial conditions are specified as:
//    s(0) = [1;0]    s(1) = [0;b/L]                 s(2) = [0;-b*f/(L*L)]
// sdot(0) = [0;0] sdot(1) = [0;-3*f*f*b/(L*L*L)] sdot(2) = [0;3*f*f*f*b/(L*L*L*L)]
//
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SinCosModel::
getExactSensSolution(int j, Scalar t) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setupInOutArgs_ must be called first!\n");
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  if (!acceptModelParams_) {
    return inArgs;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, np_ );
  Scalar exact_t = t;
  Scalar b = b_;
  Scalar f = f_;
  Scalar L = L_;
  Scalar phi = phi_;
  inArgs.set_t(exact_t);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_s = createMember(x_space_);
   { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_s_view(*exact_s);
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
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_s_dot = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_s_dot_view(*exact_s_dot);
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
  return(inArgs);
}



Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SinCosModel::
get_x_space() const
{
  return x_space_;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SinCosModel::
get_f_space() const
{
  return f_space_;
}

Thyra::ModelEvaluatorBase::InArgs<Scalar>
SinCosModel::
getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
SinCosModel::
create_W() const
{
   using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory = this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<Scalar> > multivec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<Scalar> vec_view( *vec );
        vec_view[0] = 0.0;
        vec_view[1] = 1.0;
      }
      V_V(multivec->col(0).ptr(),*vec);
      {
        Thyra::DetachedVectorView<Scalar> vec_view( *vec );
        vec_view[0] = 1.0;
        vec_view[1] = 0.0;
      }
      V_V(multivec->col(1).ptr(),*vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
    Thyra::linearOpWithSolve<Scalar>(*W_factory, matrix);
  return W;
}

Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
SinCosModel::
create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix = Thyra::createMembers(x_space_, dim_);
  return(matrix);
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
SinCosModel::
get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}


Thyra::ModelEvaluatorBase::InArgs<Scalar>
SinCosModel::
createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

Thyra::ModelEvaluatorBase::OutArgs<Scalar>
SinCosModel::
createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

void
SinCosModel::
evalModelImpl( const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
               const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Teuchos::rcp_dynamic_cast;
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setupInOutArgs_ must be called first!\n");

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view( *x_in );

  //Scalar t = inArgs.get_t();
  Scalar a = a_;
  Scalar f = f_;
  Scalar L = L_;
  if (acceptModelParams_) {
    const RCP<const VectorBase<Scalar> > p_in =
      inArgs.get_p(0).assert_not_null();
    Thyra::ConstDetachedVectorView<Scalar> p_in_view( *p_in );
    a = p_in_view[0];
    f = p_in_view[1];
    L = p_in_view[2];
  }
  RCP<const MultiVectorBase<Scalar> > DxDp_in, DxdotDp_in;
  if (acceptModelParams_) {
    if (inArgs.get_p(1) != Teuchos::null)
      DxDp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(1))->getMultiVector();
    if (inArgs.get_p(2) != Teuchos::null)
      DxdotDp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(2))->getMultiVector();
  }

  Scalar beta = inArgs.get_beta();

  const RCP<VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<Scalar> > DfDp_out;
  if (acceptModelParams_) {
    Thyra::ModelEvaluatorBase::Derivative<Scalar> DfDp = outArgs.get_DfDp(0);
    DfDp_out = DfDp.getMultiVector();
  }
  if (inArgs.get_x_dot().is_null()) {

    // Evaluate the Explicit ODE f(x,t) [= 0]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
      f_out_view[0] = x_in_view[1];
      f_out_view[1] = (f/L)*(f/L)*(a-x_in_view[0]);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view( *matrix );
      matrix_view(0,0) = 0.0;               // d(f0)/d(x0_n)
      matrix_view(0,1) = +beta;             // d(f0)/d(x1_n)
      matrix_view(1,0) = -beta*(f/L)*(f/L); // d(f1)/d(x0_n)
      matrix_view(1,1) = 0.0;               // d(f1)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
      DfDp_out_view(0,1) = 0.0;
      DfDp_out_view(0,2) = 0.0;
      DfDp_out_view(1,0) = (f/L)*(f/L);
      DfDp_out_view(1,1) = (2.0*f/(L*L))*(a-x_in_view[0]);
      DfDp_out_view(1,2) = -(2.0*f*f/(L*L*L))*(a-x_in_view[0]);
 // Compute df/dp + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(DxDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxDp( *DxDp_in );
        DfDp_out_view(0,0) +=  DxDp(1,0);
        DfDp_out_view(0,1) +=  DxDp(1,1);
        DfDp_out_view(0,2) +=  DxDp(1,2);
        DfDp_out_view(1,0) += -(f/L)*(f/L) * DxDp(0,0);
        DfDp_out_view(1,1) += -(f/L)*(f/L) * DxDp(0,1);
        DfDp_out_view(1,2) += -(f/L)*(f/L) * DxDp(0,2);
      }
    }
  } else {

    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const VectorBase<Scalar> > x_dot_in;
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view( *x_dot_in );
      f_out_view[0] = x_dot_in_view[0] - x_in_view[1];
      f_out_view[1] = x_dot_in_view[1] - (f/L)*(f/L)*(a-x_in_view[0]);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view( *matrix );
      matrix_view(0,0) = alpha;             // d(f0)/d(x0_n)
      matrix_view(0,1) = -beta;             // d(f0)/d(x1_n)
      matrix_view(1,0) = +beta*(f/L)*(f/L); // d(f1)/d(x0_n)
      matrix_view(1,1) = alpha;             // d(f1)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view( *DfDp_out );
      DfDp_out_view(0,0) = 0.0;
      DfDp_out_view(0,1) = 0.0;
      DfDp_out_view(0,2) = 0.0;
      DfDp_out_view(1,0) = -(f/L)*(f/L);
      DfDp_out_view(1,1) = -(2.0*f/(L*L))*(a-x_in_view[0]);
      DfDp_out_view(1,2) = +(2.0*f*f/(L*L*L))*(a-x_in_view[0]);

      // Compute df/dp + (df/dx_dot) * (dx_dot/dp) + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(DxdotDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxdotDp( *DxdotDp_in );
        DfDp_out_view(0,0) += DxdotDp(0,0);
        DfDp_out_view(0,1) += DxdotDp(0,1);
        DfDp_out_view(0,2) += DxdotDp(0,2);
        DfDp_out_view(1,0) += DxdotDp(1,0);
        DfDp_out_view(1,1) += DxdotDp(1,1);
        DfDp_out_view(1,2) += DxdotDp(1,2);
      }
      if (useDfDpAsTangent_ && !is_null(DxDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxDp( *DxDp_in );
        DfDp_out_view(0,0) += -DxDp(1,0);
        DfDp_out_view(0,1) += -DxDp(1,1);
        DfDp_out_view(0,2) += -DxDp(1,2);
        DfDp_out_view(1,0) +=  (f/L)*(f/L) * DxDp(0,0);
        DfDp_out_view(1,1) +=  (f/L)*(f/L) * DxDp(0,1);
        DfDp_out_view(1,2) +=  (f/L)*(f/L) * DxDp(0,2);
      }
    }
  }

  // Responses:  g = x
  if (acceptModelParams_) {
    RCP<VectorBase<Scalar> > g_out = outArgs.get_g(0);
    if (g_out != Teuchos::null)
      Thyra::assign(g_out.ptr(), *x_in);

    RCP<Thyra::MultiVectorBase<Scalar> > DgDp_out =
      outArgs.get_DgDp(0,0).getMultiVector();
    if (DgDp_out != Teuchos::null)
      Thyra::assign(DgDp_out.ptr(), Scalar(0.0));

    RCP<Thyra::MultiVectorBase<Scalar> > DgDx_out =
      outArgs.get_DgDx(0).getMultiVector();
    if (DgDx_out != Teuchos::null) {
      Thyra::DetachedMultiVectorView<Scalar> DgDx_out_view( *DgDx_out );
      DgDx_out_view(0,0) = 1.0;
      DgDx_out_view(0,1) = 0.0;
      DgDx_out_view(1,0) = 0.0;
      DgDx_out_view(1,1) = 1.0;
    }
  }
}


Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SinCosModel::
get_p_space(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
  if (l == 0)
    return p_space_;
  else if (l == 1 || l == 2)
    return DxDp_space_;
  return Teuchos::null;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
SinCosModel::
get_p_names(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
  Teuchos::RCP<Teuchos::Array<std::string> > p_strings =
    Teuchos::rcp(new Teuchos::Array<std::string>());
  if (l == 0) {
    p_strings->push_back("Model Coefficient:  a");
    p_strings->push_back("Model Coefficient:  f");
    p_strings->push_back("Model Coefficient:  L");
  }
  else if (l == 1)
    p_strings->push_back("DxDp");
  else if (l == 2)
    p_strings->push_back("Dx_dotDp");
  return p_strings;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SinCosModel::
get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
  return g_space_;
}

void
SinCosModel::
setupInOutArgs_() const
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
    inArgs.setSupports( MEB::IN_ARG_t );
    inArgs.setSupports( MEB::IN_ARG_x );
    inArgs.setSupports( MEB::IN_ARG_beta );
    inArgs.setSupports( MEB::IN_ARG_x_dot );
    inArgs.setSupports( MEB::IN_ARG_alpha );
    if (acceptModelParams_) {
      inArgs.set_Np(Np_);
    }
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( MEB::OUT_ARG_f );
    outArgs.setSupports( MEB::OUT_ARG_W_op );
    if (acceptModelParams_) {
      outArgs.set_Np_Ng(Np_,Ng_);
      outArgs.setSupports( MEB::OUT_ARG_DfDp,0,
                           MEB::DERIV_MV_JACOBIAN_FORM );
      outArgs.setSupports( MEB::OUT_ARG_DgDp,0,0,
                           MEB::DERIV_MV_JACOBIAN_FORM );
      outArgs.setSupports( MEB::OUT_ARG_DgDx,0,
                           MEB::DERIV_MV_GRADIENT_FORM );
    }
    outArgs_ = outArgs;
  }
   // Set up nominal values
  nominalValues_ = inArgs_;
  if (haveIC_)
  {
    nominalValues_.set_t(t0_ic_);
    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
      x_ic_view[0] = a_+b_*sin((f_/L_)*t0_ic_+phi_);
      x_ic_view[1] = b_*(f_/L_)*cos((f_/L_)*t0_ic_+phi_);
    }
    nominalValues_.set_x(x_ic);
    if (acceptModelParams_) {
      const RCP<Thyra::VectorBase<Scalar> > p_ic = createMember(p_space_);
      {
        Thyra::DetachedVectorView<Scalar> p_ic_view( *p_ic );
        p_ic_view[0] = a_;
        p_ic_view[1] = f_;
        p_ic_view[2] = L_;
      }
      nominalValues_.set_p(0,p_ic);
    }
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view( *x_dot_ic );
      x_dot_ic_view[0] = b_*(f_/L_)*cos((f_/L_)*t0_ic_+phi_);
      x_dot_ic_view[1] = -b_*(f_/L_)*(f_/L_)*sin((f_/L_)*t0_ic_+phi_);
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;

}

void
SinCosModel::
setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL = Teuchos::rcp(new ParameterList("SinCosModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();
  bool acceptModelParams = get<bool>(*pl,"Accept model parameters");
  bool haveIC = get<bool>(*pl,"Provide nominal values");
  bool useDfDpAsTangent = get<bool>(*pl, "Use DfDp as Tangent");
  if ( (acceptModelParams != acceptModelParams_) ||
       (haveIC != haveIC_)
     ) {
    isInitialized_ = false;
  }
  acceptModelParams_ = acceptModelParams;
  haveIC_ = haveIC;
  useDfDpAsTangent_ = useDfDpAsTangent;
  a_ = get<Scalar>(*pl,"Coeff a");
  f_ = get<Scalar>(*pl,"Coeff f");
  L_ = get<Scalar>(*pl,"Coeff L");
  x0_ic_ = get<Scalar>(*pl,"IC x0");
  x1_ic_ = get<Scalar>(*pl,"IC x1");
  t0_ic_ = get<Scalar>(*pl,"IC t0");
  calculateCoeffFromIC_();
  setupInOutArgs_();
}

Teuchos::RCP<const Teuchos::ParameterList>
SinCosModel::
getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters", false);
    pl->set("Provide nominal values", true);
    pl->set("Use DfDp as Tangent", false);
    pl->set<std::string>("Output File Name", "Tempus_BDF2_SinCos");
    Teuchos::setDoubleParameter(
        "Coeff a", 0.0, "Coefficient a in model", &*pl);
    Teuchos::setDoubleParameter(
        "Coeff f", 1.0, "Coefficient f in model", &*pl);
    Teuchos::setDoubleParameter(
        "Coeff L", 1.0, "Coefficient L in model", &*pl);
    Teuchos::setDoubleParameter(
        "IC x0", 0.0, "Initial Condition for x0", &*pl);
    Teuchos::setDoubleParameter(
        "IC x1", 1.0, "Initial Condition for x1", &*pl);
    Teuchos::setDoubleParameter(
        "IC t0", 0.0, "Initial time t0", &*pl);
    Teuchos::setIntParameter(
        "Number of Time Step Sizes", 1, "Number time step sizes for convergence study", &*pl);
    validPL = pl;
  }
  return validPL;
}

void
SinCosModel::
calculateCoeffFromIC_()
{
  phi_ = atan(((f_/L_)/x1_ic_)*(x0_ic_-a_))-(f_/L_)*t0_ic_;
  b_ = x1_ic_/((f_/L_)*cos((f_/L_)*t0_ic_+phi_));
}

Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
SinCosModelAdjoint::
create_W() const
{
  using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory = this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<Scalar> > multivec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(this->f_space_);
      {
        Thyra::DetachedVectorView<Scalar> vec_view( *vec );
        vec_view[0] = 0.0;
        vec_view[1] = 1.0;
      }
      V_V(multivec->col(0).ptr(),*vec);
      {
        Thyra::DetachedVectorView<Scalar> vec_view( *vec );
        vec_view[0] = 1.0;
        vec_view[1] = 0.0;
      }
      V_V(multivec->col(1).ptr(),*vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
    Thyra::linearOpWithSolve<Scalar>(
      *W_factory,
      matrix
      );
  return W;
}

Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
SinCosModelAdjoint::
create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix = Thyra::createMembers(this->f_space_, this->dim_);
  return(matrix);
}

Thyra::ModelEvaluatorBase::InArgs<Scalar>
SinCosModelAdjoint::
createInArgs() const
{
  // This ME should use the same InArgs as the base SinCosModel.  However
  // we can't just use it's InArgs directly because the description won't
  // match (which is checked in debug builds).  Instead create a new InArgsSetup
  // initialized by SinCosModel::createInArgs() and set the description
  // appropriately.
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs = SinCosModel::createInArgs();
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}

Thyra::ModelEvaluatorBase::OutArgs<Scalar>
SinCosModelAdjoint::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports( MEB::OUT_ARG_f ); // Apparently all models have to support f
  outArgs.setSupports( MEB::OUT_ARG_W_op );
  outArgs.set_Np_Ng(this->Np_,0);
  return outArgs;
}

void
SinCosModelAdjoint::
evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Teuchos::rcp_dynamic_cast;
  TEUCHOS_TEST_FOR_EXCEPTION( !this->isInitialized_, std::logic_error,
      "Error, setupInOutArgs_ must be called first!\n");

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view( *x_in );

  //double t = inArgs.get_t();
  //Scalar a = this->a_;
  Scalar f = this->f_;
  Scalar L = this->L_;
  if (this->acceptModelParams_) {
    const RCP<const VectorBase<Scalar> > p_in =
      inArgs.get_p(0).assert_not_null();
    Thyra::ConstDetachedVectorView<Scalar> p_in_view( *p_in );
    //a = p_in_view[0];
    f = p_in_view[1];
    L = p_in_view[2];
  }

  Scalar beta = inArgs.get_beta();

  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  if (inArgs.get_x_dot().is_null()) {

    // Evaluate the Explicit ODE f(x,t) [= 0]
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view( *matrix );
      matrix_view(0,0) = 0.0;               // d(f0)/d(x0_n)
      matrix_view(1,0) = +beta;             // d(f0)/d(x1_n)
      matrix_view(0,1) = -beta*(f/L)*(f/L); // d(f1)/d(x0_n)
      matrix_view(1,1) = 0.0;               // d(f1)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
  } else {

    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const VectorBase<Scalar> > x_dot_in;
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view( *matrix );
      matrix_view(0,0) = alpha;             // d(f0)/d(x0_n)
      matrix_view(1,0) = -beta;             // d(f0)/d(x1_n)
      matrix_view(0,1) = +beta*(f/L)*(f/L); // d(f1)/d(x0_n)
      matrix_view(1,1) = alpha;             // d(f1)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
  }
}

