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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "VanderPolModel.hpp"

#ifdef Rythmos_ENABLE_Sacado

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Sacado.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"

// Nonlinear ODE system with manufactured solution based on asymptotic solution 
// for small epsilon
//
// f[0] = x_dot[0] - x[1];
// f[1] = x_dot[1] - (eps*(1.0-x[0]*x[0])*x[1]-x[0]) - forcing_term;
//
// forcing term is defined so that exact solution is given by
//
// x_0(t) =  2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t))
// x_1(t) = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t))
//
// initial conditions for exact solution
//
// x_0(0) = 2
// x_1(0) = 0
//
// exact time derivatives
//
// x_dot_0(t) = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t))
// x_dot_1(t) = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t))
//


#ifdef VANDERPOLMODEL_DEBUG
#include <iostream>
#endif // VANDERPOLMODEL_DEBUG

namespace {
  const std::string Implicit_name = "Implicit model formulation";
  const bool Implicit_default = false;

  const std::string AcceptModelParams_name = "Accept model parameters";
  const bool AcceptModelParams_default = false;

  const std::string HaveIC_name = "Provide nominal values";
  const bool HaveIC_default = true;

  const std::string Coeff_epsilon_name = "Coeff epsilon";
  const double Coeff_epsilon_default = 0.5;

  const std::string IC_t0_name = "IC t_0";
  const double IC_t0_default = 0.0;

  const std::string IC_x0_name = "IC x_0";
  const double IC_x0_default = 2.0;

  const std::string IC_x1_name = "IC x_1";
  const double IC_x1_default = 0.0;


} // namespace


namespace Rythmos {


// Constructor
VanderPolModel::VanderPolModel()
{
  isInitialized_ = false;
  dim_ = 2;
  Np_ = 1; // Number of parameter vectors (1)
  np_ = 1; // Number of parameters in this vector (1)
  Ng_ = 0; // Number of observation functions (0)
  ng_ = 0; // Number of elements in this observation function (0)
  isImplicit_ = Implicit_default;
  acceptModelParams_ = AcceptModelParams_default;
  haveIC_ = HaveIC_default;
  epsilon_ = Coeff_epsilon_default;
  t0_ic_ = IC_t0_default;
  x0_ic_ = IC_x0_default;
  x1_ic_ = IC_x1_default;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  // Create p_space and g_space
  p_space_ = Thyra::defaultSpmdVectorSpace<double>(np_);
  g_space_ = Thyra::defaultSpmdVectorSpace<double>(ng_);
}

void VanderPolModel::setImplicitFlag(bool implicit) 
{
  if (isImplicit_ != implicit) {
    isInitialized_ = false;
  }
  isImplicit_ = implicit;
  setupInOutArgs_();
}


ModelEvaluatorBase::InArgs<double> VanderPolModel::getExactSolution(double t) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);

  double eps = epsilon_;
  RCP<VectorBase<double> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_x_view(*exact_x);

    exact_x_view[0] =  2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
    exact_x_view[1] = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  }
  inArgs.set_x(exact_x);
  if (isImplicit_) {
    RCP<VectorBase<double> > exact_x_dot = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_x_dot_view(*exact_x_dot);

      exact_x_dot_view[0] = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
      exact_x_dot_view[1] = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t));
    }
    inArgs.set_x_dot(exact_x_dot);
  }
  return(inArgs);
}


ModelEvaluatorBase::InArgs<double>
VanderPolModel::getExactSensSolution(int j, double t) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );

#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, np_ );
#endif

  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);

  RCP<VectorBase<double> > exact_s = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_s_view(*exact_s);

    exact_s_view[0] = (0.75*sin(t)-0.25*sin(3.0*t));
    exact_s_view[1] = (0.75*cos(t)-0.75*cos(3.0*t));
  }
  inArgs.set_x(exact_s);
  if (isImplicit_) {
    RCP<VectorBase<double> > exact_s_dot = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_s_dot_view(*exact_s_dot);

      exact_s_dot_view[0] = (0.75*cos(t)-0.75*cos(3.0*t));
      exact_s_dot_view[1] = (-0.75*sin(t)+3.0*0.75*sin(3.0*t));
    }
    inArgs.set_x_dot(exact_s_dot);
  }
  return(inArgs);
}


RCP<const Thyra::VectorSpaceBase<double> >
VanderPolModel::get_x_space() const
{
  return x_space_;
}


RCP<const Thyra::VectorSpaceBase<double> >
VanderPolModel::get_f_space() const
{
  return f_space_;
}


ModelEvaluatorBase::InArgs<double>
VanderPolModel::getNominalValues() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );
  return nominalValues_;
}


RCP<Thyra::LinearOpWithSolveBase<double> >
VanderPolModel::create_W() const
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
        vec_view[0] = 0.0;
        vec_view[1] = 1.0;
      }
      V_V(&*(multivec->col(0)),*vec);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        vec_view[0] = 1.0;
        vec_view[1] = 1.0;
      }
      V_V(&*(multivec->col(1)),*vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<double> > W = 
    Thyra::linearOpWithSolve<double>(
      *W_factory,
      matrix
      );
  return W;
}


RCP<Thyra::LinearOpBase<double> >
VanderPolModel::create_W_op() const
{
  RCP<Thyra::MultiVectorBase<double> > matrix = Thyra::createMembers(x_space_, dim_);
  return(matrix);
}


RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > 
VanderPolModel::get_W_factory() const
{
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  return W_factory;
}

ModelEvaluatorBase::InArgs<double>
VanderPolModel::createInArgs() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );
  return inArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class ScalarT>
Array<Sacado::Fad::DFad<ScalarT> > convertToPassiveFadArray(
  const ArrayView<const ScalarT> &av
  )
{
  using Sacado::Fad::DFad;
  Array<DFad<ScalarT> > a_fad(av.size());
  for (int i = 0; i < av.size(); ++i) {
    a_fad[i] = DFad<double>(av[i]);
  }
  return a_fad;
}


template<class ScalarT>
Array<Sacado::Fad::DFad<ScalarT> > convertToIndepVarFadArray(
  const ArrayView<const ScalarT> &av, const int num_derivs,
  const Ptr<int> &deriv_idx)
{
  using Sacado::Fad::DFad;
  Array<DFad<ScalarT> > a_fad(av.size());
  for (int i = 0; i < av.size(); ++i) {
    a_fad[i] = DFad<double>(num_derivs, (*deriv_idx)++, av[i]);
  }
  return a_fad;
}


ModelEvaluatorBase::OutArgs<double>
VanderPolModel::createOutArgsImpl() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );
  return outArgs_;
}


void VanderPolModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::optInArg;
  using Teuchos::inOutArg;
  using Sacado::Fad::DFad;

  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setParameterList must be called first!\n"
      );

  const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<double> x_in_view( *x_in ); 

  double t = inArgs.get_t();
  double eps = epsilon_;
  if (acceptModelParams_) {
    const RCP<const VectorBase<double> > p_in = inArgs.get_p(0).assert_not_null();
    Thyra::ConstDetachedVectorView<double> p_in_view( *p_in ); 
    eps = p_in_view[0];
  }

  RCP<const VectorBase<double> > x_dot_in;
  double alpha = -1.0;
  double beta = -1.0;
  if (isImplicit_) {
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    alpha = inArgs.get_alpha();
    beta = inArgs.get_beta();
  }

  const RCP<VectorBase<double> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<double> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<double> > DfDp_out; 
  if (acceptModelParams_) {
    Derivative<double> DfDp = outArgs.get_DfDp(0); 
    DfDp_out = DfDp.getMultiVector();
  }

  // Determine how many derivatives we will compute
  
  int num_derivs = 0;
  if (nonnull(W_out)) {
    num_derivs += 2;
    if (isImplicit_) {
      num_derivs += 2;
    }
  }
  if (nonnull(DfDp_out))
    num_derivs += 1;

  // Set up the FAD derivative objects

  int deriv_i = 0;

  Array<DFad<double> > x_dot_fad;
  int x_dot_idx_offset = 0;
  if (isImplicit_) {
    Thyra::ConstDetachedVectorView<double> x_dot_in_view( *x_dot_in );
    if (nonnull(W_out)) {
      x_dot_idx_offset = deriv_i;
      x_dot_fad = convertToIndepVarFadArray<double>(x_dot_in_view.sv().values()(),
        num_derivs, inOutArg(deriv_i));
    }
    else {
      x_dot_fad = convertToPassiveFadArray<double>(x_dot_in_view.sv().values()());
    }
  }

  Array<DFad<double> > x_fad;
  int x_idx_offset = 0;
  if (nonnull(W_out)) {
    x_idx_offset = deriv_i;
    x_fad = convertToIndepVarFadArray<double>(x_in_view.sv().values()(),
      num_derivs, inOutArg(deriv_i));
  }
  else {
    x_fad = convertToPassiveFadArray<double>(x_in_view.sv().values()());
  }

  DFad<double> eps_fad(eps); // Default passive
  int eps_idx_offset = 0;
  if (nonnull(DfDp_out)) {
    eps_idx_offset = deriv_i;
    eps_fad = DFad<double>(num_derivs, deriv_i++, eps);
  }
  
  // Compute the function

  Array<DFad<double> > f_fad(2);
  this->eval_f<DFad<double> >( x_dot_fad, x_fad, eps_fad, t, f_fad );

  // Extract the output

  if (nonnull(f_out)) {
    Thyra::DetachedVectorView<double> f_out_view( *f_out ); 
    for ( int i = 0; i < as<int>(f_fad.size()); ++i )
      f_out_view[i] = f_fad[i].val();
  }

  if (nonnull(W_out)) {
    const RCP<Thyra::MultiVectorBase<double> > matrix =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out, true);
    Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
    if (isImplicit_) {
      for ( int i = 0; i < matrix_view.subDim(); ++i) {
        for ( int j = 0; j < matrix_view.numSubCols(); ++j) {
          matrix_view(i, j) = alpha * f_fad[i].dx(x_dot_idx_offset+j)
            + beta * f_fad[i].dx(x_idx_offset + j);
        }
      }
    }
    else {
      for ( int i = 0; i < matrix_view.subDim(); ++i) {
        for ( int j = 0; j < matrix_view.numSubCols(); ++j) {
          matrix_view(i, j) = f_fad[i].dx(x_idx_offset + j);
        }
      }
    }
  }

  if (nonnull(DfDp_out)) {
    Thyra::DetachedMultiVectorView<double> DfDp_out_view( *DfDp_out );
    for ( int i = 0; i < DfDp_out_view.subDim(); ++i )
      DfDp_out_view(i,0) = f_fad[i].dx(eps_idx_offset);
  }

}


RCP<const Thyra::VectorSpaceBase<double> > VanderPolModel::get_p_space(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_;
}

RCP<const Teuchos::Array<std::string> > VanderPolModel::get_p_names(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  RCP<Teuchos::Array<std::string> > p_strings = 
    Teuchos::rcp(new Teuchos::Array<std::string>());
  p_strings->push_back("Model Coefficient:  epsilon");
  return p_strings;
}

RCP<const Thyra::VectorSpaceBase<double> > VanderPolModel::get_g_space(int j) const
{
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_;
}

// private

void VanderPolModel::setupInOutArgs_() 
{
  if (!isInitialized_) {
    
    {
      // Set up prototypical InArgs
      ModelEvaluatorBase::InArgsSetup<double> inArgs;
      inArgs.setModelEvalDescription(this->description());
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_t );
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x );
      if (isImplicit_) {
        inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x_dot );
        inArgs.setSupports( ModelEvaluatorBase::IN_ARG_alpha );
        inArgs.setSupports( ModelEvaluatorBase::IN_ARG_beta );
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
        x_ic_view[0] = x0_ic_;
        x_ic_view[1] = x1_ic_;
      }
      nominalValues_.set_x(x_ic);
      if (acceptModelParams_) {
        const RCP<VectorBase<double> > p_ic = createMember(p_space_);
        {
          Thyra::DetachedVectorView<double> p_ic_view( *p_ic );
          p_ic_view[0] = epsilon_;
        }
        nominalValues_.set_p(0,p_ic);
      }
      if (isImplicit_) {
        const RCP<VectorBase<double> > x_dot_ic = createMember(x_space_);
        { // scope to delete DetachedVectorView
          Thyra::DetachedVectorView<double> x_dot_ic_view( *x_dot_ic );
          x_dot_ic_view[0] = 0.0;
          x_dot_ic_view[1] = 0.0;
        }
        nominalValues_.set_x_dot(x_dot_ic);
      }
    }
    isInitialized_ = true;

  }

}

void VanderPolModel::setParameterList(RCP<ParameterList> const& paramList)
{
  using Teuchos::get;
  TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
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
  epsilon_ = get<double>(*pl,Coeff_epsilon_name);
  x0_ic_ = get<double>(*pl,IC_x0_name);
  x1_ic_ = get<double>(*pl,IC_x1_name);
  t0_ic_ = get<double>(*pl,IC_t0_name);
  setupInOutArgs_();
}

RCP<const ParameterList> VanderPolModel::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(Implicit_name,Implicit_default);
    pl->set(AcceptModelParams_name, AcceptModelParams_default);
    pl->set(HaveIC_name, HaveIC_default);
    Teuchos::setDoubleParameter(
        Coeff_epsilon_name, Coeff_epsilon_default, "Coefficient epsilon in model", 
        &*pl 
        );
    Teuchos::setDoubleParameter(
        IC_x0_name, IC_x0_default, "Initial Condition for x0", 
        &*pl 
        );
    Teuchos::setDoubleParameter(
        IC_x1_name, IC_x1_default, "Initial Condition for x1", 
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


template<class ScalarT>
void VanderPolModel::eval_f(
  const ArrayView<const ScalarT> &x_dot,
  const ArrayView<const ScalarT> &x,
  const ScalarT &eps,
  const ScalarT &t,
  const ArrayView<ScalarT> &f
  ) const
{
  
  const ScalarT x0 = 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
  const ScalarT x1 = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  
  const ScalarT x1prime = -2*cos(t) + eps*(-0.75*sin(t)+9.0/4.0*sin(3.0*t));
  const ScalarT forcing_term = x1prime-eps*(1.0-x0*x0)*x1+x0;

  if (!isImplicit_) { // EXPLICIT
    f[0] = x[1];
    f[1] = eps*(1.0-x[0]*x[0])*x[1]-x[0] + forcing_term;
  }
  else { // IMPLICIT
    f[0] = x_dot[0] - x[1];
    f[1] = x_dot[1] - (eps*(1.0-x[0]*x[0])* x[1] - x[0]) - forcing_term;
  }

}


} // namespace Rythmos


// non-member Constructors


Teuchos::RCP<Rythmos::VanderPolModel>
Rythmos::vanderPolModel() 
{
  RCP<VanderPolModel> model = rcp(new VanderPolModel);
  return(model);
}


Teuchos::RCP<Rythmos::VanderPolModel>
Rythmos::vanderPolModel(bool implicit) 
{
  RCP<VanderPolModel> model = vanderPolModel();
  model->setImplicitFlag(implicit);
  return(model);
}


Teuchos::RCP<Rythmos::VanderPolModel>
Rythmos::vanderPolModel(const RCP<ParameterList> &pl)
{
  RCP<VanderPolModel> model = vanderPolModel();
  model->setParameterList(pl);
  return(model);
}

#endif // Rythmos_ENABLE_Sacado
