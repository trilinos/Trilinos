// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Epetra_RowMatrix.h"


namespace {


const std::string StateFunctionScaling_name = "State Function Scaling";
Teuchos::RefCountPtr<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::EpetraModelEvaluator::EStateFunctionScaling
    >
  >
stateFunctionScalingValidator;
const std::string StateFunctionScaling_default = "None";

const std::string DoStateVariableScalingFirst_name = "Do State Variable Scaling First";
const bool DoStateVariableScalingFirst_default = false;

Teuchos::RefCountPtr<Epetra_RowMatrix>
get_Epetra_RowMatrix(
  const EpetraExt::ModelEvaluator::OutArgs &epetraOutArgs
  )
{
  const Teuchos::RefCountPtr<Epetra_Operator>
    eW = epetraOutArgs.get_W();
  const Teuchos::RefCountPtr<Epetra_RowMatrix>
    ermW = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(eW,false);
  TEST_FOR_EXCEPTION(
    is_null(ermW), std::logic_error,
    "Thyra::EpetraModelEvaluator::evalModel(...): Error, if\n"
    "scaling is turned on, the underlying Epetra_Operator created\n"
    "an initialized by the underlying epetra model evaluator\n"
    "\"" << epetraOutArgs.modelEvalDescription() << "\"\n"
    "must support the Epetra_RowMatrix interface through a dynamic cast.\n"
    "The concrete type " << Teuchos::typeName(*eW) << " does not support\n"
    "Epetra_RowMatrix!"
    );
  return ermW;
}


} // namespace


namespace Thyra {


// Constructors/initializers/accessors.


EpetraModelEvaluator::EpetraModelEvaluator()
  :finalPointWasSolved_(false),stateFunctionScaling_(STATE_FUNC_SCALING_NONE),
  doStateVariableScalingFirst_(false)
{}


EpetraModelEvaluator::EpetraModelEvaluator(
  const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> &epetraModel,
  const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > &W_factory
  )
  :finalPointWasSolved_(false),stateFunctionScaling_(STATE_FUNC_SCALING_NONE),
  doStateVariableScalingFirst_(false)
{
  initialize(epetraModel,W_factory);
}

void EpetraModelEvaluator::initialize(
  const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> &epetraModel,
  const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > &W_factory
  )
{
  using Teuchos::implicit_cast;
  typedef ModelEvaluatorBase MEB;
  //
  epetraModel_ = epetraModel;
  //
  W_factory_ = W_factory;
  //
  x_space_ = create_VectorSpace( x_map_ = epetraModel_->get_x_map() );
  f_space_ = create_VectorSpace( f_map_ = epetraModel_->get_f_map() );
  //
  EpetraExt::ModelEvaluator::InArgs inArgs = epetraModel_->createInArgs();
  p_map_.resize(inArgs.Np()); p_space_.resize(inArgs.Np());
  for( int l = 0; l < implicit_cast<int>(p_space_.size()); ++l )
    p_space_[l] = create_VectorSpace( p_map_[l] = epetraModel_->get_p_map(l) );
  //
  EpetraExt::ModelEvaluator::OutArgs outArgs = epetraModel_->createOutArgs();
  g_map_.resize(outArgs.Ng()); g_space_.resize(outArgs.Ng());
  for( int j = 0; j < implicit_cast<int>(g_space_.size()); ++j )
    g_space_[j] = create_VectorSpace( g_map_[j] = epetraModel_->get_g_map(j) );
  //
  initialGuess_ = this->createInArgs();
  lowerBounds_ = this->createInArgs();
  upperBounds_ = this->createInArgs();
  if(initialGuess_.supports(MEB::IN_ARG_x)) {
    initialGuess_.set_x( create_Vector( epetraModel_->get_x_init(), x_space_ ) );
    lowerBounds_.set_x( create_Vector( epetraModel_->get_x_lower_bounds(), x_space_ ) );
    upperBounds_.set_x( create_Vector( epetraModel_->get_x_upper_bounds(), x_space_ ) );
  }
  if(initialGuess_.supports(MEB::IN_ARG_x_dot)) {
    initialGuess_.set_x_dot( create_Vector( epetraModel_->get_x_dot_init(), x_space_ ) );
  }
  for( int l = 0; l < implicit_cast<int>(p_space_.size()); ++l ) {
    initialGuess_.set_p( l, create_Vector( epetraModel_->get_p_init(l), p_space_[l] ) );
    lowerBounds_.set_p( l, create_Vector( epetraModel_->get_p_lower_bounds(l), p_space_[l] ) );
    upperBounds_.set_p( l, create_Vector( epetraModel_->get_p_upper_bounds(l), p_space_[l] ) );
  }
  if(initialGuess_.supports(MEB::IN_ARG_t)) {
    initialGuess_.set_t(epetraModel_->get_t_init());
    lowerBounds_.set_t(epetraModel_->get_t_lower_bound());
    upperBounds_.set_t(epetraModel_->get_t_upper_bound());
  }
  finalPointWasSolved_ = false;
  stateFunctionScalingVec_ = Teuchos::null;
}


Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>
EpetraModelEvaluator::getEpetraModel() const
{
  return epetraModel_;
}


void EpetraModelEvaluator::setInitialGuess(
  const ModelEvaluatorBase::InArgs<double>& initialGuess
 )
{
  initialGuess_.setArgs(initialGuess);
}


void EpetraModelEvaluator::setStateVariableScalingVec(
  const Teuchos::RefCountPtr<const Epetra_Vector> &stateVariableScalingVec
  )
{

  typedef ModelEvaluatorBase MEB;

  stateVariableScalingVec_ = stateVariableScalingVec.assert_not_null();
  invStateVariableScalingVec_ = Teuchos::rcp(
    new Epetra_Vector(*x_map_)
    );
  invStateVariableScalingVec_->Reciprocal(*stateVariableScalingVec_);
  
  if (
    initialGuess_.supports(MEB::IN_ARG_x)
    && !is_null(initialGuess_.get_x())
    )
  {
    Teuchos::RefCountPtr<const VectorBase<double> >
      x_unscaled = initialGuess_.get_x();
    Teuchos::RefCountPtr<Epetra_Vector>
      e_x_scaled = Teuchos::rcp(new Epetra_Vector(*x_map_));
    e_x_scaled->Multiply( 1.0, *stateVariableScalingVec_,
      *Thyra::get_Epetra_Vector(*x_map_,x_unscaled), 0.0
      );
    initialGuess_.set_x( create_Vector( e_x_scaled, x_space_ ) );
  }
  if (
    initialGuess_.supports(MEB::IN_ARG_x_dot)
    && !is_null(initialGuess_.get_x_dot())
    )
  {
    Teuchos::RefCountPtr<const VectorBase<double> >
      x_dot_unscaled = initialGuess_.get_x_dot();
    Teuchos::RefCountPtr<Epetra_Vector>
      e_x_dot_scaled = Teuchos::rcp(new Epetra_Vector(*x_map_));
    e_x_dot_scaled->Multiply( 1.0, *stateVariableScalingVec_,
      *Thyra::get_Epetra_Vector(*x_map_,x_dot_unscaled), 0.0
      );
    initialGuess_.set_x_dot( create_Vector( e_x_dot_scaled, x_space_ ) );
  }

}


Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEvaluator::getStateVariableScalingVec() const
{
  return stateVariableScalingVec_;
}


void EpetraModelEvaluator::setStateFunctionScalingVec(
  const Teuchos::RefCountPtr<const Epetra_Vector> &stateFunctionScalingVec
  )
{
  stateFunctionScalingVec_ = stateFunctionScalingVec;
}


Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEvaluator::getStateFunctionScalingVec() const
{
  return stateFunctionScalingVec_;
}


void EpetraModelEvaluator::uninitialize(
  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> *epetraModel,
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> > *W_factory
  )
{
  if(epetraModel) *epetraModel = epetraModel_;
  if(W_factory) *W_factory = W_factory_;
  epetraModel_ = Teuchos::null;
  W_factory_ = Teuchos::null;
  stateFunctionScalingVec_ = Teuchos::null;
  stateVariableScalingVec_ = Teuchos::null;
  invStateVariableScalingVec_ = Teuchos::null;
}


const ModelEvaluatorBase::InArgs<double>&
EpetraModelEvaluator::getFinalPoint() const
{
  return finalPoint_;
}


bool EpetraModelEvaluator::finalPointWasSolved() const
{
  return finalPointWasSolved_;
}


// Overridden from Teuchos::ParameterListAcceptor


void EpetraModelEvaluator::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters(),0); // Just validate my params
  paramList_ = paramList;
  const EStateFunctionScaling stateFunctionScaling_old = stateFunctionScaling_; 
  stateFunctionScaling_ = stateFunctionScalingValidator->getIntegralValue(
    *paramList_, StateFunctionScaling_name, StateFunctionScaling_default
    );
  if( stateFunctionScaling_ != stateFunctionScaling_old )
    stateFunctionScalingVec_ = Teuchos::null;
  doStateVariableScalingFirst_ = paramList_->get(
    DoStateVariableScalingFirst_name,DoStateVariableScalingFirst_default );
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}


Teuchos::RefCountPtr<Teuchos::ParameterList>
EpetraModelEvaluator::getParameterList()
{
  return paramList_;
}


Teuchos::RefCountPtr<Teuchos::ParameterList>
EpetraModelEvaluator::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


Teuchos::RefCountPtr<const Teuchos::ParameterList>
EpetraModelEvaluator::getParameterList() const
{
  return paramList_;
}


Teuchos::RefCountPtr<const Teuchos::ParameterList>
EpetraModelEvaluator::getValidParameters() const
{
  using Teuchos::rcp;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::tuple;
  static Teuchos::RefCountPtr<const Teuchos::ParameterList> validPL;
  if(is_null(validPL)) {
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    stateFunctionScalingValidator = rcp(
      new StringToIntegralParameterEntryValidator<EStateFunctionScaling>(
        tuple<std::string>(
          "None",
          "Row Sum"
          ),
        tuple<std::string>(
          "Do not scale the state function f(...) in this class.",

          "Scale the state function f(...) and all its derivatives\n"
          "using the row sum scaling from the initial Jacobian\n"
          "W=d(f)/d(x).  Note, this only works with Epetra_CrsMatrix\n"
          "currently."
          ),
        tuple<EStateFunctionScaling>(
          STATE_FUNC_SCALING_NONE,
          STATE_FUNC_SCALING_ROW_SUM
          ),
        StateFunctionScaling_name
        )
      );
    pl->set(StateFunctionScaling_name,StateFunctionScaling_default,
      "Determines if and how the state function f(...) and all of its\n"
      "derivatives are scaled.  The scaling is done explicitly so there should\n"
      "be no impact on the meaning of inner products or tolerances for\n"
      "linear solves.",
      stateFunctionScalingValidator
      );
    pl->set(
      DoStateVariableScalingFirst_name,DoStateVariableScalingFirst_default,
      "If true, the Jacobian D=d(f)/d(x) will be scaled from the left with\n"
      "inv(S_x) before scaling is computed on the right."
      );
    validPL = pl;
  }
  return validPL;
}


// Overridden from ModelEvaulator.


int EpetraModelEvaluator::Np() const
{
  return p_space_.size();
}


int EpetraModelEvaluator::Ng() const
{
  return g_space_.size();
}


Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_x_space() const
{
  return x_space_;
}


Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_f_space() const
{
  return f_space_;
}


Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_p_space(int l) const
{
  TEST_FOR_EXCEPT( ! ( 0 <= l && l < this->Np() ) );
  return p_space_[l];
}


Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_g_space(int j) const
{
  TEST_FOR_EXCEPT( ! ( 0 <= j && j < this->Ng() ) );
  return g_space_[j];
}


ModelEvaluatorBase::InArgs<double> EpetraModelEvaluator::getNominalValues() const
{
  return initialGuess_;
}


ModelEvaluatorBase::InArgs<double> EpetraModelEvaluator::getLowerBounds() const
{
  return lowerBounds_;
}


ModelEvaluatorBase::InArgs<double> EpetraModelEvaluator::getUpperBounds() const
{
  return upperBounds_;
}


Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
EpetraModelEvaluator::create_W() const
{
  TEST_FOR_EXCEPTION(
    W_factory_.get()==NULL, std::logic_error
    ,"Thyra::EpetraModelEvaluator::create_W(): "
    "Error, the client did not set a LinearOpWithSolveFactoryBase"
    " object for W!"
    );
  W_factory_->setOStream(this->getOStream());
  W_factory_->setVerbLevel(this->getVerbLevel());
  return W_factory_->createOp();
}


Teuchos::RefCountPtr<LinearOpBase<double> >
EpetraModelEvaluator::create_W_op() const
{
  return Teuchos::rcp(new Thyra::EpetraLinearOp());
}


Teuchos::RefCountPtr<LinearOpBase<double> >
EpetraModelEvaluator::create_DfDp_op(int l) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


Teuchos::RefCountPtr<LinearOpBase<double> >
EpetraModelEvaluator::create_DgDx_op(int j) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


Teuchos::RefCountPtr<LinearOpBase<double> >
EpetraModelEvaluator::create_DgDp_op( int j, int l ) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


ModelEvaluatorBase::InArgs<double> EpetraModelEvaluator::createInArgs() const
{
  const EpetraExt::ModelEvaluator &epetraModel = *epetraModel_;
  InArgsSetup<double> inArgs;
  typedef EpetraExt::ModelEvaluator EME;
  EME::InArgs epetraInArgs = epetraModel.createInArgs();
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(epetraInArgs.Np());
  inArgs.setSupports(IN_ARG_x_dot,epetraInArgs.supports(EME::IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x,epetraInArgs.supports(EME::IN_ARG_x));
  inArgs.setSupports(IN_ARG_x_dot_poly,
         epetraInArgs.supports(EME::IN_ARG_x_dot_poly));
  inArgs.setSupports(IN_ARG_x_poly,epetraInArgs.supports(EME::IN_ARG_x_poly));
  inArgs.setSupports(IN_ARG_t,epetraInArgs.supports(EME::IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha,epetraInArgs.supports(EME::IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta,epetraInArgs.supports(EME::IN_ARG_beta));
  return inArgs;
}


ModelEvaluatorBase::OutArgs<double> EpetraModelEvaluator::createOutArgs() const
{
  const EpetraExt::ModelEvaluator &epetraModel = *epetraModel_;
  OutArgsSetup<double> outArgs;
  typedef EpetraExt::ModelEvaluator EME;
  EME::InArgs  epetraInArgs  = epetraModel.createInArgs();
  EME::OutArgs epetraOutArgs = epetraModel.createOutArgs();
  const int Np = epetraOutArgs.Np(), Ng = epetraOutArgs.Ng();
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng);
  outArgs.setSupports(OUT_ARG_f,epetraOutArgs.supports(EME::OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W,epetraOutArgs.supports(EME::OUT_ARG_W)&&W_factory_.get()!=NULL);
  outArgs.setSupports(OUT_ARG_W_op,epetraOutArgs.supports(EME::OUT_ARG_W));
  outArgs.set_W_properties(convert(epetraOutArgs.get_W_properties()));
  for(int l=0; l<Np; ++l) {
    outArgs.setSupports(OUT_ARG_DfDp,l,convert(epetraOutArgs.supports(EME::OUT_ARG_DfDp,l)));
    if(!outArgs.supports(OUT_ARG_DfDp,l).none())
      outArgs.set_DfDp_properties(l,convert(epetraOutArgs.get_DfDp_properties(l)));
  }
  for(int j=0; j<Ng; ++j) {
    outArgs.setSupports(OUT_ARG_DgDx,j,convert(epetraOutArgs.supports(EME::OUT_ARG_DgDx,j)));
    if(!outArgs.supports(OUT_ARG_DgDx,j).none())
      outArgs.set_DgDx_properties(j,convert(epetraOutArgs.get_DgDx_properties(j)));
  }
  for(int j=0; j<Ng; ++j) for(int l=0; l<Np; ++l) {
    outArgs.setSupports(OUT_ARG_DgDp,j,l,convert(epetraOutArgs.supports(EME::OUT_ARG_DgDp,j,l)));
    if(!outArgs.supports(OUT_ARG_DgDp,j,l).none())
      outArgs.set_DgDp_properties(j,l,convert(epetraOutArgs.get_DgDp_properties(j,l)));
  }
  outArgs.setSupports(OUT_ARG_f_poly,epetraOutArgs.supports(EME::OUT_ARG_f_poly));
  return outArgs;
}


void EpetraModelEvaluator::evalModel( const InArgs<double>& inArgs_in, const OutArgs<double>& outArgs ) const
{

  using Thyra::get_Epetra_Vector;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::implicit_cast;

  typedef EpetraExt::ModelEvaluator EME;

  // State function Scaling
  const bool firstTimeStateFuncScaling
    = (
      stateFunctionScaling_ != STATE_FUNC_SCALING_NONE
      && is_null(stateFunctionScalingVec_)
      );

  InArgs<double> inArgs = initialGuess_; // Make sure we grab the initial guess first!
  inArgs.setArgs(inArgs_in);

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::EpetraModelEvaluator::evalModel(...) ...\n";

  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\ninArgs =\n" << Teuchos::describe(inArgs,verbLevel)
      << "\noutArgs on input =\n" << Teuchos::describe(outArgs,Teuchos::VERB_LOW);
  
  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveFactoryBase<double> > VOTSLOWSF;
  VOTSLOWSF W_factory_outputTempState(W_factory_,out,verbLevel);
  
  // InArgs
  
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nSetting-up/creating input arguments ...\n";
  timer.start(true);

  EME::InArgs epetraInArgs = epetraModel_->createInArgs();

  RefCountPtr<const VectorBase<double> > x_dot;
  if( inArgs.supports(IN_ARG_x_dot) && (x_dot = inArgs.get_x_dot()).get() ) {
    RefCountPtr<const Epetra_Vector> e_x_dot = get_Epetra_Vector(*x_map_,x_dot);
    if(!is_null(invStateVariableScalingVec_)) {
      if(is_null(x_dot_unscaled_))
        x_dot_unscaled_ = rcp(new Epetra_Vector(*x_map_));
      x_dot_unscaled_->Multiply( 1.0, *invStateVariableScalingVec_, *e_x_dot, 0.0 );
      epetraInArgs.set_x_dot(x_dot_unscaled_);
    }
    else {
      epetraInArgs.set_x_dot(e_x_dot);
    }
  }

  RefCountPtr<const VectorBase<double> > x;
  if( inArgs.supports(IN_ARG_x) && (x = inArgs.get_x()).get() ) {
    RefCountPtr<const Epetra_Vector> e_x = get_Epetra_Vector(*x_map_,x);
    if(!is_null(invStateVariableScalingVec_)) {
      if(is_null(x_unscaled_))
        x_unscaled_ = rcp(new Epetra_Vector(*x_map_));
      x_unscaled_->Multiply( 1.0, *invStateVariableScalingVec_, *e_x, 0.0 );
      epetraInArgs.set_x(x_unscaled_);
    }
    else {
      epetraInArgs.set_x(e_x);
    }
  }

  {
    RefCountPtr<const VectorBase<double> > p_l;
    for(int l = 0;  l < outArgs.Np(); ++l ) {
      p_l = inArgs.get_p(l);
      if(p_l.get()) epetraInArgs.set_p(l,get_Epetra_Vector(*p_map_[l],p_l));
    }
  }

  RefCountPtr<const Teuchos::Polynomial< VectorBase<double> > > x_dot_poly;
  Teuchos::RefCountPtr<Epetra_Vector> epetra_ptr;
  if( inArgs.supports(IN_ARG_x_dot_poly) && \
    (x_dot_poly = inArgs.get_x_dot_poly()).get() ) {
    RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > epetra_x_dot_poly = 
      Teuchos::rcp(new Teuchos::Polynomial<Epetra_Vector>(x_dot_poly->degree()));
    for (unsigned int i=0; i<=x_dot_poly->degree(); i++) {
      epetra_ptr = 
        Teuchos::rcp_const_cast<Epetra_Vector>(
          get_Epetra_Vector(*x_map_, x_dot_poly->getCoefficient(i)));
      epetra_x_dot_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraInArgs.set_x_dot_poly(epetra_x_dot_poly);
  }
  
  RefCountPtr<const Teuchos::Polynomial< VectorBase<double> > > x_poly;
  if( inArgs.supports(IN_ARG_x_poly) && \
    (x_poly = inArgs.get_x_poly()).get() ) {
    RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > epetra_x_poly = 
      Teuchos::rcp(new Teuchos::Polynomial<Epetra_Vector>(x_poly->degree()));
    for (unsigned int i=0; i<=x_poly->degree(); i++) {
      epetra_ptr = 
        Teuchos::rcp_const_cast<Epetra_Vector>(
          get_Epetra_Vector(*x_map_, x_poly->getCoefficient(i)));
      epetra_x_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraInArgs.set_x_poly(epetra_x_poly);
  }

  if( inArgs.supports(IN_ARG_t) )
    epetraInArgs.set_t(inArgs.get_t());

  if( inArgs.supports(IN_ARG_alpha) )
    epetraInArgs.set_alpha(inArgs.get_alpha());

  if( inArgs.supports(IN_ARG_beta) )
    epetraInArgs.set_beta(inArgs.get_beta());

  timer.stop();
  if(out.get() && verbLevel >= Teuchos::VERB_LOW)
    OSTab(out).o() << "\nTime to setup InArgs = "<<timer.totalElapsedTime()<<" sec\n";

  // OutArgs
  
  if(out.get() && verbLevel >= Teuchos::VERB_LOW)
    *out << "\nSetting-up/creating output arguments ...\n";
  timer.start(true);

  EME::OutArgs epetraOutArgs = epetraModel_->createOutArgs();

  RefCountPtr<VectorBase<double> > f;
  if( outArgs.supports(OUT_ARG_f) && (f = outArgs.get_f()).get() )
    epetraOutArgs.set_f(get_Epetra_Vector(*f_map_,f));

  {
    Teuchos::RefCountPtr<VectorBase<double> > g_j;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      g_j = outArgs.get_g(j);
      if(g_j.get()) epetraOutArgs.set_g(j,get_Epetra_Vector(*g_map_[j],g_j));
    }
  }
  
  RefCountPtr<LinearOpWithSolveBase<double> > W;
  RefCountPtr<LinearOpBase<double> >          W_op;
  RefCountPtr<const LinearOpBase<double> >    fwdW;
  RefCountPtr<EpetraLinearOp>                 efwdW;
  if( outArgs.supports(OUT_ARG_W) && (W = outArgs.get_W()).get() ) {
    Thyra::uninitializeOp<double>(*W_factory_,&*W,&fwdW);
    if(fwdW.get()) {
      efwdW = rcp_const_cast<EpetraLinearOp>(
        rcp_dynamic_cast<const EpetraLinearOp>(fwdW,true));
    }
    else {
      efwdW = Teuchos::rcp(new EpetraLinearOp());
      fwdW = efwdW;
    }
  }
  if( outArgs.supports(OUT_ARG_W_op) && (W_op = outArgs.get_W_op()).get() ) {
    if( W_op.get() && !efwdW.get() )
      efwdW = rcp_const_cast<EpetraLinearOp>(
        rcp_dynamic_cast<const EpetraLinearOp>(W_op,true));
  }
  RefCountPtr<Epetra_Operator> eW;
  if(efwdW.get()) {
    eW = efwdW->epetra_op();
    if(!eW.get())
      eW = epetraModel_->create_W();
    epetraOutArgs.set_W(eW);
  }
  // NOTE: Above, if both W and W_op are set and have been through at least
  // one prior evaluation (and therefore have Epetra_Operator objects embedded
  // in them), then we will use the Epetra_Operator embedded in W to pass to
  // the EpetraExt::ModelEvaluator object and ignore the Epetra_Operator
  // object in W_op.  In the standard use case, these will be the same
  // Epetra_Operator objects.  However, it is possible that the client could
  // use this interface in such a way that these would have different
  // Epetra_Operator objects embedded in them.  In this (very unlikely) case,
  // the Epetra_Operator embedded in W_op will be discarded!  This might be
  // surprising to a client but it is very unlikely that this will ever be a
  // problem, but the issue is duly noted here!  Only dangerous programming
  // use of this interface would cause any problem.

  // Note: The following derivative objects update in place!

  {
    Derivative<double> DfDp_l;
    for(int l = 0;  l < outArgs.Np(); ++l ) {
      if( !outArgs.supports(OUT_ARG_DfDp,l).none()
        && !(DfDp_l = outArgs.get_DfDp(l)).isEmpty() )
      {
        epetraOutArgs.set_DfDp(l,convert(DfDp_l,f_map_,p_map_[l]));
      }
    }
  }

  {
    Derivative<double> DgDx_j;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      if( !outArgs.supports(OUT_ARG_DgDx,j).none()
        && !(DgDx_j = outArgs.get_DgDx(j)).isEmpty() )
      {
        epetraOutArgs.set_DgDx(j,convert(DgDx_j,g_map_[j],x_map_));
      }
    }
  }

  {
    Derivative<double> DgDp_j_l;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      for(int l = 0;  l < outArgs.Np(); ++l ) {
        if( !outArgs.supports(OUT_ARG_DgDp,j,l).none()
          && !(DgDp_j_l = outArgs.get_DgDp(j,l)).isEmpty() )
        {
          epetraOutArgs.set_DgDp(j,l,convert(DgDp_j_l,g_map_[j],p_map_[l]));
        }
      }
    }
  }

  RefCountPtr<const Teuchos::Polynomial< VectorBase<double> > > f_poly;
  if( outArgs.supports(OUT_ARG_f_poly) && (f_poly = outArgs.get_f_poly()).get() )
  {
    RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > epetra_f_poly = 
      Teuchos::rcp(new Teuchos::Polynomial<Epetra_Vector>(f_poly->degree()));
    for (unsigned int i=0; i<=f_poly->degree(); i++) {
      epetra_ptr = 
        Teuchos::rcp_const_cast<Epetra_Vector>(get_Epetra_Vector(*f_map_,
            f_poly->getCoefficient(i)));
      epetra_f_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraOutArgs.set_f_poly(epetra_f_poly);
  }

  bool createdTempEpetraW = false;
  if (
    firstTimeStateFuncScaling
    && ( stateFunctionScaling_ == STATE_FUNC_SCALING_ROW_SUM )
    && (
      epetraOutArgs.supports(EME::OUT_ARG_f) 
      && epetraOutArgs.funcOrDerivesAreSet(EME::OUT_ARG_f)
      )
    && (
      epetraOutArgs.supports(EME::OUT_ARG_W)
      && is_null(epetraOutArgs.get_W())
      )
    )
  {
    // This is the first pass through with scaling turned on and the client
    // turned on automatic scaling but did not ask for W.  We must compute W
    // in order to compute the scale factors so we must allocate a temporary W
    // just to compute the scale factors and then throw it away.  If the
    // client wants to evaluate W at the same point, then it should have
    // passed W in but that is not our problem here.  The ModelEvaluator
    // relies on the client to set up the calls to allow for efficient
    // evaluation.

    if(out.get() && verbLevel >= Teuchos::VERB_LOW)
      *out
        << "\nCreating a temporary Epetra W to compute scale factors"
        << " for f(...) ...\n";
    epetraOutArgs.set_W(epetraModel_->create_W());
    createdTempEpetraW = true; // This flag will tell us to delete W later!
  }

  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to setup OutArgs = "<<timer.totalElapsedTime()<<" sec\n";

  // Do the evaluation
  
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEvaluating the output functions ...\n";
  timer.start(true);

  epetraModel_->evalModel(epetraInArgs,epetraOutArgs);

  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to evaluate output functions = "<<timer.totalElapsedTime()<<" sec\n";

  // Postprocess arguments
  
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nPost processing the output objects ...\n";
  timer.start(true);

  if (doStateVariableScalingFirst_) {

    // Apply scaling to Epetra for scaling in x.  We must do this before the
    // right matrix scaling below in order to ensure that good scale factors are
    // computed for a good matrix scaling for the linear solve.
    if (
      !is_null(invStateVariableScalingVec_)
      && epetraOutArgs.supports(EME::OUT_ARG_W)
      && !is_null(epetraOutArgs.get_W())
      )
    {
      if(out.get() && verbLevel >= Teuchos::VERB_LOW)
        *out << "\nScaling Epetra W from the right with invStateVariableScalingVec ...\n";
      const RefCountPtr<Epetra_RowMatrix>
        ermW = get_Epetra_RowMatrix(epetraOutArgs);
      ermW->RightScale(*invStateVariableScalingVec_);
    }

  }

  // Compute the scale factors for the state function f(...)
  if (firstTimeStateFuncScaling) {
    switch(stateFunctionScaling_) {
      case STATE_FUNC_SCALING_ROW_SUM: {
        // Compute the inverse row-sum scaling from W
        const RefCountPtr<Epetra_RowMatrix>
          ermW = get_Epetra_RowMatrix(epetraOutArgs);
        // Note: Above, we get the Epetra W object directly from the Epetra
        // OutArgs object since this might be a temporary matrix just to
        // compute scaling factors.  In this case, the stack funtion variable
        // eW might be empty!
        RefCountPtr<Epetra_Vector>
          invRowSums = rcp(new Epetra_Vector(ermW->OperatorRangeMap()));
        // Above, From the documentation is seems that the RangeMap should be
        // okay but who knows for sure!
        ermW->InvRowSums(*invRowSums);
        if(out.get() && verbLevel >= Teuchos::VERB_LOW) {
          *out
            << "\nComputed inverse row sum scaling from W that"
            " will be used to scale f(...) and its derivatives:\n";
          double minVal = 0, maxVal = 0, avgVal = 0;
          invRowSums->MinValue(&minVal);
          invRowSums->MaxValue(&maxVal);
          invRowSums->MeanValue(&avgVal);
          OSTab tab(out);
          *out
            << "min(invRowSums) = " << minVal << "\n"
            << "max(invRowSums) = " << maxVal << "\n"
            << "avg(invRowSums) = " << avgVal << "\n";
        }
        stateFunctionScalingVec_ = invRowSums;
        break;
      }
      default:
        TEST_FOR_EXCEPT("Should never get here!");
    }
  }

  // Apply scaling to Epetra objects f(...) objects
  if (!is_null(stateFunctionScalingVec_)) {

    // Scale f
    RefCountPtr<Epetra_Vector> ef;
    if( epetraOutArgs.supports(EME::OUT_ARG_f) && !is_null(ef=epetraOutArgs.get_f()) )
    {
      if(out.get() && verbLevel >= Teuchos::VERB_LOW)
        *out << "\nScaling Epetra f with stateFunctionScalingVec ...\n";
      ef->Multiply(1.0, *stateFunctionScalingVec_, *ef, 0.0);
    }

    // Scale W
    if(!is_null(eW)) {
      if(out.get() && verbLevel >= Teuchos::VERB_LOW)
        *out << "\nScaling Epetra W from the left with stateFunctionScalingVec ...\n";
      const RefCountPtr<Epetra_RowMatrix>
        ermW = rcp_dynamic_cast<Epetra_RowMatrix>(eW); // Should not throw
      ermW->LeftScale(*stateFunctionScalingVec_);
      // Note: The above scaled W will be used to create the linear solver
      // below ...
      // Also note that we don't grab the Epetra W out of epetraOutArgs since
      // this may just be a temporary matrix that was just used to compute scale
      // vectors.
    }

    // Scale DfDp(l)
    for ( int l = 0; l < epetraOutArgs.Np(); ++l ) {
      EME::Derivative deriv_DfDp_l;
      if (
        !epetraOutArgs.supports(EME::OUT_ARG_DfDp,l).none()
        && !(deriv_DfDp_l=epetraOutArgs.get_DfDp(l)).isEmpty()
        )
      {
        RefCountPtr<Epetra_MultiVector>
          e_mv_DfDp_l = deriv_DfDp_l.getDerivativeMultiVector().getMultiVector();
        RefCountPtr<Epetra_Operator>
          e_op_DfDp_l = deriv_DfDp_l.getLinearOp();
        if (!is_null(e_mv_DfDp_l)) {
          TEST_FOR_EXCEPTION(
            deriv_DfDp_l.getDerivativeMultiVector().getOrientation() != EME::DERIV_MV_BY_COL,
            std::logic_error,
            "Error, the model \""<<epetraOutArgs.modelEvalDescription()<<"\" computes DfDp("<<l<<") "
            "with teh wrong orientation!"
            );
          if(out.get() && verbLevel >= Teuchos::VERB_LOW)
            *out << "\nScaling Epetra DfDp("<<l<<") multi-vector with stateFunctionScalingVec ...\n";
          e_mv_DfDp_l->Multiply(1.0, *stateFunctionScalingVec_, *e_mv_DfDp_l, 0.0);
        }
        if (!is_null(e_op_DfDp_l)) {
          TEST_FOR_EXCEPTION(true,std::logic_error,
            "Error, can't handle scaling of a Epetra_Operator version of DfDp("<<l<<") yet!"
            );
        }
      }
      // ToDo: Scale second derivatives of f???
    }

  }

  // Apply scaling to Epetra objects for scaling in x
  if(!is_null(invStateVariableScalingVec_)) {

    // Scale W
    if ( !doStateVariableScalingFirst_ && !is_null(eW) ) {
      if(out.get() && verbLevel >= Teuchos::VERB_LOW)
        *out << "\nScaling Epetra W from the right with invStateVariableScalingVec ...\n";
      const RefCountPtr<Epetra_RowMatrix>
        ermW = rcp_dynamic_cast<Epetra_RowMatrix>(eW); // Should not throw
      ermW->RightScale(*invStateVariableScalingVec_);
      // Note: The above scaled W will be used to create the linear solver
      // below ...
      // Also note that we don't grab the Epetra W out of epetraOutArgs since
      // this may just be a temporary matrix that was just used to compute scale
      // vectors.
    }

    // Scale DgDx
    {
      EME::Derivative DgDx_j;
      for(int j = 0;  j < epetraOutArgs.Ng(); ++j ) {
        if( !epetraOutArgs.supports(EME::OUT_ARG_DgDx,j).none()
          && !(DgDx_j = epetraOutArgs.get_DgDx(j)).isEmpty() )
        {
          TEST_FOR_EXCEPT("Can't scale DgDx yet!");
        }
      }
    }

    // ToDo: Scale second derivatives when added!

  }

  // ToDo: Put all of the above scaling code into separate private functions
  // to clean this up!
  
  // Set the objects that need to be set

  if(efwdW.get()) {
    efwdW->initialize(eW);  // This will directly update W_op if W.get()==NULL!
  }
  
  if( W.get() ) {
    Thyra::initializeOp<double>(*W_factory_,fwdW,&*W);
    W->setOStream(this->getOStream());
  }

  if( W_op.get() ) {
    if( W_op.shares_resource(efwdW) ) {
      // W_op was already updated above since *efwdW is the same object as *W_op
    }
    else {
      rcp_dynamic_cast<EpetraLinearOp>(W_op,true)->initialize(eW);
    }
  }

  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to process output objects = "<<timer.totalElapsedTime()<<" sec\n";

  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\noutArgs on output =\n" << Teuchos::describe(outArgs,verbLevel);

  totalTimer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal evaluation time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::EpetraModelEvaluator::evalModel(...) ...\n";
  
}


void EpetraModelEvaluator::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<double>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  finalPoint_ = this->createInArgs();
  finalPoint_.setArgs(finalPoint);
  finalPointWasSolved_ = wasSolved;
}


// Public functions overridden from Teuchos::Describable


std::string EpetraModelEvaluator::description() const
{
  std::ostringstream oss;
  oss << "Thyra::EpetraModelEvaluator{";
  oss << "epetraModel=";
  if(epetraModel_.get())
    oss << "\'"<<epetraModel_->description()<<"\'";
  else
    oss << "NULL";
  oss << ",W_factory=";
  if(W_factory_.get())
    oss << "\'"<<W_factory_->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


} // namespace Thyra


//
// Non-member utility functions
//


Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation
Thyra::convert( const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation &mvOrientation )
{
  switch(mvOrientation) {
    case EpetraExt::ModelEvaluator::DERIV_MV_BY_COL :
      return ModelEvaluatorBase::DERIV_MV_BY_COL;
    case EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW :
      return ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DERIV_MV_BY_COL; // Should never be called!
}


EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation
Thyra::convert( const ModelEvaluatorBase::EDerivativeMultiVectorOrientation &mvOrientation )
{
  switch(mvOrientation) {
    case ModelEvaluatorBase::DERIV_MV_BY_COL :
      return EpetraExt::ModelEvaluator::DERIV_MV_BY_COL;
    case ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW :
      return EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return EpetraExt::ModelEvaluator::DERIV_MV_BY_COL; // Should never be called!
}


Thyra::ModelEvaluatorBase::DerivativeProperties
Thyra::convert( const EpetraExt::ModelEvaluator::DerivativeProperties &derivativeProperties )
{
  ModelEvaluatorBase::EDerivativeLinearity linearity;
  switch(derivativeProperties.linearity) {
    case EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN;
      break;
    case EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_CONST;
      break;
    case  EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_NONCONST;
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
  ModelEvaluatorBase::ERankStatus rank;
  switch(derivativeProperties.rank) {
    case EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN:
      rank = ModelEvaluatorBase::DERIV_RANK_UNKNOWN;
      break;
    case EpetraExt::ModelEvaluator::DERIV_RANK_FULL:
      rank = ModelEvaluatorBase::DERIV_RANK_FULL;
      break;
    case EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT:
      rank = ModelEvaluatorBase::DERIV_RANK_DEFICIENT;
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DerivativeProperties(linearity,rank,derivativeProperties.supportsAdjoint);
}


Thyra::ModelEvaluatorBase::DerivativeSupport
Thyra::convert( const EpetraExt::ModelEvaluator::DerivativeSupport &derivativeSupport )
{
  ModelEvaluatorBase::DerivativeSupport ds;
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_LINEAR_OP))
    ds.plus(ModelEvaluatorBase::DERIV_LINEAR_OP);
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_MV_BY_COL))
    ds.plus(ModelEvaluatorBase::DERIV_MV_BY_COL);
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW))
    ds.plus(ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW);
  return ds;
}


EpetraExt::ModelEvaluator::Derivative
Thyra::convert(
  const ModelEvaluatorBase::Derivative<double>        &derivative
  ,const Teuchos::RefCountPtr<const Epetra_Map>       &fnc_map
  ,const Teuchos::RefCountPtr<const Epetra_Map>       &var_map
  )
{
  typedef ModelEvaluatorBase MEB;
  if(derivative.getLinearOp().get()) {
    return EpetraExt::ModelEvaluator::Derivative(
      Teuchos::rcp_const_cast<Epetra_Operator>(
        Teuchos::dyn_cast<const EpetraLinearOp>(*derivative.getLinearOp()).epetra_op()
        )
      );
  }
  else if(derivative.getDerivativeMultiVector().getMultiVector().get()) {
    return EpetraExt::ModelEvaluator::Derivative(
      EpetraExt::ModelEvaluator::DerivativeMultiVector(
        get_Epetra_MultiVector(
          ( derivative.getDerivativeMultiVector().getOrientation() == MEB::DERIV_MV_BY_COL
            ? *fnc_map
            : *var_map
            )
          ,derivative.getDerivativeMultiVector().getMultiVector()
          )
        ,convert(derivative.getDerivativeMultiVector().getOrientation())
        )
      );
  }
  return EpetraExt::ModelEvaluator::Derivative();
}
