//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef Rythmos_IMPLICIT_RK_STEPPER_DEF_H
#define Rythmos_IMPLICIT_RK_STEPPER_DEF_H

#include "Rythmos_ImplicitRKStepper_decl.hpp"

#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_ImplicitRKModelEvaluator.hpp"
#include "Rythmos_DiagonalImplicitRKModelEvaluator.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_RKButcherTableauHelpers.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_TestingTools.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {

/** \brief Nonmember constructor.
 *
 * \relates ImplicitRKStepper
 */
template<class Scalar>
RCP<ImplicitRKStepper<Scalar> >
implicitRKStepper()
{
  RCP<ImplicitRKStepper<Scalar> > stepper(new ImplicitRKStepper<Scalar>());
  return stepper;
}

template<class Scalar>
RCP<ImplicitRKStepper<Scalar> >
implicitRKStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& irk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> >& irkbt
  )
{
  RCP<ImplicitRKStepper<Scalar> > stepper(new ImplicitRKStepper<Scalar>());

  stepper->setModel(model);
  stepper->setSolver(solver);
  stepper->set_W_factory(irk_W_factory);
  stepper->setRKButcherTableau(irkbt);
  return stepper;
}


// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
ImplicitRKStepper<Scalar>::ImplicitRKStepper()
{
  this->defaultInitializeAll_();
  irkButcherTableau_ = rKButcherTableau<Scalar>();
  numSteps_ = 0;
}

template<class Scalar>
void ImplicitRKStepper<Scalar>::defaultInitializeAll_()
{
  isInitialized_ = false;
  model_ = Teuchos::null;
  solver_ = Teuchos::null;
  irk_W_factory_ = Teuchos::null;
  paramList_ = Teuchos::null;
  //basePoint_;
  x_ = Teuchos::null;
  x_old_ = Teuchos::null;
  x_dot_ = Teuchos::null;
  //timeRange_;
  irkModel_ = Teuchos::null;
  irkButcherTableau_ = Teuchos::null;
  isDirk_ = false;
  numSteps_ = -1;
  haveInitialCondition_ = false;
  x_stage_bar_ = Teuchos::null;
}

template<class Scalar>
void ImplicitRKStepper<Scalar>::set_W_factory(
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory
  )
{
  TEUCHOS_ASSERT( !is_null(irk_W_factory) ); 
  irk_W_factory_ = irk_W_factory;
}

template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > ImplicitRKStepper<Scalar>::get_W_factory() const
{
  return irk_W_factory_;
}

// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void ImplicitRKStepper<Scalar>::setSolver(
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  TEST_FOR_EXCEPT(is_null(solver));
  solver_ = solver;
}


template<class Scalar>
RCP<Thyra::NonlinearSolverBase<Scalar> >
ImplicitRKStepper<Scalar>::getNonconstSolver()
{
  return solver_;
}


template<class Scalar>
RCP<const Thyra::NonlinearSolverBase<Scalar> >
ImplicitRKStepper<Scalar>::getSolver() const
{
  return solver_;
}


// Overridden from StepperBase
 

template<class Scalar>
bool ImplicitRKStepper<Scalar>::isImplicit() const
{
  return true;
}

template<class Scalar>
bool ImplicitRKStepper<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepperBase<Scalar> >
ImplicitRKStepper<Scalar>::cloneStepperAlgorithm() const
{
  // Just use the interface to clone the algorithm in a basically
  // uninitialized state
  RCP<ImplicitRKStepper<Scalar> >
    stepper = Teuchos::rcp(new ImplicitRKStepper<Scalar>());

  if (!is_null(model_)) {
    stepper->setModel(model_); // Shallow copy is okay!
  }

  if (!is_null(irkButcherTableau_)) {
    // 06/16/09 tscoffe:  should we clone the RKBT here?
    stepper->setRKButcherTableau(irkButcherTableau_);
  }

  if (!is_null(solver_)) {
    stepper->setSolver(solver_->cloneNonlinearSolver().assert_not_null());
  }

  if (!is_null(irk_W_factory_)) {
    // 06/16/09 tscoffe:  should we clone the W_factory here?
    stepper->set_W_factory(irk_W_factory_);
  }

  if (!is_null(paramList_)) {
    stepper->setParameterList(Teuchos::parameterList(*paramList_));
  }

  return stepper;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{
  TEST_FOR_EXCEPT(is_null(model));
  assertValidModel( *this, *model );
  model_ = model;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  this->setModel(model); // TODO 09/09/09 tscoffe:  use ConstNonconstObjectContainer!
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ImplicitRKStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
ImplicitRKStepper<Scalar>::getNonconstModel() 
{
  return Teuchos::null;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
#endif

  x_ = x_init->clone_v();

  // x_dot

  x_dot_ = createMember(x_->space());

  RCP<const Thyra::VectorBase<Scalar> >
    x_dot_init = initialCondition.get_x_dot();

  if (!is_null(x_dot_init))
    assign(x_dot_.ptr(),*x_dot_init);
  else
    assign(x_dot_.ptr(),ST::zero());
  
  // t

  const Scalar t =
    (
      initialCondition.supports(MEB::IN_ARG_t)
      ? initialCondition.get_t()
      : ST::zero()
      );

  timeRange_ = timeRange(t,t);

  // x_old
  x_old_ = x_->clone_v();

  haveInitialCondition_ = true;

}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
ImplicitRKStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}


template<class Scalar>
Scalar ImplicitRKStepper<Scalar>::takeStep(Scalar dt, StepSizeType stepSizeType)
{


  using Teuchos::as;
  using Teuchos::incrVerbLevel;
  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::NonlinearSolverBase<Scalar> NSB;
  typedef Teuchos::VerboseObjectTempState<NSB> VOTSNSB;

  RCP<FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"takeStep");
  VOTSNSB solver_outputTempState(solver_,out,incrVerbLevel(verbLevel,-1));

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<ImplicitRKStepper<Scalar> >::name()
      << "::takeStep("<<dt<<","<<toString(stepSizeType)<<") ...\n"; 
  }

  if (!isInitialized_) {
    initialize_();
  }

  TEST_FOR_EXCEPT( stepSizeType != STEP_TYPE_FIXED ); // ToDo: Handle variable case later

  // A) Set up the IRK ModelEvaluator so that it can represent the time step
  // equation to be solved.

  // Set irkModel_ with x_old_, t_old_, and dt
  V_V( x_old_.ptr(), *x_ );
  Scalar current_dt = dt;
  Scalar t = timeRange_.upper();

  // B) Solve the timestep equation

  // Set the guess for the stage derivatives to zero (unless we can think of
  // something better)
  V_S( Teuchos::rcp_dynamic_cast<Thyra::VectorBase<Scalar> >(x_stage_bar_).ptr(), ST::zero() );

  if (!isDirk_) { // General Implicit RK Case:
    RCP<ImplicitRKModelEvaluator<Scalar> > firkModel_ = 
      Teuchos::rcp_dynamic_cast<ImplicitRKModelEvaluator<Scalar> >(irkModel_,true);
    firkModel_->setTimeStepPoint( x_old_, t, current_dt );

    // Solve timestep equation
    solver_->solve( &*x_stage_bar_ );

  } else { // Diagonal Implicit RK Case:

    RCP<DiagonalImplicitRKModelEvaluator<Scalar> > dirkModel_ = 
      Teuchos::rcp_dynamic_cast<DiagonalImplicitRKModelEvaluator<Scalar> >(irkModel_,true);
    dirkModel_->setTimeStepPoint( x_old_, t, current_dt );
    int numStages = irkButcherTableau_->numStages();
    for (int stage=0 ; stage < numStages ; ++stage) {
      dirkModel_->setCurrentStage(stage);
      solver_->solve( &*(x_stage_bar_->getNonconstVectorBlock(stage)) );
      dirkModel_->setStageSolution( stage, *(x_stage_bar_->getVectorBlock(stage)) );
    }

  }

  // C) Complete the step ...
  
  // Combine the stage derivatives with the Butcher tableau "b" vector to obtain the solution at the final time.
  // x_{k+1} = x_k + dt*sum_{i}^{p}(b_i*x_stage_bar_[i])

  assembleIRKSolution( irkButcherTableau_->b(), current_dt, *x_old_, *x_stage_bar_,
    outArg(*x_)
    );

  // Update time range
  timeRange_ = timeRange(t,t+current_dt);
  numSteps_++;

  return current_dt;

}


template<class Scalar>
const StepStatus<Scalar> ImplicitRKStepper<Scalar>::getStepStatus() const
{
  StepStatus<Scalar> stepStatus;

  if (!isInitialized_) {
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
    stepStatus.message = "This stepper is uninitialized.";
    return stepStatus;
  } 
  if (numSteps_ > 0) {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED;
  } 
  else {
    stepStatus.stepStatus = STEP_STATUS_UNKNOWN;
  }
  stepStatus.stepSize = timeRange_.length();
  stepStatus.order = irkButcherTableau_->order();
  stepStatus.time = timeRange_.upper();
  stepStatus.solution = x_;
  return(stepStatus);
}


// Overridden from InterpolationBufferBase


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitRKStepper<Scalar>::get_x_space() const
{
  return ( !is_null(model_) ? model_->get_x_space() : Teuchos::null );
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::addPoints(
    const Array<Scalar>& time_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    )
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
TimeRange<Scalar> ImplicitRKStepper<Scalar>::getTimeRange() const
{
  return timeRange_;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
  ,Array<ScalarMag>* accuracy_vec) const
{
  using Teuchos::constOptInArg;
  using Teuchos::null;
  TEUCHOS_ASSERT(haveInitialCondition_);
  defaultGetPoints<Scalar>(
    timeRange_.lower(), constOptInArg(*x_old_),
    Ptr<const VectorBase<Scalar> >(null), // Sun
    timeRange_.upper(), constOptInArg(*x_),
    Ptr<const VectorBase<Scalar> >(null), // Sun
    time_vec,
    ptr(x_vec), ptr(xdot_vec), ptr(accuracy_vec),
    Ptr<InterpolatorBase<Scalar> >(null) // For Sun
    );
  // 04/17/09 tscoffe:  Currently, we don't have x_dot to pass out (TODO)
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (!haveInitialCondition_) {
    return;
  }
  time_vec->push_back(timeRange_.lower());
  if (numSteps_ > 0) {
    time_vec->push_back(timeRange_.upper());
  }
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
int ImplicitRKStepper<Scalar>::getOrder() const
{
  return irkButcherTableau_->order();
}


// Overridden from Teuchos::ParameterListAcceptor


template <class Scalar>
void ImplicitRKStepper<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  paramList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}


template <class Scalar>
RCP<ParameterList>
ImplicitRKStepper<Scalar>::getNonconstParameterList()
{
  return(paramList_);
}


template <class Scalar>
RCP<ParameterList>
ImplicitRKStepper<Scalar>::unsetParameterList()
{
  RCP<ParameterList>
    temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
RCP<const ParameterList>
ImplicitRKStepper<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from Teuchos::Describable


template<class Scalar>
void ImplicitRKStepper<Scalar>::describe(
  FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  using Teuchos::as;
  if (!isInitialized_) {
    out << this->description() << " : This stepper is not initialized yet" << std::endl;
    return;
  }
  if (
    as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT)
    ||
    as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)
    )
  {
    out << this->description() << ":" << endl;
    Teuchos::OSTab tab(out);
    out << "model = " << Teuchos::describe(*model_,verbLevel);
    out << "solver = " << Teuchos::describe(*solver_,verbLevel);
    out << "irk_W_factory = " << Teuchos::describe(*irk_W_factory_,verbLevel);
  }
}


// private


template <class Scalar>
void ImplicitRKStepper<Scalar>::initialize_()
{

  typedef ScalarTraits<Scalar> ST;
  using Teuchos::rcp_dynamic_cast;

  TEST_FOR_EXCEPT(is_null(model_));
  TEST_FOR_EXCEPT(is_null(solver_));
  TEST_FOR_EXCEPT(irkButcherTableau_->numStages() == 0);
  TEUCHOS_ASSERT(haveInitialCondition_);

#ifdef RYTHMOS_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::ImplicitRKStepper::initialize_(...)",
    *x_->space(), *model_->get_x_space() );
#endif


  // Set up the IRK mdoel

  if (!isDirk_) { // General Implicit RK 
    TEST_FOR_EXCEPT(is_null(irk_W_factory_));
    irkModel_ = implicitRKModelEvaluator(
      model_,basePoint_,irk_W_factory_,irkButcherTableau_);
  } else { // Diagonal Implicit RK
    irkModel_ = diagonalImplicitRKModelEvaluator(
      model_,basePoint_,irk_W_factory_,irkButcherTableau_);
  }

  solver_->setModel(irkModel_);

  // Set up the vector of stage derivatives ...
  const int numStages = irkButcherTableau_->numStages();
  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > pvs = productVectorSpace(model_->get_x_space(),numStages);
  RCP<const Thyra::VectorSpaceBase<Scalar> > vs = rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(pvs,true);
  x_stage_bar_ = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(createMember(vs),true);
//  x_stage_bar_ = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
//    createMember(irkModel_->get_x_space())
//    );

  isInitialized_ = true;

}

template <class Scalar>
void ImplicitRKStepper<Scalar>::setRKButcherTableau( const RCP<const RKButcherTableauBase<Scalar> > &rkButcherTableau )
{
  TEUCHOS_ASSERT( !is_null(rkButcherTableau) );
  TEST_FOR_EXCEPTION( isInitialized_, std::logic_error,
      "Error!  The RK Butcher Tableau cannot be changed after internal initialization!"
      );
  validateIRKButcherTableau(*rkButcherTableau);
  irkButcherTableau_ = rkButcherTableau;
  E_RKButcherTableauTypes rkType = determineRKBTType<Scalar>(*irkButcherTableau_);
  if (
         (rkType == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK) 
      || (rkType == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK) 
      || (irkButcherTableau_->numStages() == 1)
     ) 
  {
    isDirk_ = true;
  } 
}

template <class Scalar>
RCP<const RKButcherTableauBase<Scalar> > ImplicitRKStepper<Scalar>::getRKButcherTableau() const
{
  return irkButcherTableau_;
}

template<class Scalar>
void ImplicitRKStepper<Scalar>::setDirk(bool isDirk)
{
  TEST_FOR_EXCEPTION(isInitialized_, std::logic_error,
      "Error!  Cannot change DIRK flag after internal initialization is completed\n"
      );
  if (isDirk == true) {
    E_RKButcherTableauTypes rkType = determineRKBTType<Scalar>(*irkButcherTableau_);
    bool RKBT_is_DIRK = (
         (rkType == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK) 
      || (rkType == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK) 
      || (irkButcherTableau_->numStages() == 1)
      );
    TEST_FOR_EXCEPTION( !RKBT_is_DIRK, std::logic_error,
        "Error!  Cannot set DIRK flag on a non-DIRK RK Butcher Tableau\n"
        );
  } else { // isDirk = false;
    isDirk_ = isDirk;
  }
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_IMPLICIT_RK_STEPPER_INSTANT(SCALAR) \
  \
  template class ImplicitRKStepper< SCALAR >; \
  \
  template RCP< ImplicitRKStepper< SCALAR > > \
  implicitRKStepper();  \
  \
  template RCP< ImplicitRKStepper< SCALAR > > \
  implicitRKStepper( \
    const RCP<const Thyra::ModelEvaluator< SCALAR > >& model, \
    const RCP<Thyra::NonlinearSolverBase< SCALAR > >& solver, \
    const RCP<Thyra::LinearOpWithSolveFactoryBase< SCALAR > >& irk_W_factory, \
    const RCP<const RKButcherTableauBase< SCALAR > >& irkbt \
      ); 

} // namespace Rythmos

#endif //Rythmos_IMPLICIT_RK_STEPPER_DEF_H
