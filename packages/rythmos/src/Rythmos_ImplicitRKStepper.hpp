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

#ifndef Rythmos_IMPLICIT_RK_STEPPER_H
#define Rythmos_IMPLICIT_RK_STEPPER_H

#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_ImplicitRKModelEvaluator.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_TestingTools.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief . */
template<class Scalar>
class ImplicitRKStepper : virtual public SolverAcceptingStepperBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \name Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  ImplicitRKStepper();
  
  /** \brief . */
  void initialize(
    const RCP<const Thyra::ModelEvaluator<Scalar> >  &model,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory
    );
  
  /** \brief . */
  void setInterpolator(RCP<InterpolatorBase<Scalar> > interpolator);
  
  /** \brief . */
  RCP<InterpolatorBase<Scalar> > unsetInterpolator();

  //@}

  /** \name Overridden from SolverAcceptingStepperBase */
  //@{

  /** \brief . */
  void setSolver(
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  /** \brief . */
  RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver();

  /** \brief . */
  RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{
 
  /** \brief Returns true. */
  bool supportsCloning() const;

  /** \brief . */
  RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

  /** \brief . */
  void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> > &model);
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getModel() const;

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Scalar takeStep(Scalar dt, StepSizeType flag);
  
  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;
  
  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );
  
  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
  
  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;
  
  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;
  
  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getValidParameters() const;
 
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  void describe(
    FancyOStream  &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  // ///////////////////////
  // Private date members

  bool isInitialized_;
  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > irk_W_factory_;
  RCP<ParameterList> paramList_;

  int numStages_; // Temp hack to select between first and second order!

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  RCP<Thyra::VectorBase<Scalar> > x_;
  RCP<Thyra::VectorBase<Scalar> > x_old_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_;

  TimeRange<Scalar> timeRange_;

  RCP<ImplicitRKModelEvaluator<Scalar> > irkModel_;
  RKButcherTableau<Scalar> irkButcherTableau_;

  // Cache
  RCP<Thyra::ProductVectorBase<Scalar> > x_stage_bar_;

  // //////////////////////////
  // Private member functions

  void initialize_();

};


/** \brief Nonmember constructor.
 *
 * \relates ImplicitRKStepper
 */
template<class Scalar>
RCP<ImplicitRKStepper<Scalar> >
implicitRKStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> >  &model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory
  )
{
  RCP<ImplicitRKStepper<Scalar> > stepper(new ImplicitRKStepper<Scalar>());
  stepper->initialize( model, solver, irk_W_factory );
  return stepper;
}


// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
ImplicitRKStepper<Scalar>::ImplicitRKStepper()
  : isInitialized_(false),
    numStages_(1)
{}


template<class Scalar>
void ImplicitRKStepper<Scalar>::initialize(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory
  )
{

  // ToDo: Validate input

  model_ = model;
  solver_ = solver;
  irk_W_factory_ = irk_W_factory;

}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setInterpolator(
  RCP<InterpolatorBase<Scalar> > interpolator
  )
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
RCP<InterpolatorBase<Scalar> >
ImplicitRKStepper<Scalar>::unsetInterpolator()
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void ImplicitRKStepper<Scalar>::setSolver(
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  TEST_FOR_EXCEPT(true);
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
bool ImplicitRKStepper<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepperBase<Scalar> >
ImplicitRKStepper<Scalar>::cloneStepperAlgorithm() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ImplicitRKStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  TEST_FOR_EXCEPT( is_null(model_) );

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::BackwardEulerStepper::setInitialCondition(...)",
    *x_init->space(), *model_->get_x_space() );
#endif

  x_ = x_init->clone_v();

  // x_dot

  x_dot_ = createMember(model_->get_x_space());

  RCP<const Thyra::VectorBase<Scalar> >
    x_dot_init = initialCondition.get_x_dot();

  if (!is_null(x_dot_init))
    assign(&*x_dot_,*x_dot_init);
  else
    assign(&*x_dot_,ST::zero());
  
  // t

  const Scalar t =
    (
      initialCondition.supports(MEB::IN_ARG_t)
      ? initialCondition.get_t()
      : ST::zero()
      );

  timeRange_ = timeRange(t,t);

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
  Teuchos::OSTab ostab(out,1,"BES::takeStep");
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
  V_V( &*x_old_, *x_ );
  Scalar current_dt = dt;
  Scalar t = timeRange_.upper();
  irkModel_->setTimeStepPoint( x_old_, t, current_dt );

  // B) Give the IRK ModelEvaluator to the nonlinear solver and solve the timestep equation

  // Set the guess for the stage derivatives to zero (unless we can think of
  // something better)
  V_S( &*x_stage_bar_, ST::zero() );

  // Solve timestep equation
  solver_->solve( &*x_stage_bar_ );

  // C) Complete the step ...
  
  // Combine the stage derivatives with the Butcher tableau "b" vector to obtain the solution at the final time.
  // x_{k+1} = x_k + dt*sum_{i}^{p}(b_i*x_stage_bar_[i])

  assembleIRKSolution( irkButcherTableau_.b(), current_dt, *x_old_, *x_stage_bar_,
    outArg(*x_)
    );

  // Update time range
  timeRange_ = timeRange(t,t+current_dt);

  return current_dt;

}


template<class Scalar>
const StepStatus<Scalar> ImplicitRKStepper<Scalar>::getStepStatus() const
{
  StepStatus<Scalar> stepStatus;
  TEST_FOR_EXCEPT(true);
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

  if (x_vec)
    x_vec->resize( time_vec.size() );
  if (xdot_vec)
    xdot_vec->resize( time_vec.size() );
  if (accuracy_vec)
    accuracy_vec->resize( time_vec.size() );

  // This is a temp hack!
  if (time_vec.size() == 1 && compareTimeValues(timeRange_.upper(),time_vec[0])==0 ) {
    if (x_vec) {
      (*x_vec)[0] = x_;
    }
    TEST_FOR_EXCEPT( 0 != xdot_vec ); // Can't handle xdot yet!

    return; // We are done!
  }

  TEST_FOR_EXCEPT(true); // ToDo: Implement the final version!

}


template<class Scalar>
void ImplicitRKStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
void ImplicitRKStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
int ImplicitRKStepper<Scalar>::getOrder() const
{
  TEST_FOR_EXCEPT(true);
  return -1;
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
  numStages_ = paramList_->get("Num Stages",numStages_); // Temp hack
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
    pl->set("Num Stages",1); // Temp hack
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
  Teuchos::OSTab tab(out);
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
  const Scalar one = ST::one();
  const Scalar zero = ST::zero();
  using Teuchos::rcp_dynamic_cast;

  TEST_FOR_EXCEPT(model_ == Teuchos::null);
  TEST_FOR_EXCEPT(solver_ == Teuchos::null);
  TEST_FOR_EXCEPT(irk_W_factory_ == Teuchos::null);

  if (is_null(x_)) {
    // If x has not been set then setInitialCondition(...) was not
    // called by the client so the model had better have the 
    // initial condition!
    this->setInitialCondition(model_->getNominalValues());
  }

  if (is_null(x_dot_)) {
    x_dot_ = createMember(model_->get_x_space());
    V_S(&*x_dot_,ScalarTraits<Scalar>::zero());
  }

  x_old_ = x_->clone_v();

  // Set up the butcher tableau

  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages_,numStages_);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages_);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages_);

  // For now just two simple versions!
  TEUCHOS_ASSERT( numStages_ == 1 || numStages_ == 2 );

  if (numStages_ == 1) {

    A(0,0) = one;

    b(0) = one;
    
    c(0) = one;

  }
  else if (numStages_ == 2) {

    const Scalar half = 0.5 * one;

    A(0,0) = half;   A(0,1) = -half;
    A(1,0) = half;   A(1,1) = half;

    b(0) = half;
    b(1) = half;
    
    c(0) = zero;
    c(1) = one;


  }
  else {
    TEST_FOR_EXCEPT(true); // Should never get here!
  }
  
  irkButcherTableau_ = RKButcherTableau<Scalar>(A,b,c);

  // Set up the IRK mdoel

  irkModel_ = implicitRKModelEvaluator(
    model_,basePoint_,irk_W_factory_,irkButcherTableau_);

  solver_->setModel(irkModel_);

  // Set up the vector of stage derivatives ...
  x_stage_bar_ = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
    createMember(irkModel_->get_x_space())
    );

  isInitialized_ = true;

}


} // namespace Rythmos

#endif //Rythmos_IMPLICIT_RK_STEPPER_H
