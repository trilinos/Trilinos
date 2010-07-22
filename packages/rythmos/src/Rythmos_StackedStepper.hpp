//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef RYTHMOS_STACKED_STEPPER_HPP
#define RYTHMOS_STACKED_STEPPER_HPP


#include "Rythmos_StepperBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

// Includes for ForwardSensitivityStackedStepper
#include "Rythmos_ForwardSensitivityModelEvaluatorBase.hpp"

namespace Rythmos {

template<class Scalar>
class StackedStepperStepStrategyBase
{
  public:
    virtual ~StackedStepperStepStrategyBase() {}
    virtual void setupNextStepper(
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        ) = 0;
    virtual Scalar evaluateStep(
        Scalar local_dt_taken, 
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        ) = 0;
};

template<class Scalar>
class DefaultStackedStepperStepStrategy
 : virtual public StackedStepperStepStrategyBase<Scalar>
{
  public:
    DefaultStackedStepperStepStrategy();
    virtual ~DefaultStackedStepperStepStrategy();
    void setupNextStepper(
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        );
    Scalar evaluateStep(
        Scalar local_dt_taken, 
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        );
  private:
    Scalar dt_taken_;
};

// Nonmember constructor declaration
template<class Scalar>
RCP<DefaultStackedStepperStepStrategy<Scalar> > 
defaultStackedStepperStepStrategy();

// Nonmember constructor definition
template<class Scalar>
RCP<DefaultStackedStepperStepStrategy<Scalar> > 
defaultStackedStepperStepStrategy()
{
  return Teuchos::rcp(new DefaultStackedStepperStepStrategy<Scalar>());
}

template<class Scalar>
DefaultStackedStepperStepStrategy<Scalar>::DefaultStackedStepperStepStrategy()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  dt_taken_ = ST::nan();
}


template<class Scalar>
DefaultStackedStepperStepStrategy<Scalar>::~DefaultStackedStepperStepStrategy()
{
}


template<class Scalar>
void DefaultStackedStepperStepStrategy<Scalar>::setupNextStepper(
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        )
{
  if (i > 0) {
    RCP<StepperBase<Scalar> > baseStepper = stepperArray[0];
    stepperArray[i]->setStepControlData(*baseStepper);
  }
}


template<class Scalar>
Scalar DefaultStackedStepperStepStrategy<Scalar>::evaluateStep(
        Scalar local_dt_taken, 
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        )
{
  if (i == 0) {
    dt_taken_ = local_dt_taken;
  }
  else {
    TEST_FOR_EXCEPTION( local_dt_taken != dt_taken_, std::logic_error,
        "Error!  sub-stepper["<<i<<"] did not take the same "
        "size step as the base stepper!"
        );
  }
  return dt_taken_;
}



template<class Scalar>
class ForwardSensitivityStackedStepperStepStrategy
 : virtual public StackedStepperStepStrategyBase<Scalar>
{
  public:
    ForwardSensitivityStackedStepperStepStrategy();
    virtual ~ForwardSensitivityStackedStepperStepStrategy();
    void setupNextStepper(
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        );
    Scalar evaluateStep(
        Scalar local_dt_taken, 
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        );
    void setStateBasePoint( 
        const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint
        );
  private:
    Scalar dt_taken_;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;
};

// Nonmember constructor declaration
template<class Scalar>
RCP<ForwardSensitivityStackedStepperStepStrategy<Scalar> >
forwardSensitivityStackedStepperStepStrategy();

// Nonmember constructor definition
template<class Scalar>
RCP<ForwardSensitivityStackedStepperStepStrategy<Scalar> >
forwardSensitivityStackedStepperStepStrategy()
{
  return Teuchos::rcp(new ForwardSensitivityStackedStepperStepStrategy<Scalar>());
}

template<class Scalar>
ForwardSensitivityStackedStepperStepStrategy<Scalar>::ForwardSensitivityStackedStepperStepStrategy()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  dt_taken_ = ST::nan();
}


template<class Scalar>
ForwardSensitivityStackedStepperStepStrategy<Scalar>::~ForwardSensitivityStackedStepperStepStrategy()
{
}


template<class Scalar>
void ForwardSensitivityStackedStepperStepStrategy<Scalar>::setupNextStepper(
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  if (i > 0) {
    RCP<StepperBase<Scalar> > baseStepper = stepperArray[0];
    RCP<Thyra::ModelEvaluator<Scalar> > modelI = 
      // TODO:  09/09/09 tscoffe:  This should be replaced by
      // getNonconstModel() when it is implemented correctly. 
      Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar> >(
          stepperArray[i]->getModel()
          );
    RCP<ForwardSensitivityModelEvaluatorBase<Scalar> > fwdSensME = 
      Teuchos::rcp_dynamic_cast<ForwardSensitivityModelEvaluatorBase<Scalar> > (
          modelI,false
          );
    if (Teuchos::nonnull(fwdSensME)) {
      bool forceUpToDateW = true;
      fwdSensME->initializePointState( Teuchos::inOutArg(*baseStepper), forceUpToDateW);
    }
    stepperArray[i]->setStepControlData(*baseStepper);
  }
}


template<class Scalar>
Scalar ForwardSensitivityStackedStepperStepStrategy<Scalar>::evaluateStep(
        Scalar local_dt_taken, 
        int i, 
        Array<RCP<StepperBase<Scalar> > > &stepperArray
        )
{
  if (i == 0) {
    dt_taken_ = local_dt_taken;
  }
  else {
    TEST_FOR_EXCEPTION( local_dt_taken != dt_taken_, std::logic_error,
        "Error!  sub-stepper["<<i<<"] did not take the same "
        "size step as the base stepper!"
        );
  }
  return dt_taken_;
}


template<class Scalar>
void ForwardSensitivityStackedStepperStepStrategy<Scalar>::setStateBasePoint( 
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint
    )
{
  stateBasePoint_ = stateBasePoint;
}


template<class Scalar> 
class StackedStepper
  : virtual public StepperBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Intializers, Misc. */
  //@{

  /** \brief Constructs to uninitialized. */
  StackedStepper();

  /** \brief Add a stepper to the stack
   */
  void addStepper(const RCP<StepperBase<Scalar> >& stepper);

  /** \brief Return the steppers that were passed into addStepper
   */
  ArrayView<const RCP<StepperBase<Scalar> > > getNonconstSteppers() const;

  /** \brief Add a Step Control Strategy
   */
  void setStackedStepperStepControlStrategy( 
      const RCP<StackedStepperStepStrategyBase<Scalar> >& stepControl
      );

  /** \brief Get the Step Control Strategy
   */
  RCP<const StackedStepperStepStrategyBase<Scalar> > 
    getStackedStepperStepCongrolStrategy() const;

  //@}

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{

  /** \brief Returns false. */
  bool acceptsModel() const;

  /** \brief Throws exception. */
  void setModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model
    );

  /** \brief Throws exception. */
  void setNonconstModel(
    const RCP<Thyra::ModelEvaluator<Scalar> >& model
    );

  /** \brief Returns <tt>getStateAndFwdSensModel()</tt>.
   *
   * Warning, currently the returned model does not implement evalModel(...) 
   * or define a W object.  It is just used for getting the spaces and for
   * creating an InArgs object for setting the initial condition.
   */
  RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> > getNonconstModel();

  /** \brief Sets the full initial condition for 
   * <tt>x_bar = [ x^{N} ]_{N=0..numSteppers} </tt>
   *
   * The InArgs object must be created using
   * <tt>this->getModel()->createInArgs()</tt> and then populated with the
   * initial values.  The product vectors for <tt>x_bar</tt> and
   * <tt>x_bar_dot</tt> can be created using
   * <tt>this->getModel()->create_x_bar_vec(...)</tt>.  All of
   * the input objects in <tt>state_ic</tt> will be cloned and
   * therefore no memory of the objects in <tt>state_ic</tt> will be
   * retained after calling this function.
   */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stacked_ic
    );

  /* \brief Get the initial condigion
   */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getInitialCondition() const;

  /** \brief . */
  Scalar takeStep( Scalar dt, StepSizeType stepType );

  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;

  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief Returns the space for <tt>x_bar</tt> and <tt>x_bar_dot</tt>.
   *
   * This space is a nested product vector space as described above.  Dynamic
   * casting is required to get at the <tt>ProductVectorSpaceBase</tt> and
   * <tt>ProductVectorBase</tt> intefaces.
   */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
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

private:

  mutable bool isInitialized_;

  Array<RCP<StepperBase<Scalar> > > stepperArray_;

  mutable Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > xSpaceArray_;
  //Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > fSpaceArray_;

  mutable Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > 
    x_bar_space_;
  //Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > 
  //  f_bar_space_;

  RCP<StackedStepperStepStrategyBase<Scalar> > stackedStepperStepStrategy_;
  
  // Private functions:
  void setupSpaces_() const;

};


/** \brief Nonmember constructor.
 *
 * \relates StackedStepper
 */
template<class Scalar> 
inline
RCP<StackedStepper<Scalar> >
stackedStepper()
{
  return Teuchos::rcp(new StackedStepper<Scalar>());
}


// Constructors, Intializers, Misc.


template<class Scalar> 
StackedStepper<Scalar>::StackedStepper()
  : isInitialized_(false)
{}

template<class Scalar> 
void StackedStepper<Scalar>::setupSpaces_() const
{
  using Teuchos::as;
  if (!isInitialized_) {
    for (int i=0 ; i < as<int>(stepperArray_.size()) ; ++i) {
      xSpaceArray_.push_back(stepperArray_[i]->get_x_space());
      //fSpaceArray_.push_back(stepperArray_[i]->get_f_space());
    }

    x_bar_space_ = Thyra::productVectorSpace<Scalar>(xSpaceArray_);
    //f_bar_space_ = Thyra::productVectorSpace<Scalar>(fSpaceArray_);

    isInitialized_ = true;
  }
}

template<class Scalar> 
void StackedStepper<Scalar>::addStepper(
    const RCP<StepperBase<Scalar> >& stepper
    )
{
  TEUCHOS_ASSERT(!is_null(stepper));
  stepperArray_.push_back(stepper);
  isInitialized_ = false;
}



template<class Scalar> 
ArrayView<const RCP<StepperBase<Scalar> > >
StackedStepper<Scalar>::getNonconstSteppers() const
{
  return stepperArray_();
}


template<class Scalar> 
void StackedStepper<Scalar>::setStackedStepperStepControlStrategy( 
    const RCP<StackedStepperStepStrategyBase<Scalar> >& stepControl
    )
{
  TEUCHOS_ASSERT( Teuchos::nonnull(stepControl) );
  stackedStepperStepStrategy_ = stepControl;
}


template<class Scalar> 
RCP<const StackedStepperStepStrategyBase<Scalar> > 
StackedStepper<Scalar>::getStackedStepperStepCongrolStrategy() const
{
  return stackedStepperStepStrategy_;
}


template<class Scalar> 
void StackedStepper<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const Teuchos::ParameterList>
StackedStepper<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from StepperBase

template<class Scalar> 
bool StackedStepper<Scalar>::acceptsModel() const
{
  return false;
}

template<class Scalar> 
void StackedStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{
  TEST_FOR_EXCEPT_MSG( true,
    "Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
}


template<class Scalar> 
void StackedStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  TEST_FOR_EXCEPT_MSG( true,
    "Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
}


template<class Scalar> 
RCP<const Thyra::ModelEvaluator<Scalar> >
StackedStepper<Scalar>::getModel() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar> 
RCP<Thyra::ModelEvaluator<Scalar> >
StackedStepper<Scalar>::getNonconstModel() 
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar> 
void StackedStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stacked_ic
  )
{
  TEST_FOR_EXCEPT(true);
}
 

template<class Scalar> 
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
StackedStepper<Scalar>::getInitialCondition() const
{
  TEST_FOR_EXCEPT(true);
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs;
  return inArgs;
}
 

template<class Scalar> 
Scalar
StackedStepper<Scalar>::takeStep(
  Scalar dt, StepSizeType stepType
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEST_FOR_EXCEPTION( stepperArray_.size() == 0, std::logic_error,
     "Error!  Rythmos::StackedStepper::takeStep: "
     "addStepper must be called at least once before takeStep!" 
     );
  this->setupSpaces_();

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:StackedStepper::takeStep");
#endif
  
  if (Teuchos::is_null(stackedStepperStepStrategy_)) {
    stackedStepperStepStrategy_ = defaultStackedStepperStepStrategy<Scalar>();
  }

  Scalar dt_taken = Scalar(-ST::one());
  for (int i=0 ; i<Teuchos::as<int>(stepperArray_.size()) ; ++i) {
    stackedStepperStepStrategy_->setupNextStepper(i,stepperArray_);
    Scalar local_dt_taken = stepperArray_[i]->takeStep(dt,stepType);
    dt_taken = stackedStepperStepStrategy_->evaluateStep(
          local_dt_taken,
          i,
          stepperArray_
          );
  }
  //// We cascade information from the first stepper to all the other steppers.     
  //RCP<StepperBase<Scalar> > baseStepper = stepperArray_[0];
  //Scalar dt_taken = baseStepper->takeStep(dt,stepType);
  //for (int i=1 ; i<Teuchos::as<int>(stepperArray_.size()) ; ++i) {
  //  // 07/27/09 tscoffe:  This line should be handled by a strategy object
  //  stepperArray_[i]->setStepControlData(*baseStepper);
  //  Scalar local_dt_taken = stepperArray_[i]->takeStep(dt,stepType);
  //  TEST_FOR_EXCEPTION( local_dt_taken != dt_taken, std::logic_error,
  //      "Error!  sub-stepper["<<i<<"] did not take the same "
  //      "size step as the base stepper!"
  //      );
  //}
  return dt_taken;
}


template<class Scalar> 
const StepStatus<Scalar>
StackedStepper<Scalar>::getStepStatus() const
{
  TEST_FOR_EXCEPT(true);
  const StepStatus<Scalar> stepStatus;
  return stepStatus;

}


// Overridden from InterpolationBufferBase


template<class Scalar> 
RCP<const Thyra::VectorSpaceBase<Scalar> >
StackedStepper<Scalar>::get_x_space() const
{
  this->setupSpaces_();
  TEST_FOR_EXCEPTION( is_null(x_bar_space_), std::logic_error,
      "Rythmos::StackedStepper::get_x_space():  "
      "addStepper must be called at least once before get_x_space()!"
      );
  return(x_bar_space_);
}


template<class Scalar> 
void StackedStepper<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  TEST_FOR_EXCEPT(
      "Not implemented addPoints(...) yet but we could if we wanted!"
      );
}


template<class Scalar> 
TimeRange<Scalar>
StackedStepper<Scalar>::getTimeRange() const
{
  TEST_FOR_EXCEPT(true);
  TimeRange<Scalar> tr;
  return tr;
}


template<class Scalar> 
void StackedStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_bar_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_bar_dot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  using Teuchos::as;
  this->setupSpaces_();
  TEST_FOR_EXCEPTION( stepperArray_.size() == 0, std::logic_error,
      "Rythmos::StackedStepper:getPoints:  Error!  "
      "addStepper must be called at least once before getPoints!"
      );
  int numSteppers = as<int>(stepperArray_.size());
  int numTimePoints = as<int>(time_vec.size());

  // Set up local x_bar_vec for filling:
  Array<RCP<Thyra::VectorBase<Scalar> > > x_bar_vec_local;
  if (x_bar_vec != NULL) {
    x_bar_vec_local.clear();
    x_bar_vec->clear();
    for (int n = 0 ; n < numTimePoints ; ++n) {
      x_bar_vec_local.push_back(Thyra::createMember(this->get_x_space()));
      x_bar_vec->push_back(x_bar_vec_local[n]);
    }
  }

  // Set up local x_bar_dot_vec for filling:
  Array<RCP<Thyra::VectorBase<Scalar> > > x_bar_dot_vec_local;
  if (x_bar_dot_vec != NULL) {
    x_bar_dot_vec_local.clear();
    x_bar_dot_vec->clear();
    for (int n = 0 ; n < numTimePoints ; ++n) {
      x_bar_dot_vec_local.push_back(Thyra::createMember(this->get_x_space()));
      x_bar_dot_vec->push_back(x_bar_dot_vec_local[n]);
    }
  }

  // Set up accuracy_vec
  if (accuracy_vec) {
    accuracy_vec->clear();
    accuracy_vec->resize(numTimePoints);
  }

  for (int i=0 ; i < numSteppers; ++i) {
    Array<RCP<const Thyra::VectorBase<Scalar> > > sub_x_bar_vec;
    Array<RCP<const Thyra::VectorBase<Scalar> > > sub_x_bar_dot_vec;
    Array<ScalarMag> sub_accuracy_vec;
    stepperArray_[i]->getPoints(
        time_vec,
        x_bar_vec ? &sub_x_bar_vec : 0,
        x_bar_dot_vec ? &sub_x_bar_dot_vec : 0,
        accuracy_vec ? &sub_accuracy_vec: 0
        );
    // Copy vectors into the sub-blocks of the 
    // x_bar_vec, x_bar_dot_vec, and accuracy_vec vectors.
    for (int n=0; n < numTimePoints ; ++n ) {
      if (x_bar_vec) {
        RCP<Thyra::ProductVectorBase<Scalar> > x_bar_pv = 
          Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
              x_bar_vec_local[n],
              true
              );
        RCP<Thyra::VectorBase<Scalar> > xi = 
          x_bar_pv->getNonconstVectorBlock(i);
        V_V(Teuchos::outArg(*xi),*sub_x_bar_vec[n]);
      }
      if (x_bar_dot_vec) {
        RCP<Thyra::ProductVectorBase<Scalar> > x_bar_dot_pv = 
          Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
              x_bar_dot_vec_local[n],
              true
              );
        RCP<Thyra::VectorBase<Scalar> > xdoti = 
          x_bar_dot_pv->getNonconstVectorBlock(i);
        V_V(Teuchos::outArg(*xdoti),*sub_x_bar_dot_vec[n]);
      }
      // Just copy the accuracy from the first stepper:
      if ( (i == 0) && accuracy_vec) {
          (*accuracy_vec)[n] = sub_accuracy_vec[n];
      }
    }
  }
}


template<class Scalar>
void StackedStepper<Scalar>::getNodes(
  Array<Scalar>* time_vec
  ) const
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar> 
void StackedStepper<Scalar>::removeNodes(
  Array<Scalar>& time_vec
  )
{
  TEST_FOR_EXCEPT("Not implemented yet but we can!");
}


template<class Scalar> 
int StackedStepper<Scalar>::getOrder() const
{
  TEST_FOR_EXCEPT(true);
  return -1;
}


// private


} // namespace Rythmos


#endif //RYTHMOS_STACKED_STEPPER_HPP
