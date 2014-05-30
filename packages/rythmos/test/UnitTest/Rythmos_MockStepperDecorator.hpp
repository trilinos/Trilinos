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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_MOCK_STEPPER_DECORATOR_HPP
#define Rythmos_MOCK_STEPPER_DECORATOR_HPP

#include "Rythmos_StepperBase.hpp"


namespace Rythmos {


/** \brief Wrap a given StepperBase object for unit testing purposes.
 */
template<class Scalar>
class MockStepperDecorator : virtual public StepperBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \name Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  MockStepperDecorator();
  
  /** \brief . */
  void initialize(
    const RCP<StepperBase<Scalar> >& stepper
    );

  /** \brief . */
  void setFailOnStepId(const int failOnStepId);

  /** \brief . */
  void setNumFailedSteps(const int numFailedSteps);

  //@}

  /** \name Overridden from StepperBase */
  //@{

  /** \brief . */
  bool isImplicit() const;

  /** \brief . */
  void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> >& model);

  /** \brief . */
  void setNonconstModel(const RCP<Thyra::ModelEvaluator<Scalar> >& model);
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> > getNonconstModel();

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getInitialCondition() const;

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
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;
 
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  void describe(
    Teuchos::FancyOStream  &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}


private:

  // ///////////////////////
  // Private date members

  RCP<StepperBase<Scalar> > stepper_;

  // User options
  int failOnStepId_;
  int numFailedSteps_;

  // Cache
  int currentStepId_;
  int numFailedStepCounter_;

  void resetTakeStepCache();

};


/** \brief Nonmember constructor.
 *
 * \relates MockStepperDecorator
 */
template<class Scalar>
RCP<MockStepperDecorator<Scalar> >
createMockStepperDecorator(
  const RCP<StepperBase<Scalar> >& stepper
  )
{
  const RCP<MockStepperDecorator<Scalar> > msd(new MockStepperDecorator<Scalar>);
  msd->initialize(stepper);
  return msd;
}


// //////////////////////////////////////////////////////
// Implementations


// Constructors, intializers, Misc.


template<class Scalar>
MockStepperDecorator<Scalar>::MockStepperDecorator()
  : failOnStepId_(-1),
    numFailedSteps_(1),
    currentStepId_(-1),
    numFailedStepCounter_(0)
{
  resetTakeStepCache();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::initialize(
  const RCP<StepperBase<Scalar> >& stepper
  )
{
  stepper_ = stepper;
  resetTakeStepCache();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::setFailOnStepId(const int failOnStepId)
{
  failOnStepId_ = failOnStepId;
}


template<class Scalar>
void MockStepperDecorator<Scalar>::setNumFailedSteps(
  const int numFailedSteps)
{
  numFailedSteps_ = numFailedSteps;
}


// Overridden from StepperBase


template<class Scalar>
bool MockStepperDecorator<Scalar>::isImplicit() const
{
  return stepper_->isImplicit();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  stepper_->setModel(model);
  resetTakeStepCache();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  stepper_->setNonconstModel(model);
  resetTakeStepCache();
}
  

template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
MockStepperDecorator<Scalar>::getModel() const
{
  return stepper_->getModel();
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
MockStepperDecorator<Scalar>::getNonconstModel()
{
  return stepper_->getNonconstModel();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  stepper_->setInitialCondition(initialCondition);
  resetTakeStepCache();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
MockStepperDecorator<Scalar>::getInitialCondition() const
{
  return stepper_->getInitialCondition();
}


template<class Scalar>
Scalar MockStepperDecorator<Scalar>::takeStep(Scalar dt, StepSizeType flag)
{
  if (currentStepId_ == failOnStepId_) {
    if (numFailedStepCounter_ < numFailedSteps_) {
      ++numFailedStepCounter_;
      return -ScalarTraits<Scalar>::one();
    }
  }
  ++currentStepId_;
  return stepper_->takeStep(dt, flag);
}
  

template<class Scalar>
const StepStatus<Scalar>
MockStepperDecorator<Scalar>::getStepStatus() const
{
  return stepper_->getStepStatus();
}
  

// Overridden from InterpolationBufferBase


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
MockStepperDecorator<Scalar>::get_x_space() const
{
  return stepper_->get_x_space();
}


template<class Scalar>
void MockStepperDecorator<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  return stepper_->addPoints(time_vec, x_vec, xdot_vec);
}
  

template<class Scalar>
TimeRange<Scalar>
MockStepperDecorator<Scalar>::getTimeRange() const
{
  return stepper_->getTimeRange();
}
  

template<class Scalar>
void MockStepperDecorator<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  return stepper_->getPoints(time_vec, x_vec, xdot_vec, accuracy_vec);
}


template<class Scalar>
void MockStepperDecorator<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  return stepper_->getNodes(time_vec);
}
  

template<class Scalar>
void MockStepperDecorator<Scalar>::removeNodes(Array<Scalar>& time_vec)
{
  return stepper_->removeNodes(time_vec);
}


template<class Scalar>
int MockStepperDecorator<Scalar>::getOrder() const
{
  return stepper_->getOrder();
}

  
// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void MockStepperDecorator<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList)
{
  return stepper_->setParameterList(paramList);
}
  

template<class Scalar>
RCP<Teuchos::ParameterList>
MockStepperDecorator<Scalar>::getNonconstParameterList()
{
  return stepper_->getNonconstParameterList();
}
  

template<class Scalar>
RCP<Teuchos::ParameterList>
MockStepperDecorator<Scalar>::unsetParameterList()
{
  return stepper_->unsetParameterList();
}
  

template<class Scalar>
RCP<const Teuchos::ParameterList>
MockStepperDecorator<Scalar>::getValidParameters() const
{
  return stepper_->getValidParameters();
}
 

// Overridden from Teuchos::Describable


template<class Scalar>
void MockStepperDecorator<Scalar>::describe(
  Teuchos::FancyOStream  &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  stepper_->describe(out, verbLevel);
}


// private


template<class Scalar>
void MockStepperDecorator<Scalar>::resetTakeStepCache()
{
  currentStepId_ = 0;
  numFailedStepCounter_ = 0;
}


} // namespace Rythmos


#endif //Rythmos_MOCK_STEPPER_DECORATOR_HPP
