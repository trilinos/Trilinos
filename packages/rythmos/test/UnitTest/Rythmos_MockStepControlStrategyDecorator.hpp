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


#ifndef RYTHMOS_MOCK_STEP_CONTROL_STRATEGY_HPP
#define RYTHMOS_MOCK_STEP_CONTROL_STRATEGY_HPP


#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include <vector>
#include <utility>  // for pair

namespace Rythmos {


/** \brief A mock step control strategy decortator that can trigger failures to test step control strategies
 *
 */
template<class Scalar>
class MockStepControlStrategyDecorator : virtual public Rythmos::StepControlStrategyBase<Scalar>
{
public:

  void initialize(const Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> >& scs);

  //! If number of consecutive failures is set to -1, then it will fail every time
  void addNonlinearSolverFailureOnStep(int step_number, int number_of_consecutive_failures = 1);

  /** \name Overridden from StepControlStrategyBase */
  //@{

  void initialize(const StepperBase<Scalar>& stepper);

  void setRequestedStepSize(
      const StepperBase<Scalar>& stepper
      , const Scalar& stepSize
      , const StepSizeType& stepSizeType
      );

  void nextStepSize(
      const StepperBase<Scalar>& stepper
      , Scalar* stepSize
      , StepSizeType* stepSizeType
      , int* order 
      );

  void setCorrection(
      const StepperBase<Scalar>& stepper
      , const RCP<const Thyra::VectorBase<Scalar> >& soln
      , const RCP<const Thyra::VectorBase<Scalar> >& ee
      , int solveStatus
      ); 

  bool acceptStep(
      const StepperBase<Scalar>& stepper
      ,Scalar* LETValue
      );

  void completeStep(
      const StepperBase<Scalar>& stepper
      );

  AttemptedStepStatusFlag rejectStep(
      const StepperBase<Scalar>& stepper
      );

  StepControlStrategyState getCurrentState();

  int getMaxOrder() const;

  void setStepControlData(const StepperBase<Scalar>& stepper);

  bool supportsCloning() const;

  RCP<StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const;

  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>&);

  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  //@}

private:

  Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> > scs_;

  std::vector<std::pair<int,int> > failure_points_;

  int numberOfSteps_;
};

// *******************************
// Implementation
// *******************************

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::
initialize(const Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> >& scs)
{
  TEUCHOS_ASSERT(nonnull(scs));
  scs_ = scs;
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::
addNonlinearSolverFailureOnStep(int step_number, int number_of_consecutive_failures)
{
  failure_points_.push_back(std::pair<int,int>(step_number,number_of_consecutive_failures));
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::initialize(const StepperBase<Scalar>& stepper)
{
  numberOfSteps_ = 0;
  scs_->initialize(stepper);
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::setRequestedStepSize(
      const StepperBase<Scalar>& stepper
      , const Scalar& stepSize
      , const StepSizeType& stepSizeType
      )
{
  scs_->setRequestedStepSize(stepper,stepSize,stepSizeType);
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::nextStepSize(
      const StepperBase<Scalar>& stepper
      , Scalar* stepSize
      , StepSizeType* stepSizeType
      , int* order 
      )
{
  scs_->nextStepSize(stepper,stepSize,stepSizeType,order);
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::setCorrection(
      const StepperBase<Scalar>& stepper
      , const RCP<const Thyra::VectorBase<Scalar> >& soln
      , const RCP<const Thyra::VectorBase<Scalar> >& ee
      , int solveStatus
      )
{
  typedef std::vector<std::pair<int,int> >::iterator it;
  bool fail_this_step = false;
  for (it i = failure_points_.begin(); i != failure_points_.end(); ++i) {
    if ( (i->first == numberOfSteps_ + 1) && (i->second > 0) ) {
      fail_this_step = true;
      i->second -= 1;
    }
  }

  if (fail_this_step)
    return scs_->setCorrection(stepper,soln,ee,-1);

  scs_->setCorrection(stepper,soln,ee,solveStatus);
}

template<class Scalar>
bool MockStepControlStrategyDecorator<Scalar>::acceptStep(
      const StepperBase<Scalar>& stepper
      ,Scalar* LETValue
      )
{
  return scs_->acceptStep(stepper,LETValue);
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::completeStep(
      const StepperBase<Scalar>& stepper
      )
{
  ++numberOfSteps_;
  scs_->completeStep(stepper);
}

template<class Scalar>
AttemptedStepStatusFlag MockStepControlStrategyDecorator<Scalar>::rejectStep(
      const StepperBase<Scalar>& stepper
      )
{
  return scs_->rejectStep(stepper);
}

template<class Scalar>
StepControlStrategyState MockStepControlStrategyDecorator<Scalar>::getCurrentState()
{
  return scs_->getCurrentState();
}

template<class Scalar>
int MockStepControlStrategyDecorator<Scalar>::getMaxOrder() const
{
  return scs_->getMaxOrder();
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::setStepControlData(const StepperBase<Scalar>& stepper)
{
  scs_->setStepControlData(stepper);
}

template<class Scalar>
bool MockStepControlStrategyDecorator<Scalar>::supportsCloning() const
{
  return scs_->supportsCloning();
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
MockStepControlStrategyDecorator<Scalar>::cloneStepControlStrategyAlgorithm() const
{
  return scs_->cloneStepControlStrategyAlgorithm();
}

template<class Scalar>
void MockStepControlStrategyDecorator<Scalar>::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  scs_->setParameterList(p);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> MockStepControlStrategyDecorator<Scalar>::getNonconstParameterList()
{
  return scs_->getNonconstParameterList();
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> MockStepControlStrategyDecorator<Scalar>::unsetParameterList()
{
  return scs_->unsetParameterList();
}

// *****************
// nonmember ctor
// *****************

/** \brief Nonmember ctor for creating MockStepControlStrategyDecorator objects

 \relates MockStepControlStrategyDecorator
*/
template<class Scalar>
Teuchos::RCP<MockStepControlStrategyDecorator<Scalar> > 
mockStepControlStrategyDecorator()
{
  Teuchos::RCP<MockStepControlStrategyDecorator<Scalar> > mscs = 
    Teuchos::rcp(new MockStepControlStrategyDecorator<Scalar>);
  return mscs;
}

} // namespace Rythmos


#endif // RYTHMOS_STEP_CONTROL_STRATEGY_BASE_HPP
