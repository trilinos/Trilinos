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

#ifdef HAVE_RYTHMOS_EXPERIMENTAL

#include "Rythmos_Types.hpp"

namespace RythmosCharon {

class CharonIntegrationControlAndObserver
  : virtual public Rythmos::IntegrationControlStrategyBase<double>,
    virtual public Rythmos::IntegrationObserverBase<double>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  CharonIntegrationControlAndObserver() { }
  virtual ~CharonIntegrationControlAndObserver() { }
  // Overridden from Rythmos::IntegrationControlStrategyBase
  RCP<IntegrationControlStrategyBase<double> > cloneIntegrationControlStrategy() const
  {
    return Teuchos::null;
  }
  void resetIntegrationControlStrategy(
    const TimeRange<double> &integrationTimeDomain
    )
  { }
  StepControlInfo<double> getNextStepControlInfo(
    const StepperBase<double> &stepper,
    const StepControlInfo<double> &stepCtrlInfoLast,
    const int timeStepIter
    )
  { }
  // Overridden from Rythmos::IntegrationObserverBase
  RCP<IntegrationObserverBase<double> > cloneIntegrationObserver() const
  {
    return Teuchos::null;
  }
  void resetIntegrationObserver(
    const TimeRange<double> &integrationTimeDomain
    )
  { }
  void observeCompletedTimeStep(
    const StepperBase<double> &stepper,
    const StepControlInfo<double> &stepCtrlInfo,
    const int timeStepIter
    )
  { }
  // Overridden from Teuchos::ParameterListAcceptorDefaultBase
  void setParameterList( RCP<Teuchos::ParameterList> const& paramList )
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class CharonImplicitBDFStepperErrWtVecCalc
  : virtual public Rythmos::ErrWtVecCalcBase<double>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  ImplicitBDFStepperErrWtVecCalc() { }
  virtual ~ImplicitBDFStepperErrWtVecCalc() { }
  void errWtVecSet(
      Thyra::VectorBase<double>* weight, 
      const Thyra::VectorBase<double>& vector, 
      double relTol, 
      double absTol
      ) const
  { }
  // Overridden from Teuchos::ParameterListAcceptorDefaultBase
  void setParameterList( RCP<Teuchos::ParameterList> const& paramList )
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class CharonImplicitBDFStepperStepControl
  : virtual public Rythmos::StepControlStrategyBase<double>,
    virtual public Teuchos:ParameterListAcceptorDefaultBase
{
  ImplicitBDFStepperStepControl() { }
  virtual ~ImplicitBDFStepperStepControl() { }
  void initialize(const StepperBase<double>& stepper)
  { }
  void setRequestedStepSize(
      const StepperBase<double>& stepper
      , const double& stepSize
      , const StepSizeType& stepSizeType
      )
  { }
  void nextStepSize(
      const StepperBase<double>& stepper
      , double* stepSize
      , StepSizeType* stepSizeType
      , int* order 
      )
  { }
  void setCorrection(
      const StepperBase<double>& stepper
      , const RCP<const Thyra::VectorBase<double> >& soln
      , const RCP<const Thyra::VectorBase<double> >& ee
      , int solveStatus
      )
  { }
  bool acceptStep(
      const StepperBase<double>& stepper
      ,double* LETValue
      )
  { 
    return false;
  }
  void completeStep(
      const StepperBase<double>& stepper
      )
  { }
  AttemptedStepStatusFlag rejectStep(
      const StepperBase<double>& stepper
      )
  { 
    return Rythmos::REP_ERR_FAIL;
  }
  StepControlStrategyState getCurrentState()
  {
    return Rythmos::UNINITIALIZED;
  }
  int getMaxOrder() const
  { 
    return 0;
  }
  void setStepControlData(const StepperBase<double>& stepper)
  { }
  bool supportsCloning() const
  { 
    return false;
  }
  RCP<StepControlStrategyBase<double> > cloneStepControlStrategyAlgorithm() const
  {
    return Teuchos::null;
  }
  // Overridden from Teuchos::ParameterListAcceptorDefaultBase
  void setParameterList( RCP<Teuchos::ParameterList> const& paramList )
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

} // namespace RythmosCharon

#endif // HAVE_RYTHMOS_EXPERIMENTAL

