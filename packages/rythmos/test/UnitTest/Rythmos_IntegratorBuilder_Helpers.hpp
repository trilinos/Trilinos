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

#ifndef Rythmos_INTEGRATOR_BUILDER_HELPERS_H
#define Rythmos_INTEGRATOR_BUILDER_HELPERS_H


#include "Rythmos_IntegratorBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include "Rythmos_StepControlInfo.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Rythmos_ErrWtVecCalcBase.hpp"
#include "Rythmos_InterpolatorBase.hpp"

namespace Rythmos {

// Classes for testing the IntegratorBuilder
class FoolishIntegrator : 
  virtual public IntegratorBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  FoolishIntegrator() {}
  ~FoolishIntegrator() {}
  // Overridden from IntegratorBase:
  void setStepper(
    const RCP<StepperBase<double> > &stepper,
    const double &finalTime,
    const bool landOnFinalTime = true
    ) 
  { }
  Teuchos::RCP<const StepperBase<double> > getStepper() const
  { 
    return Teuchos::null;
  }
  Teuchos::RCP<StepperBase<double> > getNonconstStepper() const
  {
    return Teuchos::null;
  }
  RCP<StepperBase<double> > unSetStepper() 
  { 
    RCP<StepperBase<double> > stepper; return stepper; 
  }
  void getFwdPoints(
    const Array<double>& time_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* xdot_vec,
    Array<double>* accuracy_vec
    ) 
  { }
  TimeRange<double> getFwdTimeRange() const
  { 
    return TimeRange<double>(); 
  }
  // Overridden from InterpolationBufferBase:
  RCP<const Thyra::VectorSpaceBase<double> >
  get_x_space() const
  {
    RCP<Thyra::VectorSpaceBase<double> > space;
    return space;
  }
  void addPoints(
    const Array<double>& time_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& xdot_vec
    )
  { }
  TimeRange<double> getTimeRange() const
  {
    return TimeRange<double>();
  }
  void getPoints(
    const Array<double>& time_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* xdot_vec,
    Array<double>* accuracy_vec
    ) const
  { }
  void getNodes(Array<double>* time_vec) const
  { }
  void removeNodes(Array<double>& time_vec)
  { }
  int getOrder() const
  { 
    return 0; 
  }
  // Overridden from ParameterListAcceptor
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class FoolishIntegrationControlStrategy : 
  virtual public IntegrationControlStrategyBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  FoolishIntegrationControlStrategy() {}
  ~FoolishIntegrationControlStrategy() {}
  // Overridden from IntegrationControlStrategyBase:
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
  { 
    return StepControlInfo<double>();
  }
  // Overridden from ParameterListAcceptorDefaultBase
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class FoolishStepControlStrategy : 
  virtual public StepControlStrategyBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  FoolishStepControlStrategy() {}
  virtual ~FoolishStepControlStrategy() {}
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
    return CONTINUE_ANYWAY;
  }
  StepControlStrategyState getCurrentState()
  {
    return UNINITIALIZED;
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
    RCP<StepControlStrategyBase<double> > scsb;
    return scsb;
  }
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class FoolishStepper :
  virtual public StepperBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  FoolishStepper() {}
  virtual ~FoolishStepper() {}
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<double> &initialCondition
    )
  { }
  Thyra::ModelEvaluatorBase::InArgs<double> getInitialCondition() const
  { Thyra::ModelEvaluatorBase::InArgs<double> inArgs; return inArgs; }
  void setModel(
    const RCP<const Thyra::ModelEvaluator<double> >& model
    )
  { }
  void setNonconstModel(
    const RCP<Thyra::ModelEvaluator<double> >& model
    )
  { }
  RCP<const Thyra::ModelEvaluator<double> > getModel() const
  {
    return Teuchos::null;
  }
  RCP<Thyra::ModelEvaluator<double> > getNonconstModel() 
  {
    return Teuchos::null;
  }
  double takeStep(double dt, StepSizeType stepType)
  {
    return 0.0;
  }
  const StepStatus<double> getStepStatus() const
  {
    StepStatus<double> ss;
    return ss;
  }
  // Overridden from InterpolationBufferBase
  RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const
  {
    return Teuchos::null;
  }
  void addPoints(
    const Array<double>& time_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& xdot_vec
    )
  { }
  TimeRange<double> getTimeRange() const
  {
    return TimeRange<double>(0.0,0.0);
  }
  void getPoints(
    const Array<double>& time_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* xdot_vec,
    Array<double>* accuracy_vec
    ) const
  { }
  void getNodes(Array<double>* time_vec) const
  { }
  void removeNodes(Array<double>& time_vec)
  { }
  int getOrder() const
  {
    return 0;
  }
  // Overridden from ParameterListAcceptorDefaultBase
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};


class FoolishInterpolationBuffer : 
  virtual public InterpolationBufferBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  FoolishInterpolationBuffer() {}
  virtual ~FoolishInterpolationBuffer() {}
  RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const
  {
    return Teuchos::null;
  }
  void addPoints(
    const Array<double>& time_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<double> > >& xdot_vec
    )
  { }
  TimeRange<double> getTimeRange() const
  {
    return TimeRange<double>();
  }
  void getPoints(
    const Array<double>& time_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<double> > >* xdot_vec,
    Array<double>* accuracy_vec
    ) const
  { }
  void getNodes(Array<double>* time_vec) const
  { }
  void removeNodes(Array<double>& time_vec)
  { }
  int getOrder() const
  {
    return 0;
  }
  // Overriden from ParameterListAcceptor
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class FoolishInterpolationBufferAppender
: virtual public InterpolationBufferAppenderBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  void append(
    const InterpolationBufferBase<double>& interpBuffSource,
    const TimeRange<double>& range,
    const Ptr<InterpolationBufferBase<double> > &interpBuffSink
    )
  { }
  // Overriden from ParameterListAcceptor
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};
 
class FoolishErrWtVecCalc
: virtual public ErrWtVecCalcBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  void errWtVecSet(
       Thyra::VectorBase<double>* weight 
      ,const Thyra::VectorBase<double>& vector
      ,double relTol
      ,double absTol
      ) const
  { }
  // Overriden from ParameterListAcceptor
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

class FoolishInterpolator
: virtual public InterpolatorBase<double>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  void setNodes(
    const RCP<const DataStore<double>::DataStoreVector_t> & nodes
    )
  { }
  void interpolate(
    const Array<double> &t_values,
    DataStore<double>::DataStoreVector_t *data_out
    ) const
  { }
  int order() const
  { 
    return 0;
  }
  // Overriden from ParameterListAcceptor
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }

};

} // namespace Rythmos 

#endif // Rythmos_INTEGRATOR_BUILDER_HELPERS_H


