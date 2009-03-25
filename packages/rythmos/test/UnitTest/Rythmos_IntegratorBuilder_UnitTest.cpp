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


#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Rythmos_IntegratorBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

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
    RCP<StepperBase<double> > stepper; return stepper; 
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
  void setModel(
    const RCP<const Thyra::ModelEvaluator<double> > &model
    )
  { }
  RCP<const Thyra::ModelEvaluator<double> > getModel() const
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
 

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, construct ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  TEST_EQUALITY_CONST( is_null(ib), false );
  TEST_NOTHROW( ib = Teuchos::null );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, constructFoolish ) {
  RCP<FoolishIntegrator> fi;
  TEST_NOTHROW( fi = rcp(new FoolishIntegrator()) );

  RCP<FoolishIntegrationControlStrategy> fics;
  TEST_NOTHROW( fics = rcp(new FoolishIntegrationControlStrategy()) );

  RCP<FoolishStepper> fs;
  TEST_NOTHROW( fs = rcp(new FoolishStepper()) );

  RCP<FoolishStepControlStrategy> fscs;
  TEST_NOTHROW( fscs = rcp(new FoolishStepControlStrategy()) );

  RCP<FoolishInterpolationBuffer> fib;
  TEST_NOTHROW( fib = rcp(new FoolishInterpolationBuffer()) );

  RCP<FoolishInterpolationBufferAppender> fiba;
  TEST_NOTHROW( fiba = rcp(new FoolishInterpolationBufferAppender()) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegratorFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setIntegratorFactory(
      Teuchos::abstractFactoryStd< IntegratorBase<double>, FoolishIntegrator >(),
      "Foolish Integrator"
      );
  ib->setIntegratorFactory(
      Teuchos::abstractFactoryStd< IntegratorBase<double>, DefaultIntegrator<double> >(),
      "Other Default Integrator"
      );
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integrator").set("Integrator Type","Foolish Integrator");
  pl->sublist("Stepper").set("Stepper Type","Explicit RK");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  TEST_NOTHROW( integrator = ib->create(model,ic,nlSolver) );
  RCP<FoolishIntegrator> fInt = Teuchos::rcp_dynamic_cast<FoolishIntegrator>(integrator, false);
  TEST_ASSERT( !is_null(fInt) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegratorFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setIntegratorFactory(
        Teuchos::abstractFactoryStd< IntegratorBase<double>, FoolishIntegrator >(),
        "Default Integrator"
        ),
      std::logic_error
      );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
      ib->setIntegratorFactory(
        Teuchos::abstractFactoryStd< IntegratorBase<double>, FoolishIntegrator >(),
        "Default Integrator"
        )
      );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegrationControlFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setIntegrationControlFactory(
      Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>, FoolishIntegrationControlStrategy >(),
      "Foolish Integration Control"
      );
  ib->setIntegrationControlFactory(
      Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>, SimpleIntegrationControlStrategy<double> >(),
      "Other Simple Integration Control"
      );
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integration Control").set("Integration Control Strategy Type","Foolish Integration Control");
  pl->sublist("Stepper").set("Stepper Type","Explicit RK");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  TEST_NOTHROW( integrator = ib->create(model,ic,nlSolver) );
  RCP<IntegrationControlStrategyAcceptingIntegratorBase<double> > specialInt = 
    Teuchos::rcp_dynamic_cast<IntegrationControlStrategyAcceptingIntegratorBase<double> >(integrator,false);
  TEST_ASSERT( !is_null(specialInt) );
  RCP<const IntegrationControlStrategyBase<double> > ibControl = specialInt->getIntegrationControlStrategy();
  TEST_ASSERT( !is_null(ibControl) );
  RCP<const FoolishIntegrationControlStrategy> specialIBControl = 
    Teuchos::rcp_dynamic_cast<const FoolishIntegrationControlStrategy>(ibControl,false);
  TEST_ASSERT( !is_null(specialIBControl) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegrationControlFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setIntegrationControlFactory(
        Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>, FoolishIntegrationControlStrategy >(),
        "Simple Integration Control Strategy"
        ),
      std::logic_error
      );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
      ib->setIntegrationControlFactory(
        Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>, FoolishIntegrationControlStrategy >(),
        "Simple Integration Control Strategy"
        )
      );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepperBuilder ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<StepperBuilder<double> > sb = stepperBuilder<double>();
  sb->setStepperFactory(
      Teuchos::abstractFactoryStd< StepperBase<double>, FoolishStepper >(),
      "Foolish Stepper"
      );
  sb->setStepperFactory(
      Teuchos::abstractFactoryStd< StepperBase<double>, BackwardEulerStepper<double> >(),
      "Other Backward Euler Stepper"
      );
  ib->setStepperBuilder(sb);

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper").set("Stepper Type","Foolish Stepper");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  TEST_ASSERT( !is_null(stepper) );
  RCP<const FoolishStepper> fStepper = Teuchos::rcp_dynamic_cast<const FoolishStepper>(stepper,false);
  TEST_ASSERT( !is_null(fStepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepControlFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

  ib->setStepControlFactory(
      Teuchos::abstractFactoryStd< StepControlStrategyBase<double>, FoolishStepControlStrategy >(),
      "Foolish Step Control"
      );
  ib->setStepControlFactory(
      Teuchos::abstractFactoryStd< StepControlStrategyBase<double>, ImplicitBDFStepperStepControl<double> >(),
      "Other Implicit BDF Step Control"
      );

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
  pl->sublist("Step Control").set("Step Control Strategy Type","Foolish Step Control");
  ib->setParameterList(pl);

  // Model:
  RCP<SinCosModel> model = sinCosModel(true);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  TEST_ASSERT( !is_null(stepper) );
  RCP<const ImplicitBDFStepper<double> > ibdfStepper = 
    Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper,false);
  TEST_ASSERT( !is_null(ibdfStepper) );
  RCP<const StepControlStrategyBase<double> > stepControl = ibdfStepper->getStepControlStrategy();
  TEST_ASSERT( !is_null(stepControl) );
  RCP<const FoolishStepControlStrategy> fStepControl = Teuchos::rcp_dynamic_cast<const FoolishStepControlStrategy>(stepControl,false);
  TEST_ASSERT( !is_null(fStepControl) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepControlFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setStepControlFactory(
        Teuchos::abstractFactoryStd< StepControlStrategyBase<double>, FoolishStepControlStrategy >(),
        "Implicit BDF Step Control Strategy"
        ),
      std::logic_error
      );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
      ib->setStepControlFactory(
        Teuchos::abstractFactoryStd< StepControlStrategyBase<double>, FoolishStepControlStrategy >(),
        "Implicit BDF Step Control Strategy"
        )
      );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolationBufferFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setInterpolationBufferFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferBase<double>, FoolishInterpolationBuffer >(),
      "Foolish InterpolationBuffer"
      );
  ib->setInterpolationBufferFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferBase<double>, InterpolationBuffer<double> >(),
      "Other InterpolationBuffer"
      );

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper").set("Stepper Type","Explicit RK");
  pl->sublist("Trailing Interpolation Buffer").set("Interpolation Buffer Type","Foolish InterpolationBuffer");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null 

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> > tibaIntegrator =
    Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<double> >(integrator,false);
  TEST_ASSERT( !is_null(tibaIntegrator) );
  RCP<const InterpolationBufferBase<double> > interpBuffer = tibaIntegrator->getTrailingInterpolationBuffer();
  TEST_ASSERT( !is_null(interpBuffer) );
  RCP<const FoolishInterpolationBuffer> fInterpBuffer = 
    Teuchos::rcp_dynamic_cast<const FoolishInterpolationBuffer>(interpBuffer,false);
  TEST_ASSERT( !is_null(fInterpBuffer) );
}
TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolationBufferFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setInterpolationBufferFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferBase<double>, FoolishInterpolationBuffer >(),
        "Interpolation Buffer"
        ),
      std::logic_error
      );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
      ib->setInterpolationBufferFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferBase<double>, FoolishInterpolationBuffer >(),
        "Interpolation Buffer"
        )
      );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolationBufferAppenderFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setInterpolationBufferAppenderFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>, FoolishInterpolationBufferAppender >(),
      "Foolish InterpolationBufferAppender"
      );
  ib->setInterpolationBufferAppenderFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>, PointwiseInterpolationBufferAppender<double> >(),
      "Other InterpolationBufferAppender"
      );

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper").set("Stepper Type","Explicit RK");
  pl->sublist("Interpolation Buffer Appender").set("Interpolation Buffer Appender Type","Foolish InterpolationBufferAppender");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null 

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<InterpolationBufferAppenderAcceptingIntegratorBase<double> > ibaaIntegrator = 
    Teuchos::rcp_dynamic_cast<InterpolationBufferAppenderAcceptingIntegratorBase<double> >(integrator,false);
  TEST_ASSERT( !is_null(ibaaIntegrator) );
  RCP<const InterpolationBufferAppenderBase<double> > interpBufferAppender = ibaaIntegrator->getInterpolationBufferAppender();
  TEST_ASSERT( !is_null(interpBufferAppender) );
  RCP<const FoolishInterpolationBufferAppender> fInterpBufferAppender = 
    Teuchos::rcp_dynamic_cast<const FoolishInterpolationBufferAppender>(interpBufferAppender,false);
  TEST_ASSERT( !is_null(fInterpBufferAppender) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolationBufferAppenderFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW(
      ib->setInterpolationBufferAppenderFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>, FoolishInterpolationBufferAppender >(),
        "Pointwise Interpolation Buffer Appender"
        ),
      std::logic_error
      );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
      ib->setInterpolationBufferAppenderFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>, FoolishInterpolationBufferAppender >(),
        "Pointwise Interpolation Buffer Appender"
        )
      );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}


TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ExplicitRK ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper").set("Stepper Type","Explicit RK");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_EQUALITY_CONST( is_null(integrator), false );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_EQUALITY_CONST( is_null(stepper_out), false );
    RCP<const ExplicitRKStepper<double> > erkStepper = Teuchos::rcp_dynamic_cast<const ExplicitRKStepper<double> >(stepper_out,false);
    TEST_EQUALITY_CONST( is_null(erkStepper), false );
    RCP<const Thyra::ModelEvaluator<double> > model_out = stepper_out->getModel();
    TEST_EQUALITY_CONST( is_null(model_out), false );
    RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_EQUALITY_CONST( is_null(sinCosModel), false );
  }
}

/*
TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ForwardEuler ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper").set("Stepper Type","Forward Euler");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_EQUALITY_CONST( is_null(integrator), false );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ForwardEulerStepper<double> > feStepper = Teuchos::rcp_dynamic_cast<const ForwardEulerStepper<double> >(stepper_out,false);
    TEST_ASSERT( !is_null(feStepper) );
  }
}
*/

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ImplicitRK ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper").set("Stepper Type","Implicit RK");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_EQUALITY_CONST( is_null(integrator), false );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_EQUALITY_CONST( is_null(stepper_out), false );
    RCP<const ImplicitRKStepper<double> > irkStepper = Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(stepper_out,false);
    TEST_EQUALITY_CONST( is_null(irkStepper), false );
    RCP<const Thyra::ModelEvaluator<double> > model_out = stepper_out->getModel();
    TEST_EQUALITY_CONST( is_null(model_out), false );
    RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_EQUALITY_CONST( is_null(sinCosModel), false );
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out = irkStepper->getSolver();
    TEST_EQUALITY_CONST( is_null(nlSolver_out), false );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_EQUALITY_CONST( is_null(tsnlSolver), false );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ImplicitBDF ) {
  {
    // Base test that all the pointers are of the right type
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_EQUALITY_CONST( is_null(integrator), false );
    
    // Stepper = Implicit BDF
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_EQUALITY_CONST( is_null(stepper_out), false );
    RCP<const ImplicitBDFStepper<double> > ibdfStepper = Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper_out,false);
    TEST_EQUALITY_CONST( is_null(ibdfStepper), false );
    // Model = SinCosModel
    RCP<const Thyra::ModelEvaluator<double> > model_out = stepper_out->getModel();
    TEST_EQUALITY_CONST( is_null(model_out), false );
    RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_EQUALITY_CONST( is_null(sinCosModel), false );
    // Solver = TimeStepNonlinearSolver
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out = ibdfStepper->getSolver();
    TEST_EQUALITY_CONST( is_null(nlSolver_out), false );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_EQUALITY_CONST( is_null(tsnlSolver), false );
    // Step Control = ImplicitBDFStepperStepControl
    RCP<const StepControlStrategyBase<double> > stepControl_out = ibdfStepper->getStepControlStrategy();
    TEST_EQUALITY_CONST( is_null(stepControl_out), false );
    RCP<const ImplicitBDFStepperStepControl<double> > ibdfStepControl = Teuchos::rcp_dynamic_cast<const ImplicitBDFStepperStepControl<double> >(stepControl_out, false);
    TEST_EQUALITY_CONST( is_null(ibdfStepControl), false );
    // Integrator = DefaultIntegrator
    RCP<DefaultIntegrator<double> > defInt = Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,false);
    TEST_EQUALITY_CONST( is_null(defInt), false );
    // Integration Control = SimpleIntegrationControlStrategy
    RCP<const IntegrationControlStrategyBase<double> > integrationControl_out = defInt->getIntegrationControlStrategy();
    TEST_EQUALITY_CONST( is_null(integrationControl_out), false );
    RCP<const SimpleIntegrationControlStrategy<double> > simpleIControl = Teuchos::rcp_dynamic_cast<const SimpleIntegrationControlStrategy<double> >(integrationControl_out, false);
    TEST_EQUALITY_CONST( is_null(simpleIControl), false );
    // TrailingInterpolationBuffer = InterpolationBuffer
    RCP<const InterpolationBufferBase<double> > interpolationBuffer = defInt->getTrailingInterpolationBuffer();
    TEST_EQUALITY_CONST( is_null(interpolationBuffer), false );
    RCP<const InterpolationBuffer<double> > defaultInterpolationBuffer = Teuchos::rcp_dynamic_cast<const InterpolationBuffer<double> >(interpolationBuffer, false);
    TEST_EQUALITY_CONST( is_null(defaultInterpolationBuffer), false );
  }
  {
    // Test that we can set Integration Control Strategy to None.
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
    pl->sublist("Integration Control").set("Integration Control Strategy Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_EQUALITY_CONST( is_null(integrator), false );
    // Integrator = DefaultIntegrator
    RCP<DefaultIntegrator<double> > defInt = Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,false);
    TEST_EQUALITY_CONST( is_null(defInt), false );
    // Integration Control = None
    RCP<const IntegrationControlStrategyBase<double> > integrationControl_out = defInt->getIntegrationControlStrategy();
    TEST_EQUALITY_CONST( is_null(integrationControl_out), true );
  }
  {
    // Test that we can set Step Control Strategy to None
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
    pl->sublist("Step Control").set("Step Control Strategy Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_EQUALITY_CONST( is_null(integrator), false );
    // Stepper = Implicit BDF
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_EQUALITY_CONST( is_null(stepper_out), false );
    RCP<const ImplicitBDFStepper<double> > ibdfStepper = Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper_out,false);
    TEST_EQUALITY_CONST( is_null(ibdfStepper), false );
    // Step Control = None
    RCP<const StepControlStrategyBase<double> > stepControl_out = ibdfStepper->getStepControlStrategy();
    TEST_EQUALITY_CONST( is_null(stepControl_out), true );
  }
  {
    // Test that we can set TrailingInterpolationBuffer to None
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
    pl->sublist("Trailing Interpolation Buffer").set("Interpolation Buffer Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_EQUALITY_CONST( is_null(integrator), false );
    // Trailing Interpolation Buffer = null
    RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> > tibaIntegrator = 
      Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<double> >(integrator,false);
    TEST_ASSERT( !is_null(tibaIntegrator) );
    RCP<const InterpolationBufferBase<double> > trailingIB = tibaIntegrator->getTrailingInterpolationBuffer();
    TEST_EQUALITY_CONST( is_null(trailingIB), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_BackwardEuler ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper").set("Stepper Type","Backward Euler");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_EQUALITY_CONST( is_null(integrator), false );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_EQUALITY_CONST( is_null(stepper_out), false );
    RCP<const BackwardEulerStepper<double> > beStepper = Teuchos::rcp_dynamic_cast<const BackwardEulerStepper<double> >(stepper_out,false);
    TEST_EQUALITY_CONST( is_null(beStepper), false );
    RCP<const Thyra::ModelEvaluator<double> > model_out = stepper_out->getModel();
    TEST_EQUALITY_CONST( is_null(model_out), false );
    RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_EQUALITY_CONST( is_null(sinCosModel), false );
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out = beStepper->getSolver();
    TEST_EQUALITY_CONST( is_null(nlSolver_out), false );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_EQUALITY_CONST( is_null(tsnlSolver), false );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_invalid ) {
  {
    // Integrator == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Integrator").set("Integrator Type","None");
    pl->sublist("Stepper").set("Stepper Type","Explicit RK");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Stepper == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper").set("Stepper Type","None");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Model == null
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<ParameterList> pl = Teuchos::parameterList();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model; 
    Thyra::ModelEvaluatorBase::InArgs<double> ic;
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Implicit Stepper and Solver == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper").set("Stepper Type","Implicit BDF");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setParameterList ) {
  // does it validate the list?
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Hello","World");
  TEST_THROW(ib->setParameterList(pl), std::logic_error);
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, getValidParameters ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<const ParameterList> pl = ib->getValidParameters();
  TEST_EQUALITY_CONST( is_null(pl), false );
  TEST_EQUALITY_CONST( 
      pl->isSublist("Integrator"), 
      true 
      );
  TEST_EQUALITY_CONST( 
      pl->sublist("Integrator").isParameter("Integrator Type"), 
      true 
      );
  TEST_EQUALITY_CONST( 
      pl->isSublist("Integration Control"), 
      true 
      );
  TEST_EQUALITY_CONST( 
      pl->sublist("Integration Control").isParameter("Integration Control Strategy Type"), 
      true 
      );
  TEST_EQUALITY_CONST( 
      pl->isSublist("Stepper"), 
      true 
      );
  TEST_EQUALITY_CONST( 
      pl->sublist("Stepper").isParameter("Stepper Type"), 
      true
      );
  TEST_EQUALITY_CONST( 
      pl->isSublist("Step Control"), 
      true 
      );
  TEST_EQUALITY_CONST(
      pl->sublist("Step Control").isParameter("Step Control Strategy Type"),
      true
      );
  TEST_EQUALITY_CONST(
      pl->sublist("Trailing Interpolation Buffer").isParameter("Interpolation Buffer Type"),
      true
      );
  TEST_EQUALITY_CONST(
      pl->sublist("Interpolation Buffer Appender").isParameter("Interpolation Buffer Appender Type"),
      true
      );
  TEST_EQUALITY_CONST(
      pl->isParameter("Final Time"),
      true
      );
  TEST_EQUALITY_CONST(
      pl->isParameter("Land On Final Time"),
      true
      );
}

} // namespace Rythmos 



