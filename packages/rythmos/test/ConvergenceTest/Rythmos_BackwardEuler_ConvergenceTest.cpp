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


#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_BackwardEuler_ConvergenceTest.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, GlobalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(true);
  RCP<SinCosModelExactSolutionObject> exactSolution =
    sinCosModelExactSolutionObject(modelFactory);
  RCP<BackwardEulerStepperFactory<double> > stepperFactory =
    backwardEulerStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(
    stepperFactory,exactSolution);

  double slope =
    computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution);
  int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
  double tol = 1.0e-1;
  TEST_FLOATING_EQUALITY( slope, 1.0*order, tol );
}

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, LocalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(true);
  RCP<SinCosModelExactSolutionObject> exactSolution =
    sinCosModelExactSolutionObject(modelFactory);
  RCP<BackwardEulerStepperFactory<double> > stepperFactory =
    backwardEulerStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(
    stepperFactory,exactSolution);

  double slope =
    computeOrderByLocalErrorConvergenceStudy(stepperFactoryAndExactSolution);
  int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
  double tol = 1.0e-1;
  int localOrder = order+1;
  TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol );
}


TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, RampingIntCtrlVariableTimeStep ) {

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

  // Set IntegratorBuilder ParameterList
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));

  pl->sublist("Integrator Settings")
         .set(    "Final Time", 1.0e-8);
  pl->sublist("Integrator Settings")
     .sublist(    "Integrator Selection")
     .sublist(    "Default Integrator")
     .sublist(    "VerboseObject")
         .set(        "Verbosity Level","none");
  pl->sublist("Integrator Settings")
     .sublist(    "Integrator Selection")
     .sublist(    "Default Integrator")
         .set(        "Max Number Time Steps", 100000);

  pl->sublist("Integration Control Strategy Selection")
         .set(    "Integration Control Strategy Type",
                  "Ramping Integration Control Strategy");
  ParameterList& RICS =
    pl->sublist("Integration Control Strategy Selection", true)
       .sublist(    "Ramping Integration Control Strategy", true);
  RICS.set("Take Variable Steps",                false);
  RICS.set("Number of Initial Constant Steps",       0);
  RICS.set("Number of Ramping Steps",                0);
  RICS.set("Initial dt",                       1.0e-11);
  RICS.set("Min dt",                           1.0e-11);
  RICS.set("Max dt",                           1.0e-08);
  RICS.set("Ramping Factor",                       1.0);

  pl->sublist("Stepper Settings")
     .sublist(    "Stepper Selection")
         .set(        "Stepper Type","Backward Euler");

  pl->sublist("Stepper Settings")
     .sublist(    "Step Control Settings")
     .sublist(        "Step Control Strategy Selection")
         .set(            "Step Control Strategy Type",
                          "Fixed Step Control Strategy");

  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
         .set(    "Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);

  // Values from ParameterList
  double finalTime =
    ib->getParameterList()->sublist("Integrator Settings")
                           .get<double>(    "Final Time");
  //double maxStepSize =
  //  ib->getParameterList()
  //    ->sublist("Integration Control Strategy Selection")
  //     .sublist(    "Ramping Integration Control Strategy")
  //     .get<double>("Max dt");

  RCP<LogTimeModel> model = logTimeModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Finish Creation of IntegratorBuilder
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

  //std::ofstream ftmp("PL.txt");
  //ib->getParameterList()
  //  ->print(ftmp,Teuchos::ParameterList::PrintOptions().showDoc(true)
  //                                                     .indent(4));
  //ftmp.close();


  Teuchos::Array<double> time_vec;

  // Constant output times.
  int numSteps = 101;
  double stepSize = finalTime/(numSteps-1);
  for (int n=0; n<numSteps; n++) {
    time_vec.push_back(stepSize*n);
  }

  Teuchos::Array<RCP<const VectorBase<double> > > x_vec;

  // Perform Time Integration
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

  // Get final time from stepper to test if we completed integration.
  double t = integrator->getStepper()->getStepStatus().time;
  TEST_FLOATING_EQUALITY( t, time_vec[time_vec.size()-1], 1.0e-14 );
  if ( abs(t-time_vec[time_vec.size()-1]) <= 1.0e-14 ) {
    std::stringstream fname;
    fname << "BackwardEuler_RampingIntCtrl_var_dt.dat";
    std::ofstream fout(fname.str().c_str());

    double temporal_error = 0.0;
    double nrm = 0.0;
    for (int n=0; n<time_vec.size(); n++) {
      RCP<const VectorBase<double> > x_exact =
        model->getExactSolution(time_vec[n]).get_x();

      RCP<VectorBase<double> > diff_vec =
        createMember(integrator->getStepper()->getModel()->get_x_space());
      Thyra::V_StVpStV(diff_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[n]));
      nrm = Thyra::norm_2(*diff_vec);
      // Simple Trapeziodal Rule
      if ( n==0 )
        temporal_error += 0.5*nrm*nrm*(time_vec[n+1]-time_vec[n]);
      else if ( n==time_vec.size()-1 )
        temporal_error += 0.5*nrm*nrm*(time_vec[n]-time_vec[n-1]);
      else
        temporal_error += 0.5*nrm*nrm*(time_vec[n+1]-time_vec[n-1]);
      fout << time_vec[n]
           << "   " << get_ele(*(x_exact),0)
           << "   " << get_ele(*(x_vec[n]),0)
           << "   " << nrm
           << std::endl;
    }
    fout.close();
    temporal_error = sqrt(temporal_error);

    //cout << "\n     Step Size        = " << maxStepSize
    //     << "\n     Error at t_final = " << nrm
    //     << "\n     Temporal Error   = " << temporal_error
    //     << "\n" << std::endl;

    TEST_FLOATING_EQUALITY( nrm, 7.85747e-05, 1.0e-06 );
    TEST_FLOATING_EQUALITY( temporal_error, 1.37496e-07, 1.0e-05 );
  }
}


TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, StepCtrlVariableTimeStep ) {

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

  // Set IntegratorBuilder ParameterList
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));

  pl->sublist("Integrator Settings")
         .set(    "Final Time", 1.0e-8);
  pl->sublist("Integrator Settings")
     .sublist(    "Integrator Selection")
     .sublist(    "Default Integrator")
     .sublist(    "VerboseObject")
         .set(        "Verbosity Level","none");
  pl->sublist("Integrator Settings")
     .sublist(    "Integrator Selection")
     .sublist(    "Default Integrator")
         .set(        "Max Number Time Steps", 100000);

  pl->sublist("Integration Control Strategy Selection")
         .set(    "Integration Control Strategy Type",
                  "Simple Integration Control Strategy");
  ParameterList& SICS =
    pl->sublist("Integration Control Strategy Selection", true)
       .sublist(    "Simple Integration Control Strategy", true);
  SICS.set("Take Variable Steps",                 true);
  SICS.set("Max dt",                           1.0e-08);
  SICS.set("Number of Time Steps",                  -1);
  SICS.set("Fixed dt",                            -1.0);

  pl->sublist("Stepper Settings")
     .sublist(    "Stepper Selection")
         .set(        "Stepper Type","Backward Euler");
  pl->sublist("Stepper Settings")
     .sublist(    "Stepper Selection")
     .sublist(    "Backward Euler")
     .sublist(    "VerboseObject")
         .set(        "Verbosity Level","none");

  pl->sublist("Stepper Settings")
     .sublist(    "Step Control Settings")
     .sublist(        "Step Control Strategy Selection")
         .set(            "Step Control Strategy Type",
                          "Simple Step Control Strategy");
  ParameterList& SCS =
    pl->sublist("Stepper Settings")
       .sublist(    "Step Control Settings")
       .sublist(        "Step Control Strategy Selection")
       .sublist(            "Simple Step Control Strategy");
  SCS.set("Initial Step Size",                       1.0e-12);
  SCS.set("Min Step Size",                           1.0e-13);
  SCS.set("Max Step Size",                           1.0e-08);
  SCS.set("Step Size Increase Factor",                  1.07);
  SCS.set("Step Size Decrease Factor",                   0.5);
  SCS.set("Maximum Number of Step Failures",              10);
  SCS.set("Solution Change Relative Tolerance",      1.0e-03);
  SCS.set("Solution Change Absolute Tolerance",      1.0e-08);

  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
         .set(    "Interpolation Buffer Type","Interpolation Buffer");

  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .sublist(    "Interpolation Buffer")
         .set(    "StorageLimit", 10000);

  pl->sublist("Interpolation Buffer Settings")
     .sublist("Interpolation Buffer Appender Selection")
         .set(    "Interpolation Buffer Appender Type",
                  "Pointwise Interpolation Buffer Appender");

  ib->setParameterList(pl);

  // Values from ParameterList
  double finalTime =
    ib->getParameterList()->sublist("Integrator Settings")
                           .get<double>(    "Final Time");

  RCP<LogTimeModel> model = logTimeModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Finish Creation of IntegratorBuilder
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

  //std::ofstream ftmp("PL.txt");
  //ib->getParameterList()
  //  ->print(ftmp,Teuchos::ParameterList::PrintOptions().showDoc(true)
  //                                                     .indent(4));
  //ftmp.close();


  Teuchos::Array<double> time_vec;

  // Constant output times.
  int numSteps = 101;
  double stepSize = finalTime/(numSteps-1);
  for (int n=0; n<numSteps; n++) {
    time_vec.push_back(stepSize*n);
  }

  Teuchos::Array<RCP<const VectorBase<double> > > x_vec;

  // Perform Time Integration
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

  // Get final time from stepper to test if we completed integration.
  double t = integrator->getStepper()->getStepStatus().time;
  TEST_FLOATING_EQUALITY( t, time_vec[time_vec.size()-1], 1.0e-14 );
  if ( abs(t-time_vec[time_vec.size()-1]) <= 1.0e-14 ) {
    std::stringstream fname;
    fname << "BackwardEuler_SimpleStepCtrl_var_dt.dat";
    std::ofstream fout(fname.str().c_str());

    RCP<DefaultIntegrator<double> > di =
      Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,true);
    RCP<const InterpolationBufferBase<double> > buf =
     di->getTrailingInterpolationBuffer();
    Teuchos::Array<double> times;
    buf->getNodes(&times);
    Teuchos::Array<RCP<const VectorBase<double> > > x;
    buf->getPoints(times, &x, NULL, NULL);

    double temporal_error = 0.0;
    double nrm = 0.0;
    for (int n=0; n<times.size(); n++) {
      RCP<const VectorBase<double> > x_exact =
        model->getExactSolution(times[n]).get_x();

      RCP<VectorBase<double> > diff_vec =
        createMember(integrator->getStepper()->getModel()->get_x_space());
      Thyra::V_StVpStV(diff_vec.ptr(), 1.0, *x_exact, -1.0, *(x[n]));
      nrm = Thyra::norm_2(*diff_vec);
      // Simple Trapeziodal Rule
      if ( n==0 )
        temporal_error += 0.5*nrm*nrm*(times[n+1]-times[n]);
      else if ( n==times.size()-1 )
        temporal_error += 0.5*nrm*nrm*(times[n]-times[n-1]);
      else
        temporal_error += 0.5*nrm*nrm*(times[n+1]-times[n-1]);
      fout << times[n]
           << "   " << get_ele(*(x_exact),0)
           << "   " << get_ele(*(x[n]),0)
           << "   " << nrm
           << "   " << ( n>0 ? times[n]-times[n-1] : 0.0 )
           << std::endl;
    }
    fout.close();

    temporal_error = sqrt(temporal_error);

    //cout << "\n     Time             = " << t
    //     << "\n     Error at t_final = " << nrm
    //     << "\n     Temporal Error   = " << temporal_error
    //     << "\n" << std::endl;

    TEST_FLOATING_EQUALITY( nrm, 0.00543419, 5.0e-06 );

    TEST_FLOATING_EQUALITY( temporal_error, 6.13193e-07, 1.0e-05 );
  }
}


TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, FirstOrderErrorStepCtrl ) {

  for (int reln=0; reln<4; reln++) {
    double relerror = pow(10.0, -2-reln);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));

    pl->sublist("Integrator Settings")
           .set(    "Final Time", 1.0e-0);
    pl->sublist("Integrator Settings")
       .sublist(    "Integrator Selection")
       .sublist(    "Default Integrator")
       .sublist(    "VerboseObject")
           .set(        "Verbosity Level","none");
    pl->sublist("Integrator Settings")
       .sublist(    "Integrator Selection")
       .sublist(    "Default Integrator")
           .set(        "Max Number Time Steps", 100000);

    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    ParameterList& RICS =
      pl->sublist("Integration Control Strategy Selection", true)
         .sublist(    "Simple Integration Control Strategy", true);
    RICS.set("Take Variable Steps",                 true);
    RICS.set("Max dt",                           1.0e-00);

    pl->sublist("Stepper Settings")
       .sublist(    "Stepper Selection")
           .set(        "Stepper Type","Backward Euler");
    pl->sublist("Stepper Settings")
       .sublist(    "Stepper Selection")
       .sublist(    "Backward Euler")
       .sublist(    "VerboseObject")
           .set(        "Verbosity Level","none");

    pl->sublist("Stepper Settings")
       .sublist(    "Step Control Settings")
       .sublist(        "Step Control Strategy Selection")
           .set(            "Step Control Strategy Type",
                            "First Order Error Step Control Strategy");
    ParameterList& SCS =
      pl->sublist("Stepper Settings")
         .sublist(    "Step Control Settings")
         .sublist(        "Step Control Strategy Selection")
         .sublist(            "First Order Error Step Control Strategy");
    SCS.set("Initial Step Size",                       1.0e-12);
    SCS.set("Min Step Size",                           1.0e-15);
    SCS.set("Max Step Size",                           1.0e-00);
    SCS.set("Max Step Size Increase Factor",               2.0);
    SCS.set("Min Step Size Decrease Factor",               0.5);
    SCS.set("Maximum Number of Step Failures",              10);
    SCS.set("Error Relative Tolerance",               relerror);
    SCS.set("Error Absolute Tolerance",                1.0e-08);

    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");

    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
       .sublist(    "Interpolation Buffer")
           .set(    "StorageLimit", 10000);

    pl->sublist("Interpolation Buffer Settings")
       .sublist("Interpolation Buffer Appender Selection")
           .set(    "Interpolation Buffer Appender Type",
                    "Pointwise Interpolation Buffer Appender");

    ib->setParameterList(pl);

    // Values from ParameterList
    double finalTime =
      ib->getParameterList()->sublist("Integrator Settings")
                             .get<double>(    "Final Time");

    RCP<LogTimeModel> model = logTimeModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream ftmp("PL.txt");
    //ib->getParameterList()
    //  ->print(ftmp,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //ftmp.close();


    Teuchos::Array<double> time_vec;

    // Constant output times.
    int numSteps = 101;
    double stepSize = finalTime/(numSteps-1);
    for (int n=0; n<numSteps; n++) {
      time_vec.push_back(stepSize*n);
    }

    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[time_vec.size()-1], 1.0e-14 );
    if ( abs(t-time_vec[time_vec.size()-1]) <= 1.0e-14 ) {
      std::stringstream fname;
      if (reln == 0)
        fname << "BackwardEuler_FirstOrderError_var_dt_RelError=1.0e-2.dat";
      else if (reln == 1)
        fname << "BackwardEuler_FirstOrderError_var_dt_RelError=1.0e-3.dat";
      else if (reln == 2)
        fname << "BackwardEuler_FirstOrderError_var_dt_RelError=1.0e-4.dat";
      else if (reln == 3)
        fname << "BackwardEuler_FirstOrderError_var_dt_RelError=1.0e-5.dat";
      std::ofstream fout(fname.str().c_str());

      RCP<DefaultIntegrator<double> > di =
        Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,true);
      RCP<const InterpolationBufferBase<double> > buf =
       di->getTrailingInterpolationBuffer();
      Teuchos::Array<double> times;
      buf->getNodes(&times);
      Teuchos::Array<RCP<const VectorBase<double> > > x;
      buf->getPoints(times, &x, NULL, NULL);

      double temporal_error = 0.0;
      double nrm = 0.0;
      for (int n=0; n<times.size(); n++) {
        RCP<const VectorBase<double> > x_exact =
          model->getExactSolution(times[n]).get_x();

        RCP<VectorBase<double> > diff_vec =
          createMember(integrator->getStepper()->getModel()->get_x_space());
        Thyra::V_StVpStV(diff_vec.ptr(), 1.0, *x_exact, -1.0, *(x[n]));
        nrm = Thyra::norm_2(*diff_vec);
        // Simple Trapeziodal Rule
        if ( n==0 )
          temporal_error += 0.5*nrm*nrm*(times[n+1]-times[n]);
        else if ( n==times.size()-1 )
          temporal_error += 0.5*nrm*nrm*(times[n]-times[n-1]);
        else
          temporal_error += 0.5*nrm*nrm*(times[n+1]-times[n-1]);
        fout << times[n]
             << "   " << get_ele(*(x_exact),0)
             << "   " << get_ele(*(x[n]),0)
             << "   " << nrm
             << "   " << ( n>0 ? times[n]-times[n-1] : 0.0 )
             << std::endl;
      }
      fout.close();

      temporal_error = sqrt(temporal_error);

      //cout << "\n     Time             = " << t
      //     << "\n     Error at t_final = " << nrm
      //     << "\n     Temporal Error   = " << temporal_error
      //     << "\n" << std::endl;

      if (reln == 0) {
        TEST_FLOATING_EQUALITY( nrm, 0.0224576, 5.0e-06 );
        TEST_FLOATING_EQUALITY( temporal_error, 0.0224098, 1.0e-05 );
      } else if (reln == 1) {
        TEST_FLOATING_EQUALITY( nrm, 0.0132634, 5.0e-06 );
        TEST_FLOATING_EQUALITY( temporal_error, 0.0132322, 1.0e-05 );
      } else if (reln == 2) {
        TEST_FLOATING_EQUALITY( nrm, 0.00482358, 5.0e-06 );
        TEST_FLOATING_EQUALITY( temporal_error, 0.00480765, 1.0e-05 );
      } else if (reln == 3) {
        TEST_FLOATING_EQUALITY( nrm, 0.00154173, 5.0e-06 );
        TEST_FLOATING_EQUALITY( temporal_error, 0.00153523, 1.0e-05 );
      }
    }
  }
}


} // namespace Rythmos
