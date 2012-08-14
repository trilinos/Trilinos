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
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "../UnitTest/Rythmos_IntegratorBuilder_Helpers.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "../UnitTest/Rythmos_UnitTestModels.hpp"

#ifdef Rythmos_ENABLE_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ForwardEuler ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ForwardEuler.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Forward Euler");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.0229174, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 1.0, 0.02 );
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, BackwardEuler ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("BackwardEuler.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Backward Euler");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.0207877, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 1.0, 0.02 );
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ImplicitBDF ) {

  int numOrders = 4;     // Should be 5, but currently there is a bug with 5.
  double error[4];
  error[0] = 0.0045542;
  error[1] = 7.36335e-05;
  error[2] = 1.33067e-06;
  error[3] = 2.58412e-08;
  for (int order=1; order<=numOrders; ++order)
  {
    // Compute Order of Accuracy
    Array<double> logStepSize;
    Array<double> logErrorNorm;
    //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

    std::stringstream fname;
    fname << "ImplicitBDF"<<order<<".dat";
    std::ofstream fout(fname.str().c_str());
    int numStepSizes = 6;
    for (int i=0; i<numStepSizes ; ++i)
    {
      double maxStepSize = 0.025/pow(2.0,i);
      double initStepSize = pow(maxStepSize, double(order))/2.0;

      RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

      // Set IntegratorBuilder ParameterList
      RCP<ParameterList> pl = Teuchos::parameterList();
      pl->setParameters(*(ib->getValidParameters()));
      //pl->sublist("Integrator Settings")
      //   .sublist("Integrator Selection")
      //   .sublist("Default Integrator")
      //   .sublist("VerboseObject")
      //       .set(    "Verbosity Level","medium");
      pl->sublist("Integrator Settings")
             .set(    "Final Time",0.5);
      pl->sublist("Stepper Settings")
         .sublist("Stepper Selection")
             .set(    "Stepper Type","Implicit BDF");
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
             .set(    "Step Control Strategy Type",
              "Implicit BDF Stepper Ramping Step Control Strategy");
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
         .sublist("Implicit BDF Stepper Ramping Step Control Strategy")
             .set(    "Initial Step Size", initStepSize);
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
         .sublist("Implicit BDF Stepper Ramping Step Control Strategy")
             .set(    "Max Step Size", maxStepSize);
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
         .sublist("Implicit BDF Stepper Ramping Step Control Strategy")
             .set(    "Max Order", order);
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
         .sublist("Implicit BDF Stepper Ramping Step Control Strategy")
             .set(    "Absolute Error Tolerance", double(1.0e-12));
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Step Control Strategy Selection")
         .sublist("Implicit BDF Stepper Ramping Step Control Strategy")
             .set(    "Relative Error Tolerance", double(1.0e-12));
      pl->sublist("Stepper Settings")
         .sublist("Step Control Settings")
         .sublist("Error Weight Vector Calculator Selection")
             .set(    "Error Weight Vector Calculator Type",
              "Implicit BDF Stepper Error Weight Vector Calculator");
      pl->sublist("Integration Control Strategy Selection")
             .set(    "Integration Control Strategy Type",
                      "Simple Integration Control Strategy");
      pl->sublist("Integration Control Strategy Selection")
         .sublist("Simple Integration Control Strategy")
             .set(    "Take Variable Steps",true);
      pl->sublist("Interpolation Buffer Settings")
         .sublist("Trailing Interpolation Buffer Selection")
             .set(    "Interpolation Buffer Type","Interpolation Buffer");
      ib->setParameterList(pl);

      double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                                .get<double>("Final Time");

      RCP<SinCosModel> model = sinCosModel(true);
      RCP<const VectorBase<double> >
        x_exact = model->getExactSolution(finalTime).get_x();
      Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
      RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
        timeStepNonlinearSolver<double>();

      // Finish Creation of IntegratorBuilder
      RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

      //std::ofstream fout("PL.txt");
      //ib->getParameterList()
      //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
      //                                                     .indent(4));
      //fout.close();

      Teuchos::Array<double> time_vec;
      Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
      time_vec.push_back(finalTime);

      // Perform Time Integration
      integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

      // Get final time from stepper to test if we completed integration.
      double t = integrator->getStepper()->getStepStatus().time;
      TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

      // Calculate the error
      RCP<VectorBase<double> > tmp_vec = 
        createMember(integrator->getStepper()->getModel()->get_x_space());
      Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
      logStepSize.push_back(log(maxStepSize));
      double nrm = Thyra::norm_inf(*tmp_vec);
      logErrorNorm.push_back(log(nrm));

      fout << maxStepSize << "       " << nrm << std::endl;
    }
    fout.close();

    if (logStepSize.size() > 1) 
    { 
      double slope = 
        computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

      // Test for shift in curve and change in slope.
      TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), error[order-1], 1.0e-5 );
      TEST_FLOATING_EQUALITY( slope, double(order), 0.05 );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_ForwardEuler ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_ForwardEuler.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type","Forward Euler");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.0229174, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 1.0, 0.02 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_4Stage ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_4Stage.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 3.8098e-07, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 4.0, 0.01 );
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_3_8_Rule ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_3_8_Rule.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type","Explicit 3/8 Rule");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 3.8098e-07, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 4.0, 0.01 );
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_4Stage3rdOrderRunge ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_4Stage3OrderRunge.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit 4 Stage 3rd order by Runge");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 1.93576e-05, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_3Stage3rdOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_3Stage3Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit 3 Stage 3rd order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 1.902e-05, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_3Stage3rdOrderHeun ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_3Stage3OrderHeun.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit 3 Stage 3rd order by Heun");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 1.902e-05, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_3Stage3rdOrderTVD ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_3Stage3OrderTVD.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit 3 Stage 3rd order TVD");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 1.902e-05, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_2Stage2ndOrderRunge ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_2Stage2OrderRunge.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit 2 Stage 2nd order by Runge");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.000758962, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 2.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, ERK_Trapezoidal ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("ERK_Trapezoidal.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Explicit Trapezoidal");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(false);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.000758962, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 2.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, IRK_BackwardEuler ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("IRK_BackwardEuler.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type","Backward Euler");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
#ifdef Rythmos_ENABLE_NOX
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      Teuchos::rcp(new Thyra::NOXNonlinearSolver);
    Teuchos::RCP<Teuchos::ParameterList> nlSolver_pl =
      Teuchos::rcp(new Teuchos::ParameterList);
    nlSolver_pl->sublist("Printing").set("Output Information", 0);
    nlSolver->setParameterList(nlSolver_pl);
#else
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
#endif // Rythmos_ENABLE_NOX

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.0207877, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 1.0, 0.02 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_2Stage2ndOrdergp5 ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_2Stage2Ordergp5.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 2 Stage 2nd order");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
       .sublist("Singly Diagonal IRK 2 Stage 2nd order")
           .set(    "gamma", 0.5);
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.0106911, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 1.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_2Stage2ndOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_2Stage2Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 2 Stage 2nd order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.000178193, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 2.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_2Stage3rdOrderAStable ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_2Stage3OrderAStable.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 2 Stage 3rd order");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
       .sublist("Singly Diagonal IRK 2 Stage 3rd order")
           .set(    "3rd Order A-stable", true);
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 3.65832e-05, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_2Stage3rdOrderLStable ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_2Stage3OrderLStable.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 2 Stage 3rd order");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
       .sublist("Singly Diagonal IRK 2 Stage 3rd order")
           .set(    "3rd Order A-stable", false);
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
       .sublist("Singly Diagonal IRK 2 Stage 3rd order")
           .set(    "2nd Order L-stable", true);
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.00482154, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 2.0, 0.04 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_5Stage5thOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_5Stage5Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 5 Stage 5th order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 5.81035e-09, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 4.58, 0.01 );  // 4.58 should be 5?
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_5Stage4thOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_5Stage4Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 5 Stage 4th order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 3.71905e-08, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 4.0, 0.01 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, SDIRK_3Stage4thOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("SDIRK_3Stage4Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Singly Diagonal IRK 3 Stage 4th order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 5.79059e-06, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 4.0, 0.02 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, DIRK_2Stage3rdOrder ) {

  // Compute Order of Accuracy
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  //std::cout << std::endl << "stepSize  errorNorm" << std::endl;

  std::ofstream fout("DIRK_2Stage3Order.dat");
  int numStepSizes = 6;
  for (int i=0; i<numStepSizes ; ++i)
  {
    double maxStepSize = 0.1/pow(2.0,i);

    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();

    // Set IntegratorBuilder ParameterList
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    //pl->sublist("Integrator Settings")
    //   .sublist("Integrator Selection")
    //   .sublist("Default Integrator")
    //   .sublist("VerboseObject")
    //       .set(    "Verbosity Level","medium");
    pl->sublist("Integrator Settings")
           .set(    "Final Time",0.5);
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
           .set(    "Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
           .set(    "Runge Kutta Butcher Tableau Type",
                    "Diagonal IRK 2 Stage 3rd order");
    pl->sublist("Integration Control Strategy Selection")
           .set(    "Integration Control Strategy Type",
                    "Simple Integration Control Strategy");
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Take Variable Steps",false);
    pl->sublist("Integration Control Strategy Selection")
       .sublist("Simple Integration Control Strategy")
           .set(    "Fixed dt", maxStepSize);
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
           .set(    "Interpolation Buffer Type","Interpolation Buffer");
    ib->setParameterList(pl);

    double finalTime = ib->getParameterList()->sublist("Integrator Settings")
                                              .get<double>("Final Time");

    RCP<SinCosModel> model = sinCosModel(true);
    RCP<const VectorBase<double> >
      x_exact = model->getExactSolution(finalTime).get_x();
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();

    // Finish Creation of IntegratorBuilder
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);

    //std::ofstream fout("PL.txt");
    //ib->getParameterList()
    //  ->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
    //                                                     .indent(4));
    //fout.close();

    Teuchos::Array<double> time_vec;
    Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
    time_vec.push_back(finalTime);

    // Perform Time Integration
    integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

    // Get final time from stepper to test if we completed integration.
    double t = integrator->getStepper()->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, time_vec[0], 1.0e-14 );

    // Calculate the error
    RCP<VectorBase<double> > tmp_vec = 
      createMember(integrator->getStepper()->getModel()->get_x_space());
    Thyra::V_StVpStV(tmp_vec.ptr(), 1.0, *x_exact, -1.0, *(x_vec[0]));
    logStepSize.push_back(log(maxStepSize));
    double nrm = Thyra::norm_inf(*tmp_vec);
    logErrorNorm.push_back(log(nrm));

    fout << maxStepSize << "       " << nrm << std::endl;
  }
  fout.close();

  if (logStepSize.size() > 1) 
  { 
    double slope = 
      computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);

    // Test for shift in curve and change in slope.
    TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 6.17632e-06, 1.0e-5 );
    TEST_FLOATING_EQUALITY( slope, 3.0, 0.01 );
  }
}

} // namespace Rythmos 
