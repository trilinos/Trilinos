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

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "../VanderPol/VanderPolModel.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_ImplicitBDFStepper, minOrder ) {
  // Model
  RCP<ParameterList> modelParamList = Teuchos::parameterList();
  sublist(modelParamList,Stratimikos_name);
  sublist(modelParamList,DiagonalTransientModel_name);
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(modelParamList);
  Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
  // Solver
  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();
  std::vector<StepSizeType> stepTypeVec(2);
  stepTypeVec[0] = STEP_TYPE_VARIABLE;
  stepTypeVec[1] = STEP_TYPE_FIXED;
  double step = -1.0;
  int maxOrder = 5;
  for (int minOrder=1 ; minOrder <= 5 ; ++minOrder ) {
    // parameter list
    RCP<ParameterList> stepperParamList = Teuchos::parameterList();
    ParameterList& pl = stepperParamList->sublist("Step Control Settings");
    pl.set("minOrder",minOrder);
    pl.set("maxOrder",maxOrder);
    ParameterList& vopl = pl.sublist("VerboseObject");
    vopl.set("Verbosity Level","none");
    for (int i=0 ; i<Teuchos::as<int>(stepTypeVec.size()) ; ++i) {
      StepSizeType stepType = stepTypeVec[i];

      // Stepper
      RCP<ImplicitBDFStepper<double> > stepper = rcp(new ImplicitBDFStepper<double>(model,nlSolver,stepperParamList));
      TEST_EQUALITY_CONST( Teuchos::is_null(stepper), false );
      stepper->setInitialCondition(model_ic);

      for (int order=1 ; order<=minOrder ; ++order) {
        step = stepper->takeStep(1.0,stepType);
        const StepStatus<double> status = stepper->getStepStatus();
        TEST_EQUALITY_CONST( status.order, order );
      }
      for (int steps=0 ; steps<4 ; ++steps) {
        step = stepper->takeStep(1.0,stepType);
        const StepStatus<double> status = stepper->getStepStatus();
        TEST_COMPARE( status.order, >=, minOrder );
        TEST_COMPARE( status.order, <=, maxOrder );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitBDFStepper, exactNumericalAnswer_BE ) {
  double a = 1.5;
  double f = 1.6;
  double L = 1.7;
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  {
    modelPL->set("Implicit model formulation",true);
    modelPL->set("Coeff a",a);
    modelPL->set("Coeff f",f);
    modelPL->set("Coeff L",L);
  }
  RCP<SinCosModel> model = sinCosModel();
  model->setParameterList(modelPL);
  Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  {
    ParameterList& pl = stepperPL->sublist("Step Control Settings");
    pl.set("minOrder",1);
    pl.set("maxOrder",1);
    ParameterList& vopl = pl.sublist("VerboseObject");
    vopl.set("Verbosity Level","none");
  }
  RCP<ImplicitBDFStepper<double> > stepper = implicitBDFStepper<double>(model,nlSolver,stepperPL);
  stepper->setInitialCondition(model_ic);
  // Take a few steps and compare the output.
  int N = 3;
  double dt = 0.1;
  double x_exact_0 = 0.0;
  double x_exact_1 = 1.0;
  for (int i=1 ; i<=N ; ++i) {
    double t = i*dt;
    double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
    TEST_ASSERT( dt_taken == dt );
    RCP<const VectorBase<double> > x;
    RCP<const VectorBase<double> > xdot;
    {
      // Get x out of stepper.
      Array<double> t_vec;
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      t_vec.resize(1); t_vec[0] = t;
      stepper->getPoints(t_vec,&x_vec,&xdot_vec,NULL);
      x = x_vec[0];
      xdot = xdot_vec[0];
    }
    // Compute exact solution:
    double c = dt/(1+f*f*dt*dt/(L*L));
    double x_e_0 = c*(x_exact_0/dt + x_exact_1 + dt*f*f/(L*L)*a);
    double x_e_1 = c*(-f*f/(L*L)*x_exact_0 + x_exact_1/dt + f*f/(L*L)*a);
    double xd_e_0 = (x_e_0-x_exact_0)/dt;
    double xd_e_1 = (x_e_1-x_exact_1)/dt;
    double tol = 1.0e-12;
    {
      Thyra::ConstDetachedVectorView<double> x_view( *x );
      TEST_FLOATING_EQUALITY( x_view[0], x_e_0, tol );
      TEST_FLOATING_EQUALITY( x_view[1], x_e_1, tol );

      Thyra::ConstDetachedVectorView<double> xdot_view( *xdot );
      TEST_FLOATING_EQUALITY( xdot_view[0], xd_e_0, tol );
      TEST_FLOATING_EQUALITY( xdot_view[1], xd_e_1, tol );
    }
    x_exact_0 = x_e_0;
    x_exact_1 = x_e_1;
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitBDFStepper, exactNumericalAnswer_BE_nonlinear ) {
  double epsilon = 0.5;
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  {
    modelPL->set("Implicit model formulation",true);
    modelPL->set("Coeff epsilon",epsilon);
  }
  RCP<VanderPolModel> model = vanderPolModel();
  model->setParameterList(modelPL);
  Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();
  {
    RCP<ParameterList> nlPL = Teuchos::parameterList();
    nlPL->set("Default Tol",1.0e-10);
    nlPL->set("Default Max Iters",20);
    nlSolver->setParameterList(nlPL);
  }
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  {
    ParameterList& pl = stepperPL->sublist("Step Control Settings");
    pl.set("minOrder",1);
    pl.set("maxOrder",1);
    ParameterList& vopl = pl.sublist("VerboseObject");
    vopl.set("Verbosity Level","none");
  }
  RCP<ImplicitBDFStepper<double> > stepper = implicitBDFStepper<double>(model,nlSolver,stepperPL);
  stepper->setInitialCondition(model_ic);
  double h = 0.1;
  std::vector<double> x_0_exact;
  std::vector<double> x_1_exact;
  std::vector<double> x_0_dot_exact;
  std::vector<double> x_1_dot_exact;
  {
    x_0_exact.push_back(2.0); // IC
    x_1_exact.push_back(0.0);

    x_0_exact.push_back(1.982896621392518e+00); // matlab 
    x_1_exact.push_back(-1.710337860748234e-01); 

    x_0_exact.push_back(1.951487185706842e+00); // matlab 
    x_1_exact.push_back(-3.140943568567556e-01); 
    
    x_0_exact.push_back(1.908249109758246e+00); // matlab 
    x_1_exact.push_back(-4.323807594859574e-01); 
    
    x_0_dot_exact.push_back(0.0);
    x_1_dot_exact.push_back(0.0);

    for ( int i=1 ; i< Teuchos::as<int>(x_0_exact.size()) ; ++i ) {
      x_0_dot_exact.push_back( (x_0_exact[i]-x_0_exact[i-1])/h );
      x_1_dot_exact.push_back( (x_1_exact[i]-x_1_exact[i-1])/h );
      //std::cout << "x_0_dot_exact["<<i<<"] = "<<x_0_dot_exact[i] << std::endl;
      //std::cout << "x_1_dot_exact["<<i<<"] = "<<x_1_dot_exact[i] << std::endl;
    }
  }
  double tol_discrete = 1.0e-12;
  double tol_continuous = 1.0e-2;
  {
    // Get IC out
    double t = 0.0;
    RCP<const VectorBase<double> > x;
    RCP<const VectorBase<double> > xdot;
    {
      // Get x out of stepper.
      Array<double> t_vec;
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      t_vec.resize(1); t_vec[0] = t;
      stepper->getPoints(t_vec,&x_vec,&xdot_vec,NULL);
      x = x_vec[0];
      xdot = xdot_vec[0];
    }
    {
      Thyra::ConstDetachedVectorView<double> x_view( *x );
      TEST_FLOATING_EQUALITY( x_view[0], x_0_exact[0], tol_discrete );
      TEST_FLOATING_EQUALITY( x_view[1], x_1_exact[0], tol_discrete );

      Thyra::ConstDetachedVectorView<double> xdot_view( *xdot );
      TEST_FLOATING_EQUALITY( xdot_view[0], x_0_dot_exact[0], tol_discrete );
      TEST_FLOATING_EQUALITY( xdot_view[1], x_1_dot_exact[0], tol_discrete );
    }
  }
  for (int i=1 ; i < Teuchos::as<int>(x_0_exact.size()); ++i) {
    // Each time step
    double t = 0.0+i*h;
    double h_taken = stepper->takeStep(h,STEP_TYPE_FIXED);
    TEST_ASSERT( h_taken == h );
    RCP<const VectorBase<double> > x;
    RCP<const VectorBase<double> > xdot;
    {
      // Get x out of stepper.
      Array<double> t_vec;
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      t_vec.resize(1); t_vec[0] = t;
      stepper->getPoints(t_vec,&x_vec,&xdot_vec,NULL);
      x = x_vec[0];
      xdot = xdot_vec[0];
    }
    {
      Thyra::ConstDetachedVectorView<double> x_view( *x );
      TEST_FLOATING_EQUALITY( x_view[0], x_0_exact[i], tol_discrete );
      TEST_FLOATING_EQUALITY( x_view[1], x_1_exact[i], tol_discrete );

      Thyra::ConstDetachedVectorView<double> xdot_view( *xdot );
      TEST_FLOATING_EQUALITY( xdot_view[0], x_0_dot_exact[i], tol_discrete );
      TEST_FLOATING_EQUALITY( xdot_view[1], x_1_dot_exact[i], tol_discrete );
    }
    // Now compare this to the continuous exact solution:
    {
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = model->getExactSolution(t);
      RCP<const VectorBase<double> > x_continuous_exact = inArgs.get_x();
      RCP<const VectorBase<double> > xdot_continuous_exact = inArgs.get_x_dot();
      {
        Thyra::ConstDetachedVectorView<double> x_view( *x );
        Thyra::ConstDetachedVectorView<double> xce_view( *x_continuous_exact );
        TEST_FLOATING_EQUALITY( x_view[0], xce_view[0], tol_continuous );
        TEST_FLOATING_EQUALITY( x_view[1], xce_view[1], tol_continuous*10 );

        Thyra::ConstDetachedVectorView<double> xdot_view( *xdot );
        Thyra::ConstDetachedVectorView<double> xdotce_view( *xdot_continuous_exact );
        TEST_FLOATING_EQUALITY( xdot_view[0], xdotce_view[0], tol_continuous*10 );
        TEST_FLOATING_EQUALITY( xdot_view[1], xdotce_view[1], tol_continuous*10 );
      }
    }
  }
}

} // namespace Rythmos

