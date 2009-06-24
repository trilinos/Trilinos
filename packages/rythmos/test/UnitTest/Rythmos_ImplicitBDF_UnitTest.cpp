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
    {
      // Get x out of stepper.
      Array<double> t_vec;
      Array<RCP<const VectorBase<double> > > x_vec;
      t_vec.resize(1); t_vec[0] = t;
      stepper->getPoints(t_vec,&x_vec,NULL,NULL);
      x = x_vec[0];
    }
    // Compute exact solution:
    double c = dt/(1+f*f*dt*dt/(L*L));
    double x_e_0 = c*(x_exact_0/dt + x_exact_1 + dt*f*f/(L*L)*a);
    double x_e_1 = c*(-f*f/(L*L)*x_exact_0 + x_exact_1/dt + f*f/(L*L)*a);
    double tol = 1.0e-12;
    {
      Thyra::ConstDetachedVectorView<double> x_view( *x );
      TEST_FLOATING_EQUALITY( x_view[0], x_e_0, tol );
      TEST_FLOATING_EQUALITY( x_view[1], x_e_1, tol );
    }
    x_exact_0 = x_e_0;
    x_exact_1 = x_e_1;
  }
}

} // namespace Rythmos

