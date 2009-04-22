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
//#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_UnitTestModels.hpp"

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

} // namespace Rythmos

