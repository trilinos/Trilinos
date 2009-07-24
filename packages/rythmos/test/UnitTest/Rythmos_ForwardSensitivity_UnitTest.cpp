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
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_UnitTestModels.hpp"

#include "Rythmos_StepperBuilder.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "../VanderPol/VanderPolModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_ForwardSensitivityStepperTester.hpp"

#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, create ) { 
  RCP<ForwardSensitivityStepper<double> > sensStepper =
    forwardSensitivityStepper<double>();
  TEST_EQUALITY_CONST( is_null(sensStepper), false );
}

// Finite Difference test for forward sensitivity calculation on
// SinCosModel with Backward Euler
TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, FDstepperTesterSinCosBE ) {

  RCP<ParameterList> ibPL = Teuchos::getParametersFromXmlString(
    "<ParameterList>"
    "  <ParameterList name=\"Stepper Settings\">"
    "    <ParameterList name=\"Stepper Selection\">"
    "      <Parameter name=\"Stepper Type\" type=\"string\" value=\"Backward Euler\"/>"
    "    </ParameterList>"
    "  </ParameterList>"
    "  <ParameterList name=\"Integration Control Strategy Selection\">"
    "    <Parameter name=\"Integration Control Strategy Type\" type=\"string\""
    "      value=\"Simple Integration Control Strategy\"/>"
    "    <ParameterList name=\"Simple Integration Control Strategy\">"
    "      <Parameter name=\"Take Variable Steps\" type=\"bool\" value=\"false\"/>"
    "      <Parameter name=\"Fixed dt\" type=\"double\" value=\"0.5\"/>"
    "    </ParameterList>"
    "  </ParameterList>"
    "</ParameterList>"
    );

  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  {
    modelPL->set("Accept model parameters",true);
    modelPL->set("Implicit model formulation",true);
    modelPL->set("Provide nominal values",true);
    modelPL->set("Coeff a", 2.0);
    modelPL->set("Coeff f", 3.0);
    modelPL->set("Coeff L", 4.0);
    modelPL->set("IC x_0", 2.0);
    modelPL->set("IC x_1", 5.0*3.0/4.0); 
  }
  stateModel->setParameterList(modelPL);
  int p_index = 0;
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();

  RCP<IntegratorBase<double> > fwdSensIntegrator = createForwardSensitivityIntegrator<double>(
      stateModel,p_index,state_ic,nonlinearSolver,ibPL);
  
  // Test the Forward Sensitivities

  RCP<ParameterList> fsstPL = Teuchos::getParametersFromXmlString(
    "<ParameterList>"
    "  <ParameterList name=\"FD Calc\">"
    "    <Parameter name=\"FD Method\" type=\"string\" value=\"order-four-central\"/>"
    "    <Parameter name=\"FD Step Length\" type=\"double\" value=\"1e-3\"/>"
    "    <Parameter name=\"FD Step Select Type\" type=\"string\" value=\"Relative\"/>"
    "  </ParameterList>"
    "  <Parameter name=\"Error Tol\" type=\"double\" value=\"1e-10\"/>"
    "</ParameterList>"
    );
  const RCP<ForwardSensitivityStepperTester<double> >
    fwdSensStepperTester = forwardSensitivityStepperTester<double>(fsstPL);

  //fwdSensStepperTester->setVerbLevel(Teuchos::VERB_EXTREME);
  fwdSensStepperTester->setVerbLevel(Teuchos::VERB_NONE);
  fwdSensStepperTester->setOStream(Teuchos::rcpFromRef(out));

  TEST_ASSERT(fwdSensStepperTester->testForwardSens(fwdSensIntegrator.ptr()));

}


// Finite Difference test for forward sensitivity calculation on
// VanderPol with Backward Euler
//TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, FDstepperTesterVanderPolBE ) {
//
//  RCP<ParameterList> ibPL = Teuchos::getParametersFromXmlString(
//    "<ParameterList>"
//    "  <ParameterList name=\"Stepper Settings\">"
//    "    <ParameterList name=\"Stepper Selection\">"
//    "      <Parameter name=\"Stepper Type\" type=\"string\" value=\"Backward Euler\"/>"
//    "    </ParameterList>"
//    "  </ParameterList>"
//    "  <ParameterList name=\"Integration Control Strategy Selection\">"
//    "    <Parameter name=\"Integration Control Strategy Type\" type=\"string\""
//    "      value=\"Simple Integration Control Strategy\"/>"
//    "    <ParameterList name=\"Simple Integration Control Strategy\">"
//    "      <Parameter name=\"Take Variable Steps\" type=\"bool\" value=\"false\"/>"
//    "      <Parameter name=\"Fixed dt\" type=\"double\" value=\"0.5\"/>"
//    "    </ParameterList>"
//    "  </ParameterList>"
//    "</ParameterList>"
//    );
//
//  RCP<VanderPolModel> stateModel = vanderPolModel();
//  RCP<ParameterList> modelPL = Teuchos::parameterList();
//  {
//    modelPL->set("Implicit model formulation",true);
//    modelPL->set("Accept model parameters",true);
//  }
//  stateModel->setParameterList(modelPL);
//  int p_index = 0;
//  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
//  typedef Thyra::ModelEvaluatorBase MEB;
//  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();
//
//  RCP<IntegratorBase<double> > fwdSensIntegrator = createForwardSensitivityIntegrator<double>(
//      stateModel,p_index,state_ic,nonlinearSolver,ibPL);
//  
//  // Test the Forward Sensitivities
//
//  RCP<ParameterList> fsstPL = Teuchos::getParametersFromXmlString(
//    "<ParameterList>"
//    "  <ParameterList name=\"FD Calc\">"
//    "    <Parameter name=\"FD Method\" type=\"string\" value=\"order-four-central\"/>"
//    "    <Parameter name=\"FD Step Length\" type=\"double\" value=\"1e-3\"/>"
//    "    <Parameter name=\"FD Step Select Type\" type=\"string\" value=\"Relative\"/>"
//    "  </ParameterList>"
//    "  <Parameter name=\"Error Tol\" type=\"double\" value=\"1e-10\"/>"
//    "</ParameterList>"
//    );
//  const RCP<ForwardSensitivityStepperTester<double> >
//    fwdSensStepperTester = forwardSensitivityStepperTester<double>(fsstPL);
//
//  fwdSensStepperTester->setVerbLevel(Teuchos::VERB_EXTREME);
//  //fwdSensStepperTester->setVerbLevel(Teuchos::VERB_NONE);
//  fwdSensStepperTester->setOStream(Teuchos::rcpFromRef(out));
//
//  TEST_ASSERT(fwdSensStepperTester->testForwardSens(fwdSensIntegrator.ptr()));
//
//}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, exactSinCosBE ) {
  //RCP<Teuchos::FancyOStream>
  //  std_out_rcp = Teuchos::VerboseObjectBase::getDefaultOStream();
  //Teuchos::FancyOStream& std_out = *std_out_rcp;
  Teuchos::FancyOStream& std_out = out; // TEUCHOS_UNIT_TEST defines "out"
  double b = 5.0;
  //double phi = 0.0;
  double a = 2.0;
  double f = 3.0;
  double L = 4.0;
  double x0 = a;
  double x1 = b*f/L;
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Accept model parameters",true);
  modelPL->set("Implicit model formulation",true);
  modelPL->set("Provide nominal values",true);
  modelPL->set("Coeff a", a);
  modelPL->set("Coeff f", f);
  modelPL->set("Coeff L", L);
  modelPL->set("IC x_0", x0);
  modelPL->set("IC x_1", x1);
  stateModel->setParameterList(modelPL);
  const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type","Backward Euler");
  builder->setParameterList(stepperPL);
  RCP<StepperBase<double> > stateStepper = builder->create();
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
  {
    // Set the nonlinear solver on the stepper.
    RCP<SolverAcceptingStepperBase<double> > SAStateStepper = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<double> >(stateStepper,true);
    SAStateStepper->setSolver(nonlinearSolver);
  }
  int p_index = 0;
  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper = 
    forwardSensitivityStepper<double>();
  stateAndSensStepper->initializeSyncedSteppers(
    stateModel, 
    p_index, 
    stateModel->getNominalValues(),
    stateStepper, 
    nonlinearSolver
    );

  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();

  RCP<Thyra::VectorBase<double> > s_bar_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  {
    // Initial conditions for Sensitivity problem.  (How do we get these correct in general?)
    RCP<Thyra::DefaultMultiVectorProductVector<double> > s_bar_mv =
      Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(s_bar_init,true);
    RCP<Thyra::VectorBase<double> > s_bar_mv_0 = s_bar_mv->getNonconstVectorBlock(0);
    Thyra::DetachedVectorView<double> s_bar_mv_0_view( *s_bar_mv_0 );
    s_bar_mv_0_view[0] = 1.0;
    s_bar_mv_0_view[1] = 0.0;
    RCP<Thyra::VectorBase<double> > s_bar_mv_1 = s_bar_mv->getNonconstVectorBlock(1);
    Thyra::DetachedVectorView<double> s_bar_mv_1_view( *s_bar_mv_1 );
    s_bar_mv_1_view[0] = 0.0;
    s_bar_mv_1_view[1] = b/L;
    RCP<Thyra::VectorBase<double> > s_bar_mv_2 = s_bar_mv->getNonconstVectorBlock(2);
    Thyra::DetachedVectorView<double> s_bar_mv_2_view( *s_bar_mv_2 );
    s_bar_mv_2_view[0] = 0.0;
    s_bar_mv_2_view[1] = -b*f/(L*L);
  }
  RCP<Thyra::VectorBase<double> > s_bar_dot_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  {
    RCP<Thyra::DefaultMultiVectorProductVector<double> > s_bar_dot_mv =
      Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(s_bar_dot_init,true);
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_0 = s_bar_dot_mv->getNonconstVectorBlock(0);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_0_view( *s_bar_dot_mv_0 );
    s_bar_dot_mv_0_view[0] = 0.0;
    s_bar_dot_mv_0_view[1] = 0.0;
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_1 = s_bar_dot_mv->getNonconstVectorBlock(1);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_1_view( *s_bar_dot_mv_1 );
    s_bar_dot_mv_1_view[0] = 0.0;
    s_bar_dot_mv_1_view[1] = -3.0*f*f*b/(L*L*L);
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_2 = s_bar_dot_mv->getNonconstVectorBlock(2);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_2_view( *s_bar_dot_mv_2 );
    s_bar_dot_mv_2_view[0] = 0.0;
    s_bar_dot_mv_2_view[1] = 3.0*f*f*f*b/(L*L*L*L);
  }

  RCP<const StateAndForwardSensitivityModelEvaluator<double> >
    stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

  MEB::InArgs<double>
    state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar = [ x; s_bar ]
  state_and_sens_ic.set_x(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
    );
  // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
  state_and_sens_ic.set_x_dot(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
    );

  stateAndSensStepper->setInitialCondition(state_and_sens_ic);

  double finalTime = 1.0e-4;
  int numTimeSteps = 2;
  RCP<IntegratorBase<double> > integrator;
  {
    RCP<ParameterList> integratorPL = Teuchos::parameterList();
    integratorPL->set( "Take Variable Steps", false );
    integratorPL->set( "Fixed dt", Teuchos::as<double>((finalTime - state_ic.get_t())/numTimeSteps) );
    RCP<IntegratorBase<double> >
      defaultIntegrator = controlledDefaultIntegrator<double>(
        simpleIntegrationControlStrategy<double>(integratorPL)
        );
    integrator = defaultIntegrator;
  }

  RCP<StepperAsModelEvaluator<double> >
    stateAndSensIntegratorAsModel = stepperAsModelEvaluator(
      Teuchos::rcp_implicit_cast<StepperBase<double> >(stateAndSensStepper),
      integrator, state_and_sens_ic
      );
  stateAndSensIntegratorAsModel->setVerbLevel(verbLevel);
    
  RCP<Thyra::VectorBase<double> > x_bar_final;
  {
  
    x_bar_final = createMember(stateAndSensIntegratorAsModel->get_g_space(0));
  
    Thyra::eval_g(
      *stateAndSensIntegratorAsModel,
      0, *state_ic.get_p(0),
      finalTime,
      0, &*x_bar_final
      );

    std_out
      << "\nx_bar_final = x_bar(p,finalTime) evaluated using stateAndSensIntegratorAsModel:\n"
      << Teuchos::describe(*x_bar_final,verbLevel);

  }

  // Now we check that the sensitivities are correct
  RCP<const Thyra::VectorBase<double> >
    DxDp_vec_final = Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(1);
  std_out << "\nDxDp_vec_final:\n"
    << Teuchos::describe(*DxDp_vec_final,verbLevel);
  RCP<const Thyra::DefaultMultiVectorProductVector<double> > DxDp_mv_final =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<double> >(DxDp_vec_final,true);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s0_final = DxDp_mv_final->getVectorBlock(0);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s1_final = DxDp_mv_final->getVectorBlock(1);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s2_final = DxDp_mv_final->getVectorBlock(2);
  

  std_out << "\nDxDp_s0_final:\n"
    << Teuchos::describe(*DxDp_s0_final,verbLevel);
  std_out << "\nDxDp_s1_final:\n"
    << Teuchos::describe(*DxDp_s1_final,verbLevel);
  std_out << "\nDxDp_s2_final:\n"
    << Teuchos::describe(*DxDp_s2_final,verbLevel);

  MEB::InArgs<double> exactSensSolution = stateModel->getExactSensSolution(0,finalTime);
  RCP<const Thyra::VectorBase<double> > ds0dp = exactSensSolution.get_x();
  exactSensSolution = stateModel->getExactSensSolution(1,finalTime);
  RCP<const Thyra::VectorBase<double> > ds1dp = exactSensSolution.get_x();
  exactSensSolution = stateModel->getExactSensSolution(2,finalTime);
  RCP<const Thyra::VectorBase<double> > ds2dp = exactSensSolution.get_x();

  std_out << "\nds0dp exact:\n"
    << Teuchos::describe(*ds0dp,verbLevel);
  std_out << "\nds1dp exact:\n"
    << Teuchos::describe(*ds1dp,verbLevel);
  std_out << "\nds2dp exact:\n"
    << Teuchos::describe(*ds2dp,verbLevel);

  /*
  // Compute finite difference Sensitivities:
  RCP<Thyra::MultiVectorBase<double> > DxDp_fd_final;
  {
    // Create (just) State integrator
    RCP<Rythmos::StepperAsModelEvaluator<double> >
      stateIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
        stateStepper, integrator, state_ic
        );
    // Create the finite difference calculator
    Thyra::DirectionalFiniteDiffCalculator<double> fdCalc;
    //fdCalc.setParameterList(sublist(paramList,FdCalc_name));
    //fdCalc.setOStream(out);
    //fdCalc.setVerbLevel(verbLevel);

    MEB::InArgs<double>
      fdBasePoint = stateIntegratorAsModel->createInArgs();
  
    fdBasePoint.set_t(finalTime);
    fdBasePoint.set_p(0,stateModel->getNominalValues().get_p(0));
  
    DxDp_fd_final = createMembers(
      stateIntegratorAsModel->get_g_space(0),
      stateIntegratorAsModel->get_p_space(0)->dim()
      );
  
    typedef Thyra::DirectionalFiniteDiffCalculatorTypes::SelectedDerivatives
      SelectedDerivatives;
  
    MEB::OutArgs<double> fdOutArgs =
      fdCalc.createOutArgs(
        *stateIntegratorAsModel,
        SelectedDerivatives().supports(MEB::OUT_ARG_DgDp,0,0)
        );
    fdOutArgs.set_DgDp(0,0,DxDp_fd_final);
  
    // Silence the model evaluators that are called.  The fdCal object
    // will show all of the inputs and outputs for each call.
    stateStepper->setVerbLevel(Teuchos::VERB_NONE);
    stateIntegratorAsModel->setVerbLevel(Teuchos::VERB_NONE);
  
    fdCalc.calcDerivatives(
      *stateIntegratorAsModel, fdBasePoint,
      stateIntegratorAsModel->createOutArgs(), // Don't bother with function value
      fdOutArgs
      );
    
    std_out
      << "\nFinite difference DxDp_fd_final = DxDp(p,finalTime): "
      << Teuchos::describe(*DxDp_fd_final,verbLevel);

  }
  RCP<const Thyra::VectorBase<double> >
    DxDp_fd_vec_final = Thyra::multiVectorProductVector(
      Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVectorSpace<double> >(
        DxDp_vec_final->range()
        ),
      DxDp_fd_final
      );
  
  verbLevel = Teuchos::VERB_EXTREME; // DEBUG

  double maxSensError = 1.0e-4;
  double s_fd_correct = Thyra::testRelNormDiffErr(
    "DxDp_vec_final", *DxDp_vec_final,
    "DxDp_fd_vec_final", *DxDp_fd_vec_final,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s_fd_correct, true );
  */

  double maxSensError = 1.0e-8;

  double s0_correct = Thyra::testRelNormDiffErr(
    "DxDp_s0_final", *DxDp_s0_final,
    "DxDp_exact_s0_final", *ds0dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s0_correct, true );

  double s1_correct = Thyra::testRelNormDiffErr(
    "DxDp_s1_final", *DxDp_s1_final,
    "DxDp_exact_s1_final", *ds1dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s1_correct, true );

  double s2_correct = Thyra::testRelNormDiffErr(
    "DxDp_s2_final", *DxDp_s2_final,
    "DxDp_exact_s2_final", *ds2dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s2_correct, true );

}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, exactSinCosFE ) {
  Teuchos::FancyOStream& std_out = out; // TEUCHOS_UNIT_TEST defines "out"
  double b = 5.0;
  //double phi = 0.0;
  double a = 2.0;
  double f = 3.0;
  double L = 4.0;
  double x0 = a;
  double x1 = b*f/L;
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Accept model parameters",true);
  modelPL->set("Implicit model formulation",false);
  modelPL->set("Provide nominal values",true);
  modelPL->set("Coeff a", a);
  modelPL->set("Coeff f", f);
  modelPL->set("Coeff L", L);
  modelPL->set("IC x_0", x0);
  modelPL->set("IC x_1", x1);
  stateModel->setParameterList(modelPL);
  const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type","Forward Euler");
  //stepperPL->set("Stepper Type","Explicit RK");
  builder->setParameterList(stepperPL);
  RCP<StepperBase<double> > stateStepper = builder->create();
  {
    // Set RKBT on ERK stepper
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Explicit 4 Stage");
    RCP<ExplicitRKStepper<double> > erkStepper = Teuchos::rcp_dynamic_cast<ExplicitRKStepper<double> >(stateStepper,false);
    if (!is_null(erkStepper)) {
      erkStepper->setRKButcherTableau(rkbt);
    }
  }
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver; 
  int p_index = 0;
  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper = 
    forwardSensitivityStepper<double>();
  stateAndSensStepper->initializeSyncedSteppers(
    stateModel, 
    p_index, 
    stateModel->getNominalValues(),
    stateStepper, 
    nonlinearSolver
    );

  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();

  RCP<Thyra::VectorBase<double> > s_bar_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  {
    // Initial conditions for Sensitivity problem.  (How do we get these correct in general?)
    RCP<Thyra::DefaultMultiVectorProductVector<double> > s_bar_mv =
      Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(s_bar_init,true);
    RCP<Thyra::VectorBase<double> > s_bar_mv_0 = s_bar_mv->getNonconstVectorBlock(0);
    Thyra::DetachedVectorView<double> s_bar_mv_0_view( *s_bar_mv_0 );
    s_bar_mv_0_view[0] = 1.0;
    s_bar_mv_0_view[1] = 0.0;
    RCP<Thyra::VectorBase<double> > s_bar_mv_1 = s_bar_mv->getNonconstVectorBlock(1);
    Thyra::DetachedVectorView<double> s_bar_mv_1_view( *s_bar_mv_1 );
    s_bar_mv_1_view[0] = 0.0;
    s_bar_mv_1_view[1] = b/L;
    RCP<Thyra::VectorBase<double> > s_bar_mv_2 = s_bar_mv->getNonconstVectorBlock(2);
    Thyra::DetachedVectorView<double> s_bar_mv_2_view( *s_bar_mv_2 );
    s_bar_mv_2_view[0] = 0.0;
    s_bar_mv_2_view[1] = -b*f/(L*L);
  }
  RCP<Thyra::VectorBase<double> > s_bar_dot_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  {
    RCP<Thyra::DefaultMultiVectorProductVector<double> > s_bar_dot_mv =
      Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(s_bar_dot_init,true);
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_0 = s_bar_dot_mv->getNonconstVectorBlock(0);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_0_view( *s_bar_dot_mv_0 );
    s_bar_dot_mv_0_view[0] = 0.0;
    s_bar_dot_mv_0_view[1] = 0.0;
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_1 = s_bar_dot_mv->getNonconstVectorBlock(1);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_1_view( *s_bar_dot_mv_1 );
    s_bar_dot_mv_1_view[0] = 0.0;
    s_bar_dot_mv_1_view[1] = -3.0*f*f*b/(L*L*L);
    RCP<Thyra::VectorBase<double> > s_bar_dot_mv_2 = s_bar_dot_mv->getNonconstVectorBlock(2);
    Thyra::DetachedVectorView<double> s_bar_dot_mv_2_view( *s_bar_dot_mv_2 );
    s_bar_dot_mv_2_view[0] = 0.0;
    s_bar_dot_mv_2_view[1] = 3.0*f*f*f*b/(L*L*L*L);
  }

  RCP<const StateAndForwardSensitivityModelEvaluator<double> >
    stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

  MEB::InArgs<double>
    state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar = [ x; s_bar ]
  state_and_sens_ic.set_x(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
    );
  // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
  //state_and_sens_ic.set_x_dot(
  //  stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
  //  );

  stateAndSensStepper->setInitialCondition(state_and_sens_ic);

  double finalTime = 1.0e-4;
  int numTimeSteps = 2;
  RCP<IntegratorBase<double> > integrator;
  {
    RCP<ParameterList> integratorPL = Teuchos::parameterList();
    integratorPL->set( "Take Variable Steps", false );
    integratorPL->set( "Fixed dt", Teuchos::as<double>((finalTime - state_ic.get_t())/numTimeSteps) );
    RCP<IntegratorBase<double> >
      defaultIntegrator = controlledDefaultIntegrator<double>(
        simpleIntegrationControlStrategy<double>(integratorPL)
        );
    integrator = defaultIntegrator;
  }

  RCP<StepperAsModelEvaluator<double> >
    stateAndSensIntegratorAsModel = stepperAsModelEvaluator(
      Teuchos::rcp_implicit_cast<StepperBase<double> >(stateAndSensStepper),
      integrator, state_and_sens_ic
      );
  stateAndSensIntegratorAsModel->setVerbLevel(verbLevel);
    
  RCP<Thyra::VectorBase<double> > x_bar_final;
  {
  
    x_bar_final = createMember(stateAndSensIntegratorAsModel->get_g_space(0));
  
    Thyra::eval_g(
      *stateAndSensIntegratorAsModel,
      0, *state_ic.get_p(0),
      finalTime,
      0, &*x_bar_final
      );

    std_out
      << "\nx_bar_final = x_bar(p,finalTime) evaluated using stateAndSensIntegratorAsModel:\n"
      << Teuchos::describe(*x_bar_final,verbLevel);

  }

  // Now we check that the sensitivities are correct
  RCP<const Thyra::VectorBase<double> >
    DxDp_vec_final = Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(1);
  std_out << "\nDxDp_vec_final:\n"
    << Teuchos::describe(*DxDp_vec_final,verbLevel);
  RCP<const Thyra::DefaultMultiVectorProductVector<double> > DxDp_mv_final =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<double> >(DxDp_vec_final,true);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s0_final = DxDp_mv_final->getVectorBlock(0);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s1_final = DxDp_mv_final->getVectorBlock(1);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s2_final = DxDp_mv_final->getVectorBlock(2);
  

  std_out << "\nDxDp_s0_final:\n"
    << Teuchos::describe(*DxDp_s0_final,verbLevel);
  std_out << "\nDxDp_s1_final:\n"
    << Teuchos::describe(*DxDp_s1_final,verbLevel);
  std_out << "\nDxDp_s2_final:\n"
    << Teuchos::describe(*DxDp_s2_final,verbLevel);

  MEB::InArgs<double> exactSensSolution;
  exactSensSolution = stateModel->getExactSensSolution(0,finalTime);
  RCP<const Thyra::VectorBase<double> > ds0dp = exactSensSolution.get_x();
  exactSensSolution = stateModel->getExactSensSolution(1,finalTime);
  RCP<const Thyra::VectorBase<double> > ds1dp = exactSensSolution.get_x();
  exactSensSolution = stateModel->getExactSensSolution(2,finalTime);
  RCP<const Thyra::VectorBase<double> > ds2dp = exactSensSolution.get_x();

  std_out << "\nds0dp exact:\n"
    << Teuchos::describe(*ds0dp,verbLevel);
  std_out << "\nds1dp exact:\n"
    << Teuchos::describe(*ds1dp,verbLevel);
  std_out << "\nds2dp exact:\n"
    << Teuchos::describe(*ds2dp,verbLevel);

  double maxSensError = 1.0e-8;

  verbLevel = Teuchos::VERB_EXTREME;
  double s0_correct = Thyra::testRelNormDiffErr(
    "DxDp_s0_final", *DxDp_s0_final,
    "DxDp_exact_s0_final", *ds0dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s0_correct, true );

  double s1_correct = Thyra::testRelNormDiffErr(
    "DxDp_s1_final", *DxDp_s1_final,
    "DxDp_exact_s1_final", *ds1dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s1_correct, true );

  double s2_correct = Thyra::testRelNormDiffErr(
    "DxDp_s2_final", *DxDp_s2_final,
    "DxDp_exact_s2_final", *ds2dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s2_correct, true );

}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, getNodes ) {
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Accept model parameters",true);
  modelPL->set("Implicit model formulation",true);
  modelPL->set("Provide nominal values",true);
  stateModel->setParameterList(modelPL);
  const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type","Backward Euler");
  builder->setParameterList(stepperPL);
  RCP<StepperBase<double> > stateStepper = builder->create();
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
  {
    // Set the nonlinear solver on the stepper.
    RCP<SolverAcceptingStepperBase<double> > SAStateStepper = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<double> >(stateStepper,true);
    SAStateStepper->setSolver(nonlinearSolver);
  }
  int p_index = 0;
  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper = 
    forwardSensitivityStepper<double>();
  {
    // Verify no nodes are returned when uninitialized
    Array<double> nodes;
    stateAndSensStepper->getNodes(&nodes);
    TEST_ASSERT( nodes.size() == 0 );
  }
  stateAndSensStepper->initializeSyncedSteppers(
    stateModel, 
    p_index, 
    stateModel->getNominalValues(),
    stateStepper, 
    nonlinearSolver
    );
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();
  RCP<Thyra::VectorBase<double> > s_bar_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  V_S(outArg(*s_bar_init),0.0);
  RCP<Thyra::VectorBase<double> > s_bar_dot_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  V_S(outArg(*s_bar_dot_init),0.0);
  RCP<const StateAndForwardSensitivityModelEvaluator<double> >
    stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

  MEB::InArgs<double>
    state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();
  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar = [ x; s_bar ]
  state_and_sens_ic.set_x(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
    );
  // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
  state_and_sens_ic.set_x_dot(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
    );
  stateAndSensStepper->setInitialCondition(state_and_sens_ic);

  {
    // Verify there is only one node when initialized
    Array<double> nodes;
    stateAndSensStepper->getNodes(&nodes);
    TEST_ASSERT( nodes.size() == 1 );
  }

  // Take a step.
  double dt = 0.1;
  stateAndSensStepper->takeStep(dt,STEP_TYPE_FIXED);
  {
    // Verify there are two nodes after a step
    Array<double> nodes;
    stateAndSensStepper->getNodes(&nodes);
    TEST_ASSERT( nodes.size() == 2 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, getPoints ) {
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Accept model parameters",true);
  modelPL->set("Implicit model formulation",true);
  modelPL->set("Provide nominal values",true);
  stateModel->setParameterList(modelPL);
  const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type","Backward Euler");
  builder->setParameterList(stepperPL);
  RCP<StepperBase<double> > stateStepper = builder->create();
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
  {
    // Set the nonlinear solver on the stepper.
    RCP<SolverAcceptingStepperBase<double> > SAStateStepper = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<double> >(stateStepper,true);
    SAStateStepper->setSolver(nonlinearSolver);
  }
  int p_index = 0;
  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper = 
    forwardSensitivityStepper<double>();

  stateAndSensStepper->initializeSyncedSteppers(
    stateModel, 
    p_index, 
    stateModel->getNominalValues(),
    stateStepper, 
    nonlinearSolver
    );
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();
  RCP<Thyra::VectorBase<double> > s_bar_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  V_S(outArg(*s_bar_init),0.0);
  RCP<Thyra::VectorBase<double> > s_bar_dot_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  V_S(outArg(*s_bar_dot_init),0.0);
  RCP<const StateAndForwardSensitivityModelEvaluator<double> >
    stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

  MEB::InArgs<double>
    state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();
  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar_ic = [ x; s_bar ]
  RCP<const VectorBase<double> > x_bar_ic = stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init);
  state_and_sens_ic.set_x(x_bar_ic);
  // Set initial condition for x_bar_dot_ic = [ x_dot; s_bar_dot ]
  RCP<const VectorBase<double> > x_bar_dot_ic = stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init);
  state_and_sens_ic.set_x_dot(x_bar_dot_ic);
  stateAndSensStepper->setInitialCondition(state_and_sens_ic);

  {
    // Verify we can get the initial condition when stepper is initialized but hasn't taken a step.
    typedef ScalarTraits<double> SMT;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    time_vec.push_back(state_and_sens_ic.get_t());
    stateAndSensStepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    TEST_ASSERT(
      Thyra::testRelNormDiffErr(
        "x_bar_ic", *x_bar_ic, "x_bar", *(x_vec[0]),
        "epsilon", SMT::eps(),
        "epsilon", SMT::eps(),
        0
        )
      );
    TEST_ASSERT(
      Thyra::testRelNormDiffErr(
        "x_bar_dot_ic", *x_bar_dot_ic, "x_bar_dot", *(xdot_vec[0]),
        "epsilon", SMT::eps(),
        "epsilon", SMT::eps(),
        0
        )
      );
  }

  // Take a step.
  double dt = 0.1;
  stateAndSensStepper->takeStep(dt,STEP_TYPE_FIXED);
  {
    // Verify we can still get the initial condition after a step.
    typedef ScalarTraits<double> SMT;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    time_vec.push_back(state_and_sens_ic.get_t());
    stateAndSensStepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    TEST_ASSERT(
      Thyra::testRelNormDiffErr(
        "x_bar_ic", *x_bar_ic, "x_bar", *(x_vec[0]),
        "epsilon", SMT::eps(),
        "epsilon", SMT::eps(),
        0
        )
      );
    TEST_ASSERT(
      Thyra::testRelNormDiffErr(
        "x_bar_dot_ic", *x_bar_dot_ic, "x_bar_dot", *(xdot_vec[0]),
        "epsilon", SMT::eps(),
        "epsilon", SMT::eps(),
        0
        )
      );
  }
}

} // namespace Rythmos

