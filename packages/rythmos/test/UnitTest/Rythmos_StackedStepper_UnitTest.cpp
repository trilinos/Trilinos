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
#include "Rythmos_UnitTestModels.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_StackedStepper.hpp"
#include "Rythmos_ForwardSensitivityExplicitModelEvaluator.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "../ConvergenceTest/Rythmos_ConvergenceTestHelpers.hpp"

#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, create ) { 
  RCP<StackedStepper<double> > sStepper = stackedStepper<double>();
  TEST_ASSERT( !is_null(sStepper) );
}


// StackedStepper will throw on takeStep if you haven't called addStepper yet.
TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, invalidTakeStep ) {
  RCP<StackedStepper<double> > sStepper = stackedStepper<double>();
  TEST_THROW( sStepper->takeStep(0.1,STEP_TYPE_FIXED), std::logic_error );
}


TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, addgetStepper ) {
  RCP<StackedStepper<double> > sStepper = stackedStepper<double>();
  {
    ArrayView<const RCP<StepperBase<double> > > stepperArray; 
    stepperArray = sStepper->getNonconstSteppers();
    TEST_ASSERT( stepperArray.size() == 0 );
  }
  {
    RCP<StepperBase<double> > stepper;
    const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Stepper Type","Explicit RK");
    builder->setParameterList(pl);
    stepper = builder->create();
    sStepper->addStepper(stepper);
    pl->set("Stepper Type","Implicit RK");
    sStepper->addStepper(builder->create());
    pl->set("Stepper Type","Implicit BDF");
    sStepper->addStepper(builder->create());
  }
  ArrayView<const RCP<StepperBase<double> > > stepperArray;
  stepperArray = sStepper->getNonconstSteppers();
  TEST_ASSERT(stepperArray.size() == 3);

  RCP<const StepperBase<double> > stepper;
  stepper = stepperArray[0];
  TEST_ASSERT( !is_null(stepper) );
  {
    RCP<const ExplicitRKStepper<double> > erkStepper = 
      Teuchos::rcp_dynamic_cast<const ExplicitRKStepper<double> >(
          stepper,
          false
          );
    TEST_ASSERT( !is_null(erkStepper) );
  }
  stepper = stepperArray[1];
  TEST_ASSERT( !is_null(stepper) );
  {
    RCP<const ImplicitRKStepper<double> > irkStepper = 
      Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(
          stepper,
          false
          );
    TEST_ASSERT( !is_null(irkStepper) );
  }
  stepper = stepperArray[2];
  TEST_ASSERT( !is_null(stepper) );
  {
    RCP<const ImplicitBDFStepper<double> > ibdfStepper = 
      Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(
          stepper,
          false
          );
    TEST_ASSERT( !is_null(ibdfStepper) );
  }
}


//TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, spaces ) {
//  // Verify get_x_space is the right size and type
//  // Verify get_x_space throws before addStepper is called
//  TEST_ASSERT( false );
//}
//
//
//TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, getPoints ) {
//  // Verify getPoints does the right thing
//  // Verify getPoints throws before addStepper is called
//  TEST_ASSERT( false );
//}


TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, sameStepper ) {
  typedef Thyra::ModelEvaluatorBase MEB;
  RCP<StepperBase<double> > fwdStepperA;
  RCP<StepperBase<double> > fwdStepperB;
  {
    const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
    RCP<ParameterList> stepperPL = Teuchos::parameterList();
    stepperPL->set("Stepper Type","Explicit RK");
    builder->setParameterList(stepperPL);
    fwdStepperA = builder->create();
    fwdStepperB = builder->create();
    RCP<ExplicitRKStepper<double> > erkFwdStepperA = 
      Teuchos::rcp_dynamic_cast<ExplicitRKStepper<double> >(
          fwdStepperA,
          true
          );
    RCP<ExplicitRKStepper<double> > erkFwdStepperB = 
      Teuchos::rcp_dynamic_cast<ExplicitRKStepper<double> >(
          fwdStepperB,
          true
          );
    RCP<RKButcherTableauBase<double> > rkbt = 
      createRKBT<double>("Forward Euler");
    erkFwdStepperA->setRKButcherTableau(rkbt);
    erkFwdStepperB->setRKButcherTableau(rkbt);
  }
  RCP<SinCosModel> fwdModelA = sinCosModel(false);
  RCP<SinCosModel> fwdModelB = sinCosModel(false);
  const MEB::InArgs<double> fwdModelA_ic = fwdModelA->getNominalValues();
  const MEB::InArgs<double> fwdModelB_ic = fwdModelB->getNominalValues();
  fwdStepperA->setModel(fwdModelA);
  fwdStepperA->setInitialCondition(fwdModelA_ic);
  fwdStepperB->setModel(fwdModelB);
  fwdStepperB->setInitialCondition(fwdModelB_ic);

  RCP<StackedStepper<double> > sStepper = stackedStepper<double>();
  sStepper->addStepper(fwdStepperA);
  sStepper->addStepper(fwdStepperB);

  double finalTime = 1.0e-4;
  int numTimeSteps = 2;
  double dt = finalTime/numTimeSteps; // Assume t_0 = 0.0;
  for (int i=0 ; i < numTimeSteps ; ++i ) {
    double dt_taken = sStepper->takeStep(dt,STEP_TYPE_FIXED);
    TEST_ASSERT( dt_taken == dt );
  }
  RCP<VectorBase<double> > x_bar_final = 
    Thyra::createMember(sStepper->get_x_space());
  {
    Array<double> t_vec;
    Array<RCP<const VectorBase<double> > > x_vec;

    t_vec.push_back(finalTime);
    sStepper->getPoints(
        t_vec,
        &x_vec,
        NULL,
        NULL
        );
    V_V(Teuchos::outArg(*x_bar_final),*x_vec[0]);
  }
  RCP<const Thyra::VectorBase<double> >
    solA = Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(0);
  RCP<const Thyra::VectorBase<double> >
    solB = Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(1);

  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  typedef Teuchos::ScalarTraits<double> ST;
  double maxSensError = ST::eps();
  double sol_correct = Thyra::testRelNormDiffErr(
    "sol A", *solA,
    "sol B", *solB,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &out, verbLevel
    );
  TEST_ASSERT( sol_correct );

  // 07/23/09 tscoffe: 
  // Ideally, I should check that the answer is correct also...


}


double computeForwardSensitivityErrorStackedStepperSinCosFE(
    int numTimeSteps,
    Array<RCP<const VectorBase<double> > >& computedSol,
    Array<RCP<const VectorBase<double> > >& exactSol
    )
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::ModelEvaluatorBase MEB;
  // Forward ODE Model:
  RCP<SinCosModel> fwdModel = sinCosModel();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",false);
    pl->set("Provide nominal values",true);
    double b = 5.0;
    //double phi = 0.0;
    double a = 2.0;
    double f = 3.0;
    double L = 4.0;
    double x0 = a;
    double x1 = b*f/L;
    pl->set("Coeff a", a);
    pl->set("Coeff f", f);
    pl->set("Coeff L", L);
    pl->set("IC x_0", x0);
    pl->set("IC x_1", x1);
    fwdModel->setParameterList(pl);
  }
  RCP<StepperBase<double> > fwdStepper;
  RCP<StepperBase<double> > fsStepper;
  {
    const RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
    RCP<ParameterList> stepperPL = Teuchos::parameterList();
    stepperPL->set("Stepper Type","Forward Euler");
    builder->setParameterList(stepperPL);
    fwdStepper = builder->create();
    fsStepper = builder->create();
  }
  // Forward Sensitivity Model:
  RCP<ForwardSensitivityExplicitModelEvaluator<double> > fsModel = 
    forwardSensitivityExplicitModelEvaluator<double>();
  int p_index = 0;
  fsModel->initializeStructure(fwdModel,p_index);

  const MEB::InArgs<double> fwdModel_ic = fwdModel->getNominalValues();
  fwdStepper->setModel(fwdModel);
  fwdStepper->setInitialCondition(fwdModel_ic);
  fsModel->initializePointState(Teuchos::inOutArg(*fwdStepper),false);

  MEB::InArgs<double> fsModel_ic = fsModel->getNominalValues();
  {
    // Set up sensitivity initial conditions so they match the initial
    // conditions in getExactSensSolution
    RCP<Thyra::VectorBase<double> > s_bar_init
      = createMember(fsModel->get_x_space());
    RCP<Thyra::DefaultMultiVectorProductVector<double> > s_bar_mv =
      rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
          s_bar_init,
          true
          );
    int np = 3; // SinCos problem number of elements in parameter vector.
    for (int j=0 ; j < np ; ++j) {
      MEB::InArgs<double> sens_ic = fwdModel->getExactSensSolution(j,0.0);
      V_V(outArg(*(s_bar_mv->getNonconstVectorBlock(j))),
          *(sens_ic.get_x())
          );
    }
    fsModel_ic.set_x(s_bar_init);
  }
  fsStepper->setModel(fsModel);
  fsStepper->setInitialCondition(fsModel_ic);

  RCP<StackedStepper<double> > sStepper = stackedStepper<double>();
  sStepper->addStepper(fwdStepper);
  sStepper->addStepper(fsStepper);
  {
    // Set up Forward Sensitivities step strategy
    RCP<ForwardSensitivityStackedStepperStepStrategy<double> > stepStrategy = 
      forwardSensitivityStackedStepperStepStrategy<double>();
    sStepper->setStackedStepperStepControlStrategy(stepStrategy);
  }

  double finalTime = 1.0e-4;
  double dt = finalTime/numTimeSteps; // Assume t_0 = 0.0;
  for (int i=0 ; i < numTimeSteps ; ++i ) {
    double dt_taken = sStepper->takeStep(dt,STEP_TYPE_FIXED);
    TEUCHOS_ASSERT( dt_taken == dt );
  }
  RCP<VectorBase<double> > x_bar_final = 
    Thyra::createMember(sStepper->get_x_space());
  {
    Array<double> t_vec;
    Array<RCP<const VectorBase<double> > > x_vec;

    t_vec.push_back(finalTime);
    sStepper->getPoints(
        t_vec,
        &x_vec,
        NULL,
        NULL
        );
    V_V(Teuchos::outArg(*x_bar_final),*x_vec[0]);
  }

  // Now we check that the sensitivities are correct
  RCP<const Thyra::VectorBase<double> > DxDp_vec_final = 
    Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(1);
  RCP<const Thyra::DefaultMultiVectorProductVector<double> > DxDp_mv_final =
    rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<double> >(
        DxDp_vec_final,
        true
        );
  RCP<const Thyra::VectorBase<double> >
    DxDp_s0_final = DxDp_mv_final->getVectorBlock(0);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s1_final = DxDp_mv_final->getVectorBlock(1);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s2_final = DxDp_mv_final->getVectorBlock(2);

  computedSol.clear();
  computedSol.push_back(DxDp_s0_final);
  computedSol.push_back(DxDp_s1_final);
  computedSol.push_back(DxDp_s2_final);
  
  MEB::InArgs<double> exactSensSolution;
  exactSensSolution = fwdModel->getExactSensSolution(0,finalTime);
  RCP<const Thyra::VectorBase<double> > ds0dp = exactSensSolution.get_x();
  exactSensSolution = fwdModel->getExactSensSolution(1,finalTime);
  RCP<const Thyra::VectorBase<double> > ds1dp = exactSensSolution.get_x();
  exactSensSolution = fwdModel->getExactSensSolution(2,finalTime);
  RCP<const Thyra::VectorBase<double> > ds2dp = exactSensSolution.get_x();

  exactSol.clear();
  exactSol.push_back(ds0dp);
  exactSol.push_back(ds1dp);
  exactSol.push_back(ds2dp);

  return dt;
}


TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, forwardSensitivities ) {
  using Teuchos::rcp_dynamic_cast;
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  Array<RCP<const VectorBase<double> > > computedSol;
  Array<RCP<const VectorBase<double> > > exactSol;

  computeForwardSensitivityErrorStackedStepperSinCosFE(
      1,
      computedSol,
      exactSol
      );

  out << "\nDxDp_s0_final:\n"
    << Teuchos::describe(*computedSol[0],verbLevel);
  out << "\nDxDp_s1_final:\n"
    << Teuchos::describe(*computedSol[1],verbLevel);
  out << "\nDxDp_s2_final:\n"
    << Teuchos::describe(*computedSol[2],verbLevel);


  out << "\nds0dp exact:\n"
    << Teuchos::describe(*exactSol[0],verbLevel);
  out << "\nds1dp exact:\n"
    << Teuchos::describe(*exactSol[1],verbLevel);
  out << "\nds2dp exact:\n"
    << Teuchos::describe(*exactSol[2],verbLevel);

  double maxSensError = 1.0e-8;

  verbLevel = Teuchos::VERB_EXTREME;
  double s0_correct = Thyra::testRelNormDiffErr(
    "DxDp_s0_final", *computedSol[0],
    "DxDp_exact_s0_final", *exactSol[0],
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &out, verbLevel
    );
  TEST_EQUALITY_CONST( s0_correct, true );

  double s1_correct = Thyra::testRelNormDiffErr(
    "DxDp_s1_final", *computedSol[1],
    "DxDp_exact_s1_final", *exactSol[1],
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &out, verbLevel
    );
  TEST_EQUALITY_CONST( s1_correct, true );

  double s2_correct = Thyra::testRelNormDiffErr(
    "DxDp_s2_final", *computedSol[2],
    "DxDp_exact_s2_final", *exactSol[2],
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &out, verbLevel
    );
  TEST_EQUALITY_CONST( s2_correct, true );

}

double local2NormError(
    const Array<RCP<const VectorBase<double> > >& vecA,
    const Array<RCP<const VectorBase<double> > >& vecB
    )
{
  TEUCHOS_ASSERT( vecA.size() == vecB.size() );
  double error = 0.0;
  for (int i = 0 ; i < Teuchos::as<int>(vecA.size()) ; ++i) {
    RCP<VectorBase<double> > errVec = Thyra::createMember(vecA[i]->space());
    V_StVpStV(Teuchos::outArg(*errVec),1.0,*vecA[i],-1.0,*vecB[i]);
    double tmpError = Thyra::norm(*errVec);
    error += tmpError*tmpError;
  }
  return pow(error,0.5);
}

TEUCHOS_UNIT_TEST( Rythmos_StackedStepper, forwardSensitivitiesConvergence ) {
  Array<double> logHVec;
  Array<double> logEVec;
  for (int i=1 ; i<= 4 ; i*=2)  {
    Array<RCP<const VectorBase<double> > > computedSol;
    Array<RCP<const VectorBase<double> > > exactSol;
    double h = computeForwardSensitivityErrorStackedStepperSinCosFE( 
        i,
        computedSol,
        exactSol
        );
    double logH = log(h);
    double logE = log( local2NormError(computedSol,exactSol) );
    logHVec.push_back( logH );
    logEVec.push_back( logE );
  }
  double slope = computeLinearRegressionSlope<double>(logHVec,logEVec);
  double tol = 1.0e-6;
  TEST_FLOATING_EQUALITY( slope, 1.0, tol ); }

} // namespace Rythmos

