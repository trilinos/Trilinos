// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "PitchforkProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 1.0;
    double beta = 0.0;
    double lambda = -2.25;
    double pi = 4.0*atan(1.0);
    double h = 2.0/(n-1);
    int maxNewtonIters = 20;

    // Create asymmetry vector
    Teuchos::RCP<NOX::LAPACK::Vector> asymLapackVec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));  // length n
    for (int i=0; i<n; i++)
      (*asymLapackVec)(i) = sin( pi/2.0 * (-1.0 + h*i) );
    Teuchos::RCP<NOX::Abstract::Vector> asymVec = asymLapackVec;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create initial values for a and b for minimally augmented method
    Teuchos::RCP<NOX::Abstract::Vector> a_vec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    *a_vec = *asymVec;

    Teuchos::RCP<NOX::Abstract::Vector> b_vec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    *b_vec = *asymVec;

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");         // Default
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", -0.01);
    stepperList.set("Max Steps", 20);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList =
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Pitchfork");
    bifurcationList.set("Bifurcation Parameter", "lambda");       // Must set
    bifurcationList.set("Antisymmetric Vector", asymVec);         // Must set

//     // For Moore-Spence Formulation
//     bifurcationList.set("Formulation", "Moore-Spence");           // Default
//     //bifurcationList.set("Solver Method", "Salinger Bordering"); // Default
//     bifurcationList.set("Solver Method", "Phipps Bordering");
//     bifurcationList.set("Bordered Solver Method",
//                       "LAPACK Direct Solve");  // For Phipps Bordering
//     bifurcationList.set("Length Normalization Vector", asymVec);  // Must set
//     bifurcationList.set("Initial Null Vector", asymVec);          // Must set

    // For minimally augmented formulation
    bifurcationList.set("Formulation", "Minimally Augmented");
    bifurcationList.set("Initial A Vector", a_vec);                // Must set
    bifurcationList.set("Initial B Vector", b_vec);                // Must set

    // For minimally augmented method, should set these for good performance
    // Direct solve of bordered equations
    bifurcationList.set("Bordered Solver Method",  "LAPACK Direct Solve");
    // Combine arc-length and turning point bordered rows & columns
    stepperList.set("Bordered Solver Method", "Nested");
    Teuchos::ParameterList& nestedList =
      stepperList.sublist("Nested Bordered Solver");
    // Direct solve of combined bordered system
    nestedList.set("Bordered Solver Method", "LAPACK Direct Solve");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");                       // Default

    // Should use for Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& firstStepPredictor
    //  = predictorList.sublist("First Step Predictor");
    //firstStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Epsilon", 1.0e-3);

    // Should use for Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& lastStepPredictor
    //  = predictorList.sublist("Last Step Predictor");
    //lastStepPredictor.set("Method", "Random");
    //lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");                      // Default
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 1.0);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information",
              NOX::Utils::OuterIteration +
              NOX::Utils::OuterIterationStatusTest +
              NOX::Utils::InnerIteration +
              NOX::Utils::Details +
              NOX::Utils::Warning +
              NOX::Utils::StepperIteration +
              NOX::Utils::StepperDetails +
              NOX::Utils::StepperParameters);

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    PitchforkProblemInterface pf(n, alpha, beta, lambda);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("lambda",lambda);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp =
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, pf));

    grp->setParams(p);

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> statusTestA =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10,
                          NOX::StatusTest::NormF::Scaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                          statusTestA, statusTestB));

    // Create the stepper
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status == LOCA::Abstract::Iterator::Finished)
      std::cout << "All examples passed" << std::endl;
    else {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out()
      << "Stepper failed to converge!" << std::endl;
    }

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
    << std::endl << "Final Parameters" << std::endl
    << "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

  return 0;
}
