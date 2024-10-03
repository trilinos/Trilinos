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
#include "ChanProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 4.0;
    double beta = 0.0;
    double scale = 1.0;
    int maxNewtonIters = 10;

    // Create output file to save solutions
    std::ofstream outFile("ChanTPContinuation.dat");
    outFile.setf(std::ios::scientific, std::ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << std::endl;

    // Create initial guess for the null vector of jacobian
    Teuchos::RCP<NOX::Abstract::Vector> nullVec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    nullVec->init(1.0);               // initial value 1.0

    // Create initial values for a and b for minimally augmented method
    Teuchos::RCP<NOX::Abstract::Vector> a_vec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    a_vec->init(1.0);

    Teuchos::RCP<NOX::Abstract::Vector> b_vec =
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    b_vec->init(1.0);

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");  // Default
    stepperList.set("Continuation Parameter", "beta");     // Must set
    stepperList.set("Initial Value", beta);                // Must set
    stepperList.set("Max Value", 1.0);                     // Must set
    stepperList.set("Min Value", 0.0);                     // Must set
    stepperList.set("Max Steps", 20);                      // Should set
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters); // Should set

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList =
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Turning Point");          // For turning point
    bifurcationList.set("Bifurcation Parameter", "alpha"); // Must set

    // For Moore-Spence formulation w/bordering
    //bifurcationList.set("Formulation", "Moore-Spence");          // Default
    //bifurcationList.set("Solver Method", "Salinger Bordering");  // Default
    //bifurcationList.set("Solver Method", "Phipps Bordering");
    //bifurcationList.set("Bordered Solver Method",
    //                    "LAPACK Direct Solve");   // For Phipps Bordering
    //bifurcationList.set("Length Normalization Vector", nullVec); // Must set
    //bifurcationList.set("Initial Null Vector", nullVec);         // Must set

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
    predictorList.set("Method", "Secant");                         // Default

    // Should use for Moore-Spence w/Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& firstStepPredictor
    //  = predictorList.sublist("First Step Predictor");
    //firstStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Epsilon", 1.0e-3);

    // Should use for Moore-Spence w/Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& lastStepPredictor
    //  = predictorList.sublist("Last Step Predictor");
    //lastStepPredictor.set("Method", "Random");
    //lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");                      // Default
    stepSizeList.set("Initial Step Size", 0.1);                  // Should set
    stepSizeList.set("Min Step Size", 1.0e-3);                   // Should set
    stepSizeList.set("Max Step Size", 1.0);                      // Should set

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
              NOX::Utils::StepperParameters);            // Should set

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp =
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));

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

    outFile.close();

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
