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
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;

  alpha = alpha / scale;

  try {

    // Create output file to save solutions
    std::ofstream outFile("ChanContinuation.dat");
    outFile.setf(std::ios::scientific, std::ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << std::endl;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");// Default
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Parameter", "alpha");  // Must set
    stepperList.set("Initial Value", alpha);             // Must set
    stepperList.set("Max Value", 5.0/scale);             // Must set
    stepperList.set("Min Value", 0.0/scale);             // Must set
    stepperList.set("Max Steps", 50);                    // Should set
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters); // Should set
    stepperList.set("Compute Eigenvalues",false);        // Default

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList =
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");                 // Default

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");               // Default
    //predictorList.set("Method", "Constant");
    //predictorList.set("Method", "Tangent");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");             // Default
    stepSizeList.set("Initial Step Size", 0.1/scale);   // Should set
    stepSizeList.set("Min Step Size", 1.0e-3/scale);    // Should set
    stepSizeList.set("Max Step Size", 10.0/scale);      // Should set

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information",
              NOX::Utils::Details +
              NOX::Utils::OuterIteration +
              NOX::Utils::InnerIteration +
              NOX::Utils::Warning +
              NOX::Utils::StepperIteration +
              NOX::Utils::StepperDetails +
              NOX::Utils::StepperParameters);  // Should set

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
    Teuchos::RCP<NOX::StatusTest::NormF> normF =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                          normF,
                          maxIters));

    // Create the stepper
    LOCA::Stepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
    if (status == LOCA::Abstract::Iterator::Finished)
      std::cout << "All examples passed" << std::endl;
    else {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out()
      << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::LAPACK::Group> finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution =
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

    // Print final solution
    globalData->locaUtils->out()
                << std::endl << "Final solution is " << std::endl;
    finalGroup->printSolution(finalSolution,
                      finalGroup->getParam("alpha"));

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
