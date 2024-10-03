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
#include "BrusselatorProblemInterface.H"

int main()
{
  int n = 100;
  double alpha = 0.25;
  double beta = 0.0;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  int maxNewtonIters = 10;

  try {

    // Create output file to save solutions
    std::ofstream outFile("BrusselatorContinuation.dat");
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
    stepperList.set("Continuation Method", "Arc Length");    // Default
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 2.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Use LAPACK solver for arclength bordered solves for efficiency
    stepperList.set("Bordered Solver Method", "LAPACK Direct Solve");

    // Create Eigensolver list using LAPACK eigensolver
    stepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& eigenList = stepperList.sublist("Eigensolver");
    eigenList.set("Method", "DGGEV");     // LAPACK eigensolver
    eigenList.set("NEV", 10);             // Show 10 eigenvalues
    eigenList.set("Sorting Order", "LR"); // Sort by largest real component

    // Create predictor sublist.  We use the constant predictor because the
    // solution is linear in the continuation parameter.  Using Secant or
    // Tangent is not very interesting, but are perfectly valid
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Constant");
    //predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");    // Default
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 10.0);
    stepSizeList.set("Aggressiveness", 0.1);

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
              NOX::Utils::StepperParameters);

    // Create LAPACK Factory (necessary for LAPACK eigensolver)
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    BrusselatorProblemInterface brus(globalData, n, alpha, beta, D1, D2,
                     outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("D1",D1);
    p.addParameter("D2",D2);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp =
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, brus));

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
