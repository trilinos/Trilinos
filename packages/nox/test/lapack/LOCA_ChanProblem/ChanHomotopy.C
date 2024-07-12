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
#include "NOX_TestCompare.H"

int main(int argc, char *argv[])
{
  int n = 100;
  double alpha = 10.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;
  int ierr = 0;

  alpha = alpha / scale;

  try {

    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v')
    verbose = true;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
       nlPrintParams.set("Output Information",
             NOX::Utils::Error +
             NOX::Utils::Details +
             NOX::Utils::OuterIterationStatusTest +
             NOX::Utils::OuterIteration +
             NOX::Utils::InnerIteration +
             NOX::Utils::Warning +
             NOX::Utils::TestDetails +
             NOX::Utils::StepperIteration +
             NOX::Utils::StepperDetails +
             NOX::Utils::StepperParameters);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::LAPACK::Group> grp =
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    grp->setParams(p);

    // Set up the status tests
    grp->computeF();
    double size = sqrt(static_cast<double>(grp->getX().length()));
    double normF_0 = grp->getNormF() / size;
    double tolerance = 1.0e-8 / normF_0 / normF_0;
    Teuchos::RCP<NOX::StatusTest::NormF> statusTestA =
      Teuchos::rcp(new NOX::StatusTest::NormF(*grp, tolerance));
    Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                           statusTestA, statusTestB));

    // Create the homotopy group
    Teuchos::RCP<LOCA::Homotopy::Group> hGrp =
      Teuchos::rcp(new LOCA::Homotopy::Group(locaParamsList, globalData, grp));

    // Override default parameters
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Create the stepper
    LOCA::Stepper stepper(globalData, hGrp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out() << "Stepper failed to converge!"
                     << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::LAPACK::Group> finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution =
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
    << std::endl << "Final Parameters" << std::endl
    << "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Check some statistics on the solution
    NOX::Utils utils(nlPrintParams);
    NOX::TestCompare testCompare(std::cout, utils);

    if (utils.isPrintType(NOX::Utils::TestDetails))
      std::cout << std::endl << "***** Checking solutions statistics *****" << std::endl;

    // Check number of steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 14;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
                  "number of continuation steps",
                  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 3;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected,
                  0.0, "number of failed continuation steps",
                  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double alpha_final =
      finalGroup->getParam("Homotopy Continuation Parameter");
    double alpha_expected = 1.0;
    ierr += testCompare.testValue(alpha_final, alpha_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected =  456.0303417;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-7,
                  "norm of final solution",
                  NOX::TestCompare::Relative);

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  if (ierr == 0)
    std::cout << "All tests passed!" << std::endl;
  else
    std::cout << ierr << " test(s) failed!" << std::endl;

  return 0;
}
