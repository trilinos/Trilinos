// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_Stepper_impl_hpp
#define Tempus_Stepper_impl_hpp


namespace Tempus {


void getValidParametersBasic(
  Teuchos::RCP<Teuchos::ParameterList> pl, std::string stepperType)
{
  pl->setName("Default Stepper - " + stepperType);
  pl->set<std::string>("Stepper Type", stepperType);

  pl->set<bool>("Use FSAL", false,
    "The First-Step-As-Last (FSAL) principle is the situation where the\n"
    "last function evaluation, f(x^{n-1},t^{n-1}) [a.k.a. xDot^{n-1}],\n"
    "can be used for the first function evaluation, f(x^n,t^n)\n"
    "[a.k.a. xDot^n].  For RK methods, this applies to the stages.\n"
    "\n"
    "Often the FSAL priniciple can be used to save an evaluation.\n"
    "However there are cases when it cannot be used, e.g., operator\n"
    "splitting where other steppers/operators have modified the solution,\n"
    "x^*, and thus require the function evaluation, f(x^*, t^{n-1}).\n"
    "\n"
    "It should be noted that when the FSAL priniciple can be used\n"
    "(can set useFSAL=true), setting useFSAL=false will give the\n"
    "same solution but at additional expense.  However, the reverse\n"
    "is not true.  When the FSAL priniciple can not be used\n"
    "(need to set useFSAL=false), setting useFSAL=true will produce\n"
    "incorrect solutions.\n"
    "\n"
    "Default in general for explicit and implicit steppers is false,\n"
    "but individual steppers can override this default.");

  pl->set<std::string>("Initial Condition Consistency", "None",
    "This indicates which type of consistency should be applied to\n"
    "the initial conditions (ICs):\n"
    "\n"
    "  'None'   - Do nothing to the ICs provided in the SolutionHistory.\n"
    "  'Zero'   - Set the derivative of the SolutionState to zero in the\n"
    "             SolutionHistory provided, e.g., xDot^0 = 0, or \n"
    "             xDotDot^0 = 0.\n"
    "  'App'    - Use the application's ICs, e.g., getNominalValues().\n"
    "  'Consistent' - Make the initial conditions for x and xDot\n"
    "             consistent with the governing equations, e.g.,\n"
    "             xDot = f(x,t), and f(x, xDot, t) = 0.  For implicit\n"
    "             ODEs, this requires a solve of f(x, xDot, t) = 0 for\n"
    "             xDot, and another Jacobian and residual may be\n"
    "             needed, e.g., boundary conditions on xDot may need\n"
    "             to replace boundary conditions on x.\n"
    "\n"
    "In general for explicit steppers, the default is 'Consistent',\n"
    "because it is fairly cheap with just one residual evaluation.\n"
    "In general for implicit steppers, the default is 'None', because\n"
    "the application often knows its IC and can set it the initial\n"
    "SolutionState.  Also, as noted above, 'Consistent' may require\n"
    "another Jacobian from the application.  Individual steppers may\n"
    "override these defaults.");

  pl->set<bool>("Initial Condition Consistency Check", true,
    "Check if the initial condition, x and xDot, is consistent with the\n"
    "governing equations, xDot = f(x,t), or f(x, xDot, t) = 0.\n"
    "\n"
    "In general for explicit and implicit steppers, the default is true,\n"
    "because it is fairly cheap with just one residual evaluation.\n"
    "Individual steppers may override this default.");
}


Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters()
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // NOX Solver ParameterList
  RCP<ParameterList> noxPL = Teuchos::parameterList();

    // Direction ParameterList
    RCP<ParameterList> directionPL = Teuchos::parameterList();
    directionPL->set<std::string>("Method", "Newton");
      RCP<ParameterList> newtonPL = Teuchos::parameterList();
      newtonPL->set<std::string>("Forcing Term Method", "Constant");
      newtonPL->set<bool>       ("Rescue Bad Newton Solve", 1);
      directionPL->set("Newton", *newtonPL);
    noxPL->set("Direction", *directionPL);

    // Line Search ParameterList
    RCP<ParameterList> lineSearchPL = Teuchos::parameterList();
    lineSearchPL->set<std::string>("Method", "Full Step");
      RCP<ParameterList> fullStepPL = Teuchos::parameterList();
      fullStepPL->set<double>("Full Step", 1);
      lineSearchPL->set("Full Step", *fullStepPL);
    noxPL->set("Line Search", *lineSearchPL);

    noxPL->set<std::string>("Nonlinear Solver", "Line Search Based");

    // Printing ParameterList
    RCP<ParameterList> printingPL = Teuchos::parameterList();
    printingPL->set<int>("Output Precision", 3);
    printingPL->set<int>("Output Processor", 0);
      RCP<ParameterList> outputPL = Teuchos::parameterList();
      outputPL->set<bool>("Error", 1);
      outputPL->set<bool>("Warning", 1);
      outputPL->set<bool>("Outer Iteration", 0);
      outputPL->set<bool>("Parameters", 0);
      outputPL->set<bool>("Details", 0);
      outputPL->set<bool>("Linear Solver Details", 1);
      outputPL->set<bool>("Stepper Iteration", 1);
      outputPL->set<bool>("Stepper Details", 1);
      outputPL->set<bool>("Stepper Parameters", 1);
      printingPL->set("Output Information", *outputPL);
    noxPL->set("Printing", *printingPL);

    // Solver Options ParameterList
    RCP<ParameterList> solverOptionsPL = Teuchos::parameterList();
    solverOptionsPL->set<std::string>("Status Test Check Type", "Minimal");
    noxPL->set("Solver Options", *solverOptionsPL);

    // Status Tests ParameterList
    RCP<ParameterList> statusTestsPL = Teuchos::parameterList();
    statusTestsPL->set<std::string>("Test Type", "Combo");
    statusTestsPL->set<std::string>("Combo Type", "OR");
    statusTestsPL->set<int>("Number of Tests", 2);
      RCP<ParameterList> test0PL = Teuchos::parameterList();
      test0PL->set<std::string>("Test Type", "NormF");
      test0PL->set<double>("Tolerance", 1e-08);
      statusTestsPL->set("Test 0", *test0PL);
      RCP<ParameterList> test1PL = Teuchos::parameterList();
      test1PL->set<std::string>("Test Type", "MaxIters");
      test1PL->set<int>("Maximum Iterations", 10);
      statusTestsPL->set("Test 1", *test1PL);
    noxPL->set("Status Tests", *statusTestsPL);

  // Solver ParameterList
  RCP<ParameterList> solverPL = Teuchos::parameterList();
  solverPL->set("NOX", *noxPL);

  return solverPL;
}


template<class Scalar>
void Stepper<Scalar>::initialize()
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();

  bool isValidSetup = this->isValidSetup(*out);

  if (isValidSetup)
    this->isInitialized_ = true;   // Only place it is set to true.
  else
    this->describe(*out, Teuchos::VERB_MEDIUM);
}


template<class Scalar>
void Stepper<Scalar>::checkInitialized()
{
#ifdef TEMPUS_HIDE_DEPRECATED_CODE
  if ( !this->isInitialized() ) {
    this->describe( *(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION( !this->isInitialized(), std::logic_error,
      "Error - " << this->description() << " is not initialized!");
  }
#endif // TEMPUS_HIDE_DEPRECATED_CODE
}


template<class Scalar>
void Stepper<Scalar>::describe(Teuchos::FancyOStream        & out,
                               const Teuchos::EVerbosityLevel verbLevel) const
{
  out << "--- Stepper ---\n"
      << "  isInitialized_      = " << Teuchos::toString(isInitialized_) << std::endl
      << "  stepperType_        = " << stepperType_ << std::endl
      << "  useFSAL_            = " << Teuchos::toString(useFSAL_) << std::endl
      << "  ICConsistency_      = " << ICConsistency_ << std::endl
      << "  ICConsistencyCheck_ = " << Teuchos::toString(ICConsistencyCheck_) << std::endl;
}


template<class Scalar>
bool Stepper<Scalar>::isValidSetup(
  Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !(ICConsistency_ == "None" || ICConsistency_ == "Zero" ||
         ICConsistency_ == "App"  || ICConsistency_ == "Consistent") ) {
    isValidSetup = false;
    out << "The IC consistency does not have a valid value!\n"
        << "('None', 'Zero', 'App' or 'Consistent')\n"
        << "  ICConsistency  = " << ICConsistency_ << "\n";
  }

  return isValidSetup;
}


} // namespace Tempus
#endif // Tempus_Stepper_impl_hpp
