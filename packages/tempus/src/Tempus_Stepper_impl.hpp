//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Stepper_impl_hpp
#define Tempus_Stepper_impl_hpp

#include "NOX_Thyra.H"

namespace Tempus {

template <class Scalar>
void Stepper<Scalar>::initialize()
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->setOutputToRootOnly(0);

  bool isValidSetup = this->isValidSetup(*out);

  if (isValidSetup)
    this->isInitialized_ = true;  // Only place it is set to true.
  else
    this->describe(*out, Teuchos::VERB_MEDIUM);
}

template <class Scalar>
void Stepper<Scalar>::checkInitialized()
{
  if (!this->isInitialized()) {
    this->describe(*(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION(
        !this->isInitialized(), std::logic_error,
        "Error - " << this->description() << " is not initialized!");
  }
}

template <class Scalar>
void Stepper<Scalar>::setUseFSALTrueOnly(bool a)
{
  if (a == false) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "Stepper::setUseFSALTrueOnly()");
    *out << "Warning -- useFSAL for '" << this->getStepperType() << "'\n"
         << "can only be set to true.  Leaving set to true." << std::endl;
  }
  useFSAL_ = true;
}

template <class Scalar>
void Stepper<Scalar>::setUseFSALFalseOnly(bool a)
{
  if (a == true) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "Stepper::setUseFSALFalseOnly()");
    *out << "Warning -- useFSAL for '" << this->getStepperType() << "'\n"
         << "can only be set to false.  Leaving set to false." << std::endl;
  }
  useFSAL_ = false;
}

template <class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> > Stepper<Scalar>::getStepperX()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperX_ == Teuchos::null, std::logic_error,
      "Error - stepperX_ has not been set in setInitialConditions() or\n"
      "        can not be set from the state!\n");

  return stepperX_;
}

template <class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> > Stepper<Scalar>::getStepperXDot()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperXDot_ == Teuchos::null, std::logic_error,
      "Error - stepperXDot_ has not been set in setInitialConditions() or\n"
      "        can not be set from the state!\n");

  return stepperXDot_;
}

template <class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> > Stepper<Scalar>::getStepperXDotDot()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperXDotDot_ == Teuchos::null, std::logic_error,
      "Error - stepperXDotDot_ has not been set in setInitialConditions() or\n"
      "        can not be set from the state!\n");

  return stepperXDotDot_;
}

template <class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> > Stepper<Scalar>::getStepperXDotDot(
    Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDotDot() != Teuchos::null)
    stepperXDotDot_ = state->getXDotDot();
  // Else use temporary storage stepperXDotDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperXDotDot_ == Teuchos::null, std::logic_error,
      "Error - stepperXDotDot_ has not been set in setInitialConditions() or\n"
      "        can not be set from the state!\n");

  return stepperXDotDot_;
}

template <class Scalar>
void Stepper<Scalar>::describe(Teuchos::FancyOStream& out,
                               const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "--- Stepper ---\n"
         << "  isInitialized_      = " << Teuchos::toString(isInitialized_)
         << std::endl
         << "  stepperType_        = " << stepperType_ << std::endl
         << "  useFSAL_            = " << Teuchos::toString(useFSAL_)
         << std::endl
         << "  ICConsistency_      = " << ICConsistency_ << std::endl
         << "  ICConsistencyCheck_ = " << Teuchos::toString(ICConsistencyCheck_)
         << std::endl
         << "  stepperX_           = " << stepperX_ << std::endl
         << "  stepperXDot_        = " << stepperXDot_ << std::endl
         << "  stepperXDotDot_     = " << stepperXDotDot_ << std::endl;
}

template <class Scalar>
bool Stepper<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (!(ICConsistency_ == "None" || ICConsistency_ == "Zero" ||
        ICConsistency_ == "App" || ICConsistency_ == "Consistent")) {
    isValidSetup = false;
    auto l_out   = Teuchos::fancyOStream(out.getOStream());
    l_out->setOutputToRootOnly(0);
    *l_out << "The IC consistency does not have a valid value!\n"
           << "('None', 'Zero', 'App' or 'Consistent')\n"
           << "  ICConsistency  = " << ICConsistency_ << "\n";
  }

  return isValidSetup;
}

template <class Scalar>
void Stepper<Scalar>::setStepperValues(Teuchos::RCP<Teuchos::ParameterList> pl)
{
  if (pl != Teuchos::null) {
    this->setStepperName(pl->name());
    auto stepperType =
        pl->get<std::string>("Stepper Type", this->getStepperType());
    TEUCHOS_TEST_FOR_EXCEPTION(
        stepperType != this->getStepperType(), std::runtime_error,
        "  ParameterList 'Stepper Type' (='"
            << stepperType
            << "')\n  does not match type for this Stepper (='"
            << this->getStepperType() << "').");
    this->setStepperType(stepperType);

    this->setUseFSAL(pl->get<bool>("Use FSAL", this->getUseFSAL()));
    this->setICConsistency(pl->get<std::string>("Initial Condition Consistency",
                                                this->getICConsistency()));
    this->setICConsistencyCheck(pl->get<bool>(
        "Initial Condition Consistency Check", this->getICConsistencyCheck()));
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList> Stepper<Scalar>::getValidParameters()
    const
{
  return this->getValidParametersBasic();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> Stepper<Scalar>::getValidParametersBasic()
    const
{
  auto pl = Teuchos::parameterList(this->getStepperName());
  pl->template set<std::string>("Stepper Type", this->getStepperType());

  pl->template set<bool>(
      "Use FSAL", this->getUseFSAL(),
      "The First-Same-As-Last (FSAL) principle is the situation where the\n"
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

  pl->template set<std::string>(
      "Initial Condition Consistency", this->getICConsistency(),
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

  pl->template set<bool>(
      "Initial Condition Consistency Check", this->getICConsistencyCheck(),
      "Check if the initial condition, x and xDot, is consistent with the\n"
      "governing equations, xDot = f(x,t), or f(x, xDot, t) = 0.\n"
      "\n"
      "In general for explicit and implicit steppers, the default is true,\n"
      "because it is fairly cheap with just one residual evaluation.\n"
      "Individual steppers may override this default.");

  return pl;
}

// Nonmember Helper Functions.
// ------------------------------------------------------------------------

template <class Scalar>
void validExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(model));
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<Scalar> inArgs   = model->createInArgs();
  const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
  const bool supports =
      inArgs.supports(MEB::IN_ARG_x) && outArgs.supports(MEB::OUT_ARG_f);

  TEUCHOS_TEST_FOR_EXCEPTION(
      supports == false, std::logic_error,
      model->description()
          << " can not support an explicit ODE with\n"
          << "  IN_ARG_x  = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
          << "  OUT_ARG_f = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
          << "Explicit ODE requires:\n"
          << "  IN_ARG_x  = true\n"
          << "  OUT_ARG_f = true\n"
          << "\n"
          << "NOTE: Currently the convention to evaluate f(x,t) is to set\n"
          << "xdot=null!  There is no InArgs support to test if xdot is null,\n"
          << "so we set xdot=null and hope the ModelEvaluator can handle "
             "it.\n");

  return;
}

template <class Scalar>
void validSecondOrderExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(model));
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<Scalar> inArgs   = model->createInArgs();
  const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
  const bool supports                = inArgs.supports(MEB::IN_ARG_x) &&
                        inArgs.supports(MEB::IN_ARG_x_dot) &&
                        outArgs.supports(MEB::OUT_ARG_f);

  TEUCHOS_TEST_FOR_EXCEPTION(
      supports == false, std::logic_error,
      model->description()
          << "can not support an explicit ODE with\n"
          << "  IN_ARG_x  = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
          << "  IN_ARG_x_dot  = " << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
          << "  OUT_ARG_f = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
          << "Explicit ODE requires:\n"
          << "  IN_ARG_x  = true\n"
          << "  IN_ARG_x_dot  = true\n"
          << "  OUT_ARG_f = true\n"
          << "\n"
          << "NOTE: Currently the convention to evaluate f(x, xdot, t) is to\n"
          << "set xdotdot=null!  There is no InArgs support to test if "
             "xdotdot\n"
          << "is null, so we set xdotdot=null and hope the ModelEvaluator can\n"
          << "handle it.\n");

  return;
}

template <class Scalar>
void validImplicitODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(model));
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<Scalar> inArgs   = model->createInArgs();
  const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
  const bool supports =
      inArgs.supports(MEB::IN_ARG_x) && inArgs.supports(MEB::IN_ARG_x_dot) &&
      inArgs.supports(MEB::IN_ARG_alpha) && inArgs.supports(MEB::IN_ARG_beta) &&
      !inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) &&
      outArgs.supports(MEB::OUT_ARG_f) && outArgs.supports(MEB::OUT_ARG_W);

  TEUCHOS_TEST_FOR_EXCEPTION(
      supports == false, std::logic_error,
      model->description() << " can not support an implicit ODE with\n"
                           << "  IN_ARG_x                 = "
                           << inArgs.supports(MEB::IN_ARG_x) << "\n"
                           << "  IN_ARG_x_dot             = "
                           << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
                           << "  IN_ARG_alpha             = "
                           << inArgs.supports(MEB::IN_ARG_alpha) << "\n"
                           << "  IN_ARG_beta              = "
                           << inArgs.supports(MEB::IN_ARG_beta) << "\n"
                           << "  IN_ARG_W_x_dot_dot_coeff = "
                           << inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff)
                           << "\n"
                           << "  OUT_ARG_f                = "
                           << outArgs.supports(MEB::OUT_ARG_f) << "\n"
                           << "  OUT_ARG_W                = "
                           << outArgs.supports(MEB::OUT_ARG_W) << "\n"
                           << "Implicit ODE requires:\n"
                           << "  IN_ARG_x                 = true\n"
                           << "  IN_ARG_x_dot             = true\n"
                           << "  IN_ARG_alpha             = true\n"
                           << "  IN_ARG_beta              = true\n"
                           << "  IN_ARG_W_x_dot_dot_coeff = false\n"
                           << "  OUT_ARG_f                = true\n"
                           << "  OUT_ARG_W                = true\n");

  return;
}

template <class Scalar>
void validSecondOrderODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(model));
  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<Scalar> inArgs   = model->createInArgs();
  const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
  const bool supports =
      inArgs.supports(MEB::IN_ARG_x) && inArgs.supports(MEB::IN_ARG_x_dot) &&
      inArgs.supports(MEB::IN_ARG_x_dot_dot) &&
      inArgs.supports(MEB::IN_ARG_alpha) && inArgs.supports(MEB::IN_ARG_beta) &&
      inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) &&
      outArgs.supports(MEB::OUT_ARG_f) && outArgs.supports(MEB::OUT_ARG_W);

  TEUCHOS_TEST_FOR_EXCEPTION(
      supports == false, std::logic_error,
      model->description() << " can not support an implicit ODE with\n"
                           << "  IN_ARG_x                 = "
                           << inArgs.supports(MEB::IN_ARG_x) << "\n"
                           << "  IN_ARG_x_dot             = "
                           << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
                           << "  IN_ARG_x_dot_dot         = "
                           << inArgs.supports(MEB::IN_ARG_x_dot_dot) << "\n"
                           << "  IN_ARG_alpha             = "
                           << inArgs.supports(MEB::IN_ARG_alpha) << "\n"
                           << "  IN_ARG_beta              = "
                           << inArgs.supports(MEB::IN_ARG_beta) << "\n"
                           << "  IN_ARG_W_x_dot_dot_coeff = "
                           << inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff)
                           << "\n"
                           << "  OUT_ARG_f                = "
                           << outArgs.supports(MEB::OUT_ARG_f) << "\n"
                           << "  OUT_ARG_W                = "
                           << outArgs.supports(MEB::OUT_ARG_W) << "\n"
                           << "Implicit Second Order ODE requires:\n"
                           << "  IN_ARG_x                 = true\n"
                           << "  IN_ARG_x_dot             = true\n"
                           << "  IN_ARG_x_dot_dot         = true\n"
                           << "  IN_ARG_alpha             = true\n"
                           << "  IN_ARG_beta              = true\n"
                           << "  IN_ARG_W_x_dot_dot_coeff = true\n"
                           << "  OUT_ARG_f                = true\n"
                           << "  OUT_ARG_W                = true\n");

  return;
}

Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters()
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  // NOX Solver ParameterList
  RCP<ParameterList> noxPL = Teuchos::parameterList();

  // Direction ParameterList
  RCP<ParameterList> directionPL = Teuchos::parameterList();
  directionPL->set<std::string>("Method", "Newton");
  RCP<ParameterList> newtonPL = Teuchos::parameterList();
  newtonPL->set<std::string>("Forcing Term Method", "Constant");
  newtonPL->set<bool>("Rescue Bad Newton Solve", 1);
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
  RCP<ParameterList> solverPL = Teuchos::parameterList("Default Solver");
  solverPL->set("NOX", *noxPL);

  return solverPL;
}

}  // namespace Tempus
#endif  // Tempus_Stepper_impl_hpp
