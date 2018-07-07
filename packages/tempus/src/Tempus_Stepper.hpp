// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_Stepper_hpp
#define Tempus_Stepper_hpp

//Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Thyra
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperObserver.hpp"


namespace Tempus {


/** \brief Thyra Base interface for time steppers.
 *
 * <b>Design Considerations</b>
 *   - Time steppers are designed to take a single time step.
 *     - a single implicit solve for a time step
 *     - a single solve for a IMEX time step
 *   - Multiple time steps should be managed by Integrators.
 *   - Steppers can be built from other Sub-Steppers.
 *     - An operator-split Stepper is possible with interoperable Steppers.
 *   - For explicit steppers, only one ModelEvaluator and one solution
 *     vector are required.
 *   - For implicit steppers, only one ModelEvaluator, one solution
 *     vector, and one solver are required.
 *   - Steppers will PASS/FAIL the time step based on Solver, error and
 *     order requirements, and not adjust the time step size.
 *   - Steppers can provide a suggested time step size for the next time step.
 *   - For more complex steppers, multiple ModelEvaluators, solution
 *     vectors, and solvers are possible when a common single time-integration
 *     method is desired for all solutions. Examples:
 *     - Solution A with ModelEvaluator A and Solution B with ModelEvaluator B
 *       using the same solver
 *     - Solution A with ModelEvaluator A using Solver A and Solution B with
 *       ModelEvaluator B using Solver B
 *     - Solution A with ModelEvaluator A using Solver A and Solutions A and B
 *       with ModelEvaluator C using Solver B
 *   - Steppers may maintain their own time history of the solution, e.g.,
 *     BDF steppers.
 *
 * <b> CS Design Considerations</b>
 *   - All input parameters (i.e., ParameterList) can be set by public methods.
 *   - The Stepper ParameterList must be consistent.
 *     - The "set" methods which update parameters in the ParameterList
 *       must update the Stepper ParameterList.
 */
template<class Scalar>
class Stepper
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Stepper<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel) = 0;
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel) = 0;
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() = 0;

    /// Set solver via ParameterList solver name.
    virtual void setSolver(std::string solverName) = 0;
    /// Set solver via solver ParameterList.
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null) = 0;
    /// Set solver.
    virtual void setSolver(
        Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver) = 0;
    /// Get solver
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >
      getSolver() const = 0;

    /// Set Observer
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null) = 0;

    /// Initialize during construction and after changing input parameters.
    virtual void initialize() = 0;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    virtual std::string getStepperType() const = 0;
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >
      getDefaultStepperState() = 0;
    virtual Scalar getOrder() const = 0;
    virtual Scalar getOrderMin() const = 0;
    virtual Scalar getOrderMax() const = 0;
    virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const = 0;
    virtual Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const=0;

    virtual bool isExplicit() const = 0;
    virtual bool isImplicit() const = 0;
    virtual bool isExplicitImplicit() const = 0;

    virtual bool isOneStepMethod() const = 0;
    virtual bool isMultiStepMethod() const = 0;
  //@}

  /// \name Functions for Steppers with subSteppers (e.g., OperatorSplit)
  //@{
    virtual void createSubSteppers(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models){}
  //@}

  /// \name Helper functions
  //@{
    /// Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
    /** Currently the convention to evaluate f(x,t) is to set xdot=null!
     *  There is no InArgs support to test if xdot is null, so we set
     *  xdot=null and hopefully the ModelEvaluator can handle it.
     */
    void validExplicitODE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            outArgs.supports(MEB::OUT_ARG_f);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an explicit ODE with\n"
        << "  IN_ARG_x  = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
        << "  OUT_ARG_f = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
        << "Explicit ODE requires:\n"
        << "  IN_ARG_x  = true\n"
        << "  OUT_ARG_f = true\n"
        << "\n"
        << "NOTE: Currently the convention to evaluate f(x,t) is to set\n"
        << "xdot=null!  There is no InArgs support to test if xdot is null,\n"
        << "so we set xdot=null and hope the ModelEvaluator can handle it.\n");

      return;
    }

    /// Validate that the model supports explicit second order ODE evaluation, f(x,xdot,t) [=xdotdot]
    /** Currently the convention to evaluate f(x,xdot,t) is to set xdotdot=null!
     *  There is no InArgs support to test if xdotdot is null, so we set
     *  xdotdot=null and hopefully the ModelEvaluator can handle it.
     */
    void validSecondOrderExplicitODE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            inArgs.supports(MEB::IN_ARG_x_dot) and
                            outArgs.supports(MEB::OUT_ARG_f);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an explicit ODE with\n"
        << "  IN_ARG_x  = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
        << "  IN_ARG_x_dot  = " << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
        << "  OUT_ARG_f = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
        << "Explicit ODE requires:\n"
        << "  IN_ARG_x  = true\n"
        << "  IN_ARG_x_dot  = true\n"
        << "  OUT_ARG_f = true\n"
        << "\n"
        << "NOTE: Currently the convention to evaluate f(x, xdot, t) is to\n"
        << "set xdotdot=null!  There is no InArgs support to test if xdotdot\n"
        << "is null, so we set xdotdot=null and hope the ModelEvaluator can\n"
        << "handle it.\n");

      return;
    }


    /// Validate ME supports implicit ODE/DAE evaluation, f(xdot,x,t) [= 0]
    void validImplicitODE_DAE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            inArgs.supports(MEB::IN_ARG_x_dot) and
                            inArgs.supports(MEB::IN_ARG_alpha) and
                            inArgs.supports(MEB::IN_ARG_beta) and
                           !inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) and
                            outArgs.supports(MEB::OUT_ARG_f) and
                            outArgs.supports(MEB::OUT_ARG_W);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an implicit ODE with\n"
        << "  IN_ARG_x                 = "
        << inArgs.supports(MEB::IN_ARG_x) << "\n"
        << "  IN_ARG_x_dot             = "
        << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
        << "  IN_ARG_alpha             = "
        << inArgs.supports(MEB::IN_ARG_alpha) << "\n"
        << "  IN_ARG_beta              = "
        << inArgs.supports(MEB::IN_ARG_beta) << "\n"
        << "  IN_ARG_W_x_dot_dot_coeff = "
        << inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) << "\n"
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

    /// Validate ME supports 2nd order implicit ODE/DAE evaluation, f(xdotdot,xdot,x,t) [= 0]
    void validSecondOrderODE_DAE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            inArgs.supports(MEB::IN_ARG_x_dot) and
                            inArgs.supports(MEB::IN_ARG_x_dot_dot) and
                            inArgs.supports(MEB::IN_ARG_alpha) and
                            inArgs.supports(MEB::IN_ARG_beta) and
                            inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) and
                            outArgs.supports(MEB::OUT_ARG_f) and
                            outArgs.supports(MEB::OUT_ARG_W);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an implicit ODE with\n"
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
        << inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff) << "\n"
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

    Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters() const
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
  //@}

};
} // namespace Tempus
#endif // Tempus_Stepper_hpp
