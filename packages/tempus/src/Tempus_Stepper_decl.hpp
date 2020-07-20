// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_Stepper_decl_hpp
#define Tempus_Stepper_decl_hpp

//Teuchos
#include "Teuchos_TimeMonitor.hpp"
//
// Thyra
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperObserver.hpp"


namespace Tempus {

enum OrderODE {
  FIRST_ORDER_ODE  = 1,  ///< Stepper integrates first-order ODEs
  SECOND_ORDER_ODE = 2,  ///< Stepper integrates second-order ODEs
};

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
 */
template<class Scalar>
class Stepper
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Stepper<Scalar> >
{
public:

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel) {}

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& /* appModel */){}

#endif // TEMPUS_HIDE_DEPRECATED_CODE
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
    { return Teuchos::null; }

    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver) {}

    /// Get solver
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
    { return Teuchos::null; }

    /// Set Observer
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null){}

    /// Get Observer
    virtual Teuchos::RCP<StepperObserver<Scalar> >  getObserver() const
    { return Teuchos::null; }

    /// Initialize after construction and changing input parameters.
    virtual void initialize();

    /// True if stepper's member data is initialized.
    virtual bool isInitialized() { return isInitialized_; }

    /// Check initialization, and error out on failure.
    virtual void checkInitialized();

    /// Set initial conditions, make them consistent, and set stepper memory.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    /// Pass initial guess to Newton solver (for implicit schemes)
    virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> >
      initialGuess = Teuchos::null) = 0;

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >
      getDefaultStepperState() = 0;
    virtual Scalar getOrder() const = 0;
    virtual Scalar getOrderMin() const = 0;
    virtual Scalar getOrderMax() const = 0;
    virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const = 0;

    virtual bool isExplicit() const = 0;
    virtual bool isImplicit() const = 0;
    virtual bool isExplicitImplicit() const = 0;

    virtual bool isOneStepMethod() const = 0;
    virtual bool isMultiStepMethod() const = 0;

    void setStepperType(std::string s) { stepperType_ = s;
      isInitialized_ = false; }
    std::string getStepperType() const { return stepperType_; }

    void setUseFSAL(bool a) { useFSAL_ = a; isInitialized_ = false; }
    bool getUseFSAL() const { return useFSAL_; }
    virtual bool getUseFSALDefault() const { return false; }

    void setICConsistency(std::string s) { ICConsistency_ = s;
      isInitialized_ = false; }
    std::string getICConsistency() const { return ICConsistency_; }
    virtual std::string getICConsistencyDefault() const { return "None"; }

    void setICConsistencyCheck(bool c) {ICConsistencyCheck_ = c;
      isInitialized_ = false; }
    bool getICConsistencyCheck() const { return ICConsistencyCheck_; }
    virtual bool getICConsistencyCheckDefault() const { return false; }

    virtual OrderODE getOrderODE() const = 0;

    /// Get x from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperX(
      Teuchos::RCP<SolutionState<Scalar> > state);

    /// Get xDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDot(
      Teuchos::RCP<SolutionState<Scalar> > state);

    /// Get xDotDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDotDot(
      Teuchos::RCP<SolutionState<Scalar> > state);
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const { return stepperType_; }
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;

  /// \name Functions for Steppers with subSteppers (e.g., OperatorSplit)
  //@{
    virtual void createSubSteppers(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > /* models */){}
  //@}

  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const = 0;

private:

  std::string stepperType_;        ///< Name of stepper type
  bool useFSAL_ = false;           ///< Use First-Step-As-Last (FSAL) principle
  std::string ICConsistency_ = std::string("None");  ///< Type of consistency to apply to ICs.
  bool ICConsistencyCheck_ = true; ///< Check if the initial condition is consistent

  // RCP to SolutionState memory or Stepper temporary memory (if needed).
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperX_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperXDotDot_;

protected:

  /// Set x for Stepper storage.
  virtual void setStepperX(Teuchos::RCP<Thyra::VectorBase<Scalar> > x)
  { stepperX_ = x; }

  /// Set xDot for Stepper storage.
  virtual void setStepperXDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot)
  { stepperXDot_ = xDot; }

  /// Set x for Stepper storage.
  virtual void setStepperXDotDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot)
  { stepperXDotDot_ = xDotDot; }

  bool isInitialized_ = false; ///< True if stepper's member data is initialized.
};


/// \name Helper functions
//@{
  /// Provide basic parameters to Steppers.
  void getValidParametersBasic(
    Teuchos::RCP<Teuchos::ParameterList> pl, std::string stepperType);

  /// Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
  /** Currently the convention to evaluate f(x,t) is to set xdot=null!
   *  There is no InArgs support to test if xdot is null, so we set
   *  xdot=null and hopefully the ModelEvaluator can handle it.
   */
  template<class Scalar>
  void validExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
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
  template<class Scalar>
  void validSecondOrderExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
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
  template<class Scalar>
  void validImplicitODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
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
  template<class Scalar>
  void validSecondOrderODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
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


  /// Returns the default solver ParameterList for implicit Steppers.
  Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters();
//@}


} // namespace Tempus
#endif // Tempus_Stepper_decl_hpp
