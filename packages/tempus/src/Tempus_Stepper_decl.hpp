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

    /// Set initial conditions, make them consistent, and set stepper memory.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;


    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    /// Pass initial guess to Newton solver (for implicit schemes)
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess = Teuchos::null) = 0;

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

    virtual OrderODE getOrderODE() const = 0;

    virtual void setUseFSAL(bool a) = 0;
    virtual bool getUseFSAL() const = 0;

    virtual void setICConsistency(std::string s) = 0;
    virtual std::string getICConsistency() const = 0;

    virtual void setICConsistencyCheck(bool c) = 0;
    virtual bool getICConsistencyCheck() const = 0;

    void getValidParametersBasic(Teuchos::RCP<Teuchos::ParameterList> pl) const;
  //@}

  virtual void modelWarning() const;

  /// \name Functions for Steppers with subSteppers (e.g., OperatorSplit)
  //@{
    virtual void createSubSteppers(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > /* models */){}
  //@}

  /// \name Helper functions
  //@{
    /// Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
    /** Currently the convention to evaluate f(x,t) is to set xdot=null!
     *  There is no InArgs support to test if xdot is null, so we set
     *  xdot=null and hopefully the ModelEvaluator can handle it.
     */
    void validExplicitODE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const;

    /// Validate that the model supports explicit second order ODE evaluation, f(x,xdot,t) [=xdotdot]
    /** Currently the convention to evaluate f(x,xdot,t) is to set xdotdot=null!
     *  There is no InArgs support to test if xdotdot is null, so we set
     *  xdotdot=null and hopefully the ModelEvaluator can handle it.
     */
    void validSecondOrderExplicitODE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const;

    /// Validate ME supports implicit ODE/DAE evaluation, f(xdot,x,t) [= 0]
    void validImplicitODE_DAE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const;

    /// Validate ME supports 2nd order implicit ODE/DAE evaluation, f(xdotdot,xdot,x,t) [= 0]
    void validSecondOrderODE_DAE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const;

    Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters() const;
  //@}

};
} // namespace Tempus
#endif // Tempus_Stepper_decl_hpp
