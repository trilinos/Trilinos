//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Stepper_decl_hpp
#define Tempus_Stepper_decl_hpp

#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

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
template <class Scalar>
class Stepper : virtual public Teuchos::Describable,
                virtual public Teuchos::VerboseObject<Stepper<Scalar> > {
 public:
  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
  {
  }

  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const
  {
    return Teuchos::null;
  }

  /// Set solver.
  virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
  {
  }

  /// Get solver
  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
  {
    return Teuchos::null;
  }

  /// Initialize after construction and changing input parameters.
  virtual void initialize();

  /// True if stepper's member data is initialized.
  virtual bool isInitialized() { return isInitialized_; }

  /// Check initialization, and error out on failure.
  virtual void checkInitialized();

  /// Set initial conditions, make them consistent, and set stepper memory.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

  /// Pass initial guess to Newton solver (for implicit schemes)
  virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> >
                                   initialGuess = Teuchos::null) = 0;

  virtual Teuchos::RCP<Tempus::StepperState<Scalar> >
  getDefaultStepperState()           = 0;
  virtual Scalar getOrder() const    = 0;
  virtual Scalar getOrderMin() const = 0;
  virtual Scalar getOrderMax() const = 0;
  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const = 0;

  virtual bool isExplicit() const         = 0;
  virtual bool isImplicit() const         = 0;
  virtual bool isExplicitImplicit() const = 0;

  virtual bool isOneStepMethod() const   = 0;
  virtual bool isMultiStepMethod() const = 0;

  /// Set the stepper name.
  void setStepperName(std::string s)
  {
    stepperName_   = s;
    isInitialized_ = false;
  }

  /** \brief Get the stepper name.
   *
   * The stepper name is just a name used to distinguish it during
   * I/O and in ParameterLists, and can be anything the user would
   * like.  One example is when two steppers of the same type
   * (see getStepperType()) are being used during the same
   * simulation.  The user can name one as "Stepper with settings 1"
   * and the other as "Stepper with settings 2".  The default
   * name is the stepper type (e.g., "BDF2" or "Bogacki-Shampine 3(2) Pair").
   */
  std::string getStepperName() const { return stepperName_; }

 protected:
  /// Set the stepper type.
  void setStepperType(std::string s)
  {
    stepperType_   = s;
    isInitialized_ = false;
  }

 public:
  /** \brief Get the stepper type.
   *  The stepper type is used as an identifier for the stepper,
   *  and can only be set by the derived Stepper class.
   */
  std::string getStepperType() const { return stepperType_; }

  virtual void setUseFSAL(bool a) { setUseFSALFalseOnly(a); }
  void setUseFSALTrueOnly(bool a);
  void setUseFSALFalseOnly(bool a);
  bool getUseFSAL() const { return useFSAL_; }

  void setICConsistency(std::string s)
  {
    ICConsistency_ = s;
    isInitialized_ = false;
  }
  std::string getICConsistency() const { return ICConsistency_; }

  void setICConsistencyCheck(bool c)
  {
    ICConsistencyCheck_ = c;
    isInitialized_      = false;
  }
  bool getICConsistencyCheck() const { return ICConsistencyCheck_; }

  virtual OrderODE getOrderODE() const = 0;

  /// Get Stepper x.
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperX();

  /// Get Stepper xDot.
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDot();

  /// Get Stepper xDotDot.
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDotDot();

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
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

  /// Set Stepper member data from ParameterList.
  void setStepperValues(const Teuchos::RCP<Teuchos::ParameterList> pl);

  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// Add basic parameters to Steppers ParameterList.
  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasic() const;

 private:
  std::string stepperName_;  ///< Name used for output and ParameterLists
  std::string stepperType_;  ///< Name of stepper type
  std::string ICConsistency_ =
      std::string("None");  ///< Type of consistency to apply to ICs.
  bool ICConsistencyCheck_ =
      false;  ///< Check if the initial condition is consistent

  // RCP to SolutionState memory or Stepper temporary memory (if needed).
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperX_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > stepperXDotDot_;

 protected:
  /// Set x for Stepper storage.
  virtual void setStepperX(Teuchos::RCP<Thyra::VectorBase<Scalar> > x)
  {
    stepperX_ = x;
  }

  /// Set xDot for Stepper storage.
  virtual void setStepperXDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot)
  {
    stepperXDot_ = xDot;
  }

  /// Set x for Stepper storage.
  virtual void setStepperXDotDot(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot)
  {
    stepperXDotDot_ = xDotDot;
  }

  bool useFSAL_ = false;  ///< Use First-Same-As-Last (FSAL) principle
  bool isInitialized_ =
      false;  ///< True if stepper's member data is initialized.
};

/// \name Helper functions
//@{
/// Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
/** Currently the convention to evaluate f(x,t) is to set xdot=null!
 *  There is no InArgs support to test if xdot is null, so we set
 *  xdot=null and hopefully the ModelEvaluator can handle it.
 */
template <class Scalar>
void validExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

/// Validate that the model supports explicit second order ODE evaluation,
/// f(x,xdot,t) [=xdotdot]
/** Currently the convention to evaluate f(x,xdot,t) is to set xdotdot=null!
 *  There is no InArgs support to test if xdotdot is null, so we set
 *  xdotdot=null and hopefully the ModelEvaluator can handle it.
 */
template <class Scalar>
void validSecondOrderExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

/// Validate ME supports implicit ODE/DAE evaluation, f(xdot,x,t) [= 0]
template <class Scalar>
void validImplicitODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

/// Validate ME supports 2nd order implicit ODE/DAE evaluation,
/// f(xdotdot,xdot,x,t) [= 0]
template <class Scalar>
void validSecondOrderODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

/// Returns the default solver ParameterList for implicit Steppers.
Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters();
//@}

}  // namespace Tempus
#endif  // Tempus_Stepper_decl_hpp
