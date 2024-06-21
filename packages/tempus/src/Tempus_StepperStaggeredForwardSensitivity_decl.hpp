//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperStaggeredForwardSensitivity_decl_hpp
#define Tempus_StepperStaggeredForwardSensitivity_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"

namespace Tempus {

enum class SensitivityStepMode { Forward,
                                 Sensitivity,
                                 Combined,
                                 Adjoint };

/** \brief A stepper implementing staggered forward sensitivity analysis.
 */
/**
 * It constructs two internal steppers, one for the state equations as usual
 * and one for the sensitivity equations using
 * Tempus::StaggeredForwardSensitivityModelEvaluator.  It's implementation
 * of takeStep() first takes a step using the state stepper, updates the
 * sensitivity model evaluator with the compute state solution and time
 * derivatives, and then takes a step using the sensitivity stepper. It
 * optionally can reuse the state solver for the sensitivity equations as well.
 */
template <class Scalar>
class StepperStaggeredForwardSensitivity
  : virtual public Tempus::Stepper<Scalar>,
    virtual public Teuchos::ParameterListAcceptor {
 public:
  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
   */
  StepperStaggeredForwardSensitivity();

  /// Constructor
  /*!
   * The first parameter list argument supplies supplies regular stepper
   * options, while the second provides sensitivity specific options:
   * <ul>
   *   <li> "Reuse State Linear Solver" (default: false) Whether to reuse the
   *        model's W matrix, solver, and preconditioner when solving the
   *        sensitivity equations.  If they can be reused, substantial savings
   *        in compute time are possible.
   *   <li> "Force W Update" (default: false) When reusing the solver as above
   *        whether to force recomputation of W.  This can be necessary when
   *        the solver overwrites it during the solve-phase (e.g., by a
   *        factorization).
   *   <li> "Use DfDp as Tangent" (default:  false) Reinterpret the df/dp
   *        out-arg as the tangent vector (df/dx)(x,p) * dx/dp + df/dp(x,p)
   *        as described in the Tempus::CombinedForwardSensitivityModelEvaluator
   *        documentation.
   *   <li> "Sensitivity Parameter Index" (default: 0) Model evaluator
   *        parameter index for which sensitivities will be computed.
   *   <li> "Sensitivity X Tangent Index" (default: 1) If "Use DfDp as Tangent"
   *        is true, the model evaluator parameter index for passing dx/dp
   *        as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot Tangent Index" (default: 2) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot/dp as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot-Dot Tangent Index" (default: 3) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot_dot/dp as a
   *        Thyra::DefaultMultiVectorProductVector (if the model supports
   *        x_dot_dot).
   * </ul>
   */
  StepperStaggeredForwardSensitivity(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_residual_model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_solve_model,
      const Teuchos::RCP<Teuchos::ParameterList>& pList      = Teuchos::null,
      const Teuchos::RCP<Teuchos::ParameterList>& sens_pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
  {
    setModel(appModel, appModel, appModel);
  }
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_residual_model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_solve_model);
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

  virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver = Teuchos::null);
  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
  {
    return stateStepper_->getSolver();
  }

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */)
  {
  }

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const { return stateStepper_->getOrder(); }
  virtual Scalar getOrderMin() const { return stateStepper_->getOrderMin(); }
  virtual Scalar getOrderMax() const { return stateStepper_->getOrderMax(); }
  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
  {
    return Scalar(1.0e+99);
  }

  virtual bool isExplicit() const
  {
    return stateStepper_->isExplicit() || sensitivityStepper_->isExplicit();
  }
  virtual bool isImplicit() const
  {
    return stateStepper_->isImplicit() || sensitivityStepper_->isImplicit();
  }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }

  virtual bool isOneStepMethod() const
  {
    return stateStepper_->isOneStepMethod() &&
           sensitivityStepper_->isOneStepMethod();
  }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }

  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }

  virtual void setUseFSAL(bool a) { stepperPL_->set<bool>("Use FSAL", a); }
  virtual bool getUseFSAL() const
  {
    return stepperPL_->get<bool>("Use FSAL", false);
  }

  virtual void setICConsistency(std::string s)
  {
    stepperPL_->set<std::string>("Initial Condition Consistency", s);
  }
  virtual std::string getICConsistency() const
  {
    return stepperPL_->get<std::string>("Initial Condition Consistency",
                                        "None");
  }

  virtual void setICConsistencyCheck(bool c)
  {
    stepperPL_->set<bool>("Initial Condition Consistency Check", c);
  }
  virtual bool getICConsistencyCheck() const
  {
    return stepperPL_->get<bool>("Initial Condition Consistency Check", false);
  }
  //@}

  /// Pass initial guess to Newton solver
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > /* initial_guess */)
  {
  }

  /// \name ParameterList methods
  //@{
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& pl);
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const
  {
    return "StepperStaggeredForwardSensitivity";
  }
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const
  {
    return stepperPL_;
  }

  //! What mode the current time integration step is in
  SensitivityStepMode getStepMode() const { return stepMode_; }

 private:
  void setParams(const Teuchos::RCP<Teuchos::ParameterList>& pl,
                 const Teuchos::RCP<Teuchos::ParameterList>& spl);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> stepperPL_;
  Teuchos::RCP<Teuchos::ParameterList> sensPL_;
  Teuchos::RCP<Stepper<Scalar> > stateStepper_;
  Teuchos::RCP<Stepper<Scalar> > sensitivityStepper_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> > combined_fsa_model_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> > fsa_model_;
  Teuchos::RCP<SolutionHistory<Scalar> > stateSolutionHistory_;
  Teuchos::RCP<SolutionHistory<Scalar> > sensSolutionHistory_;
  bool reuse_solver_;
  bool force_W_update_;
  SensitivityStepMode stepMode_;
};

}  // namespace Tempus

#endif  // Tempus_StepperStaggeredForwardSensitivity_decl_hpp
