// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperStaggeredForwardSensitivity_decl_hpp
#define Tempus_StepperStaggeredForwardSensitivity_decl_hpp

#include "Tempus_Stepper.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"

namespace Tempus {

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
template<class Scalar>
class StepperStaggeredForwardSensitivity :
    virtual public Tempus::Stepper<Scalar>
{
public:

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
    const Teuchos::RCP<Teuchos::ParameterList>& pList = Teuchos::null,
    const Teuchos::RCP<Teuchos::ParameterList>& sens_pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel();

    virtual void setSolver(std::string solverName);
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
    { return stateStepper_->getSolver(); }

    /// Set Observer
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null){}

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual std::string getStepperType() const
     { return stepperPL_->get<std::string>("Stepper Type"); }

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >
      getDefaultStepperState();
    virtual Scalar getOrder() const {return stateStepper_->getOrder();}
    virtual Scalar getOrderMin() const {return stateStepper_->getOrderMin();}
    virtual Scalar getOrderMax() const {return stateStepper_->getOrderMax();}
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const
      {return std::numeric_limits<Scalar>::max();}

    virtual bool isExplicit()         const
      {return stateStepper_->isExplicit() or sensitivityStepper_->isExplicit();}
    virtual bool isImplicit()         const
      {return stateStepper_->isImplicit() or sensitivityStepper_->isImplicit();}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}

    virtual bool isOneStepMethod()   const
      {return stateStepper_->isOneStepMethod() and
              sensitivityStepper_->isOneStepMethod();}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
    
  //@}
    
  /// Pass initial guess to Newton solver (only relevant for explicit schemes)  
  virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
     {initial_guess_ = initial_guess;}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

private:

  /// Default Constructor -- not allowed
  StepperStaggeredForwardSensitivity();

  void setParams(const Teuchos::RCP<Teuchos::ParameterList> & pl,
                 const Teuchos::RCP<Teuchos::ParameterList> & spl);

private:

  Teuchos::RCP<Teuchos::ParameterList>               stepperPL_;
  Teuchos::RCP<Teuchos::ParameterList>               sensPL_;
  Teuchos::RCP<Stepper<Scalar> >                     stateStepper_;
  Teuchos::RCP<Stepper<Scalar> >                     sensitivityStepper_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> > combined_fsa_model_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> > fsa_model_;
  Teuchos::RCP<SolutionHistory<Scalar> > stateSolutionHistory_;
  Teuchos::RCP<SolutionHistory<Scalar> > sensSolutionHistory_;
  bool reuse_solver_;
  bool force_W_update_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;  
};

} // namespace Tempus

#endif // Tempus_StepperStaggeredForwardSensitivity_decl_hpp
