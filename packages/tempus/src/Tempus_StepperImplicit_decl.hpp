// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperImplicit_decl_hpp
#define Tempus_StepperImplicit_decl_hpp

// Tempus
#include "Tempus_Stepper.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"


namespace Tempus {


template<class Scalar>
class ImplicitODEParameters
{
  public:
    /// Constructor
    ImplicitODEParameters()
      : timeDer_(Teuchos::null), timeStepSize_(Scalar(0.0)), stageNumber_(0),
        alpha_(Scalar(0.0)), beta_(Scalar(0.0)), evaluationType_(SOLVE_FOR_X)
    {}
    /// Constructor
    ImplicitODEParameters(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
                          Scalar timeStepSize, Scalar alpha, Scalar beta,
                          EVALUATION_TYPE evaluationType = SOLVE_FOR_X)
      : timeDer_(timeDer), timeStepSize_(timeStepSize), stageNumber_(0),
        alpha_(alpha), beta_(beta), evaluationType_(evaluationType)
    {}

    Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
    Scalar                                timeStepSize_;
    int                                   stageNumber_;
    Scalar                                alpha_;
    Scalar                                beta_;
    EVALUATION_TYPE                       evaluationType_;
};

/** \brief Thyra Base interface for implicit time steppers.
 *
 */
template<class Scalar>
class StepperImplicit : virtual public Tempus::Stepper<Scalar>
{
public:

  /// \name Basic implicit stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return wrapperModel_->getAppModel();}

    /// Set solver via ParameterList solver name.
    virtual void setSolver(std::string solverName);
    /// Set solver via solver ParameterList.
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return solver_; }

    virtual std::string getStepperType() const
     { return stepperPL_->get<std::string>("Stepper Type"); }

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Return alpha = d(xDot)/dx.
    virtual Scalar getAlpha(const Scalar dt) const = 0;
    /// Return beta  = d(x)/dx.
    virtual Scalar getBeta (const Scalar dt) const = 0;

    /// Solve problem using x in-place.  (Needs to be deprecated!)
    const Thyra::SolveStatus<Scalar> solveImplicitODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x);

    /// Solve implicit ODE, f(x, xDot, t, p) = 0.
    const Thyra::SolveStatus<Scalar> solveImplicitODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
      const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p );

    /// Evaluate implicit ODE, f(x, xDot, t, p), residual.
    void evaluateImplicitODE(
            Teuchos::RCP<Thyra::VectorBase<Scalar> > & f,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
      const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p );

    /// Pass initial guess to Newton solver (only relevant for implicit solvers)
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
      {initial_guess_ = initial_guess;}

    /// Set parameter so that the initial guess is set to zero (=True) or use last timestep (=False).
    virtual void setZeroInitialGuess(bool zIG)
      { stepperPL_->set<bool>("Zero Initial Guess", zIG); }
    virtual bool getZeroInitialGuess() const
      { return stepperPL_->get<bool>("Zero Initial Guess", false); }
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return Scalar(1.0e+99);}

    virtual bool getEmbedded() const
      { return stepperPL_->get<bool>("Use Embedded", false); }

    virtual void setUseFSAL(bool a) {stepperPL_->set<bool>("Use FSAL", a);}
    virtual bool getUseFSAL() const
      {
        bool defaultUseFSAL =
          this->getDefaultParameters()->template get<bool>("Use FSAL");
        return stepperPL_->get<bool>("Use FSAL", defaultUseFSAL);
      }

    virtual void setICConsistency(std::string s)
      {stepperPL_->set<std::string>("Initial Condition Consistency", s);}
    virtual std::string getICConsistency() const
      {
        std::string defaultICConsistency = this->getDefaultParameters()->
          template get<std::string>("Initial Condition Consistency");
        return stepperPL_->get<std::string>("Initial Condition Consistency",
                                            defaultICConsistency);
      }

    virtual void setICConsistencyCheck(bool c)
      {stepperPL_->set<bool>("Initial Condition Consistency Check", c);}
    virtual bool getICConsistencyCheck() const
      {
        bool defaultICConsistencyCheck = this->getDefaultParameters()->
          template get<bool>("Initial Condition Consistency Check");
        return stepperPL_->get<bool>("Initial Condition Consistency Check",
                                     defaultICConsistencyCheck);
      }

    /// Set xDot for Stepper storage.
    virtual void setStepperXDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot)
      { stepperXDot_ = xDot; }

    /// Get xDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDot(
      Teuchos::RCP<SolutionState<Scalar> > state);

    /// Get xDotDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDotDot(
      Teuchos::RCP<SolutionState<Scalar> > state);
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>                stepperPL_;
  Teuchos::RCP<WrapperModelEvaluator<Scalar> >        wrapperModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >   solver_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;

  Teuchos::RCP<StepperObserver<Scalar> >              stepperObserver_;

  // RCP to state or temporary storage if needed.
  Teuchos::RCP<Thyra::VectorBase<Scalar> >            stepperXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >            stepperXDotDot_;
};

} // namespace Tempus
#endif // Tempus_StepperImplicit_decl_hpp
