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
      : timeDer_(Teuchos::null), timeStepSize_(Scalar(0.0)),
        alpha_(Scalar(0.0)), beta_(Scalar(0.0)), evaluationType_(SOLVE_FOR_X),
        stageNumber_(0)
    {}
    /// Constructor
    ImplicitODEParameters(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
                          Scalar timeStepSize, Scalar alpha, Scalar beta,
                          EVALUATION_TYPE evaluationType = SOLVE_FOR_X,
                          int stageNumber = 0)
      : timeDer_(timeDer), timeStepSize_(timeStepSize),
        alpha_(alpha), beta_(beta), evaluationType_(evaluationType),
        stageNumber_(stageNumber)
    {}

    Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
    Scalar                                timeStepSize_;
    Scalar                                alpha_;
    Scalar                                beta_;
    EVALUATION_TYPE                       evaluationType_;
    int                                   stageNumber_;
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
    virtual Teuchos::RCP<const WrapperModelEvaluator<Scalar> >
      getWrapperModel(){return wrapperModel_;}

    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver = Teuchos::null);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return solver_; }

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

    /// Evaluate implicit ODE residual, f(x, xDot, t, p).
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
    virtual void setZeroInitialGuess(bool zIG) { zeroInitialGuess_ = zIG; }
    virtual bool getZeroInitialGuess() const { return zeroInitialGuess_; }

    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return Scalar(1.0e+99);}

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

  Teuchos::RCP<WrapperModelEvaluator<Scalar> >        wrapperModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >   solver_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;
  bool zeroInitialGuess_;

  Teuchos::RCP<StepperObserver<Scalar> >              stepperObserver_;

  // RCP to state or temporary storage if needed.
  Teuchos::RCP<Thyra::VectorBase<Scalar> >            stepperXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >            stepperXDotDot_;
};

} // namespace Tempus
#endif // Tempus_StepperImplicit_decl_hpp
