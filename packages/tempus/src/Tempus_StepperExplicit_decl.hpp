// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicit_decl_hpp
#define Tempus_StepperExplicit_decl_hpp

// Tempus
#include "Tempus_Stepper.hpp"


template<class Scalar>
class ExplicitODEParameters
{
  public:
    /// Constructor
    ExplicitODEParameters()
      : timeStepSize_(Scalar(0.0)), stageNumber_(0)
    {}
    /// Constructor
    ExplicitODEParameters(Scalar timeStepSize, int stageNumber = 0)
      : timeStepSize_(timeStepSize), stageNumber_(stageNumber)
    {}

    Scalar                                timeStepSize_;
    int                                   stageNumber_;
};


namespace Tempus {


/** \brief Thyra Base interface for implicit time steppers.
 *
 */
template<class Scalar>
class StepperExplicit : virtual public Tempus::Stepper<Scalar>
{
public:

  /// \name Basic explicit stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return appModel_;}

    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return std::numeric_limits<Scalar>::max();}

    /// Set the initial conditions, make them consistent, and set needed memory.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);

    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return Teuchos::null; }

    /// Pass initial guess to Newton solver (only relevant for implicit solvers)
    //  thus a no-op for explicit steppers.
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > /* initial_guess */){}

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()    const {return true;}
    virtual bool isMultiStepMethod()  const {return !isOneStepMethod();}

    /// Evaluate xDot = f(x,t).
    virtual void evaluateExplicitODE(
      Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      const Scalar time,
      const Teuchos::RCP<ExplicitODEParameters<Scalar> > & p );

    /// Evaluate xDotDot = f(x, xDot, t).
    virtual void evaluateExplicitODE(
      Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDot,
      const Scalar time,
      const Teuchos::RCP<ExplicitODEParameters<Scalar> > & p );
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;

protected:

  /// Explicit ODE ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>          inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         outArgs_;

  Teuchos::RCP<StepperObserver<Scalar> >             stepperObserver_;
 };

} // namespace Tempus
#endif // Tempus_StepperExplicit_decl_hpp
