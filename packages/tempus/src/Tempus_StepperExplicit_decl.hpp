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
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return appModel_;}

    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return std::numeric_limits<Scalar>::max();}

    /// Set the initial conditions, make them consistent, and set needed memory.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Set solver via ParameterList solver name.
    virtual void setSolver(std::string solverName);
    /// Set solver via solver ParameterList.
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return Teuchos::null; }

    virtual std::string getStepperType() const
     { return stepperPL_->get<std::string>("Stepper Type"); }

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

    /// Set x for Stepper storage.
    virtual void setStepperX(Teuchos::RCP<Thyra::VectorBase<Scalar> > x)
      { stepperX_ = x; }
    /// Set xDot for Stepper storage.
    virtual void setStepperXDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot)
      { stepperXDot_ = xDot; }
    /// Set x for Stepper storage.
    virtual void setStepperXDotDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot)
      { stepperXDotDot_ = xDotDot; }

    /// Get x from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperX(
      Teuchos::RCP<SolutionState<Scalar> > state);
    /// Get xDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDot(
      Teuchos::RCP<SolutionState<Scalar> > state);
    /// Get xDotDot from SolutionState or Stepper storage.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getStepperXDotDot(
      Teuchos::RCP<SolutionState<Scalar> > state);

    /// Evaluate xDot = f(x,t).
    virtual void evaluateExplicitODE(
      Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      const Scalar time);

    /// Evaluate xDotDot = f(x, xDot, t).
    virtual void evaluateExplicitODE(
      Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDot,
      const Scalar time);
  //@}


protected:

  Teuchos::RCP<Teuchos::ParameterList>               stepperPL_;

  /// Explicit ODE ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>          inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         outArgs_;

  Teuchos::RCP<StepperObserver<Scalar> >             stepperObserver_;

  // RCP to state or temporary storage if needed.
  Teuchos::RCP<Thyra::VectorBase<Scalar> >           stepperX_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >           stepperXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >           stepperXDotDot_;
};

} // namespace Tempus
#endif // Tempus_StepperExplicit_decl_hpp
