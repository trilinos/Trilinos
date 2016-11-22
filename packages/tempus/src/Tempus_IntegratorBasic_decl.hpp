#ifndef Tempus_IntegratorBasic_decl_hpp
#define Tempus_IntegratorBasic_decl_hpp

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"
// Tempus
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
// Tempus
#include "Tempus_Integrator.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_IntegratorObserver.hpp"

#include <string>

namespace Tempus {


/** \brief Basic time integrator
 */
template<class Scalar>
class IntegratorBasic : virtual public Tempus::Integrator<Scalar>
{
public:

  /** \brief Constructor with ParameterList, model and optional observer. */
  IntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /// Destructor
  virtual ~IntegratorBasic() {}

  /// \name Basic integrator methods
  //@{
    /// Advance the solution to timeMax, and return true if successful.
    virtual bool advanceTime();
    /// Advance the solution to timeFinal, and return true if successful.
    virtual bool advanceTime(const Scalar timeFinal);
    /// Perform tasks before start of integrator.
    virtual void startIntegrator();
    /// Start time step.
    virtual void startTimeStep();
    /// Only accept step after meeting time step criteria.
    virtual void acceptTimeStep();
    /// Perform tasks after end of integrator.
    virtual void endIntegrator();
  //@}

  /// \name Accessor methods
  //@{
    /// Get current time
    virtual Scalar getTime() const {return solutionHistory_->getCurrentTime();}
    /// Get current index
    virtual Scalar getIndex() const {return solutionHistory_->getCurrentIndex();}
    /// Get the Stepper
    virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const
      {return stepper_;}
    /// Set the Stepper
    virtual void setStepper(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);
    /// Set the Stepper
    virtual void setStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
    /// Get the SolutionHistory
    virtual Teuchos::RCP<SolutionHistory<Scalar> > getSolutionHistory()
      {return solutionHistory_;}
    /// Set the SolutionHistory
    virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar> > sh = Teuchos::null);
    /// Get the TimeStepControl
    virtual Teuchos::RCP<TimeStepControl<Scalar> > getTimeStepControl()
      {return timeStepControl_;}
    /// Set the TimeStepControl
    virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar> > tsc = Teuchos::null);
    /// Get the Observer
    virtual Teuchos::RCP<IntegratorObserver<Scalar> > getObserver()
      {return integratorObserver_;}
    /// Set the Observer
    virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
    /// Initializes the Integrator after set* function calls
    virtual void initialize();


    /// Get current the solution, x
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getX() const
      {return solutionHistory_->getCurrentState()->getX();}
    /// Get current the time derivative of the solution, xdot
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getXdot() const
      {return solutionHistory_->getCurrentState()->getXDot();}
    /// Get current the second time derivative of the solution, xdotdot
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getXdotdot() const
      {return solutionHistory_->getCurrentState()->getXDotDot();}

    /// Get current state
    virtual Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
      {return solutionHistory_->getCurrentState();}
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const;
    void describe(Teuchos::FancyOStream        & out,
                  const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>      tempusPL_;
  Teuchos::RCP<Teuchos::ParameterList>      integratorPL_;
  Teuchos::RCP<SolutionHistory<Scalar> >    solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> >    timeStepControl_;
  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver_;
  Teuchos::RCP<Stepper<Scalar> >            stepper_;

  Teuchos::RCP<Teuchos::Time>  integratorTimer_;
  Teuchos::RCP<Teuchos::Time>  stepperTimer_;

  std::vector<int>    outputScreenIndices_;///< Vector of screen output indices.

  /** The integratorStatus is primarily in the WORKING Status, and
   *  PASSED/FAILED are noted at the end of the run.  A FAILED value
   *  is used to jump out of the time-integration loop.
   */
  Status integratorStatus_;
};

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model);

} // namespace Tempus

#endif // Tempus_IntegratorBasic_decl_hpp
