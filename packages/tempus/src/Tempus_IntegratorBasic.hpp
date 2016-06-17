#ifndef TEMPUS_INTEGRATORBASIC_HPP
#define TEMPUS_INTEGRATORBASIC_HPP

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
#include "Tempus_Stepper.hpp"
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
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<IntegratorObserver<Scalar> >&    observer=Teuchos::null);

  /// Destructor
  virtual ~IntegratorBasic() {}

  /// \name Basic integrator methods
  //@{
    /// Perform tasks before start of integrator.
    void startIntegrator();
    /// Advance the solution to timeMax, and return true if successful.
    bool advanceTime();
    /// Advance the solution to timeFinal, and return true if successful.
    bool advanceTime(const Scalar timeFinal);
    /// Only accept step after meeting time step criteria.
    void acceptTimeStep();
    /// Perform tasks after end of integrator.
    void endIntegrator();
  //@}

  /// \name Accessor methods
  //@{
    /// Get current time
    Scalar getTime() const {return solutionHistory_->getCurrentTime();}
    /// Get current index
    Scalar getIndex() const {return solutionHistory_->getCurrentIndex();}
    /// Get current the solution, x
    Teuchos::RCP<Thyra::VectorBase<double> > getX() const
      {return solutionHistory_->getCurrentState()->getX();}
    /// Get current the time derivative of the solution, xdot
    Teuchos::RCP<Thyra::VectorBase<double> > getXdot() const
      {return solutionHistory_->getCurrentState()->getXdot();}
    /// Get current the second time derivative of the solution, xdotdot
    Teuchos::RCP<Thyra::VectorBase<double> > getXdotdot() const
      {return solutionHistory_->getCurrentState()->getXdotdot();}

    /// Get SolutionHistory
    Teuchos::RCP<SolutionHistory<Scalar> > getSolutionHistory()
    { return solutionHistory_; }
    /// Get current state
    Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
    { return solutionHistory_->getCurrentState(); }
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

  Teuchos::RCP<Teuchos::ParameterList>      pList_;
  Teuchos::RCP<SolutionHistory<Scalar> >    solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> >    timeStepControl_;
  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver_;
  Teuchos::RCP<Stepper<Scalar> >            stepper_;

  Teuchos::RCP<Teuchos::Time>  integratorTimer;
  Teuchos::RCP<Teuchos::Time>  stepperTimer;

  std::vector<int>    outputScreenIndices; ///< Vector of screen output indices.

  /** The integratorStatus is primarily in the WORKING Status, and
   *  PASSED/FAILED are noted at the end of the run.  A FAILED value
   *  is used to jump out of the time-integration loop.
   */
  Status integratorStatus;
};
} // namespace Tempus

#include "Tempus_IntegratorBasic_impl.hpp"

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  const Teuchos::RCP<Tempus::IntegratorObserver<Scalar> >& ob=Teuchos::null)
{
  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>(pList, model, ob));
  return(integrator);
}


#endif // TEMPUS_INTEGRATORBASIC_HPP
