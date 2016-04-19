#ifndef TEMPUS_INTEGRATORBASIC_HPP
#define TEMPUS_INTEGRATORBASIC_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
// Tempus
#include "Thyra_ModelEvaluator.hpp"
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

  /** \brief Constructor with ParameterList, models, initial conditions
   *  and optional solvers. */
  IntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                pList_,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >&     x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >&     xdot=Teuchos::null,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >&     xdotdot=Teuchos::null );

  /// Destructor
  virtual ~IntegratorBasic() {}

  //! Unique name for this integrator.
  virtual std::string name() const = 0;

  /// \name Basic integrator methods
  //@{
    /// Advance the solution to time, and return true if successful.
    virtual bool advanceTime(const Scalar time) = 0;
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pl);
    virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}
  /// \name Accessor methods
  //@{
    /// Get time
    virtual Scalar getTime() const{return workingState->getTime();}
    /// Get index
    virtual Scalar getIndex() const{return workingState->getIndex();}
  //@}

  /// \name Undo type capabilities
  //@{
    /// Only accept step after meeting time step criteria.
    virtual bool acceptTimeStep();
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>      pList;
  Teuchos::RCP<SolutionHistory<Scalar> >    solutionHistory;
  Teuchos::RCP<TimeStepControl<Scalar> >    timeStepControl;
  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver;
  Teuchos::RCP<Stepper<Scalar> >            stepper;

  Teuchos::RCP<SolutionState<Scalar> >      currentState; ///< The last accepted state
  Teuchos::RCP<SolutionState<Scalar> >      workingState; ///< The state being worked on

};
} // namespace Tempus
#endif // TEMPUS_INTEGRATORBASIC_HPP
