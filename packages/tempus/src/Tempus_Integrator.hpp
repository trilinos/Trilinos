#ifndef TEMPUS_INTEGRATOR_HPP
#define TEMPUS_INTEGRATOR_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_IntegratorObserver.hpp"
#include <string>

namespace tempus {

  /** \brief Thyra Base interface for time integrators.
   *  Time integrators are designed to advance the solution from an initial
   *  time, \f$t_0\f$, to a final time, \f$t_f\f$.
   *
   * <b>Design Considerations</b>
   *   - Integrators manage multiple time steps
   *   - Integrators have a single Stepper
   *   - Time-step ramping and startup are handled by TimeStepControl.
   *   - Solution output, e.g., solution plotting and check pointing,
   *     is coordinated in the Integrator.
   *   - Solution stability is handled in the timeStepControl, e.g., CFL
   *     constraint.
   *   - Error control over multiple time steps is handled in the Integrators,
   *     while error control over a single time step is handled in the Steppers.
   *   - Integrators will collect error control information from the Stepper
   *     and determine the next time step size and order.
   *   - Integrator maintains its own copy of the time history in the
   *     SolutionHistory, which may be just a single time step up to
   *     the entire solution history.
   *   - Integrators should compute the next time, and after accepting
   *     the time step advance the solution.  This allows a simple undo
   *     capability, if a solution is not acceptable.
   */
  template<class Scalar>
  class Integrator : public Teuchos::Describable,
                     public Teuchos::VerboseObject<tempus::Integrator<Scalar> >,
                     public Teuchos::ParameterListAcceptor {
  public:

    /// Destructor
    virtual ~Integrator();

    //! Unique name for this integrator.
    virtual std::string name() const = 0;

    /// \name Basic integrator methods
    //@{
    /// Advance the solution to time, and return true if successful.
    virtual bool advanceTime(const Scalar time_final) = 0;

    /// Get the current solution, x.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXs() = 0;

    /// Get the current time derivative of the solution, xdot.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXDots();

    /// Get the current second time derivative of the solution, xdotdot.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXDotDots();
    //@}

    /// \name Overridden from Teuchos::ParameterListAcceptor
    //@{
    virtual void setParameterList(RCP<ParameterList> const& pl);
    virtual RCP<const ParameterList> getValidParameters() const;
    virtual RCP<const ParameterList> getParameterList() const;
    virtual RCP<ParameterList> getNonconstParameterList();
    virtual RCP<ParameterList> unsetParameterList();
    //@}

    /// \name Accessor methods
    //@{
    virtual std::string description() const;
    virtual void describe( Teuchos::FancyOStream        & out,
                           const Teuchos::EVerbosityLevel verbLevel) const;
    /// Get time
    virtual Scalar getTime() const{return currentState->getTime();}
    /// Get index
    virtual Scalar getIndex() const{return currentState->getIndex();}
    //@}

    /// \name Undo type capabilities
    //@{
    /// Only accept step after meeting time step criteria.
    virtual bool acceptStep() {return false;}
    //@}

  protected:

    RCP<ParameterList>               parameterList;
    RCP<SolutionHistory<Scalar> >    solutionHistory;
    RCP<SolutionState<Scalar> >      currentState;
    RCP<TimeStepControl<Scalar> >    timeStepControl;
    RCP<IntegratorObserver<Scalar> > integratorObserver;
    RCP<Stepper<Scalar> >            stepper;

};
} // namespace tempus
#endif // TEMPUS_INTEGRATOR_HPP
