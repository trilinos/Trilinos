#ifndef TEMPUS_SOLUTIONSTATE_HPP
#define TEMPUS_SOLUTIONSTATE_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
// Thrya
#include "Thyra_VectorBase.hpp"
// Tempus
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_StepperState.hpp"

namespace Tempus {

/** \brief Solution state for integrators and steppers.
 *  SolutionState contains the metadata for solutions and the solutions
 *  themselves.
 *
 *  For simple time integration, the SolutionState is sufficient for
 *  checkpointing, restart and undo operations (i.e., it is the Memento
 *  object).
 *
 *  For more complex time integration where the physics has additional
 *  state information or the time integrator is not a one-step method
 *  (i.e., can not accurately start from a single time step), this class
 *  can be inherited and the physics state or additional time-integration
 *  parameters can be managed.
 *
 *  SolutionStates can be interpolated to generate solutions at various
 *  times (see SolutionHistory).  However not all metadata or state
 *  information can be interpolated.  Thus interpolated solutions may not
 *  be suitable for checkpointing, restart and undo operations, but may
 *  be useful for adjoint sensitivities.
 */
template<class Scalar>
class SolutionState :
  public Teuchos::Describable,
  public Teuchos::VerboseObject<Tempus::SolutionState<Scalar> >
{
public:

  /// Destructor
  virtual ~SolutionState() {};

  /** \brief Default constructor. */
  SolutionState();

  /** \brief. */
  SolutionState(
    const Teuchos::RCP<SolutionStateMetaData<Scalar> > ssmd,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState);

  SolutionState(
    const Scalar time,
    const Scalar dt,
    const int    iStep,
    const Scalar errorAbs,
    const Scalar errorRel,
    const int    order,
    const int    nFailures,
    const int    nConsecutiveFailures,
    const Status solutionStatus,
    const bool   output,
    const bool   isRestartable,
    const bool   isInterpolated,
    const Scalar accuracy,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState);

  /// This is a shallow copy constructor, use clone for a deep copy
  SolutionState(const SolutionState<Scalar>& ss_);

  /// This is a deep copy and copies the underlying vectors
  virtual Teuchos::RCP<SolutionState<Scalar> > clone() const;

  /// \name Accessor methods
  //@{
    /// Get time
    virtual Scalar getTime() const {return metaData->time;}

    /// Get index
    virtual Scalar getIndex() const {return metaData->iStep;}

    /// Get time step
    virtual Scalar getTimeStep() const {return metaData->dt;}

    /// Return the Solution status
    virtual Status getSolutionStatus() const
      {return metaData->solutionStatus;};

    /// Return the Stepper status
    virtual Status getStepperStatus() const
      {return stepperState->stepperStatus;};

    /// Get the current solution, x.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getX() {return x;}

    /// Get the current time derivative of the solution, xdot.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDot() {return xdot;}

    /// Get the current time second derivative of the solution, xdotdot.
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDotDot()
      {return xdotdot;}
  //@}


  /// \name Comparison methods
  //@{
    /// Less than comparison for sorting based on time
    bool operator< (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator<= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator< (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator<= (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator> (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator>= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator> (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator>= (const Scalar& t) const;

    /// Equality comparison for matching
    bool operator== (const SolutionState<Scalar>& ss) const;

    /// Equality comparison for matching
    bool operator== (const Scalar& t) const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream          &out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  // Member Data

  /// Meta Data for the solution state
  Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData;

  /// Solution
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x;

  /// Time derivative of the solution
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot;

  /// Second time derivative of the solution
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot;

  /// StepperState for this SolutionState
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState;

};
} // namespace Tempus

#include "Tempus_SolutionState_impl.hpp"

#endif // TEMPUS_SOLUTIONSTATE_HPP
