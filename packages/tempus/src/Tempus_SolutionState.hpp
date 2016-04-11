#ifndef TEMPUS_SOLUTIONSTATE_HPP
#define TEMPUS_SOLUTIONSTATE_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

namespace tempus {

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
  public Teuchos::VerboseObject<tempus::SolutionState<Scalar> >
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /// Destructor
    virtual ~SolutionState() {};

    /** \brief Default constructor. */
    SolutionState();

    /** \brief. */
    SolutionState(
      const Teuchos::RCP<SolutionStateMetaData<Scalar> > ssmd,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot);

    SolutionState(
      const Scalar time,
      const Scalar dt,
      const int    iStep,
      const Scalar errorAbs,
      const Scalar errorRel,
      const int    order,
      const int    nFailures,
      const int    nConsecutiveFailures,
      const SolutionStatus status,
      const bool   output,
      const bool   isAccepted,
      const bool   isRestartable,
      const bool   isInterpolated,
      const Scalar accuracy,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot);

    /** \brief. */
    // This is a shallow copy constructor, use clone for a deep copy
    SolutionState(const SolutionState<Scalar>& ss_);

    /** \brief. */
    // This is a deep clone and copies the underlying vectors
    virtual RCP<SolutionState<Scalar> > clone() const;

    /// Meta Data for the solution state
    Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData;

    /// Solution
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x;

    /// Time derivative of the solution
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot;

    /// Second time derivative of the solution
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot;

    /// Get time
    virtual Scalar getTime() const{return metaData->time;}

    /// Get time
    virtual Scalar getIndex() const{return metaData->iStep;}

    /// Less than comparison for sorting based on time:
    bool operator< (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time:
    bool operator<= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time:
    bool operator< (const Scalar& t) const;

    /// Less than comparison for sorting based on time:
    bool operator<= (const Scalar& t) const;

    /// Less than comparison for sorting based on time:
    bool operator> (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time:
    bool operator>= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time:
    bool operator> (const Scalar& t) const;

    /// Less than comparison for sorting based on time:
    bool operator>= (const Scalar& t) const;

    /// Equality comparison for matching:
    bool operator== (const SolutionState<Scalar>& ss) const;

    /// Equality comparison for matching:
    bool operator== (const Scalar& t) const;

    /// Inherited from Describable:
    /** \brief . */
    virtual std::string description() const;

    /** \brief . */
    /** \brief . */
    virtual void describe( Teuchos::FancyOStream          &out,
                           const Teuchos::EVerbosityLevel verbLevel) const;

};
} // namespace tempus
#endif // TEMPUS_SOLUTIONSTATE_HPP
