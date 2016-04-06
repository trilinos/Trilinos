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
    // This is a shallow copy constructor, use clone for a deep copy
    SolutionState(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x_,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot_,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot_,
      const Scalar time_,
      const Scalar dt_,
      const Scalar dtMin_,
      const Scalar dtMax_,
      const int    iStep_,
      const int    order_,
      const Scalar error_,
      const bool   isInterpolated_,
      const Scalar accuracy_);

    /** \brief. */
    // This is a shallow copy constructor, use clone for a deep copy
    SolutionState(const SolutionState<Scalar>& ss_);

    /** \brief. */
    // This is a deep clone and copies the underlying vectors
    virtual RCP<SolutionState<Scalar> > clone() const;

    Scalar time;            ///< Time of solution
    Scalar dt;              ///< Time step for this solution
    Scalar dtMin;           ///< Minimum allowed time step
    Scalar dtMax;           ///< Maximum allowed time step
    int    iStep;           ///< Time step index for this solution
    int    order;           ///< Order of this solution
    Scalar error;           ///< Local truncation error of this solution
    bool   isInterpolated;  ///< F - soln is time integrated; T - soln is interpolated
    bool   isRestartable;   ///< T - soln can be used as a restart
    Scalar accuracy;        ///< Interpolation accuracy of solution

    /// Solution at above time:
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x;

    /// Time derivative of the solution at above time:
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot;

    /// Second time derivative of the solution at above time:
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot;

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
