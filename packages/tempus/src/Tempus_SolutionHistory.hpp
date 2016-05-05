#ifndef TEMPUS_SOLUTIONHISTORY_HPP
#define TEMPUS_SOLUTIONHISTORY_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_Array.hpp"
// Tempus
#include "Tempus_SolutionState.hpp"


namespace Tempus {

enum HistoryPolicy {
  HISTORY_POLICY_INVALID     = 0,  ///< Invalid policy
  HISTORY_POLICY_KEEP_NEWEST = 1,  ///< Keep the single newest state
  HISTORY_POLICY_UNDO        = 2,  ///< Keep the 2 newest states for undo
  HISTORY_POLICY_STATIC      = 3,  ///< Keep a fix number of states
  HISTORY_POLICY_UNLIMITED   = 4,  ///< Grow the history as needed
};


/** \brief SolutionHistory is bascially a container of SolutionStates.
 *  SolutionHistory maintains a collection of SolutionStates for later
 *  retrival and reuse, such as checkpointing, restart, and undo
 *  operations.
 *
 *  The actual storage of the SolutionStates may take several forms:
 *   - in memory
 *   - on disk
 *   - combination
 *   - use ATDM DataWareHouse
 *  but the interface should be unchanged and very similar to other
 *  containers.
 *
 *  SolutionHistory can fill in SolutionStates between other SolutionStates
 *  by either
 *   - Integrating from one SolutionState to the desired time
 *     - This might include methods like Griewank's algorithm.
 *   - Interpolating between SolutionStates
 *     - Interpolated SolutionStates may not be suitable for adjoint
 *       solutions, restart, or undo operations (see SolutionState).
 */
template<class Scalar>
class SolutionHistory
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<SolutionHistory<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
//, virtual public InterpolatorAcceptingObjectBase<Scalar>
{
public:

  /// Contructor
  SolutionHistory(Teuchos::RCP<Teuchos::ParameterList> pList_ = Teuchos::null);

  /// Destructor
  ~SolutionHistory() {};

  /// \name Basic SolutionHistory Methods
  //@{
    /// Add solution state to history
    void addState(const Teuchos::RCP<SolutionState<Scalar> >& state);

    /// Remove solution state
    void removeState(const Teuchos::RCP<SolutionState<Scalar> >& state);

    /// Remove solution state based on time
    void removeState(const Scalar time);

    /// Find solution state at requested time (no interpolation)
    Teuchos::RCP<SolutionState<Scalar> > findState(const Scalar time) const;

    /// Generate and interpolate a new solution state at requested time
    Teuchos::RCP<SolutionState<Scalar> > interpolateState(const Scalar time) const;

    /// Return the current state, i.e., the last accepted state
    Teuchos::RCP<SolutionState<Scalar> > getCurrentState() const;

    /// Return the working state
    Teuchos::RCP<SolutionState<Scalar> > getWorkingState() const;

    /// Initialize the working state
    Teuchos::RCP<SolutionState<Scalar> > initWorkingState(const bool stepperStatus);

    /// Promote the working state to current state
    void promoteWorkingState();

  //@}

  /// \name Accessor methods
  //@{
    /// Set the maximum storage of this history
    void setStorage(int storage);

    /// Get the maximum storage of this history
    int getStorage() const;

    HistoryPolicy getHistoryPolicy();

    /// Return the current minimum time of the SolutionStates
    Scalar minTime() const;

    /// Return the current maximum time of the SolutionStates
    Scalar maxTime() const;
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

//  /// \name Interpolation Methods
//  //@{
//    /// Set the interpolator for this history
//    void setInterpolator(const Teuchos::RCP<InterpolatorBase<Scalar> >& interpolator);
//    Teuchos::RCP<InterpolatorBase<Scalar> > getNonconstInterpolator();
//    Teuchos::RCP<const InterpolatorBase<Scalar> > getInterpolator() const;
//    /// Unset the interpolator for this history
//    Teuchos::RCP<InterpolatorBase<Scalar> > unSetInterpolator();
//  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>                  pList;
  Teuchos::RCP<Teuchos::Array<Teuchos::RCP<SolutionState<Scalar> > > > history;
//  Teuchos::RCP<InterpolatorBase<Scalar> >      interpolator;
  HistoryPolicy                       historyPolicy;
  int                                 storage_limit;

  Teuchos::RCP<SolutionState<Scalar> > currentState;   ///< The last accepted state
  Teuchos::RCP<SolutionState<Scalar> > workingState;   ///< The state being worked on
};


/** \brief Nonmember constructor.
 *
 * \relates SolutionHistory.
 */
template<class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory(
  Teuchos::RCP<Teuchos::ParameterList> pList_ = Teuchos::null )
{
  Teuchos::RCP<SolutionHistory<Scalar> > sh=rcp(new SolutionHistory<Scalar>(pList_));
  return sh;
}

} // namespace Tempus

#include "Tempus_SolutionHistory_impl.hpp"

#endif // TEMPUS_SOLUTIONHISTORY_HPP
