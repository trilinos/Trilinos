#ifndef TEMPUS_SOLUTIONHISTORY_HPP
#define TEMPUS_SOLUTIONHISTORY_HPP

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
// Tempus
#include "Tempus_SolutionState.hpp"


namespace tempus {

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
    virtual public Teuchos::ParameterListAcceptor,
    virtual public InterpolatorAcceptingObjectBase<Scalar>
{
public:

  /// Contructor
  SolutionHistory( RCP<ParameterList> pList_ = Teuchos::null );

  /// Destructor
  ~SolutionHistory() {};

  /// \name Basic SolutionHistory Methods
  //@{
    /// Add solution state to history
    void addState( const RCP<SolutionState<Scalar> >& state );

    /// Find solution state at requested time (no interpolation)
    RCP<SolutionState<Scalar> > findState(const Scalar time) const;

    /// Generate and interpolate a new solution state at requested time
    RCP<SolutionState<Scalar> > interpolateState(const Scalar time) const;
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
    virtual void setParameterList(RCP<ParameterList> const& pl);
    virtual RCP<const ParameterList> getValidParameters() const;
    virtual RCP<const ParameterList> getParameterList() const;
    virtual RCP<ParameterList> getNonconstParameterList();
    virtual RCP<ParameterList> unsetParameterList();
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Interpolation Methods
  //@{
    /// Set the interpolator for this history
    void setInterpolator(const RCP<InterpolatorBase<Scalar> >& interpolator);
    RCP<InterpolatorBase<Scalar> > getNonconstInterpolator();
    RCP<const InterpolatorBase<Scalar> > getInterpolator() const;
    /// Unset the interpolator for this history
    RCP<InterpolatorBase<Scalar> > unSetInterpolator();
  //@}

private:

  RCP<ParameterList> pList;
  RCP<Array<SolutionState<Scalar> > > history;
  RCP<InterpolatorBase<Scalar> > interpolator;
  HistoryPolicy historyPolicy;
  int storage_limit;

};


/** \brief Nonmember constructor.
 *
 * \relates SolutionHistory.
 */
template<class Scalar>
RCP<SolutionHistory<Scalar> > solutionHistory(
  RCP<ParameterList> pList_ = Teuchos::null )
{
  RCP<SolutionHistory<Scalar> > sh=rcp(new SolutionHistory<Scalar>(pList_));
  return sh;
}

} // namespace tempus
#endif // TEMPUS_SOLUTIONHISTORY_HPP
