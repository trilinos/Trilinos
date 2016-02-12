#ifndef TEMPUS_SOLUTIONHISTORY_HPP
#define TEMPUS_SOLUTIONHISTORY_HPP

#include "Tempus_SolutionState.hpp"


namespace tempus {

enum HistoryPolicy {
  HISTORY_POLICY_INVALID     = 0,  ///< Invalid policy
  HISTORY_POLICY_KEEP_NEWEST = 1,  ///< Keep the single newest state
  HISTORY_POLICY_UNDO        = 2,  ///< Keep the 2 newest states for undo
  HISTORY_POLICY_STATIC      = 3,  ///< Keep a fix number of states
  HISTORY_POLICY_UNLIMITED   = 4,  ///< Grow the history as needed
};


/** \brief Caretaker class for the history of solution states.
 *
 */
template<class Scalar>
class SolutionHistory :
  virtual public InterpolatorAcceptingObjectBase<Scalar>
{
public:

  /** \brief. */
  SolutionHistory( RCP<Teuchos::ParameterList> paramList_ = Teuchos::null );

  /// Set the interpolator for this history
  void setInterpolator(const RCP<InterpolatorBase<Scalar> >& interpolator);

  /** \brief . */
  RCP<InterpolatorBase<Scalar> > getNonconstInterpolator();

  /** \brief . */
  RCP<const InterpolatorBase<Scalar> > getInterpolator() const;

  /// Unset the interpolator for this history
  RCP<InterpolatorBase<Scalar> > unSetInterpolator();


  /// Set the maximum storage of this history
  void setStorage( int storage );

  /// Get the maximum storage of this history
  int getStorage() const;

  /** \brief . */
  HistoryPolicy getHistoryPolicy();

  /// Destructor
  ~SolutionHistory() {};

  /// Add solution state to history
  void addState( const RCP<SolutionState<Scalar> >& state );

  /// Get solution state from history
  RCP<SolutionState<Scalar> > getState(const Scalar time) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /// Get solution state times in history
  void getTimes(Array<Scalar>* time_vec) const;

  /// Get order of interpolation
  int getInterpolationOrder() const;

  /// Remove solution state
  void removeStates(Array<Scalar>& time_vec);

  /// Redefined from Teuchos::Describable
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(Teuchos::FancyOStream          &out,
                const Teuchos::EVerbosityLevel verbLevel) const;

  /// Redefined from Teuchos::ParameterListAcceptor
  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();

  RCP<const Teuchos::ParameterList> getValidParameters() const;

private:

  RCP<Teuchos::ParameterList> paramList;
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
  RCP<Teuchos::ParameterList> paramList_ = Teuchos::null )
{
  RCP<SolutionHistory<Scalar> > sh=rcp(new SolutionHistory<Scalar>(paramList_));
  return sh;
}

} // namespace tempus
#endif // TEMPUS_SOLUTIONHISTORY_HPP
