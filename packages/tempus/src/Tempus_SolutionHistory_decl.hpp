//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_SolutionHistory_decl_hpp
#define Tempus_SolutionHistory_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_Interpolator.hpp"

namespace Tempus {

enum StorageType {
  STORAGE_TYPE_INVALID     = 0,  ///< Invalid storage type
  STORAGE_TYPE_KEEP_NEWEST = 1,  ///< Keep the single newest state
  STORAGE_TYPE_UNDO        = 2,  ///< Keep the 2 newest states for undo
  STORAGE_TYPE_STATIC      = 3,  ///< Keep a fix number of states
  STORAGE_TYPE_UNLIMITED   = 4,  ///< Grow the history as needed
};

/** \brief SolutionHistory is basically a container of SolutionStates.
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
 *
 *  There are two sets of indices associated with SolutionHistory:
 *  timestep index and the SolutionHistory index.  Users are interested
 *  in the timestep index as it indicates the solution increment, e.g.,
 *  \f$[x^0, ... , x^{n-2}, x^{n-1}, x^{n}]\f$, where \f$n\f$ is the timestep
 *  index.  While developers can also be interested in the SolutionHistory
 *  index to access the correct SolutionState in the SolutionHistory, e.g.,
 *  \f$[S^0, ... , S^{M-1}]\f$, where \f$M\f$ is the number of SolutionStates
 *  in the SolutionHistory (1 \f$\le M \le\f$ StorageLimit_).
 *
 *  SolutionHistory will hold upto the StorageLimit_ number of states.
 *  Rules of thumb for the <b>minimum</b> storage limit:
 *   - Explicit one-step methods can update the solution state in-place
 *     --> StorageLimit_ = 1.
 *   - Implicit one-step methods need two solution states --> StorageLimit_ = 2.
 *   - MultiStep methods require k past solution states --> StorageLimit_ = k+1.
 *
 *  The states contained in the SolutionHistory will often be the last
 *  \f$M\f$ states for forward time integration, e.g.,
 *  \f$[S^0(x^{n-M+1}), ... , S^{M-3}(x^{n-2}), S^{M-2}(x^{n-1}),
 *      S^{M-1}(x^{n})]\f$, but this is not guaranteed.
 *  For example, during transient adjoint sensitivity calculations, not
 *  all the forward states can be stored and therefore will be recalculated.
 *  Thus intermediate states are kept to reduce recomputation costs, i.e.,
 *  so one does not recalculate from the initial conditions to the desired time.
 *  This means the states in the SolutionHistory may not have consecutive
 *  timestep indices, e.g., \f$[..., S^{M-4}(x^{n-200}), S^{M-3}(x^{n-100}),
 *  S^{M-2}(x^{n-1}), S^{M-1}(x^{n})]\f$.
 *
 *  The SolutionHistory is kept in chronological order, starting with the
 *  oldest state and ending with the latest.
 *
 *  The "current" state is the latest solution that has been successfully
 *  solved.  The solution from the last successful time step or the initial
 *  conditions.
 *
 *  The "working" state is the state which is being worked on.  It is valid
 *  from initialization (e.g., copied from the "current" state), during the
 *  the timestep, and until it is accepted and promoted to the new "current"
 *  state.  Between the promotion and the initialization, the workingState_
 *  is invalid (i.e., Teuchos::null).  If the timestep fails, the
 *  workingState_ is maintained, and the timestep is retried until the
 *  Integrator declares a successful timestep or the time integration is a
 *  failure.  The SolutionHistory indices associated with the currentState_
 *  and the workingState_ vary during the time step loop, and are
 *  <table>
 *  <tr><th> Loop Portion <th> currentState_ <th> workingState_
 *  <tr><td> Before initializing working state
 *      <td> PASS -> \f$M\f$-1 <br>
 *           FAIL -> \f$M\f$-2
 *      <td> PASS -> Invalid <br>
 *           FAIL -> \f$M\f$-1
 *  <tr><td> After initializing working state <br>
 *           During timestep <br>
 *           Before accepting timestep
 *      <td> \f$M\f$-2
 *      <td> \f$M\f$-1
 *  <tr><td> After accepting timestep
 *      <td> PASS -> \f$M\f$-1 <br>
 *           FAIL -> \f$M\f$-2
 *      <td> PASS -> Invalid <br>
 *           FAIL -> \f$M\f$-1
 *  </table>
 *  Initial conditions are considered PASSing, which sets up the loop.
 *  The difference between the currentState_ timestep index and the
 *  workingState_ timestep index is guaranteed to be one (except for when
 *  StorageLimit_ = 1, i.e., explicit one-step methods that can update
 *  the solution in-place).
 */
template <class Scalar>
class SolutionHistory
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<SolutionHistory<Scalar> > {
 public:
  /// Default Contructor
  SolutionHistory();

  /// Contructor
  SolutionHistory(
      std::string name,
      Teuchos::RCP<std::vector<Teuchos::RCP<SolutionState<Scalar> > > > history,
      Teuchos::RCP<Interpolator<Scalar> > interpolator, StorageType storageType,
      int storageLimit);

  /// Destructor
  ~SolutionHistory() {}

  /// \name Basic SolutionHistory Methods
  //@{
  /// Add solution state to history
  void addState(const Teuchos::RCP<SolutionState<Scalar> >& state,
                bool doChecks = true);  // <-- Perform internal checks.

  /// Add a working solution state to history
  void addWorkingState(const Teuchos::RCP<SolutionState<Scalar> >& state,
                       const bool updateTime = true);

  /// Remove solution state
  void removeState(const Teuchos::RCP<SolutionState<Scalar> >& state);

  /// Remove solution state based on time
  void removeState(const Scalar time);

  /// Find solution state at requested time (no interpolation)
  Teuchos::RCP<SolutionState<Scalar> > findState(const Scalar time) const;

  /// Generate and interpolate a new solution state at requested time
  Teuchos::RCP<SolutionState<Scalar> > interpolateState(
      const Scalar time) const;

  /// Interpolate solution state at requested time and store in supplied state
  void interpolateState(const Scalar time,
                        SolutionState<Scalar>* state_out) const;

  /// Initialize the working state
  void initWorkingState();

  /// Promote the working state to current state
  void promoteWorkingState();

  /// Clear the history.
  void clear()
  {
    history_->clear();
    isInitialized_ = false;
  }

  /// Make a shallow copy of SolutionHistory (i.e., only RCPs to states and
  /// interpolator).
  void copy(Teuchos::RCP<const SolutionHistory<Scalar> > sh);
  //@}

  /// \name Accessor methods
  //@{
  /// Get this SolutionHistory's name
  std::string getName() const { return name_; }

  /// Set this SolutionHistory's name
  void setName(std::string name) { name_ = name; }

  /// Get underlining history
  Teuchos::RCP<std::vector<Teuchos::RCP<SolutionState<Scalar> > > > getHistory()
      const
  {
    return history_;
  }

  /// Set underlining history
  void setHistory(
      Teuchos::RCP<std::vector<Teuchos::RCP<SolutionState<Scalar> > > > h)
  {
    history_       = h;
    isInitialized_ = false;
  }

  /// Subscript operator
  Teuchos::RCP<SolutionState<Scalar> > operator[](const int i)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !((0 <= i) && (i < (int)history_->size())), std::out_of_range,
        "Error - SolutionHistory index is out of range.\n"
            << "    [Min, Max] = [ 0, " << history_->size() << "]\n"
            << "    index = " << i << "\n");
    return (*history_)[i];
  }

  /// Subscript operator (const version)
  Teuchos::RCP<const SolutionState<Scalar> > operator[](const int i) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        !((0 <= i) && (i < (int)history_->size())), std::out_of_range,
        "Error - SolutionHistory index is out of range.\n"
            << "    [Min, Max] = [ 0, " << history_->size() << "]\n"
            << "    index = " << i << "\n");
    return (*history_)[i];
  }

  /// Return the current state, i.e., the last accepted state
  Teuchos::RCP<SolutionState<Scalar> > getCurrentState() const
  {
    const int m = history_->size();
    Teuchos::RCP<SolutionState<Scalar> > state;
    if (m == 0)
      state = Teuchos::null;
    else if (m == 1 || workingState_ == Teuchos::null)
      state = (*history_)[m - 1];
    else if (m > 1)
      state = (*history_)[m - 2];
    return state;
  }

  /// Return the working state
  Teuchos::RCP<SolutionState<Scalar> > getWorkingState(bool warn = true) const
  {
    if (workingState_ == Teuchos::null && warn) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::getWorkingState()");
      *out << "Warning - WorkingState is null and likely has been promoted.  "
           << "You might want to call getCurrentState() instead.\n"
           << std::endl;
    }
    return workingState_;
  }

  /// Get the number of states
  int getNumStates() const { return history_->size(); }

  /// Get the current time
  Scalar getCurrentTime() const { return getCurrentState()->getTime(); }

  /// Get the current timestep index
  int getCurrentIndex() const { return getCurrentState()->getIndex(); }

  /// Set the maximum storage of this history
  void setStorageLimit(int storage_limit);

  /// Get the maximum storage of this history
  int getStorageLimit() const { return storageLimit_; }

  /// Set the storage type via enum.
  void setStorageType(StorageType st);

  /// Get the enum storage type.
  StorageType getStorageType() const { return storageType_; }

  /// Set the storage type via string.
  void setStorageTypeString(std::string st);

  /// Set the string storage type.
  std::string getStorageTypeString() const;

  /// Return the current minimum time of the SolutionStates
  Scalar minTime() const { return (history_->front())->getTime(); }

  /// Return the current maximum time of the SolutionStates
  Scalar maxTime() const { return (history_->back())->getTime(); }

  /** \brief Get the state with timestep index equal to n
   *
   *  If not available return Teuchos::null;
   */
  Teuchos::RCP<SolutionState<Scalar> > getStateTimeIndexN(
      bool warn = true) const;

  /** \brief Get the state with timestep index equal to n-1
   *
   *  If not available return Teuchos::null;
   */
  Teuchos::RCP<SolutionState<Scalar> > getStateTimeIndexNM1(
      bool warn = true) const;

  /** \brief Get the state with timestep index equal to n-2
   *
   *  If not available return Teuchos::null;
   */
  Teuchos::RCP<SolutionState<Scalar> > getStateTimeIndexNM2(
      bool warn = true) const;

  /** \brief Get the state with timestep index equal to "index"
   *
   *  If not available return Teuchos::null;
   */
  Teuchos::RCP<SolutionState<Scalar> > getStateTimeIndex(
      int index, bool warn = true) const;
  //@}

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// Return a valid non-const ParameterList with current settings.
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const;
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Interpolation Methods
  //@{
  /// Set the interpolator for this history
  void setInterpolator(const Teuchos::RCP<Interpolator<Scalar> >& interpolator);
  Teuchos::RCP<Interpolator<Scalar> > getNonconstInterpolator();
  Teuchos::RCP<const Interpolator<Scalar> > getInterpolator() const;
  /// Unset the interpolator for this history
  Teuchos::RCP<Interpolator<Scalar> > unSetInterpolator();
  //@}

  /// Print information on States in the SolutionHistory.
  void printHistory(std::string verb = "low") const;

  /** \brief Initialize SolutionHistory
   *
   *  This function will check if all member data is initialized
   *  and is consistent.  This function does not make member data
   *  consistent, but just checks it.  This ensures it is inexpensive.
   */
  void initialize() const;

  /// Return if SolutionHistory is initialized.
  bool isInitialized() { return isInitialized_; }

 protected:
  std::string name_;
  Teuchos::RCP<std::vector<Teuchos::RCP<SolutionState<Scalar> > > > history_;
  Teuchos::RCP<Interpolator<Scalar> > interpolator_;
  StorageType storageType_;
  int storageLimit_;

  Teuchos::RCP<SolutionState<Scalar> >
      workingState_;  ///< The state being worked on

  mutable bool isInitialized_;  ///< Bool if SolutionHistory is initialized.
};

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistory();

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryPL(
    Teuchos::RCP<Teuchos::ParameterList> pList);

/// Nonmember contructor from a SolutionState.
template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryState(
    const Teuchos::RCP<SolutionState<Scalar> >& state);

/// Nonmember contructor from a Thyra ModelEvaluator.
template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryME(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

}  // namespace Tempus

#endif  // Tempus_SolutionHistory_decl_hpp
