//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeStepControl_decl_hpp
#define Tempus_TimeStepControl_decl_hpp

#include <iostream>
#include <iterator>
#include <sstream>

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_TimeEventComposite.hpp"

namespace Tempus {

/** \brief TimeStepControl manages the time step size.
 *  There several mechanisms that effect the time step size and
 *  handled with this class:
 *   - Maximum and minimum time
 *   - Maximum and minimum time index
 *   - Maximum and minimum time step size
 *   - Maximum and minimum error
 *   - Startup considerations (e.g., ramping)
 *   - Solution and/or diagnostic output
 *  Additional step control can be added through the step control observer,
 *  or inheriting from this class.
 *   - Stability limits (e.g., CFL number)
 *
 *  Using TimeStepControlStrategy allows applications to define their
 *  very own strategy used to determine the next time step size
 *  (`setNextTimeStep()`).  Applications can define multiple strategies
 *  and add it to a vector of strategies TimeStepControlStrategyComposite
 *  using setTimeStepControlStrategy().  TimeStepControlStrategyComposite
 *  iterates over the list of strategies to determine the "optimal"
 *  next time step size.
 *
 */
template <class Scalar>
class TimeStepControl
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Tempus::TimeStepControl<Scalar>> {
 public:
  /// Default Constructor
  TimeStepControl();

  /// Constructor
  TimeStepControl(
      Scalar initTime, Scalar finalTime, Scalar minTimeStep,
      Scalar initTimeStep, Scalar maxTimeStep, int initIndex, int finalIndex,
      Scalar maxAbsError, Scalar maxRelError, int maxFailures,
      int maxConsecFailures, int numTimeSteps, bool printDtChanges,
      bool outputExactly, std::vector<int> outputIndices,
      std::vector<Scalar> outputTimes, int outputIndexInterval,
      Scalar outputTimeInterval,
      Teuchos::RCP<TimeEventComposite<Scalar>> timeEvent,
      Teuchos::RCP<TimeStepControlStrategy<Scalar>> stepControlStrategy);

  /// Destructor
  virtual ~TimeStepControl() {}

  /** \brief Determine the time step size.*/
  virtual void setNextTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory,
      Status& integratorStatus);

  /** \brief Check if time is within minimum and maximum time. */
  virtual bool timeInRange(const Scalar time) const;

  /** \brief Check if time step index is within minimum and maximum index. */
  virtual bool indexInRange(const int iStep) const;

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const;
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Get accessors
  //@{
  virtual std::string getStepType() const;
  virtual Scalar getInitTime() const { return initTime_; }
  virtual Scalar getFinalTime() const { return finalTime_; }
  virtual Scalar getMinTimeStep() const { return minTimeStep_; }
  virtual Scalar getInitTimeStep() const { return initTimeStep_; }
  virtual Scalar getMaxTimeStep() const { return maxTimeStep_; }
  virtual int getInitIndex() const { return initIndex_; }
  virtual int getFinalIndex() const { return finalIndex_; }
  virtual Scalar getMaxAbsError() const { return maxAbsError_; }
  virtual Scalar getMaxRelError() const { return maxRelError_; }

  /** \brief Return if the output needs to exactly happen on output time.
   *
   *  Output TimeEvents, "Output Time Interval" and "Output Time List",
   *  will always have the same OutputExactly value, so just
   *  find one or the other.  Otherwise return the default of true.
   *
   *  \return True if output needs to happen exactly on output times.
   */
  virtual bool getOutputExactly() const;

  virtual std::vector<int> getOutputIndices() const;
  virtual std::vector<Scalar> getOutputTimes() const;
  virtual int getOutputIndexInterval() const;
  virtual Scalar getOutputTimeInterval() const;
  virtual int getMaxFailures() const { return maxFailures_; }
  virtual int getMaxConsecFailures() const { return maxConsecFailures_; }
  virtual bool getPrintDtChanges() const { return printDtChanges_; }
  virtual int getNumTimeSteps() const { return numTimeSteps_; }
  virtual Teuchos::RCP<TimeEventComposite<Scalar>> getTimeEvents() const
  {
    return timeEvent_;
  }
  virtual Teuchos::RCP<TimeStepControlStrategy<Scalar>>
  getTimeStepControlStrategy() const
  {
    return stepControlStrategy_;
  }
  //@}

  /// \name Set accessors
  //@{
  virtual void setInitTime(Scalar t)
  {
    initTime_      = t;
    isInitialized_ = false;
  }
  virtual void setFinalTime(Scalar t)
  {
    finalTime_     = t;
    isInitialized_ = false;
  }
  virtual void setMinTimeStep(Scalar t)
  {
    minTimeStep_   = t;
    isInitialized_ = false;
  }
  virtual void setInitTimeStep(Scalar t)
  {
    initTimeStep_  = t;
    isInitialized_ = false;
  }
  virtual void setMaxTimeStep(Scalar t)
  {
    maxTimeStep_   = t;
    isInitialized_ = false;
  }
  virtual void setInitIndex(int i)
  {
    initIndex_     = i;
    isInitialized_ = false;
  }
  virtual void setFinalIndex(int i)
  {
    finalIndex_    = i;
    isInitialized_ = false;
  }
  virtual void setMaxAbsError(Scalar e)
  {
    maxAbsError_   = e;
    isInitialized_ = false;
  }
  virtual void setMaxRelError(Scalar e)
  {
    maxRelError_   = e;
    isInitialized_ = false;
  }
  virtual void setMaxFailures(int i)
  {
    maxFailures_   = i;
    isInitialized_ = false;
  }
  virtual void setMaxConsecFailures(int i)
  {
    maxConsecFailures_ = i;
    isInitialized_     = false;
  }
  virtual void setPrintDtChanges(bool b)
  {
    printDtChanges_ = b;
    isInitialized_  = false;
  }
  virtual void setNumTimeSteps(int numTimeSteps);
  virtual void setOutputExactly(bool b);
  virtual void setOutputIndices(std::vector<int> v);
  virtual void setOutputTimes(std::vector<Scalar> outputTimes);
  virtual void setOutputIndexInterval(int i);
  virtual void setOutputTimeInterval(Scalar t);
  virtual void setTimeEvents(
      Teuchos::RCP<TimeEventComposite<Scalar>> teb = Teuchos::null);
  virtual void setTimeStepControlStrategy(
      Teuchos::RCP<TimeStepControlStrategy<Scalar>> tscs = Teuchos::null);
  //@}

  virtual void printDtChanges(int istep, Scalar dt_old, Scalar dt_new,
                              std::string reason) const;

  virtual void initialize() const;
  virtual bool isInitialized() { return isInitialized_; }
  virtual void checkInitialized();

 protected:
  mutable bool isInitialized_;  ///< Bool if TimeStepControl is initialized.
  Scalar initTime_;             ///< Initial Time
  Scalar finalTime_;            ///< Final Time
  Scalar minTimeStep_;          ///< Minimum Time Step
  Scalar initTimeStep_;         ///< Initial Time Step
  Scalar maxTimeStep_;          ///< Maximum Time Step
  int initIndex_;               ///< Initial Time Index
  int finalIndex_;              ///< Final Time Index
  Scalar maxAbsError_;          ///< Maximum Absolute Error
  Scalar maxRelError_;          ///< Maximum Relative Error
  int maxFailures_;             ///< Maximum Number of Stepper Failures
  int maxConsecFailures_;       ///< Maximum Number of Consecutive Stepper Failures
  int numTimeSteps_;            ///< Number of time steps for Constant time step
  bool printDtChanges_;         ///< Print timestep size when it changes

  Teuchos::RCP<TimeEventComposite<Scalar>> timeEvent_;

  bool teAdjustedDt_;        ///< Flag indicating that dt was adjusted for time event.
  Scalar dtAfterTimeEvent_;  ///< dt to reinstate after TimeEvent step.

  Teuchos::RCP<TimeStepControlStrategy<Scalar>> stepControlStrategy_;
};

/// Nonmember constructor from ParameterList.
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar>> createTimeStepControl(
    Teuchos::RCP<Teuchos::ParameterList> const& pList,
    bool runInitialize = true);

/// Nonmember function to return ParameterList with default values.
template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> getTimeStepControlPL()
{
  auto tsc = rcp(new Tempus::TimeStepControl<Scalar>());
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      tsc->getValidParameters());
}

}  // namespace Tempus

#endif  // Tempus_TimeStepControl_decl_hpp
