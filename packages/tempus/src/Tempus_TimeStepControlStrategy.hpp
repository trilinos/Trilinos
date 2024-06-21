//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeStepControlStrategy_hpp
#define Tempus_TimeStepControlStrategy_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

template <class Scalar>
class TimeStepControl;

/** \brief TimeStepControlStrategy class for TimeStepControl
 *
 *  This is the base class for TimeStepControlStrategies.
 *  The primary function required from derived classes is setNextTimeStep(),
 *  which will
 *   - determine the next step from information in the TimeStepControl
 *     and SolutionHistory (i.e., SolutionStates)
 *   - set the next time step on the workingState in the SolutionHistory
 *  If a valid timestep can not be determined the Status is set to FAILED.
 */
template <class Scalar>
class TimeStepControlStrategy : virtual public Teuchos::Describable,
                                virtual public Teuchos::VerboseObject<
                                    Tempus::TimeStepControlStrategy<Scalar> > {
 public:
  /// Constructor
  TimeStepControlStrategy()
    : strategyType_("Base Strategy"),
      stepType_("Constant"),
      name_("Base Strategy"),
      isInitialized_(false)
  {
  }

  /// Destructor
  virtual ~TimeStepControlStrategy() {}

  /// Set the time step size.
  virtual void setNextTimeStep(const TimeStepControl<Scalar>& /* tsc */,
                               Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
                               Status& /* integratorStatus */)
  {
  }

  virtual void initialize() const { isInitialized_ = true; }
  virtual bool isInitialized() { return isInitialized_; }
  virtual void checkInitialized()
  {
    if (!isInitialized_) {
      this->describe(*(this->getOStream()), Teuchos::VERB_MEDIUM);
      TEUCHOS_TEST_FOR_EXCEPTION(
          !isInitialized_, std::logic_error,
          "Error - " << this->description() << " is not initialized!");
    }
  }

  virtual void setName(std::string s) { name_ = s; }

  virtual std::string getStrategyType() const { return strategyType_; }
  virtual std::string getStepType() const { return stepType_; }
  virtual std::string getName() const { return name_; }

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }

 protected:
  virtual void setStrategyType(std::string s) { strategyType_ = s; }
  virtual void setStepType(std::string s) { stepType_ = s; }

  std::string strategyType_;    ///< Strategy type
  std::string stepType_;        ///< Step Type - "Constant" or "Variable"
  std::string name_;            ///< Name of strategy.
  mutable bool isInitialized_;  ///< Bool if strategy is initialized.
};

}  // namespace Tempus
#endif  // Tempus_TimeStepControlStrategy_hpp
