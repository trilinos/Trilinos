//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeStepControlStrategyComposite_hpp
#define Tempus_TimeStepControlStrategyComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"
#include "Tempus_TimeStepControlStrategyIntegralController.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

/** \brief TimeStepControlStrategyComposite loops over a vector of
 * TimeStepControlStrategies.
 *
 *
 * Essentially, this is an <b>and</b> case if each strategies do a `min`
 * \f$ \Delta t = \min_{i \leq N} \{ \Delta t_i \}\f$
 *
 * The assumption is that each strategy will simply
 * update (or override) the step size `dt` with `metadata->setDt(dt)`
 * sequentially.
 *
 *  Examples of TimeStepControlStrategy:
 *   - TimeStepControlStrategyConstant
 *   - TimeStepControlStrategyBasicVS
 *   - TimeStepControlStrategyIntegralController
 *
 * <b>Note:</b> The ordering in the TimeStepControlStrategyComposite
 * list is very important.  The final TimeStepControlStrategy from
 * the composite could negate all previous step size updates.
 */
template <class Scalar>
class TimeStepControlStrategyComposite
  : virtual public TimeStepControlStrategy<Scalar> {
 public:
  /// Constructor
  TimeStepControlStrategyComposite()
  {
    this->setStrategyType("Composite");
    this->setStepType("Variable");
    this->isInitialized_ = false;
  }

  /// Destructor
  virtual ~TimeStepControlStrategyComposite() {}

  /** \brief Determine the time step size.*/
  virtual void setNextTimeStep(const TimeStepControl<Scalar>& tsc,
                               Teuchos::RCP<SolutionHistory<Scalar>> sh,
                               Status& integratorStatus) override
  {
    for (auto& s : strategies_) s->setNextTimeStep(tsc, sh, integratorStatus);
  }

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override
  {
    return "Tempus::TimeStepControlComposite";
  }

  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel) const override
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, this->description());
    l_out->setOutputToRootOnly(0);

    *l_out << "\n--- " << this->description() << " ---" << std::endl;

    if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
      *l_out << "  Strategy Type = " << this->getStrategyType() << std::endl
             << "  Step Type     = " << this->getStepType() << std::endl;

      std::stringstream sList;
      for (std::size_t i = 0; i < strategies_.size(); ++i) {
        sList << strategies_[i]->getStrategyType();
        if (i < strategies_.size() - 1) sList << ", ";
      }
      *l_out << "  Strategy List = " << sList.str() << std::endl;

      for (auto& s : strategies_) s->describe(*l_out, verbLevel);

      *l_out << std::string(this->description().length() + 8, '-') << std::endl;
    }
  }
  //@}

  /** \brief Append strategy to the composite list.*/
  void addStrategy(
      const Teuchos::RCP<TimeStepControlStrategy<Scalar>>& strategy)
  {
    if (Teuchos::nonnull(strategy)) {
      if (this->size() == 0) this->setStepType(strategy->getStepType());

      TEUCHOS_TEST_FOR_EXCEPTION(
          this->getStepType() != strategy->getStepType(), std::logic_error,
          "Error - Cannot mix 'Constant' and 'Variable' step types.\n"
          "strategies in composite!  Need at least one.\n");

      strategies_.push_back(strategy);
    }
  }

  /** \brief Size of composite list.*/
  virtual int size() const { return strategies_.size(); }

  /** \brief Return composite list.*/
  virtual std::vector<Teuchos::RCP<TimeStepControlStrategy<Scalar>>>
  getStrategies() const
  {
    return strategies_;
  }

  /** \brief Clear the composite list.*/
  void clearStrategies() { strategies_.clear(); }

  virtual void initialize() const override
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        strategies_.size() == 0, std::logic_error,
        "Error - No strategies in composite!  Need at least one.\n");

    for (auto& s : strategies_) s->initialize();

    auto strategy0 = strategies_[0];
    for (auto& s : strategies_) {
      TEUCHOS_TEST_FOR_EXCEPTION(s->isInitialized() == false, std::logic_error,
                                 "Error - Composite strategy, "
                                     << s->getName()
                                     << " is not initialized!\n");

      if (strategy0->getStepType() != s->getStepType()) {
        std::ostringstream msg;
        msg << "Error - All the Strategy Step Types must match.\n";
        for (std::size_t i = 0; i < strategies_.size(); ++i) {
          msg << "  Strategy[" << i << "] = " << strategies_[i]->getStepType()
              << "  (" << strategies_[i]->getName() << ")\n";
        }
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
      }
    }

    this->isInitialized_ = true;  // Only place where this is set to true!
  }

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::parameterList("Time Step Control Strategy");

    pl->set<std::string>("Strategy Type", this->getStrategyType(), "Composite");

    std::stringstream sList;
    for (std::size_t i = 0; i < strategies_.size(); ++i) {
      sList << strategies_[i]->getStrategyType();
      if (i < strategies_.size() - 1) sList << ", ";
    }
    pl->set<std::string>("Strategy List", sList.str());

    for (auto& s : strategies_) pl->set(s->getName(), *s->getValidParameters());

    return pl;
  }

 private:
  std::vector<Teuchos::RCP<TimeStepControlStrategy<Scalar>>> strategies_;
};

// Nonmember constructor - ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<TimeStepControlStrategyComposite<Scalar>>
createTimeStepControlStrategyComposite(
    Teuchos::RCP<Teuchos::ParameterList> const& pList,
    std::string name = "Composite")
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  std::vector<std::string> tscsList;

  TEUCHOS_TEST_FOR_EXCEPTION(
      pList->get<std::string>("Strategy Type") != "Composite", std::logic_error,
      "Error - Strategy Type != 'Composite'.  (='" +
          pList->get<std::string>("Strategy Type") + "')\n");

  // string tokenizer
  tscsList.clear();
  std::string str = pList->get<std::string>("Strategy List");
  std::string delimiters(",");
  const char* WhiteSpace = " \t\v\r\n";
  // Skip delimiters at the beginning
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find the first delimiter
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    // Found a token, add it to the vector
    std::string token = str.substr(lastPos, pos - lastPos);

    std::size_t start = token.find_first_not_of(WhiteSpace);
    std::size_t end   = token.find_last_not_of(WhiteSpace);
    token =
        (start == end ? std::string() : token.substr(start, end - start + 1));

    tscsList.push_back(token);
    if (pos == std::string::npos) break;

    lastPos = str.find_first_not_of(delimiters, pos);  // Skip delimiters
    pos     = str.find_first_of(delimiters, lastPos);  // Find next delimiter
  }

  auto tscsc = Teuchos::rcp(new TimeStepControlStrategyComposite<Scalar>());

  // For each sublist name tokenized, add the TSCS
  for (auto tscsName : tscsList) {
    RCP<ParameterList> pl =
        Teuchos::rcp(new ParameterList(pList->sublist(tscsName, true)));

    auto strategyType = pl->get<std::string>("Strategy Type", "Unknown");
    if (strategyType == "Constant") {
      tscsc->addStrategy(
          createTimeStepControlStrategyConstant<Scalar>(pl, tscsName));
    }
    else if (strategyType == "Basic VS") {
      tscsc->addStrategy(
          createTimeStepControlStrategyBasicVS<Scalar>(pl, tscsName));
    }
    else if (strategyType == "Integral Controller") {
      tscsc->addStrategy(
          createTimeStepControlStrategyIntegralController<Scalar>(pl,
                                                                  tscsName));
    }
    else if (strategyType == "Composite") {
      tscsc->addStrategy(
          createTimeStepControlStrategyComposite<Scalar>(pl, tscsName));
    }
    else {
      RCP<Teuchos::FancyOStream> out =
          Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(out, 1, "createTimeStepControlStrategyComposite()");
      *out << "Warning -- Unknown strategy type!\n"
           << "'Strategy Type' = '" << strategyType << "'\n"
           << "Should call addStrategy() with this\n"
           << "(app-specific?) strategy, and initialize().\n"
           << std::endl;
    }
  }

  tscsc->setName(name);

  if (tscsc->size() == 0) {
    RCP<Teuchos::FancyOStream> out =
        Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "createTimeStepControlStrategyComposite()");
    *out << "Warning -- Did not find a Tempus strategy to create!\n"
         << "Should call addStrategy() with (app-specific?) strategy(ies),\n"
         << "and initialize().\n"
         << std::endl;
  }
  else {
    tscsc->initialize();
  }

  return tscsc;
}

/// Nonmember function to return ParameterList with default values.
template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> getTimeStepControlStrategyCompositePL()
{
  auto t    = rcp(new Tempus::TimeStepControlStrategyComposite<Scalar>());
  auto tscs = rcp(new Tempus::TimeStepControlStrategyConstant<Scalar>());
  t->addStrategy(tscs);
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      t->getValidParameters());
}

}  // namespace Tempus
#endif  // Tempus_TimeStepControlStrategy_hpp
