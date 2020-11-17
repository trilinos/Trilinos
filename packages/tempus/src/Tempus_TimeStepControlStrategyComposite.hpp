// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategyComposite_hpp
#define Tempus_TimeStepControlStrategyComposite_hpp

#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

template<class Scalar> class TimeStepControl;

/** \brief TimeStepControlStrategyComposite loops over a vector of TimeStepControlStrategy.
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
 * <b>Note:<b> The ordering in the TimeStepControlStrategyComposite list is very important.
 * The final TimeStepControlStrategy from the composite could negate all previous step size updates.
 */
template<class Scalar>
class TimeStepControlStrategyComposite : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Constructor
  TimeStepControlStrategyComposite(){}

  /// Destructor
  virtual ~TimeStepControlStrategyComposite(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(const TimeStepControl<Scalar> tsc, Teuchos::RCP<SolutionHistory<Scalar> > sh ,
        Status & integratorStatus) override {
     for(auto& s : strategies_)
        s->getNextTimeStep(tsc, sh, integratorStatus);
  }

  /** \brief Append strategy to the composite list.*/
  void addStrategy(const Teuchos::RCP<TimeStepControlStrategy<Scalar> > &strategy){
     if (Teuchos::nonnull(strategy))
        strategies_.push_back(strategy);
  }

  /** \brief Clear the composite list.*/
  void clearObservers(){
     strategies_.clear();
  }

private:
  std::vector<Teuchos::RCP<TimeStepControlStrategy<Scalar > > > strategies_;

};
} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_hpp
