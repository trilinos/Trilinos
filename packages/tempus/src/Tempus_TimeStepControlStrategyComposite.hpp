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

/** \brief StepControlStrategy class for TimeStepControl
 *
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

  // add strategy to the composite strategy list
  void addStrategy(const Teuchos::RCP<TimeStepControlStrategy<Scalar> > &strategy){
     if (Teuchos::nonnull(strategy))
        strategies_.push_back(strategy);
  }

  void clearObservers(){
     strategies_.clear();
  }

private:
  std::vector<Teuchos::RCP<TimeStepControlStrategy<Scalar > > > strategies_;

};
} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_hpp
