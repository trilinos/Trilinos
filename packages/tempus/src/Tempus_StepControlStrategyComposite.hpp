// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepControlStrategyComposite_hpp
#define Tempus_StepControlStrategyComposite_hpp

#include "Tempus_StepControlStrategy.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

template<class Scalar> class TimeStepControl;

/** \brief StepControlStrategy class for TimeStepControl
 *
 */
template<class Scalar>
class StepControlStrategyComposite 
{
public:

  /// Constructor
  StepControlStrategyComposite(){}

  /// Destructor
  virtual ~StepControlStrategyComposite(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void getNextTimeStep(TimeStepControl<Scalar> tsc, Teuchos::RCP<SolutionHistory<Scalar> > sh ){
     for(auto& s : strategies_)  
        o->getNextTimeStep(tsc, sh);
  }

  // add observer to the composite observer list
  void addObserver(const Teuchos::RCP<StepControlStrategy<Scalar> > &strategy){
  strategies_.push_back(strategy);
  }

  void clearObservers(){
     strategies_.clear();
  }

private:
  std::vector<Teuchos::RCP<StepControlStrategy<Scalar > > > strategies_;

};
} // namespace Tempus
#endif // Tempus_StepControlStrategy_hpp
