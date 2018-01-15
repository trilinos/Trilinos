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
class StepControlStrategyComposite : virtual public StepControlStrategy<Scalar>
{
public:

  /// Constructor
  StepControlStrategyComposite(){}

  /// Destructor
  virtual ~StepControlStrategyComposite(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(TimeStepControl<Scalar> tsc, Teuchos::RCP<SolutionHistory<Scalar> > sh ,
        Status & integratorStatus) override {
     for(auto& s : strategies_)  
        s->getNextTimeStep(tsc, sh, integratorStatus);
  }

  // add strategy to the composite strategy list
  void addStrategy(const Teuchos::RCP<StepControlStrategy<Scalar> > &strategy){
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
