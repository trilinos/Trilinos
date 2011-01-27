#ifndef PANZER_BCSTRATEGY_BASE_HPP
#define PANZER_BCSTRATEGY_BASE_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"

namespace PHX {
  template<typename T> class FieldManager;
}

namespace panzer {

  class PhysicsBlock;

  //! Non-templated empty base class for EquationSet objects
  class BCStrategyBase {
    
  public:
    
    BCStrategyBase() {}
    
    virtual ~BCStrategyBase() {}
    
    virtual void 
      buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				 const panzer::PhysicsBlock& pb) const = 0;

    virtual void 
      buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
  				              const panzer::PhysicsBlock& pb,
                                              const LinearObjFactory<panzer::Traits> & lof) const = 0;

  };
  
}

#endif
