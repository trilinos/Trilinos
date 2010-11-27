#ifndef PANZER_BCSTRATEGY_HPP
#define PANZER_BCSTRATEGY_HPP

#include "Panzer_BCStrategy_Base.hpp"
#include "Panzer_BC.hpp"

namespace panzer {

  template <typename EvalT>
    class BCStrategy : public panzer::BCStrategyBase {
    
  public:    
    
    BCStrategy(const panzer::BC& bc);
    
    virtual void 
      buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				 const panzer::PhysicsBlock& pb) const = 0;

  protected:

    const panzer::BC m_bc;

  };
  
}

#include "Panzer_BCStrategyT.hpp"

#endif
