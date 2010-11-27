
#ifndef PANZER_BCSTRATEGY_DIRICHLET_CONSTANT_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_CONSTANT_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

  template <typename EvalT>
    class BCStrategy_Dirichlet_Constant : public panzer::BCStrategy<EvalT> {

  public:    

    BCStrategy_Dirichlet_Constant(const panzer::BC& bc);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb) const;

  };

}

#include "user_app_BCStrategy_Dirichlet_ConstantT.hpp"

#endif
