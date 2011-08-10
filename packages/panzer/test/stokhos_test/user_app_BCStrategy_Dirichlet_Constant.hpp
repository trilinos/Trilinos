
#ifndef PANZER_BCSTRATEGY_DIRICHLET_CONSTANT_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_CONSTANT_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Basis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

  template <typename EvalT>
  class BCStrategy_Dirichlet_Constant : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
    
  public:    
    
    BCStrategy_Dirichlet_Constant(const panzer::BC& bc);
    
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& user_data);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				    const Teuchos::ParameterList& models,
				    const Teuchos::ParameterList& user_data) const;

    std::string residual_name;
    Teuchos::RCP<panzer::Basis> basis;

  };

}

#include "user_app_BCStrategy_Dirichlet_ConstantT.hpp"

#endif
