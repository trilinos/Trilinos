
#ifndef PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_BCStrategy.hpp"
#include "Panzer_Traits.hpp"

#include "Phalanx_FieldManager.hpp"

namespace panzer {
  
  template <typename EvalT>
    class BCStrategy_Dirichlet_DefaultImpl : public panzer::BCStrategy<EvalT> {

  public:    

    BCStrategy_Dirichlet_DefaultImpl(const panzer::BC& bc);
    
    virtual ~BCStrategy_Dirichlet_DefaultImpl();
    
    virtual void setup(const panzer::PhysicsBlock& side_pb, const Teuchos::ParameterList& user_data) = 0;
      
    virtual void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					    const panzer::PhysicsBlock& pb,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					    const Teuchos::ParameterList& models,
					    const Teuchos::ParameterList& user_data) const = 0;

    void 
    virtual buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
						    const panzer::PhysicsBlock& pb,
						    const panzer::LinearObjFactory<panzer::Traits> & lof,
						    const Teuchos::ParameterList& user_data) const;

  protected:

      std::vector<std::string> required_dof_names;

      std::map<std::string,std::string> residual_to_dof_names_map;
  };

}

#include "Panzer_BCStrategy_Dirichlet_DefaultImplT.hpp"

#endif
