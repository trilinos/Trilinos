#ifndef PHX_FIELD_EVALUATOR_FACTORY_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_HPP

#include <map>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_FieldEvaluator_TemplateManager.hpp"

namespace PHX {

  template<typename Traits, typename FactoryTraits>
  class FieldEvaluatorFactory {
    
  public:
    typename Teuchos::RCP< std::vector< Teuchos::RCP<PHX::FieldEvaluator_TemplateManager<Traits> > > > 
    buildFieldEvaluators(const std::map<std::string, 
			 Teuchos::RCP<Teuchos::ParameterList> >& data);
    
  };
} 

#include "Phalanx_FieldEvaluator_Factory_Def.hpp"

#endif 
