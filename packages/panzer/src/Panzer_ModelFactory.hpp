#ifndef PANZER_MODEL_FACTORY_HPP
#define PANZER_MODEL_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_InputEquationSet.hpp"

namespace panzer {

  template<typename EvalT>
  class ModelFactory : public panzer::Base {

  public:

    ModelFactory() {}
    
    virtual ~ModelFactory() {}
    
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
      virtual  buildModels(const panzer::InputEquationSet& set,
			   const std::vector<Teuchos::ParameterList>& models,
			   const Teuchos::ParameterList& default_params) const = 0;

  };
  
}

#endif
