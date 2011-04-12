#ifndef PANZER_CLOSURE_MODEL_FACTORY_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_InputEquationSet.hpp"
#include <string>
#include <vector>

namespace panzer {

  template<typename EvalT>
  class ClosureModelFactory : public panzer::Base {

  public:

    ClosureModelFactory() {}
    
    virtual ~ClosureModelFactory() {}
    
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    virtual  buildClosureModels(const std::string& model_id,
				const panzer::InputEquationSet& set,
				const Teuchos::ParameterList& models,
				const Teuchos::ParameterList& default_params,
				const Teuchos::ParameterList& user_data) const = 0;

  };
  
}

#endif
