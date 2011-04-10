#ifndef USER_APP_STK_CLOSURE_MODEL_FACTORY_HPP
#define USER_APP_STK_CLOSURE_MODEL_FACTORY_HPP

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer {
  class InputEquationSet;
}

namespace user_app {

  template<typename EvalT>
  class STKModelFactory : public panzer::ClosureModelFactory<EvalT> {

  public:

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
			 const panzer::InputEquationSet& set,
			 const Teuchos::ParameterList& models,
			 const Teuchos::ParameterList& default_params) const;

  };

}

#include "user_app_STKClosureModel_FactoryT.hpp"

#endif
