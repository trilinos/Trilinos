#ifndef USER_APP_MODEL_FACTORY_HPP
#define USER_APP_MODEL_FACTORY_HPP

#include "Panzer_ModelFactory.hpp"

namespace panzer {
  class InputEquationSet;
}

namespace user_app {

  template<typename EvalT>
  class MyModelFactory : public panzer::ModelFactory<EvalT> {

  public:

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
      buildModels(const panzer::InputEquationSet& set,
		  const std::vector<Teuchos::ParameterList>& models,
		  const Teuchos::ParameterList& default_params) const;

  };

}

#include "user_app_ModelFactoryT.hpp"

#endif
