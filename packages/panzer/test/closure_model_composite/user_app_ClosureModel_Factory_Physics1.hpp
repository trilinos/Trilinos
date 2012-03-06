#ifndef USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_1_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_1_HPP

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer {
  class InputEquationSet;
}

namespace user_app {

  template<typename EvalT>
  class MyModelFactory_Physics1 : public panzer::ClosureModelFactory<EvalT> {

  public:

    MyModelFactory_Physics1(bool throw_if_model_not_found);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const panzer::InputEquationSet& set,
		       const Teuchos::ParameterList& models,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:

    
    bool m_throw_if_model_not_found;

  };

}

#include "user_app_ClosureModel_Factory_Physics1_impl.hpp"

#endif
