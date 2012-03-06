#ifndef USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_2_TEMPLATE_BUILDER_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_2_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "user_app_ClosureModel_Factory_Physics2.hpp"

namespace user_app {

  class MyModelFactory_Physics2_TemplateBuilder {

    bool m_throw_if_model_not_found;

  public:

    MyModelFactory_Physics2_TemplateBuilder(bool throw_if_model_not_found) :
      m_throw_if_model_not_found(throw_if_model_not_found) {}

    template <typename EvalT>
    Teuchos::RCP<panzer::Base> build() const {
      return Teuchos::rcp( static_cast<panzer::Base*>
			   (new user_app::MyModelFactory_Physics2<EvalT>(m_throw_if_model_not_found)) );
    }
    
  };
  
}

#endif 
