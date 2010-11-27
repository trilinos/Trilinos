
#ifndef USER_APP_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define USER_APP_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "user_app_ModelFactory.hpp"

namespace user_app {

  class MyModelFactory_TemplateBuilder {

  public:
    
    template <typename EvalT>
    Teuchos::RCP<panzer::Base> build() const {
      return Teuchos::rcp( static_cast<panzer::Base*>
			   (new user_app::MyModelFactory<EvalT>) );
    }
    
  };
  
}

#endif 
