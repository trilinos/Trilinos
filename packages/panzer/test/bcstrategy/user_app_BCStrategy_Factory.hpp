#ifndef USER_APP_BCSTRATEGY_FACTORY_HPP
#define USER_APP_BCSTRATEGY_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"

// Add my bcstrategies here
#include "user_app_BCStrategy_Dirichlet_Constant.hpp"

namespace user_app {
  
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER("Constant", 
					     user_app::BCStrategy_Dirichlet_Constant,
					     BCStrategy_Dirichlet_Constant)

  struct MyFactory : public panzer::BCStrategyFactory {

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant", 
				      user_app::BCStrategy_Dirichlet_Constant,
				      BCStrategy_Dirichlet_Constant)

      TEST_FOR_EXCEPTION(!found, std::logic_error, 
			 "Error - the BC Strategy called \"" << bc.strategy() <<
			 "\" is not a valid identifier in the BCStrategyFactory.  Either add a valid implementation to your factory or fix your input file.  The relevant boundary condition is:\n\n" << bc << std::endl);
      
      return bcs_tm;

    }

  };
  
}

#endif
