
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "user_app_EquationSet_Energy.hpp"

namespace user_app {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("Energy", user_app::EquationSet_Energy,
					EquationSet_Energy)

  class MyFactory : public panzer::EquationSetFactory {

  public:

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const panzer::InputEquationSet& ies,
		     const panzer::CellData& cell_data) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
	Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
      
      bool found = false;
      
      PANZER_BUILD_EQSET_OBJECTS("Energy", my_app::EquationSet_Energy,
				 EquationSet_Energy)
      
      if (!found) {
	std::string msg = "Error - the \"Equation Set\" called \"" + ies.name +
	  "\" is not a valid equation set identifier. Please supply the correct factory.\n";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      
      return eq_set;
    }
    
  };

}

