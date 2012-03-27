#include "Panzer_EquationSet_Factory_Composite.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"

namespace panzer {
  
  EquationSet_FactoryComposite::EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories) :
    m_factories(factories)
  { }
  
  Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
  EquationSet_FactoryComposite::buildEquationSet(const panzer::InputEquationSet& ies,
						 const panzer::CellData& cell_data,
						 const Teuchos::RCP<panzer::GlobalData>& global_data,
						 const bool build_transient_support) const
  {
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
    
    
    for (std::vector<Teuchos::RCP<panzer::EquationSetFactory> >::const_iterator factory = m_factories.begin();
	 factory != m_factories.end(); ++factory) {
      eq_set = (*factory)->buildEquationSet(ies,cell_data,global_data,build_transient_support);

      if (nonnull(eq_set))
	break;
    }
    
    if (is_null(eq_set)) {
      std::string msg = "Error - the \"Equation Set\" called \"" + ies.name +
	"\" is not a valid equation set identifier. Please supply the correct factory.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(eq_set), std::logic_error, msg);
    }
    
    return eq_set;
  }
  
}

