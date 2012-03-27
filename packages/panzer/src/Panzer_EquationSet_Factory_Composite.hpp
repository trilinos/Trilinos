#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"

namespace panzer {

  class EquationSet_FactoryComposite : public panzer::EquationSetFactory {

    std::vector<Teuchos::RCP<panzer::EquationSetFactory> > m_factories;

  public:
    
    EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories);

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const panzer::InputEquationSet& ies,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const;
    
  };

}

