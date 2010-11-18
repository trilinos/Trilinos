#ifndef PANZER_EQUATION_SET_FACTORY_H
#define PANZER_EQUATION_SET_FACTORY_H

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {
  
  template<typename T> class EquationSet_TemplateManager;

  struct EquationSetFactory {

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const panzer::InputEquationSet& ies,
		     const panzer::CellData& cell_data) const;

  };
  
}

#endif
