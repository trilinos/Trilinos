#ifndef PANZER_BCSTRATEGY_FACTORY_HPP
#define PANZER_BCSTRATEGY_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {
  
  class BC;
  
  template<typename T> class BCStrategy_TemplateManager;

  struct BCStrategyFactory {

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc) const;

  };
  
}

#endif
