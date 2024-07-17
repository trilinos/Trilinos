// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_FACTORY_COMPOSITE_DECL_HPP
#define PANZER_BCSTRATEGY_FACTORY_COMPOSITE_DECL_HPP

#include "Panzer_BCStrategy_Factory.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>

namespace panzer {
  
  struct BCFactoryComposite : public panzer::BCStrategyFactory {

  public:

    BCFactoryComposite(const std::vector<Teuchos::RCP<panzer::BCStrategyFactory> >& factories);

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const;

  private:

    std::vector<Teuchos::RCP<panzer::BCStrategyFactory> > m_bc_strategy_factories;
    
  };
  
}

#endif
