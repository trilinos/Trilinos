// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_BCStrategy_Factory_Composite.hpp"
#include "Panzer_BC.hpp"

namespace panzer {
  
  BCFactoryComposite::BCFactoryComposite(const std::vector<Teuchos::RCP<panzer::BCStrategyFactory> >& factories) :
    m_bc_strategy_factories(factories)
  {
    
  }
  
  Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
  BCFactoryComposite::
  buildBCStrategy(const panzer::BC& bc, 
		  const Teuchos::RCP<panzer::GlobalData>& global_data) const
  {    
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm;
    
    bool found = false;

    for (std::vector<Teuchos::RCP<panzer::BCStrategyFactory> >::const_iterator factory = m_bc_strategy_factories.begin(); 
         factory != m_bc_strategy_factories.end(); ++factory) {

      bcs_tm = (*factory)->buildBCStrategy(bc,global_data);
      
      if (nonnull(bcs_tm)) {
        found = true;
        break;
      }

    }
        
    TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
			       "Error - the BC Strategy called \"" << bc.strategy() <<
			       "\" is not a valid identifier in the BCStrategyFactory.  Either add " <<
                               "a valid implementation to the factory or fix the input file.  The " <<
                               "relevant boundary condition is:\n\n" << bc << std::endl);
    
    return bcs_tm;
    
  }
  
}
