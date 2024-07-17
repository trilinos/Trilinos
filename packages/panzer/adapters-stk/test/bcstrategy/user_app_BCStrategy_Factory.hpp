// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_BCSTRATEGY_FACTORY_HPP
#define USER_APP_BCSTRATEGY_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_GlobalData.hpp"

// Add my bcstrategies here
#include "user_app_BCStrategy_Dirichlet_Constant.hpp"
#include "user_app_BCStrategy_Neumann_Constant.hpp"

namespace user_app {
  
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    user_app::BCStrategy_Dirichlet_Constant, BCStrategy_Dirichlet_Constant)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    user_app::BCStrategy_Neumann_Constant, BCStrategy_Neumann_Constant)

  struct BCFactory : public panzer::BCStrategyFactory {

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant",
        BCStrategy_Dirichlet_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Constant",
        BCStrategy_Neumann_Constant)

      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
			 "Error - the BC Strategy called \"" << bc.strategy() <<
			 "\" is not a valid identifier in the BCStrategyFactory.  Either add a valid implementation to your factory or fix your input file.  The relevant boundary condition is:\n\n" << bc << std::endl);
      
      return bcs_tm;

    }

  };
  
}

#endif
