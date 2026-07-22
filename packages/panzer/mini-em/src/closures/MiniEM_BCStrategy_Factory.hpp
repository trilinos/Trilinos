// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_BCStrategy_Factory_hpp_
#define _MiniEM_BCStrategy_Factory_hpp_

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"

#include "MiniEM_BCStrategy_Dirichlet_Constant.hpp"
#include "MiniEM_BCStrategy_Dirichlet_AuxConstant.hpp"

namespace mini_em {
  
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(mini_em::BCStrategy_Dirichlet_Constant,
					     BCStrategy_Dirichlet_Constant)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(mini_em::BCStrategy_Dirichlet_AuxConstant,
					     BCStrategy_Dirichlet_AuxConstant)

  struct BCFactory : public panzer::BCStrategyFactory {

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
	Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant", BCStrategy_Dirichlet_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("AuxConstant", BCStrategy_Dirichlet_AuxConstant)

      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
			 "Error - the BC Strategy called \"" << bc.strategy() <<
			 "\" is not a valid identifier in the BCStrategyFactory.  Either add a valid implementation to your factory or fix your input file.  The relevant boundary condition is:\n\n" << bc << std::endl);
      
      return bcs_tm;

    }

  };
  
}

#endif
