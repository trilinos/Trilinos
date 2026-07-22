// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Assert.hpp>

#include "Panzer_Traits.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy.hpp"
#include <iostream>

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_GlobalData.hpp"

#include "user_app_BCStrategy_Factory_Physics1.hpp"
#include "user_app_BCStrategy_Factory_Physics2.hpp"

#include "Panzer_BCStrategy_Factory_Composite.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"

#include "Epetra_MpiComm.h"

namespace panzer {

  TEUCHOS_UNIT_TEST(bcstrategy, basic_construction)
  {

    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "4";
    std::string element_block_id = "fluid";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		  strategy, p);

    panzer::BC bc1(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		   "Constant 1", p);

    panzer::BC bc2(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		   "Constant 2", p);
    

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    // Build factory 1 and test that it builds bc 1 and throws on
    // failure to build bc 2
    {
      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs;
      user_app::BCFactoryPhysics1 my_factory1(true);
      TEST_NOTHROW(bcs = my_factory1.buildBCStrategy(bc1,gd));
      TEST_ASSERT(nonnull(bcs));
      TEST_THROW(bcs = my_factory1.buildBCStrategy(bc2,gd),std::logic_error);
    }

    // Build factory 2 and test that it builds bc 2 and throws on
    // failure to build bc 1
    {
      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs;
      user_app::BCFactoryPhysics2 my_factory2(true);
      TEST_NOTHROW(bcs = my_factory2.buildBCStrategy(bc2,gd));
      TEST_ASSERT(nonnull(bcs));
      TEST_THROW(bcs = my_factory2.buildBCStrategy(bc1,gd),std::logic_error);
    }
    
    // Build composite and test that it can build both bc 1 and bc 2
    {
      std::vector<Teuchos::RCP<panzer::BCStrategyFactory> > factories;
      Teuchos::RCP<user_app::BCFactoryPhysics1> factory1 = Teuchos::rcp(new user_app::BCFactoryPhysics1(false));
      Teuchos::RCP<user_app::BCFactoryPhysics2> factory2 = Teuchos::rcp(new user_app::BCFactoryPhysics2(false));
      factories.push_back(factory1);
      factories.push_back(factory2);
      panzer::BCFactoryComposite composite(factories);

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs;
      TEST_NOTHROW(bcs = composite.buildBCStrategy(bc1,gd));
      TEST_ASSERT(nonnull(bcs));
      bcs = Teuchos::null;
      TEST_NOTHROW(bcs = composite.buildBCStrategy(bc2,gd));
      TEST_ASSERT(nonnull(bcs));
    }

  }

}
