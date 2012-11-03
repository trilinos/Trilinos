// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
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
