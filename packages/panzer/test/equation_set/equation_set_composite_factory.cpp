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

#include "Panzer_CellData.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "user_app_EquationSetFactory_Physics1.hpp"
#include "user_app_EquationSetFactory_Physics2.hpp"
#include "Panzer_EquationSet_Factory_Composite.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(equation_set, composite_factory)
  {    
    const int num_cells = 20;
    const int cell_dim = 2;

    panzer::CellData cell_data(num_cells, cell_dim,
			       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();

    panzer::InputEquationSet ies1;
    ies1.name = "Energy 1";
    ies1.basis = "Q2";
    ies1.integration_order = 1;
    ies1.model_id = "solid";
    ies1.prefix = "ION_";

    panzer::InputEquationSet ies2;
    ies2.name = "Energy 2";
    ies2.basis = "Q2";
    ies2.integration_order = 1;
    ies2.model_id = "solid";
    ies2.prefix = "ION_";

    panzer::InputEquationSet ies_bad;
    ies_bad.name = "UNKNOWN";
    ies_bad.basis = "Q2";
    ies_bad.integration_order = 1;
    ies_bad.model_id = "solid";
    ies_bad.prefix = "ION_";    

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    // Check that eqf1 can build Energy 1 and fails to build Energy 2 
    {
      user_app::MyFactory1 my_factory1(true);
      TEST_NOTHROW(eq_set = my_factory1.buildEquationSet(ies1, cell_data, global_data, false));
      TEST_ASSERT(nonnull(eq_set));
      TEST_THROW(eq_set = my_factory1.buildEquationSet(ies2, cell_data, global_data, false),std::logic_error);
    }

    // Check that eqf2 can build Energy 2 and fails to build Energy 1 
    {
      eq_set = Teuchos::null;
      user_app::MyFactory2 my_factory2(true);
      TEST_NOTHROW(eq_set = my_factory2.buildEquationSet(ies2, cell_data, global_data, false));
      TEST_ASSERT(nonnull(eq_set));
      TEST_THROW(eq_set = my_factory2.buildEquationSet(ies1, cell_data, global_data, false),std::logic_error);
    }

    Teuchos::RCP<panzer::EquationSetFactory> eqsf1 = Teuchos::rcp(new user_app::MyFactory1(false));
    Teuchos::RCP<panzer::EquationSetFactory> eqsf2 = Teuchos::rcp(new user_app::MyFactory2(false));

    std::vector<Teuchos::RCP<panzer::EquationSetFactory> > factories;
    factories.push_back(eqsf1);
    factories.push_back(eqsf2);

    Teuchos::RCP<panzer::EquationSet_FactoryComposite> composite_factory = 
      Teuchos::rcp(new panzer::EquationSet_FactoryComposite(factories));

    // Build each physics
    eq_set = Teuchos::null;
    TEST_NOTHROW(eq_set = composite_factory->buildEquationSet(ies1, cell_data, global_data, false));
    TEST_ASSERT(nonnull(eq_set));

    eq_set = Teuchos::null;
    TEST_NOTHROW(eq_set = composite_factory->buildEquationSet(ies2, cell_data, global_data, false));
    TEST_ASSERT(nonnull(eq_set));
    
    // try bad eq set
    eq_set = Teuchos::null;
    TEST_THROW(eq_set = composite_factory->buildEquationSet(ies_bad, cell_data, global_data, false),std::logic_error);
    TEST_ASSERT(is_null(eq_set));
    
  }

}
