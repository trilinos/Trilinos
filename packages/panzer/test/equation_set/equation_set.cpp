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
#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(equation_set, steady_state)
  {
    panzer::InputEquationSet ies;
    ies.name = "Energy";
    ies.basis = "Q2";
    ies.integration_order = 1;
    ies.model_id = "solid";
    ies.prefix = "ION_";
    
    int num_cells = 20;
    int cell_dim = 2;
    panzer::CellData cell_data(num_cells, cell_dim,
                        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    user_app::MyFactory my_factory;
    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();
    eq_set = my_factory.buildEquationSet(ies, cell_data, global_data, false);
  }

  TEUCHOS_UNIT_TEST(equation_set, transient)
  {
    panzer::InputEquationSet ies;
    ies.name = "Energy";
    ies.basis = "Q2";
    ies.integration_order = 1;
    ies.model_id = "solid";
    ies.prefix = "ION_";
    
    int num_cells = 20;
    int cell_dim = 2;
    panzer::CellData cell_data(num_cells, cell_dim,
                        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    user_app::MyFactory my_factory;
    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();
    eq_set = my_factory.buildEquationSet(ies, cell_data, global_data, true);
  }

}
