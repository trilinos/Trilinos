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
#include <Teuchos_DefaultComm.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"
#include "Phalanx_FieldManager.hpp"

// for composite closure model test
#include "Panzer_ClosureModel_Factory_Composite.hpp"
#include "Panzer_ClosureModel_Factory_Composite_TemplateBuilder.hpp"
#include "user_app_ClosureModel_Factory_Physics1.hpp"
#include "user_app_ClosureModel_Factory_Physics1_TemplateBuilder.hpp"
#include "user_app_ClosureModel_Factory_Physics2.hpp"
#include "user_app_ClosureModel_Factory_Physics2_TemplateBuilder.hpp"


#include <iostream>
#include <vector>

namespace panzer {

  TEUCHOS_UNIT_TEST(closure_model_factory, composite_factory)
  {
    panzer::InputEquationSet ies;
    {
      ies.name = "Momentum";
      ies.basis = "Q2";
      ies.integration_order = 1;
      ies.model_id = "fluid model";
      ies.prefix = "";
      ies.params.set<int>("junk", 1);
    }

    Teuchos::ParameterList default_params; 
    {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
    
      const int num_cells = 20;
      const int base_cell_dimension = 3;
      const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
      const int cubature_degree = 2;      
      Teuchos::RCP<panzer::IntegrationRule> ir = 
	Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
      default_params.set("IR",ir);
      Teuchos::RCP<panzer::BasisIRLayout> basis = 
	Teuchos::rcp(new panzer::BasisIRLayout("Q1", *ir));
      default_params.set("Basis",basis);
      
    }

    Teuchos::ParameterList p("Closure Models");
    {
      p.sublist("fluid model").sublist("Density").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Viscosity").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Heat Capacity").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Thermal Conductivity").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Global Statistics").set("Value","UX");
    }

    // default data
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators;
    Teuchos::ParameterList user_data("User Data");
    user_data.set("Comm",Teuchos::DefaultComm<int>::getComm());
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    
    // Prove Physics 1 builds all constant evaluators
    {
      user_app::MyModelFactory_Physics1<panzer::Traits::Residual> p1(false);
      PHX::FieldManager<panzer::Traits> fm;
      evaluators = p1.buildClosureModels(ies.model_id, ies, p, default_params, user_data, gd, fm);
      TEST_EQUALITY(evaluators->size(), 8);
    }

    // Prove Physics 2 builds the global statistics evaluator
    {
      user_app::MyModelFactory_Physics2<panzer::Traits::Residual> p2(false);
      PHX::FieldManager<panzer::Traits> fm;
      evaluators = p2.buildClosureModels(ies.model_id, ies, p, default_params, user_data, gd, fm);
      TEST_EQUALITY(evaluators->size(), 1);
    }
    
    PHX::FieldManager<panzer::Traits> fm;


    user_app::MyModelFactory_Physics1_TemplateBuilder builder1(false);
    Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > model_factory1 = 
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
    model_factory1->buildObjects(builder1);

    user_app::MyModelFactory_Physics2_TemplateBuilder builder2(false);
    Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > model_factory2 = 
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
    model_factory2->buildObjects(builder2);

    panzer::ClosureModelFactoryComposite_TemplateBuilder builder_composite;
    builder_composite.addFactory(model_factory1);
    builder_composite.addFactory(model_factory2);
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> model_factory_composite;
    model_factory_composite.buildObjects(builder_composite);

    evaluators = model_factory_composite.getAsObject<panzer::Traits::Residual>()->buildClosureModels(ies.model_id, ies, p, default_params, user_data, gd, fm);
    
    TEST_EQUALITY(evaluators->size(), 9);
    
    // Add an unsupported type 

    // RPP: :disabling for now due to issues with global statistics.
    // The jacobian was returning a null pointer that translated intot
    // not found and throws an error.  Neeed to better design closure
    // model factory interface to return whether model was found.
    // This can be caught in field manager.

    //p.sublist("fluid model").sublist("garbage").set<std::string>("Value","help!");
    
    //TEST_THROW(model_factory_composite.getAsObject<panzer::Traits::Residual>()->buildClosureModels(ies.model_id, ies, p, default_params, user_data, gd, fm), std::logic_error);

  }


}
