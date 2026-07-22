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

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {

TEUCHOS_UNIT_TEST(tEpetraScatter, constructor)
{

   using Teuchos::RCP;
   using Teuchos::rcp;

   typedef panzer::Traits::Residual Residual;
   typedef panzer::Traits::Jacobian Jacobian;

   Teuchos::RCP<shards::CellTopology> topo
     = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

   // auxiliary information needed to construct basis object
   std::size_t numCells = 10;
   std::string basisType = "Q1";
   panzer::CellData cellData(numCells,topo);

   // build basis
   RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,1,cellData));

   std::string scatterName = "Residual_NS";

   // build DOF names
   RCP<std::vector<std::string> > evaluatedNames = rcp(new std::vector<std::string>); 
   evaluatedNames->push_back("Residual_ux"); // in practice these probably would not be scattered together!
   evaluatedNames->push_back("Residual_p");

   // build evaluated map
   RCP<std::map<std::string,std::string> > evaluatedMap = rcp(new std::map<std::string,std::string>); 
   evaluatedMap->insert(std::make_pair("Residual_ux","ux")); // in practice these probably would not be scattered together!
   evaluatedMap->insert(std::make_pair("Residual_p","p"));

   // build scatter parameter list
   Teuchos::ParameterList scatterParams;
   scatterParams.set<std::string>("Scatter Name",scatterName);
   scatterParams.set<RCP<std::vector<std::string> > >("Dependent Names",evaluatedNames);
   scatterParams.set<RCP<std::map<std::string,std::string> > >("Dependent Map",evaluatedMap);
   scatterParams.set<RCP<panzer::PureBasis> >("Basis",basis);
   scatterParams.set<int>("Side Subcell Dimension",1);
   scatterParams.set<int>("Local Side ID",2);
   scatterParams.set<bool>("Check Apply BC",false);

   // test residual scatter evaluator
   {
      panzer::ScatterDirichletResidual_Epetra<Residual,panzer::Traits,int,int> scatterResidual(Teuchos::null,Teuchos::null,scatterParams);

      const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual.evaluatedFields();
      TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager

      const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual.dependentFields();
      TEST_EQUALITY(fields.size(),2); // these store the residual values

      TEST_EQUALITY(fields[0]->name(),"Residual_ux");
      TEST_EQUALITY(fields[1]->name(),"Residual_p");

      TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   }

   // test jacobian scatter evaluator
   {
      panzer::ScatterDirichletResidual_Epetra<Jacobian,panzer::Traits,int,int> scatterJacobian(Teuchos::null,Teuchos::null,scatterParams);

      const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterJacobian.evaluatedFields();
      TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager

      const std::vector<RCP<PHX::FieldTag> > & fields = scatterJacobian.dependentFields();
      TEST_EQUALITY(fields.size(),2); // these store the residual values

      TEST_EQUALITY(fields[0]->name(),"Residual_ux");
      TEST_EQUALITY(fields[1]->name(),"Residual_p");

      TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   }
}

}
