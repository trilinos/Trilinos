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

#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {

TEUCHOS_UNIT_TEST(tEpetraGather, constructor)
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

   // build DOF names
   RCP<std::vector<std::string> > dofNames = rcp(new std::vector<std::string>); 
   dofNames->push_back("ux"); // in practice these probably would not be gathered together!
   dofNames->push_back("p");

   // build basis
   RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,1,cellData));

   // build gather parameter list
   Teuchos::ParameterList gatherParams;
   gatherParams.set<RCP<std::vector<std::string> > >("DOF Names",dofNames);
   gatherParams.set<RCP<std::vector<std::string> > >("Indexer Names",dofNames);
   gatherParams.set<RCP<panzer::PureBasis> >("Basis",basis);

   // test residual gather evaluator
   {
      panzer::GatherSolution_Epetra<Residual,panzer::Traits,int,int> gatherResidual(Teuchos::null,gatherParams);

      const std::vector<RCP<PHX::FieldTag> > & fields = gatherResidual.evaluatedFields();
      TEST_EQUALITY(fields.size(),2);
 
      TEST_EQUALITY(fields[0]->name(),"ux");
      TEST_EQUALITY(fields[1]->name(),"p");

      TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   }

   // test jacobian gather evaluator
   {
      panzer::GatherSolution_Epetra<Jacobian,panzer::Traits,int,int> gatherJacobian(Teuchos::null,gatherParams);

      const std::vector<RCP<PHX::FieldTag> > & fields = gatherJacobian.evaluatedFields();
      TEST_EQUALITY(fields.size(),2);
 
      TEST_EQUALITY(fields[0]->name(),"ux");
      TEST_EQUALITY(fields[1]->name(),"p");

      TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   }
}

}
