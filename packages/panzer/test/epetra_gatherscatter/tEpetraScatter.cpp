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

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "Panzer_ScatterResidual_Epetra.hpp"

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
   int baseCellDim = 2;
   int cubatureDegree = 2;
   std::string basisType = "Q1";
   panzer::CellData cellData(numCells,baseCellDim,topo);

   // build basis
   RCP<const panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,cellData));

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
   scatterParams.set("Basis",basis);

   // test residual scatter evaluator
   {
      panzer::ScatterResidual_Epetra<Residual,panzer::Traits,int,int> scatterResidual(Teuchos::null,scatterParams);

      const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual.evaluatedFields();
      TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager

      const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual.dependentFields();
      TEST_EQUALITY(fields.size(),2); // these store the residual values

      TEST_EQUALITY(fields[0]->name(),"Residual_ux");
      TEST_EQUALITY(fields[1]->name(),"Residual_p");

      TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   }

   // test jacobian scatter evaluator
   {
      panzer::ScatterResidual_Epetra<Jacobian,panzer::Traits,int,int> scatterJacobian(Teuchos::null,scatterParams);

      const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterJacobian.evaluatedFields();
      TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager

      const std::vector<RCP<PHX::FieldTag> > & fields = scatterJacobian.dependentFields();
      TEST_EQUALITY(fields.size(),2); // these store the residual values

      TEST_EQUALITY(fields[0]->name(),"Residual_ux");
      TEST_EQUALITY(fields[1]->name(),"Residual_p");

      TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   }
}

}
