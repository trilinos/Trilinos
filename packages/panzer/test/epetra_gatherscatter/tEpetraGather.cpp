#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>
#include <vector>

#include "Panzer_GatherSolution_Epetra.hpp"

#include "Panzer_Basis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {

TEUCHOS_UNIT_TEST(tEpetraGather, constructor)
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   typedef panzer::Traits::Residual Residual;
   typedef panzer::Traits::Jacobian Jacobian;

   // auxiliary information needed to construct basis object
   std::size_t numCells = 10;
   int baseCellDim = 2;
   int cubatureDegree = 2;
   std::string basisType = "Q1";
   panzer::CellData cellData(numCells,baseCellDim);
   panzer::IntegrationRule intRule(cubatureDegree,cellData);

   // build DOF names
   RCP<std::vector<std::string> > dofNames = rcp(new std::vector<std::string>); 
   dofNames->push_back("ux"); // in practice these probably would not be gathered together!
   dofNames->push_back("p");

   // build basis
   RCP<panzer::Basis> basis = rcp(new panzer::Basis(basisType,intRule));

   // build gather parameter list
   Teuchos::ParameterList gatherParams;
   gatherParams.set<RCP<std::vector<std::string> > >("DOF Names",dofNames);
   gatherParams.set<RCP<panzer::Basis> >("Basis",basis);

   // test residual gather evaluator
   {
      panzer::GatherSolution_Epetra<Residual,panzer::Traits> gatherResidual(gatherParams);

      const std::vector<RCP<PHX::FieldTag> > & fields = gatherResidual.evaluatedFields();
      TEST_EQUALITY(fields.size(),2);
 
      TEST_EQUALITY(fields[0]->name(),"ux");
      TEST_EQUALITY(fields[1]->name(),"p");

      TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   }

   // test jacobian gather evaluator
   {
      panzer::GatherSolution_Epetra<Jacobian,panzer::Traits> gatherJacobian(gatherParams);

      const std::vector<RCP<PHX::FieldTag> > & fields = gatherJacobian.evaluatedFields();
      TEST_EQUALITY(fields.size(),2);
 
      TEST_EQUALITY(fields[0]->name(),"ux");
      TEST_EQUALITY(fields[1]->name(),"p");

      TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

      TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
      TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   }
}

}
