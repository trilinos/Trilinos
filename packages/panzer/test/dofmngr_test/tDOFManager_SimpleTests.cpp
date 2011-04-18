#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "dofmngr/Panzer_DOFManager.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"

// 3D basis 
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

TEUCHOS_UNIT_TEST(tDOFManager_SimpleTests,validFieldOrder)
{
   DOFManager<int,int> dofManager; 

   std::set<std::string> validFields;
   validFields.insert("horse");
   validFields.insert("cat");
   validFields.insert("monkey");
   validFields.insert("dog");

   {
      std::vector<std::string> order;
      order.push_back("cat");
      order.push_back("dog");
      order.push_back("horse");
      order.push_back("monkey");
 
      TEST_ASSERT(dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::string> order;
      order.push_back("cat");
      order.push_back("horse");
      order.push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::string> order;
      order.push_back("cat");
      order.push_back("dog");
      order.push_back("horse");
      order.push_back("monkey");
      order.push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::string> order;
      order.push_back("cat");
      order.push_back("dog");
      order.push_back("horse");
      order.push_back("tank");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }
}


}
