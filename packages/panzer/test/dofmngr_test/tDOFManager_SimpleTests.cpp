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
#include <set>

#include "Panzer_DOFManager.hpp"

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
