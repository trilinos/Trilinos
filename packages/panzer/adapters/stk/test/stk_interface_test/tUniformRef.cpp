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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

/*
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_percept/PerceptMesh.hpp"
#include "stk_adapt/UniformRefiner.hpp"
#include "stk_adapt/UniformRefinerPattern.hpp"
*/
#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace panzer_stk {

/*
TEUCHOS_UNIT_TEST(tUniformRef, stk_fixture)
{
   typedef stk::mesh::Field<int> ProcIdFieldType;

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",3);

   SquareQuadMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

   // Step 1: Build Percept data structures and pattern
   stk::mesh::fem::FEMMetaData * metaData = &*mesh->getMetaData();

   mesh->instantiateBulkData(MPI_COMM_WORLD);

   stk::mesh::BulkData * bulkData = &*mesh->getBulkData();

   stk::percept::PerceptMesh perceptMesh(metaData,bulkData,false);
    
   const std::string refine="DEFAULT";
   const std::string enrich="";
   const std::string convert="";
   
   stk::adapt::BlockNamesType block_names;
   
   Teuchos::RCP<stk::adapt::UniformRefinerPatternBase> uniformRefinePattern
      = stk::adapt::UniformRefinerPatternBase::createPattern(refine, 
                                                             enrich, 
                                                             convert, 
                                                             perceptMesh, 
                                                             block_names);

   // end Step 1

   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   mesh->writeToExodus("oldmesh.exo");

   // Step 2: Do uniform refinement

   stk::adapt::UniformRefiner breaker(perceptMesh, *uniformRefinePattern, mesh->getProcessorIdField());
   breaker.setDoProgressMeter(true);

   int num_uniform_refines = 2;
   for (int iBreak = 0; iBreak < num_uniform_refines; iBreak++) {
     out << "  ref level " << iBreak +1 << std::endl;
     breaker.doBreak();   
     stk::adapt::RefinementInfoByType::printTable(
       out, 
       breaker.getRefinementInfoByType(), 
       iBreak , 
       true);
   }

   // end Step 2

   mesh->writeToExodus("newmesh.exo");
}
*/

TEUCHOS_UNIT_TEST(tUniformRef, stk_exodus)
{
   typedef stk::mesh::Field<int> ProcIdFieldType;

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/basic.gen");

   STK_ExodusReaderFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

   // Step 1: Build Percept data structures and pattern
   stk::mesh::fem::FEMMetaData * metaData = &*mesh->getMetaData();

   mesh->instantiateBulkData(MPI_COMM_WORLD);

   stk::mesh::BulkData * bulkData = &*mesh->getBulkData();

   stk::percept::PerceptMesh perceptMesh(metaData,bulkData,false);
    
   const std::string refine="DEFAULT";
   const std::string enrich="";
   const std::string convert="";
   
   stk::adapt::BlockNamesType block_names;
   
   Teuchos::RCP<stk::adapt::UniformRefinerPatternBase> uniformRefinePattern
      = stk::adapt::UniformRefinerPatternBase::createPattern(refine, 
                                                             enrich, 
                                                             convert, 
                                                             perceptMesh, 
                                                             block_names);

   // end Step 1

   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   mesh->writeToExodus("oldmesh.exo");

   // Step 2: Do uniform refinement

   stk::adapt::UniformRefiner breaker(perceptMesh, *uniformRefinePattern, mesh->getProcessorIdField());
   breaker.setDoProgressMeter(true);

   int num_uniform_refines = 2;
   for (int iBreak = 0; iBreak < num_uniform_refines; iBreak++) {
     out << "  ref level " << iBreak +1 << std::endl;
     breaker.doBreak();   
     stk::adapt::RefinementInfoByType::printTable(
       out, 
       breaker.getRefinementInfoByType(), 
       iBreak , 
       true);
   }

   // end Step 2

   mesh->writeToExodus("newmesh.exo");
}

}
