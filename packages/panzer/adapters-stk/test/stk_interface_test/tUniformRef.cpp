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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
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
   stk::mesh::MetaData * metaData = &*mesh->getMetaData();

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
   stk::mesh::MetaData * metaData = &*mesh->getMetaData();

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
