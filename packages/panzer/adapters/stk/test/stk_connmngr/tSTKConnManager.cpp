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
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Shards_BasicTopologies.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

using Teuchos::RCP;
using Teuchos::rcp;

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer_stk {

typedef shards::Quadrilateral<4> QuadTopo;

Teuchos::RCP<STK_Interface> build2DMesh(int xElements,int yElements,int xBlocks,int yBlocks)
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",xBlocks);
   pl->set("Y Blocks",yBlocks);
   pl->set("X Elements",xElements);
   pl->set("Y Elements",yElements);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   
   Teuchos::RCP<STK_Interface> meshPtr = factory.buildMesh(MPI_COMM_WORLD);

   return meshPtr;
}

template <typename IntrepidType>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
   return pattern;
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKConnManager, 2_blocks)
{
   using Teuchos::RCP;

   int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   TEUCHOS_ASSERT(numProcs<=2);

   RCP<STK_Interface> mesh = build2DMesh(2,1,2,1);
   TEST_ASSERT(mesh!=Teuchos::null);

   RCP<const panzer::FieldPattern> fp 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer> >();

   STKConnManager connMngr(mesh);
   connMngr.buildConnectivity(*fp);

   // did we get the element block correct?
   /////////////////////////////////////////////////////////////

   TEST_EQUALITY(connMngr.numElementBlocks(),2);
   TEST_EQUALITY(connMngr.getBlockId(0),"eblock-0_0");

   // check that each element is correct size
   std::vector<std::string> elementBlockIds;
   connMngr.getElementBlockIds(elementBlockIds);
   for(std::size_t blk=0;blk<connMngr.numElementBlocks();++blk) {
      std::string blockId = elementBlockIds[blk];
      const std::vector<int> & elementBlock = connMngr.getElementBlock(blockId);
      for(std::size_t elmt=0;elmt<elementBlock.size();++elmt)
         TEST_EQUALITY(connMngr.getConnectivitySize(elementBlock[elmt]),9); 
   }

   STKConnManager::GlobalOrdinal maxEdgeId = mesh->getMaxEntityId(mesh->getEdgeRank());
   STKConnManager::GlobalOrdinal nodeCount = mesh->getEntityCounts(mesh->getNodeRank());

   if(numProcs==1) {
      const int * conn1 = connMngr.getConnectivity(1);
      const int * conn2 = connMngr.getConnectivity(2);
      TEST_EQUALITY(conn1[0],1);  
      TEST_EQUALITY(conn1[1],2);  
      TEST_EQUALITY(conn1[2],7);  
      TEST_EQUALITY(conn1[3],6);  
   
      TEST_EQUALITY(conn2[0],2);  
      TEST_EQUALITY(conn2[1],3);  
      TEST_EQUALITY(conn2[2],8);  
      TEST_EQUALITY(conn2[3],7);  
   
      TEST_EQUALITY(conn1[5],conn2[7]);
 
      TEST_EQUALITY(conn1[8],nodeCount+(maxEdgeId+1)+2);
      TEST_EQUALITY(conn2[8],nodeCount+(maxEdgeId+1)+3);

   }
   else {
      const int * conn0 = connMngr.getConnectivity(0);
      const int * conn1 = connMngr.getConnectivity(1);
      TEST_EQUALITY(conn0[0],0+myRank);  
      TEST_EQUALITY(conn0[1],1+myRank);  
      TEST_EQUALITY(conn0[2],6+myRank);  
      TEST_EQUALITY(conn0[3],5+myRank);  

      TEST_EQUALITY(conn1[0],2+myRank);  
      TEST_EQUALITY(conn1[1],3+myRank);  
      TEST_EQUALITY(conn1[2],8+myRank);  
      TEST_EQUALITY(conn1[3],7+myRank);  

      TEST_EQUALITY(conn0[8],nodeCount+(maxEdgeId+1)+1+myRank);
      TEST_EQUALITY(conn1[8],nodeCount+(maxEdgeId+1)+3+myRank);
   }
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKConnManager, single_block_2d)
{
   using Teuchos::RCP;

   int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   TEUCHOS_ASSERT(numProcs<=2);

   RCP<STK_Interface> mesh = build2DMesh(5,5,1,1);
   TEST_ASSERT(mesh!=Teuchos::null);

   RCP<const panzer::FieldPattern> fp 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   STKConnManager connMngr(mesh);
   connMngr.buildConnectivity(*fp);

   // did we get the element block correct?
   /////////////////////////////////////////////////////////////

   TEST_EQUALITY(connMngr.numElementBlocks(),1);

   const std::vector<int> & elementBlock = connMngr.getElementBlock("eblock-0_0");
   std::vector<int> nc_elementBlock = elementBlock;
   if(numProcs==1)                   { TEST_EQUALITY(elementBlock.size(),25); }
   else if(numProcs==2 && myRank==0) { TEST_EQUALITY(elementBlock.size(),15); }
   else if(numProcs==2 && myRank==1) { TEST_EQUALITY(elementBlock.size(),10); }
   else                              { TEST_ASSERT(false); }

   // check that the local elements are correctly numbered
   std::sort(nc_elementBlock.begin(),nc_elementBlock.end());
   bool check_local_blocks_passed = true;
   for(std::size_t i=0;i<elementBlock.size();i++)
      check_local_blocks_passed &= (nc_elementBlock[i]==(int) i);
   TEST_ASSERT(check_local_blocks_passed);

   TEST_EQUALITY(connMngr.getBlockId(9),"eblock-0_0");

   // test connectivities
   /////////////////////////////////////////////////////////////
   TEST_EQUALITY(connMngr.getConnectivitySize(9),4); 
   TEST_EQUALITY(connMngr.getConnectivitySize(8),4); 

  
   std::size_t localId;
   if(myRank==0)
      localId = mesh->elementLocalId(17);
   else
      localId = mesh->elementLocalId(20);
   
   {
      int conn_true[4]; 
      if(numProcs==1) {
         conn_true[0] = 20;
         conn_true[1] = 21;
         conn_true[2] = 27;
         conn_true[3] = 26;
      }
      else if(numProcs==2 && myRank==0) {
         conn_true[0] = 20;
         conn_true[1] = 21;
         conn_true[2] = 27;
         conn_true[3] = 26;
      }
      else if(numProcs==2 && myRank==1) {
         conn_true[0] = 23;
         conn_true[1] = 24;
         conn_true[2] = 30;
         conn_true[3] = 29;
      }
      else {
         TEST_ASSERT(false); 
      }

      const int * conn = connMngr.getConnectivity(localId);
      for(std::size_t i=0;(int) i<connMngr.getConnectivitySize(localId);i++)
         TEST_EQUALITY(conn[i],conn_true[i]-1);
   }
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKConnManager, four_block_2d)
{
   using Teuchos::RCP;

   int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
   // int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   RCP<STK_Interface> mesh = build2DMesh(2,2,2,2); // 4x4 elements
   TEST_ASSERT(mesh!=Teuchos::null);

   RCP<const panzer::FieldPattern> fp 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   STKConnManager connMngr(mesh);
   connMngr.buildConnectivity(*fp);

   TEUCHOS_ASSERT(numProcs<=2);

   // did we get the element block correct?
   /////////////////////////////////////////////////////////////

   TEST_EQUALITY(connMngr.numElementBlocks(),4);

   std::vector<std::string> elementBlockIds;
   connMngr.getElementBlockIds(elementBlockIds);
   for(std::size_t blk=0;blk<connMngr.numElementBlocks();blk++) {
      std::string blockId = elementBlockIds[blk];
      const std::vector<int> & elementBlock = connMngr.getElementBlock(blockId);
      std::vector<int> nc_elementBlock = elementBlock;
      if(numProcs==1)      { TEST_EQUALITY(elementBlock.size(),4); }
      else if(numProcs==2) { TEST_EQUALITY(elementBlock.size(),2); }
   
      bool check_blockid_lookup = true;
      for(std::size_t i=0;i<elementBlock.size();i++)
         check_blockid_lookup &= (connMngr.getBlockId(elementBlock[i])==blockId);
      TEST_ASSERT(check_blockid_lookup); 
   }

   // 
   /////////////////////////////////////////////////////////////
}

}
