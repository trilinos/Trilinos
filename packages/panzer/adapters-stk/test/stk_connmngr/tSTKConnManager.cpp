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
#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

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

template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

RCP<const panzer::FieldPattern> buildConstantFieldPattern(const shards::CellTopology & ct)
{
   typedef Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double> Intrepid2Type;

   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type(ct));
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
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
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double> >();

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

   if(numProcs==1) {
       TEST_EQUALITY(connMngr.getNeighborElementBlock("eblock-0_0").size(),0);
       TEST_EQUALITY(connMngr.getNeighborElementBlock("eblock-1_0").size(),0);
   }
   else {
       TEST_EQUALITY(connMngr.getNeighborElementBlock("eblock-0_0").size(),1);
       TEST_EQUALITY(connMngr.getNeighborElementBlock("eblock-1_0").size(),1);

       for(std::size_t blk=0;blk<connMngr.numElementBlocks();++blk) {
         const std::vector<int> & elementBlock = connMngr.getNeighborElementBlock(elementBlockIds[blk]);
         for(std::size_t elmt=0;elmt<elementBlock.size();++elmt)
           TEST_EQUALITY(connMngr.getConnectivitySize(elementBlock[elmt]),9);
       }
   }

   STKConnManager::GlobalOrdinal maxEdgeId = mesh->getMaxEntityId(mesh->getEdgeRank());
   STKConnManager::GlobalOrdinal nodeCount = mesh->getEntityCounts(mesh->getNodeRank());

   if(numProcs==1) {
      const auto * conn1 = connMngr.getConnectivity(1);
      const auto * conn2 = connMngr.getConnectivity(2);
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
      const auto * conn0 = connMngr.getConnectivity(0);
      const auto * conn1 = connMngr.getConnectivity(1);

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

      const auto * conn2 = connMngr.getConnectivity(2); // this is the "neighbor element"
      const auto * conn3 = connMngr.getConnectivity(3); // this is the "neighbor element"

      int otherRank = myRank==0 ? 1 : 0;

      TEST_EQUALITY(conn2[0],0+otherRank);
      TEST_EQUALITY(conn2[1],1+otherRank);
      TEST_EQUALITY(conn2[2],6+otherRank);
      TEST_EQUALITY(conn2[3],5+otherRank);

      TEST_EQUALITY(conn3[0],2+otherRank);
      TEST_EQUALITY(conn3[1],3+otherRank);
      TEST_EQUALITY(conn3[2],8+otherRank);
      TEST_EQUALITY(conn3[3],7+otherRank);

      TEST_EQUALITY(conn2[8],nodeCount+(maxEdgeId+1)+1+otherRank);
      TEST_EQUALITY(conn3[8],nodeCount+(maxEdgeId+1)+3+otherRank);
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
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

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

      const auto * conn = connMngr.getConnectivity(localId);
      for(std::size_t i=0;(int) i<connMngr.getConnectivitySize(localId);i++)
         TEST_EQUALITY(conn[i],conn_true[i]-1);
   }
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKConnManager, noConnectivityClone)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   TEUCHOS_ASSERT(numProcs<=2);

   RCP<STK_Interface> mesh = build2DMesh(5,5,1,1);
   TEST_ASSERT(mesh!=Teuchos::null);

   RCP<shards::CellTopology> ct = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<QuadTopo>()));

   RCP<const panzer::FieldPattern> fp_hgrad
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
   RCP<const panzer::FieldPattern> fp_const
         = buildConstantFieldPattern(*ct);

   STKConnManager connMngr_const(mesh);
   connMngr_const.buildConnectivity(*fp_const);

   RCP<STKConnManager> connMngr_hgrad = rcp_dynamic_cast<STKConnManager>(connMngr_const.noConnectivityClone());
   TEST_ASSERT(connMngr_hgrad!=Teuchos::null);
   connMngr_hgrad->buildConnectivity(*fp_hgrad);

   // test constant conn manager
   {
      // did we get the element block correct?
      /////////////////////////////////////////////////////////////

      TEST_EQUALITY(connMngr_const.numElementBlocks(),1);

      const std::vector<int> & elementBlock = connMngr_const.getElementBlock("eblock-0_0");
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

      TEST_EQUALITY(connMngr_const.getBlockId(9),"eblock-0_0");

      // test connectivities
      /////////////////////////////////////////////////////////////
      TEST_EQUALITY(connMngr_const.getConnectivitySize(9),1);
      TEST_EQUALITY(connMngr_const.getConnectivitySize(8),1);


      std::size_t localId;
      if(myRank==0)
         localId = mesh->elementLocalId(17);
      else
         localId = mesh->elementLocalId(20);

      {
         int conn_true[1];
         if(numProcs==1) {
            conn_true[0] = 16;
         }
         else if(numProcs==2 && myRank==0) {
            conn_true[0] = 16;
         }
         else if(numProcs==2 && myRank==1) {
            conn_true[0] = 19;
         }
         else {
            TEST_ASSERT(false);
         }

         const auto * conn = connMngr_const.getConnectivity(localId);
         for(std::size_t i=0;(int) i<connMngr_const.getConnectivitySize(localId);i++)
            TEST_EQUALITY(conn[i],conn_true[i]);
      }
   }

   // test hgrad conn manager
   {
      // did we get the element block correct?
      /////////////////////////////////////////////////////////////

      TEST_EQUALITY(connMngr_hgrad->numElementBlocks(),1);

      const std::vector<int> & elementBlock = connMngr_hgrad->getElementBlock("eblock-0_0");
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

      TEST_EQUALITY(connMngr_hgrad->getBlockId(9),"eblock-0_0");

      // test connectivities
      /////////////////////////////////////////////////////////////
      TEST_EQUALITY(connMngr_hgrad->getConnectivitySize(9),4);
      TEST_EQUALITY(connMngr_hgrad->getConnectivitySize(8),4);


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

         const auto * conn = connMngr_hgrad->getConnectivity(localId);
         for(std::size_t i=0;(int) i<connMngr_hgrad->getConnectivitySize(localId);i++)
            TEST_EQUALITY(conn[i],conn_true[i]-1);
      }
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
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

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

namespace {
void testAssociatedNeighbors(const STKConnManager& connMngr,
  const std::vector<std::vector<int> > vals, Teuchos::FancyOStream& out,
  bool& success)
{
  for (int i = 0; i < static_cast<int>(vals.size()); ++i)
  {
    const int sz(connMngr.getAssociatedNeighbors(vals[i][0]).size());
    TEST_EQUALITY(sz, vals[i][1]);
    if (sz)
      TEST_EQUALITY(connMngr.getAssociatedNeighbors(vals[i][0])[0], vals[i][2]);
  }
}
}

TEUCHOS_UNIT_TEST(tSTKConnManager, 2_blocks_interface)
{
   using Teuchos::RCP;

   const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
   const int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   const RCP<STK_Interface> mesh = build2DMesh(2,1,2,1);
   TEST_ASSERT( ! mesh.is_null());

   RCP<const panzer::FieldPattern>
     fp = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double> >();

   STKConnManager connMngr(mesh);
   connMngr.associateElementsInSideset("vertical_0");
   connMngr.associateElementsInSideset("left");
   connMngr.buildConnectivity(*fp);
   {
     Teuchos::RCP<Teuchos::Comm<int> >
       comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD));
     std::vector<std::string> ss = connMngr.checkAssociateElementsInSidesets(*comm);
     TEST_EQUALITY(ss.size(), 1);
     TEST_EQUALITY(ss[0], "left");
   }

   if (numProcs == 1) {
     const std::vector<std::vector<int> > vals{{0, 0, 0}, {1, 1, 2}, {2, 1, 1}, {3, 0, 0}};
     testAssociatedNeighbors(connMngr, vals, out, success);
   } else if (numProcs == 2 && myRank == 0) {
     const std::vector<std::vector<int> > vals{{0, 0, 0}, {1, 1, 2}};
     testAssociatedNeighbors(connMngr, vals, out, success);
   } else if (numProcs == 2 && myRank == 1) {
     const std::vector<std::vector<int> > vals{{0, 1, 3}, {1, 0, 0}};
     testAssociatedNeighbors(connMngr, vals, out, success);
   }
   else {
     // We'll not check any other cases here. Interface connection
     // correctness is tested in much greater detail in
     //     adapters/stk/example/PoissonInterfaceExample/main.cpp.
   }
}

}
