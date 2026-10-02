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
#include "Panzer_ElemFieldPattern.hpp"
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
    RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
    RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
    return pattern;
  }

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

    STKConnManager::cacheConnectivity();
    STKConnManager connMngr(mesh);
    {
      connMngr.buildConnectivity(*fp);
      // Build with a new field pattern to test caching
      connMngr.buildConnectivity(*fp); // reuse cached
      auto fp_dummy = Teuchos::make_rcp<panzer::ElemFieldPattern>(fp->getCellTopology());
      connMngr.buildConnectivity(*fp_dummy);
      connMngr.buildConnectivity(*fp); // reuse cached
      connMngr.buildConnectivity(*fp); // reuse cached
    }

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
    STKConnManager::clearCachedConnectivityData();
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),3);
  }

  TEUCHOS_UNIT_TEST(tSTKConnManager, cache_key_mesh_and_sidesets)
  {
    using Teuchos::RCP;

    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(numProcs<=2);

    RCP<STK_Interface> meshA = build2DMesh(2,1,2,1);
    RCP<STK_Interface> meshB = build2DMesh(4,1,2,1);

    RCP<const panzer::FieldPattern> fp
      = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double> >();

    STKConnManager::cacheConnectivity();
    const int startCount = STKConnManager::getCachedReuseCount();

    // A different mesh with the same field pattern must not reuse cached data
    STKConnManager cmA(meshA);
    cmA.buildConnectivity(*fp);
    STKConnManager cmB(meshB);
    cmB.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount);

    // Same mesh and field pattern reuses cached data
    STKConnManager cmB2(meshB);
    cmB2.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+1);
    {
      std::vector<stk::mesh::Entity> myElementsA, myElementsB;
      meshA->getMyElements(myElementsA);
      meshB->getMyElements(myElementsB);
      TEST_EQUALITY(cmA.getOwnedElementCount(),myElementsA.size());
      TEST_EQUALITY(cmB.getOwnedElementCount(),myElementsB.size());
      TEST_EQUALITY(cmB2.getOwnedElementCount(),myElementsB.size());
    }

    // A different sideset association list must not reuse cached data
    STKConnManager cmC(meshA);
    cmC.associateElementsInSideset("left");
    cmC.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+1);

    // Same mesh, field pattern and sideset list reuses cached data
    STKConnManager cmD(meshA);
    cmD.associateElementsInSideset("left");
    cmD.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+2);

    // Clearing meshA releases its two cache entries and leaves meshB cached
    const int strongCountABefore = meshA.strong_count();
    const int strongCountBBefore = meshB.strong_count();
    STKConnManager::clearCachedConnectivityData(meshA);
    TEST_EQUALITY(meshA.strong_count(),strongCountABefore-2);
    TEST_EQUALITY(meshB.strong_count(),strongCountBBefore);

    STKConnManager cmE(meshA);
    cmE.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+2);

    STKConnManager cmF(meshB);
    cmF.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+3);

    // Clearing meshB releases its one cache entry and leaves meshA's rebuilt entry cached
    {
      const int countA = meshA.strong_count();
      const int countB = meshB.strong_count();
      STKConnManager::clearCachedConnectivityData(meshB);
      TEST_EQUALITY(meshA.strong_count(),countA);
      TEST_EQUALITY(meshB.strong_count(),countB-1);

      // Clearing a mesh with no cache entries is a no-op
      STKConnManager::clearCachedConnectivityData(meshB);
      TEST_EQUALITY(meshA.strong_count(),countA);
      TEST_EQUALITY(meshB.strong_count(),countB-1);
    }

    STKConnManager cmG(meshB);
    cmG.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+3);

    STKConnManager cmH(meshA);
    cmH.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+4);

    STKConnManager::clearCachedConnectivityData();
  }
}
