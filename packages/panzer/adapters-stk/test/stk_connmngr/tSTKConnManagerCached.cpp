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

#include <numeric>
#include <type_traits>

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
    RCP<Intrepid2::Basis<PHX::Device,double,double> > basis = rcp(new Intrepid2Type);
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
      = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::Device,double,double> >();

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
      = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::Device,double,double> >();

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

  // The flat connectivity array must hold exactly the IDs of the current build
  bool connectivityViewIsConsistent(const STKConnManager& cm)
  {
    const auto sizes = cm.getConnectivitySizeView();
    const auto offsets = cm.getElementLidToConnView();
    const auto conn = cm.getConnectivityView();
    STKConnManager::LocalOrdinal total = 0;
    for (std::size_t e=0; e<sizes.extent(0); ++e) {
      if (offsets(e) != total) return false;
      total += sizes(e);
    }
    return static_cast<std::size_t>(total) == conn.extent(0);
  }

  TEUCHOS_UNIT_TEST(tSTKConnManager, cache_shares_data)
  {
    using Teuchos::RCP;
    using GO = STKConnManager::GlobalOrdinal;

    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(numProcs<=2);

    // Shared data must only be exposed read-only
    static_assert(std::is_const<std::remove_pointer_t<decltype(std::declval<STKConnManager>().getConnectivityView().data())>>::value,
                  "getConnectivityView() must return a const view");
    static_assert(std::is_const<std::remove_pointer_t<decltype(std::declval<STKConnManager>().getConnectivitySizeView().data())>>::value,
                  "getConnectivitySizeView() must return a const view");
    static_assert(std::is_const<std::remove_pointer_t<decltype(std::declval<STKConnManager>().getElementLidToConnView().data())>>::value,
                  "getElementLidToConnView() must return a const view");

    RCP<STK_Interface> mesh = build2DMesh(2,1,2,1);
    RCP<const panzer::FieldPattern> fp
      = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::Device,double,double> >();
    auto fp_dummy = Teuchos::make_rcp<panzer::ElemFieldPattern>(fp->getCellTopology());

    STKConnManager::cacheConnectivity();
    const int startCount = STKConnManager::getCachedReuseCount();

    STKConnManager cm1(mesh);
    cm1.buildConnectivity(*fp);
    TEST_ASSERT(connectivityViewIsConsistent(cm1));
    const auto view1 = cm1.getConnectivityView();
    const std::vector<GO> reference(view1.data(),view1.data()+view1.extent(0));

    // A cache hit shares storage instead of copying it
    STKConnManager cm2(mesh);
    cm2.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+1);
    TEST_EQUALITY(cm2.getConnectivity(0),cm1.getConnectivity(0));
    TEST_EQUALITY(cm2.getConnectivityView().data(),cm1.getConnectivityView().data());
    TEST_EQUALITY(cm2.getConnectivitySizeView().data(),cm1.getConnectivitySizeView().data());
    TEST_EQUALITY(cm2.getElementLidToConnView().data(),cm1.getElementLidToConnView().data());

    // Rebuilding after a hit allocates new storage and leaves shared data untouched
    cm2.buildConnectivity(*fp_dummy);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+1);
    TEST_ASSERT(cm2.getConnectivityView().data() != cm1.getConnectivityView().data());
    TEST_ASSERT(connectivityViewIsConsistent(cm2));
    {
      const auto view = cm1.getConnectivityView();
      TEST_ASSERT(std::vector<GO>(view.data(),view.data()+view.extent(0)) == reference);
    }

    STKConnManager cm3(mesh);
    cm3.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+2);
    {
      const auto view = cm3.getConnectivityView();
      TEST_ASSERT(std::vector<GO>(view.data(),view.data()+view.extent(0)) == reference);
    }

    // Repeated builds on one manager must not accumulate stale IDs
    RCP<STK_Interface> mesh2 = build2DMesh(2,1,2,1);
    STKConnManager cm4(mesh2);
    cm4.buildConnectivity(*fp);
    cm4.buildConnectivity(*fp_dummy);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+2);
    TEST_ASSERT(connectivityViewIsConsistent(cm4));

    STKConnManager::clearCachedConnectivityData();
  }

  TEUCHOS_UNIT_TEST(tSTKConnManager, cache_key_sideset_order)
  {
    using Teuchos::RCP;

    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(numProcs<=2);

    RCP<STK_Interface> mesh = build2DMesh(2,1,2,1);
    RCP<const panzer::FieldPattern> fp
      = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::Device,double,double> >();

    STKConnManager::cacheConnectivity();
    const int startCount = STKConnManager::getCachedReuseCount();

    STKConnManager cmA(mesh);
    cmA.associateElementsInSideset("left");
    cmA.associateElementsInSideset("right");
    cmA.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount);

    // The same sidesets in a different order are a distinct cache key
    STKConnManager cmB(mesh);
    cmB.associateElementsInSideset("right");
    cmB.associateElementsInSideset("left");
    cmB.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount);

    // A list with a subset of the cached sidesets is also a distinct cache key
    STKConnManager cmC(mesh);
    cmC.associateElementsInSideset("left");
    cmC.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount);

    // Each order reuses its own cached entry
    STKConnManager cmD(mesh);
    cmD.associateElementsInSideset("left");
    cmD.associateElementsInSideset("right");
    cmD.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+1);
    TEST_EQUALITY(cmD.getConnectivity(0),cmA.getConnectivity(0));

    STKConnManager cmE(mesh);
    cmE.associateElementsInSideset("right");
    cmE.associateElementsInSideset("left");
    cmE.buildConnectivity(*fp);
    TEST_EQUALITY(STKConnManager::getCachedReuseCount(),startCount+2);
    TEST_EQUALITY(cmE.getConnectivity(0),cmB.getConnectivity(0));

    STKConnManager::clearCachedConnectivityData();
  }
}
