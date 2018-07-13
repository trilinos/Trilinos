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
#include <Teuchos_DefaultMpiComm.hpp>

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Kokkos_DynRankView.hpp"

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "CartesianConnManager.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace unit_test {

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;
template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
  // build a geometric pattern from a single basis
  RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2Type);
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}

// test that you can correctly compute a rank index from a global processor id
TEUCHOS_UNIT_TEST(tCartesianTop, computeMyRankTriplet)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<int> Triplet;

  // test 2D
  {
    auto t0 = CCM::computeMyRankTriplet(3,2,Triplet(4,3,1));
    auto t1 = CCM::computeMyRankTriplet(6,2,Triplet(4,3,1));
    auto t2 = CCM::computeMyRankTriplet(9,2,Triplet(4,3,1));

    TEST_EQUALITY(t0.x,3); TEST_EQUALITY(t0.y,0); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,2); TEST_EQUALITY(t1.y,1); TEST_EQUALITY(t1.z,0);
    TEST_EQUALITY(t2.x,1); TEST_EQUALITY(t2.y,2); TEST_EQUALITY(t2.z,0);
  }

  // test 3D
  {
    auto t0 = CCM::computeMyRankTriplet( 3,3,Triplet(2,3,3));
    auto t1 = CCM::computeMyRankTriplet(10,3,Triplet(2,3,3));
    auto t2 = CCM::computeMyRankTriplet(13,3,Triplet(2,3,3));

    TEST_EQUALITY(t0.x,1); TEST_EQUALITY(t0.y,1); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,0); TEST_EQUALITY(t1.y,2); TEST_EQUALITY(t1.z,1);
    TEST_EQUALITY(t2.x,1); TEST_EQUALITY(t2.y,0); TEST_EQUALITY(t2.z,2);
  }
}

// test that you can correctly compute a rank index from a grobal processor id
TEUCHOS_UNIT_TEST(tCartesianTop, computeLocalElementGlobalTriplet)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> Triplet;

  // test 2D
  {
    auto t0 = CCM::computeLocalElementGlobalTriplet( 1,Triplet(4,3,1),Triplet(7,2,0));
    auto t1 = CCM::computeLocalElementGlobalTriplet( 4,Triplet(4,3,1),Triplet(7,2,0));
    auto t2 = CCM::computeLocalElementGlobalTriplet(11,Triplet(4,3,1),Triplet(7,2,0));

    TEST_EQUALITY(t0.x,1+7); TEST_EQUALITY(t0.y,0+2); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,0+7); TEST_EQUALITY(t1.y,1+2); TEST_EQUALITY(t1.z,0);
    TEST_EQUALITY(t2.x,3+7); TEST_EQUALITY(t2.y,2+2); TEST_EQUALITY(t2.z,0);
  }

  // test 3D
  {
    auto t0 = CCM::computeLocalElementGlobalTriplet( 1+2*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t1 = CCM::computeLocalElementGlobalTriplet( 4+1*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t2 = CCM::computeLocalElementGlobalTriplet(11+4*12,Triplet(4,3,5),Triplet(7,2,3));

    TEST_EQUALITY(t0.x,1+7); TEST_EQUALITY(t0.y,0+2); TEST_EQUALITY(t0.z,2+3);
    TEST_EQUALITY(t1.x,0+7); TEST_EQUALITY(t1.y,1+2); TEST_EQUALITY(t1.z,1+3);
    TEST_EQUALITY(t2.x,3+7); TEST_EQUALITY(t2.y,2+2); TEST_EQUALITY(t2.z,4+3);
  }
}

TEUCHOS_UNIT_TEST(tCartesianTop, computeLocalElementIndex)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> Triplet;

  // test 2D
  {
    auto t0 = CCM::computeLocalElementGlobalTriplet( 1,Triplet(4,3,1),Triplet(7,2,0));
    auto t1 = CCM::computeLocalElementGlobalTriplet( 4,Triplet(4,3,1),Triplet(7,2,0));
    auto t2 = CCM::computeLocalElementGlobalTriplet(11,Triplet(4,3,1),Triplet(7,2,0));

    TEST_EQUALITY(CCM::computeLocalElementIndex(t0,Triplet(4,3,1),Triplet(7,2,0)), 1);
    TEST_EQUALITY(CCM::computeLocalElementIndex(t1,Triplet(4,3,1),Triplet(7,2,0)), 4);
    TEST_EQUALITY(CCM::computeLocalElementIndex(t2,Triplet(4,3,1),Triplet(7,2,0)),11);
  }

  // test 3D
  {
    auto t0 = CCM::computeLocalElementGlobalTriplet( 1+2*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t1 = CCM::computeLocalElementGlobalTriplet( 4+1*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t2 = CCM::computeLocalElementGlobalTriplet(11+4*12,Triplet(4,3,5),Triplet(7,2,3));

    TEST_EQUALITY(CCM::computeLocalElementIndex(t0,Triplet(4,3,5),Triplet(7,2,3)), 1+2*12);
    TEST_EQUALITY(CCM::computeLocalElementIndex(t1,Triplet(4,3,5),Triplet(7,2,3)), 4+1*12);
    TEST_EQUALITY(CCM::computeLocalElementIndex(t2,Triplet(4,3,5),Triplet(7,2,3)),11+4*12);
  }
}

TEUCHOS_UNIT_TEST(tCartesianTop, computeGlobalElementIndex)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> Triplet;

  // test 2D
  {
    Ordinal64 t0 = CCM::computeGlobalElementIndex(Triplet(1,1,0),Triplet(4,3,1));
    Ordinal64 t1 = CCM::computeGlobalElementIndex(Triplet(3,2,0),Triplet(4,3,1));

    TEST_EQUALITY(t0, 5);
    TEST_EQUALITY(t1,11);
  }

  // test 3D
  {
    Ordinal64 t0 = CCM::computeGlobalElementIndex(Triplet( 1,1,3),Triplet(4,3,6));
    Ordinal64 t1 = CCM::computeGlobalElementIndex(Triplet( 3,2,5),Triplet(4,3,6));

    TEST_EQUALITY(t0, 5+3*12);
    TEST_EQUALITY(t1,11+5*12);
  }
}

// This test checks functions used to generate the topology are correct
TEUCHOS_UNIT_TEST(tCartesianTop, connmanager_2d_1dpart_helpers)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> TripletGO;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  int rank = comm.getRank();

  // field pattern for basis required
  RCP<const panzer::FieldPattern> fp
        = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer> >();

  // mesh description
  Ordinal64 nx = 10, ny = 7;
  int px = np, py = 1;
  int bx =  1, by = 2;

  RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
  connManager->initialize(comm,nx,ny,px,py,bx,by);

  // test element blocks are computed properly and sized appropriately
  {
    TEST_EQUALITY(Teuchos::as<int>(connManager->numElementBlocks()),bx*by);

    std::vector<std::string> eBlocks;
    connManager->getElementBlockIds(eBlocks);
    TEST_EQUALITY(eBlocks.size(),connManager->numElementBlocks());
    for(std::size_t i=1;i<eBlocks.size();i++) {
      out << "compare \"" << eBlocks[i-1] << "\" < \"" << eBlocks[i] << "\"" << std::endl;
      TEST_ASSERT(eBlocks[i-1]<eBlocks[i]);
    }
  }

  // test that owned and offset elements are correct
  { 
    auto myElements = connManager->getMyElementsTriplet();
    auto myOffset = connManager->getMyOffsetTriplet();

    Ordinal64 n = nx / px;
    Ordinal64 r = nx - n * px;

    TEST_EQUALITY(myElements.x,n + (r>rank ? 1 : 0));
    TEST_EQUALITY(myOffset.x,n*rank+std::min(Teuchos::as<int>(r),rank));

    TEST_EQUALITY(myElements.y,by*ny);
    TEST_EQUALITY(myOffset.y,0);

    TEST_EQUALITY(myElements.z,1);
    TEST_EQUALITY(myOffset.z,0);
  }

  // test that the elements are in the right blocks
  {
    Ordinal64 blk0[4],blk1[4];
    blk0[0] = connManager->computeLocalElementIndex(TripletGO(2,2,0));
    blk0[1] = connManager->computeLocalElementIndex(TripletGO(5,2,0));
    blk0[2] = connManager->computeLocalElementIndex(TripletGO(7,2,0));
    blk0[3] = connManager->computeLocalElementIndex(TripletGO(9,2,0));

    blk1[0] = connManager->computeLocalElementIndex(TripletGO(2,12,0));
    blk1[1] = connManager->computeLocalElementIndex(TripletGO(5,12,0));
    blk1[2] = connManager->computeLocalElementIndex(TripletGO(7,12,0));
    blk1[3] = connManager->computeLocalElementIndex(TripletGO(9,12,0));

    bool found = false;
    for(int i=0;i<4;i++) {
      if(blk0[i]!=-1) {
        TEST_EQUALITY("eblock-0_0",connManager->getBlockId(blk0[i]));
        found = true;
      }
    }
    TEST_ASSERT(found); // every processor must find at least one

    found = false;
    for(int i=0;i<4;i++) {
      if(blk1[i]!=-1) {
        TEST_EQUALITY("eblock-0_1",connManager->getBlockId(blk1[i]));
        found = true;
      }
    }
    TEST_ASSERT(found); // every processor must find at least one
  }

  {
    // check that all elements are in the right block
    const std::vector<int> & elmts0 = connManager->getElementBlock("eblock-0_0");
    for(std::size_t i=0;i<elmts0.size();i++) {
      TEST_EQUALITY(connManager->getBlockId(elmts0[i]),"eblock-0_0");
    }

    // check that all elements are accounted for
    int totalCount = 0;
    int count = Teuchos::as<int>(elmts0.size());
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&count,&totalCount);
    TEST_EQUALITY(totalCount,nx*ny);

    // check that all elements are in the right block
    const std::vector<int> & elmts1 = connManager->getElementBlock("eblock-0_1");
    for(std::size_t i=0;i<elmts1.size();i++) {
      TEST_EQUALITY(connManager->getBlockId(elmts1[i]),"eblock-0_1");
    }

    // check that all elements are accounted for
    totalCount = 0;
    count = Teuchos::as<int>(elmts1.size());
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&count,&totalCount);
    TEST_EQUALITY(totalCount,nx*ny);
  }
}

TEUCHOS_UNIT_TEST(tCartesianTop, connmanager_2d_1dpart)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  // int rank = comm.getRank();

  // mesh description
  Ordinal64 nx = 10, ny = 7;
  int px = np, py = 1;
  int bx =  1, by = 2;

  // test 2D nodal discretization
  {
    // field pattern for basis required
    RCP<const panzer::FieldPattern> fp = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

    // build the topology
    RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
    connManager->initialize(comm,nx,ny,px,py,bx,by);
    connManager->buildConnectivity(*fp);

    // test all the elements and functions
    std::string blocks[] = {"eblock-0_0","eblock-0_1"};
    for(int b=0;b<2;b++) {
      const std::vector<int> & elmts0 = connManager->getElementBlock(blocks[b]);
  
      for(std::size_t i=0;i<elmts0.size();i++) {
        TEST_EQUALITY(connManager->getConnectivitySize(elmts0[i]),4);
  
        auto global = CCM::computeLocalElementGlobalTriplet(elmts0[i],connManager->getMyElementsTriplet(),
                                                                      connManager->getMyOffsetTriplet());
        auto * conn = connManager->getConnectivity(elmts0[i]);
        TEST_ASSERT(conn!=0);
  
        TEST_EQUALITY(conn[0], global.x + (nx*bx+1)*global.y + 0 + (   0   ));
        TEST_EQUALITY(conn[1], global.x + (nx*bx+1)*global.y + 1 + (   0   ));
        TEST_EQUALITY(conn[2], global.x + (nx*bx+1)*global.y + 1 + (nx*bx+1));
        TEST_EQUALITY(conn[3], global.x + (nx*bx+1)*global.y + 0 + (nx*bx+1));
      }
    }
  }

  // test 2D Q2 discretization
  {
    panzer::Ordinal64 totalNodes = (bx*nx+1)*(by*ny+1);
    panzer::Ordinal64 totalEdges = (bx*nx+1)*(by*ny)+(bx*nx)*(by*ny+1);

    // field pattern for basis required
    RCP<const panzer::FieldPattern> fp = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer> >();

    // build the topology
    RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
    connManager->initialize(comm,nx,ny,px,py,bx,by);
    connManager->buildConnectivity(*fp);

    // test all the elements and functions
    std::string blocks[] = {"eblock-0_0","eblock-0_1"};
    for(int b=0;b<2;b++) {
      const std::vector<int> & elmts0 = connManager->getElementBlock(blocks[b]);
  
      for(std::size_t i=0;i<elmts0.size();i++) {
        TEST_EQUALITY(connManager->getConnectivitySize(elmts0[i]),9);
  
        auto global = CCM::computeLocalElementGlobalTriplet(elmts0[i],connManager->getMyElementsTriplet(),
                                                                      connManager->getMyOffsetTriplet());
        auto * conn = connManager->getConnectivity(elmts0[i]);
        TEST_ASSERT(conn!=0);
  
        // nodes
        TEST_EQUALITY(conn[0], global.x + (nx*bx+1)*global.y + 0 + (   0   ));
        TEST_EQUALITY(conn[1], global.x + (nx*bx+1)*global.y + 1 + (   0   ));
        TEST_EQUALITY(conn[2], global.x + (nx*bx+1)*global.y + 1 + (nx*bx+1));
        TEST_EQUALITY(conn[3], global.x + (nx*bx+1)*global.y + 0 + (nx*bx+1));

        // edges
        TEST_EQUALITY(conn[4], totalNodes + global.x + (2*nx*bx+1)*global.y + 0);
        TEST_EQUALITY(conn[5], totalNodes + global.x + (2*nx*bx+1)*global.y + 1 + nx*bx);
        TEST_EQUALITY(conn[6], totalNodes + global.x + (2*nx*bx+1)*global.y + (2*nx*bx+1));
        TEST_EQUALITY(conn[7], totalNodes + global.x + (2*nx*bx+1)*global.y + 0 + nx*bx);

        // cells
        TEST_EQUALITY(conn[8], totalNodes + totalEdges + global.x + nx*bx*global.y);
      }
    }
  }
}

// This test checks functions used to generate the topology are correct
TEUCHOS_UNIT_TEST(tCartesianTop, connmanager_3d_1dpart_helpers)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> TripletGO;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  int rank = comm.getRank();

  // field pattern for basis required
  RCP<const panzer::FieldPattern> fp
        = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer> >();

  // mesh description
  Ordinal64 nx = 10, ny = 7, nz = 4;
  int px = np, py = 1, pz = 1;
  int bx =  1, by = 2, bz = 1;

  RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
  connManager->initialize(comm,nx,ny,nz,px,py,pz,bx,by,bz);

  // test element blocks are computed properly and sized appropriately
  {
    TEST_EQUALITY(Teuchos::as<int>(connManager->numElementBlocks()),bx*by*bz);

    std::vector<std::string> eBlocks;
    connManager->getElementBlockIds(eBlocks);
    TEST_EQUALITY(eBlocks.size(),connManager->numElementBlocks());
    for(std::size_t i=1;i<eBlocks.size();i++) {
      out << "compare \"" << eBlocks[i-1] << "\" < \"" << eBlocks[i] << "\"" << std::endl;
      TEST_ASSERT(eBlocks[i-1]<eBlocks[i]);
    }
  }

  // test that owned and offset elements are correct
  { 
    auto myElements = connManager->getMyElementsTriplet();
    auto myOffset = connManager->getMyOffsetTriplet();

    Ordinal64 n = nx / px;
    Ordinal64 r = nx - n * px;

    TEST_EQUALITY(myElements.x,n + (r>rank ? 1 : 0));
    TEST_EQUALITY(myOffset.x,n*rank+std::min(Teuchos::as<int>(r),rank));

    TEST_EQUALITY(myElements.y,by*ny);
    TEST_EQUALITY(myOffset.y,0);

    TEST_EQUALITY(myElements.z,bz*nz);
    TEST_EQUALITY(myOffset.z,0);
  }

  // test that the elements are in the right blocks
  {
    Ordinal64 blk0[4],blk1[4];
    blk0[0] = connManager->computeLocalElementIndex(TripletGO(2,2,3));
    blk0[1] = connManager->computeLocalElementIndex(TripletGO(5,2,3));
    blk0[2] = connManager->computeLocalElementIndex(TripletGO(7,2,3));
    blk0[3] = connManager->computeLocalElementIndex(TripletGO(9,2,3));

    blk1[0] = connManager->computeLocalElementIndex(TripletGO(2,12,1));
    blk1[1] = connManager->computeLocalElementIndex(TripletGO(5,12,1));
    blk1[2] = connManager->computeLocalElementIndex(TripletGO(7,12,1));
    blk1[3] = connManager->computeLocalElementIndex(TripletGO(9,12,1));

    bool found = false;
    for(int i=0;i<4;i++) {
      if(blk0[i]!=-1) {
        TEST_EQUALITY("eblock-0_0_0",connManager->getBlockId(blk0[i]));
        found = true;
      }
    }
    TEST_ASSERT(found); // every processor must find at least one

    found = false;
    for(int i=0;i<4;i++) {
      if(blk1[i]!=-1) {
        TEST_EQUALITY("eblock-0_1_0",connManager->getBlockId(blk1[i]));
        found = true;
      }
    }
    TEST_ASSERT(found); // every processor must find at least one
  }

  {
    // check that all elements are in the right block
    const std::vector<int> & elmts0 = connManager->getElementBlock("eblock-0_0_0");
    for(std::size_t i=0;i<elmts0.size();i++) {
      TEST_EQUALITY(connManager->getBlockId(elmts0[i]),"eblock-0_0_0");
    }

    // check that all elements are accounted for
    int totalCount = 0;
    int count = Teuchos::as<int>(elmts0.size());
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&count,&totalCount);
    TEST_EQUALITY(totalCount,nx*ny*nz);

    // check that all elements are in the right block
    const std::vector<int> & elmts1 = connManager->getElementBlock("eblock-0_1_0");
    for(std::size_t i=0;i<elmts1.size();i++) {
      TEST_EQUALITY(connManager->getBlockId(elmts1[i]),"eblock-0_1_0");
    }

    // check that all elements are accounted for
    totalCount = 0;
    count = Teuchos::as<int>(elmts1.size());
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&count,&totalCount);
    TEST_EQUALITY(totalCount,nx*ny*nz);
  }
}

TEUCHOS_UNIT_TEST(tCartesianTop, connmanager_3d_1dpart)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  // int rank = comm.getRank();

  // mesh description
  Ordinal64 nx = 10, ny = 7, nz = 4;
  int px = np, py = 1, pz = 1;
  int bx =  1, by = 2, bz = 1;
/*
  Ordinal64 nx = 4, ny = 1, nz = 2;
  int px = np, py = 1, pz = 1;
  int bx =  1, by = 1, bz = 1;
*/

  // test 3D nodal discretization
  {
    // field pattern for basis required
    RCP<const panzer::FieldPattern> fp = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer> >();

    // build the topology
    RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
    connManager->initialize(comm,nx,ny,nz,px,py,pz,bx,by,bz);
    connManager->buildConnectivity(*fp);

    // test all the elements and functions
    std::string blocks[] = {"eblock-0_0_0","eblock-0_1_0"};
    for(int b=0;b<2;b++) {
      const std::vector<int> & elmts0 = connManager->getElementBlock(blocks[b]);
  
      for(std::size_t i=0;i<elmts0.size();i++) {
        TEST_EQUALITY(connManager->getConnectivitySize(elmts0[i]),8);
  
        auto global = CCM::computeLocalElementGlobalTriplet(elmts0[i],connManager->getMyElementsTriplet(),
                                                                      connManager->getMyOffsetTriplet());
        auto * conn = connManager->getConnectivity(elmts0[i]);
        TEST_ASSERT(conn!=0);
  
        auto basePoint = global.x + (nx*bx+1)*global.y + (nx*bx+1)*(ny*by+1)*global.z;
        TEST_EQUALITY(conn[0], basePoint + 0 + (   0   ) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[1], basePoint + 1 + (   0   ) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[2], basePoint + 1 + (nx*bx+1) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[3], basePoint + 0 + (nx*bx+1) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[4], basePoint + 0 + (   0   ) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[5], basePoint + 1 + (   0   ) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[6], basePoint + 1 + (nx*bx+1) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[7], basePoint + 0 + (nx*bx+1) + (nx*bx+1)*(ny*by+1));
      }
    }
  }

  // test 3D Q2 discretization
  {
    panzer::Ordinal64 totalNodes = (bx*nx+1)*(by*ny+1)*(bz*nz+1);
    panzer::Ordinal64 totalEdges = (bx*nx+1)*(by*ny)*(bz*nz+1)+(bx*nx)*(by*ny+1)*(bz*nz+1)+(bx*nx+1)*(by*ny+1)*(bz*nz);
    panzer::Ordinal64 totalFaces = (bx*nx+1)*(by*ny)*(bz*nz)+(bx*nx)*(by*ny+1)*(bz*nz)+(bx*nx)*(by*ny)*(bz*nz+1);

    // field pattern for basis required
    RCP<const panzer::FieldPattern> fp = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer> >();

    // build the topology
    RCP<CartesianConnManager<int,Ordinal64> > connManager = rcp(new CartesianConnManager<int,Ordinal64>);
    connManager->initialize(comm,nx,ny,nz,px,py,pz,bx,by,bz);
    connManager->buildConnectivity(*fp);

    // test all the elements and functions
    std::string blocks[] = {"eblock-0_0_0","eblock-0_1_0"};
    for(int b=0;b<2;b++) {
      const std::vector<int> & elmts0 = connManager->getElementBlock(blocks[b]);
  
      for(std::size_t i=0;i<elmts0.size();i++) {
        TEST_EQUALITY(connManager->getConnectivitySize(elmts0[i]),27);
  
        auto global = CCM::computeLocalElementGlobalTriplet(elmts0[i],connManager->getMyElementsTriplet(),
                                                                      connManager->getMyOffsetTriplet());
        out << "Element Triplet: " << global.x << ", " << global.y << ", " << global.z << std::endl;
        auto * conn = connManager->getConnectivity(elmts0[i]);
        TEST_ASSERT(conn!=0);
  
        // nodes
        auto nodeBasePoint = global.x + (nx*bx+1)*global.y + (nx*bx+1)*(ny*by+1)*global.z;
        TEST_EQUALITY(conn[0], nodeBasePoint + 0 + (   0   ) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[1], nodeBasePoint + 1 + (   0   ) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[2], nodeBasePoint + 1 + (nx*bx+1) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[3], nodeBasePoint + 0 + (nx*bx+1) + (   0   )*(   0   ));
        TEST_EQUALITY(conn[4], nodeBasePoint + 0 + (   0   ) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[5], nodeBasePoint + 1 + (   0   ) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[6], nodeBasePoint + 1 + (nx*bx+1) + (nx*bx+1)*(ny*by+1));
        TEST_EQUALITY(conn[7], nodeBasePoint + 0 + (nx*bx+1) + (nx*bx+1)*(ny*by+1));

        // edges
        auto e_ks = (bx*nx+1)*by*ny + bx*nx*(by*ny+1) + (nx*bx+1)*(ny*by+1);
        auto e_kp = (bx*nx+1)*by*ny + bx*nx*(by*ny+1);
        auto edgeBasePoint = totalNodes + global.x + global.y *(2*bx*nx+1) + global.z*e_ks;

        // horizontal edges: bottom 
        TEST_EQUALITY(conn[ 8], edgeBasePoint +           0);
        TEST_EQUALITY(conn[ 9], edgeBasePoint +   bx*nx + 1);
        TEST_EQUALITY(conn[10], edgeBasePoint + 2*bx*nx + 1);
        TEST_EQUALITY(conn[11], edgeBasePoint +   bx*nx + 0);

        // horizontal edges: top
        TEST_EQUALITY(conn[12], edgeBasePoint + e_ks +           0);
        TEST_EQUALITY(conn[13], edgeBasePoint + e_ks +   bx*nx + 1);
        TEST_EQUALITY(conn[14], edgeBasePoint + e_ks + 2*bx*nx + 1);
        TEST_EQUALITY(conn[15], edgeBasePoint + e_ks +   bx*nx + 0);

        // vertical edges
        TEST_EQUALITY(conn[16], edgeBasePoint + e_kp - global.y*bx*nx);
        TEST_EQUALITY(conn[17], edgeBasePoint + e_kp - global.y*bx*nx     + 1);
        TEST_EQUALITY(conn[18], edgeBasePoint + e_kp - (global.y-1)*bx*nx + 2);
        TEST_EQUALITY(conn[19], edgeBasePoint + e_kp - (global.y-1)*bx*nx + 1);

        // cells
        TEST_EQUALITY(conn[26], totalNodes + totalEdges + totalFaces + global.x + nx*bx*global.y + nx*bx*ny*by*global.z);

        // faces
        auto f_ks = nx*bx*ny*by + (nx*bx+1)*ny*by + nx*bx*(ny*by+1);
        auto f_kp = nx*bx*ny*by;
        auto faceBasePoint = totalNodes + totalEdges + global.x + global.y*nx*bx + global.z*f_ks;
        TEST_EQUALITY(conn[20], faceBasePoint + f_kp + (global.y+1)*(nx*bx+1) - nx*bx -1);
        TEST_EQUALITY(conn[21], faceBasePoint + f_kp + (global.y+1)*(nx*bx+1));
        TEST_EQUALITY(conn[22], faceBasePoint + f_kp + (global.y+1)*(nx*bx+1) + nx*bx);
        TEST_EQUALITY(conn[23], faceBasePoint + f_kp + (global.y+1)*(nx*bx+1) - 1);
        TEST_EQUALITY(conn[24], faceBasePoint);
        TEST_EQUALITY(conn[25], faceBasePoint + f_ks);
      }
    }
  }
}

} // end unit test
} // end panzer

