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

#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_FieldContainer.hpp"

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

typedef Intrepid2::FieldContainer<double> FieldContainer;
template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
  // build a geometric pattern from a single basis
  RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2Type);
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}

// test that you can correctly compute a rank index from a grobal processor id
TEUCHOS_UNIT_TEST(tCartesianTop, computeMyRankIndex)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<int> Triplet;

  // test 2D
  {
    auto t0 = CCM::computeMyRankIndex(3,2,Triplet(4,3,1));
    auto t1 = CCM::computeMyRankIndex(6,2,Triplet(4,3,1));
    auto t2 = CCM::computeMyRankIndex(9,2,Triplet(4,3,1));

    TEST_EQUALITY(t0.x,3); TEST_EQUALITY(t0.y,0); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,2); TEST_EQUALITY(t1.y,1); TEST_EQUALITY(t1.z,0);
    TEST_EQUALITY(t2.x,1); TEST_EQUALITY(t2.y,2); TEST_EQUALITY(t2.z,0);
  }

  // test 3D
  {
    auto t0 = CCM::computeMyRankIndex( 3,3,Triplet(2,3,3));
    auto t1 = CCM::computeMyRankIndex(10,3,Triplet(2,3,3));
    auto t2 = CCM::computeMyRankIndex(13,3,Triplet(2,3,3));

    TEST_EQUALITY(t0.x,1); TEST_EQUALITY(t0.y,1); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,0); TEST_EQUALITY(t1.y,2); TEST_EQUALITY(t1.z,1);
    TEST_EQUALITY(t2.x,1); TEST_EQUALITY(t2.y,0); TEST_EQUALITY(t2.z,2);
  }
}

// test that you can correctly compute a rank index from a grobal processor id
TEUCHOS_UNIT_TEST(tCartesianTop, computeGlobalElementTriplet)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef CCM::Triplet<Ordinal64> Triplet;

  // test 2D
  {
    auto t0 = CCM::computeGlobalElementTriplet( 1,Triplet(4,3,1),Triplet(7,2,0));
    auto t1 = CCM::computeGlobalElementTriplet( 4,Triplet(4,3,1),Triplet(7,2,0));
    auto t2 = CCM::computeGlobalElementTriplet(11,Triplet(4,3,1),Triplet(7,2,0));

    TEST_EQUALITY(t0.x,1+7); TEST_EQUALITY(t0.y,0+2); TEST_EQUALITY(t0.z,0);
    TEST_EQUALITY(t1.x,0+7); TEST_EQUALITY(t1.y,1+2); TEST_EQUALITY(t1.z,0);
    TEST_EQUALITY(t2.x,3+7); TEST_EQUALITY(t2.y,2+2); TEST_EQUALITY(t2.z,0);
  }

  // test 3D
  {
    auto t0 = CCM::computeGlobalElementTriplet( 1+2*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t1 = CCM::computeGlobalElementTriplet( 4+1*12,Triplet(4,3,5),Triplet(7,2,3));
    auto t2 = CCM::computeGlobalElementTriplet(11+4*12,Triplet(4,3,5),Triplet(7,2,3));

    TEST_EQUALITY(t0.x,1+7); TEST_EQUALITY(t0.y,0+2); TEST_EQUALITY(t0.z,2+3);
    TEST_EQUALITY(t1.x,0+7); TEST_EQUALITY(t1.y,1+2); TEST_EQUALITY(t1.z,1+3);
    TEST_EQUALITY(t2.x,3+7); TEST_EQUALITY(t2.y,2+2); TEST_EQUALITY(t2.z,4+3);
  }
}

// This test checks the the topology generated is what is expected
TEUCHOS_UNIT_TEST(tCartesianTop, connmanager_2d_1dpart)
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
  }
}

} // end unit test
} // end panzer

