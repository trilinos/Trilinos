// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>

#ifdef STK_HAVE_INTREPID2

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Intrepid2_CellTools.hpp>
#include <array>

using intrepid2ct = Intrepid2::CellTools<Kokkos::DefaultHostExecutionSpace>;

TEST(SearchUtil, Intrepid2_hasReferenceCell)
{
  //volume-elements; i.e., not beams/shells, not 2D tris/quads etc.
  std::array<stk::topology,14> stkVolumeElemTopos =
  {
    stk::topology::TET_4,
    stk::topology::TET_8,
    stk::topology::TET_10,
    stk::topology::TET_11,
    stk::topology::PYRAMID_5,
    stk::topology::PYRAMID_13,
    stk::topology::PYRAMID_14,
    stk::topology::WEDGE_6,
    stk::topology::WEDGE_12,
    stk::topology::WEDGE_15,
    stk::topology::WEDGE_18,
    stk::topology::HEX_8,
    stk::topology::HEX_20,
    stk::topology::HEX_27
  };

  for(stk::topology stkElem : stkVolumeElemTopos) {
    shards::CellTopology shardsCellTopo = stk::mesh::get_cell_topology(stkElem);
    if (shardsCellTopo.isValid()) {
      const bool hasReferenceCell = intrepid2ct::hasReferenceCell(shardsCellTopo);
      std::cout<<"stk '"<<stkElem<<"' -> shards '"<<shardsCellTopo.getName()
               <<"' (base "<<shardsCellTopo.getBaseName()
               <<"), hasReferenceCell="<<std::boolalpha<<hasReferenceCell<<std::endl;
    }
    else {
      std::cout<<" !!! stk-topo '"<<stkElem<<"' doesn't map to a valid shards CellTopology. !!!"<<std::endl;
    }
  }
}

TEST(SearchUtil, Intrepid2_checkPointwiseInclusion)
{
  using Intrepid2FloatView = Kokkos::DynRankView<float,Kokkos::DefaultHostExecutionSpace>;
  using Intrepid2IntView = Kokkos::DynRankView<int,Kokkos::DefaultHostExecutionSpace>;

  constexpr unsigned spatialDim=3, numNodes=4, numCells=1, numPoints=2;
  Intrepid2FloatView tet4coords("Tet4Cell", numCells, numNodes, spatialDim);
  Intrepid2FloatView points("Point", numCells, numPoints, spatialDim);

  tet4coords(0, 0, 0) = 0.0;
  tet4coords(0, 0, 1) = 0.0;
  tet4coords(0, 0, 2) = 0.0;
  tet4coords(0, 1, 0) = 1.0;
  tet4coords(0, 1, 1) = 0.0;
  tet4coords(0, 1, 2) = 0.0;
  tet4coords(0, 2, 0) = 0.0;
  tet4coords(0, 2, 1) = 1.0;
  tet4coords(0, 2, 2) = 0.0;
  tet4coords(0, 3, 0) = 0.0;
  tet4coords(0, 3, 1) = 0.0;
  tet4coords(0, 3, 2) = 1.0;

  points(0, 0, 0) = 0.2; //this point is inside the tet4
  points(0, 0, 1) = 0.2;
  points(0, 0, 2) = 0.2;

  points(0, 1, 0) = 1.2; //this point is NOT inside the tet4
  points(0, 1, 1) = 1.2;
  points(0, 1, 2) = 1.2;

  shards::CellTopology tet4 = stk::mesh::get_cell_topology(stk::topology::TET_4);
  ASSERT_TRUE(tet4.isValid());

  Intrepid2IntView inCell("inCell", numCells, numPoints);

  intrepid2ct::checkPointwiseInclusion(inCell, points, tet4coords, tet4);

  EXPECT_EQ(1, inCell(0,0)) << "got wrong answer: point 0 should be in Tet4!";
  EXPECT_EQ(0, inCell(0,1)) << "got wrong answer: point 1 should not be in Tet4!";
}

TEST(SearchUtil, Intrepid2_mapToReferenceFrame_checkPointwiseInclusion)
{
  using Intrepid2FloatView = Kokkos::DynRankView<float,Kokkos::DefaultHostExecutionSpace>;
  using Intrepid2IntView = Kokkos::DynRankView<int,Kokkos::DefaultHostExecutionSpace>;

  constexpr unsigned spatialDim=3, numNodes=4, numCells=1, numPoints=2;
  Intrepid2FloatView tet4coords("Tet4Cell", numCells, numNodes, spatialDim);
  Intrepid2FloatView points("Point", numCells, numPoints, spatialDim);
  Intrepid2FloatView refPoints("Point", numCells, numPoints, spatialDim);

  tet4coords(0, 0, 0) = 0.0;
  tet4coords(0, 0, 1) = 0.0;
  tet4coords(0, 0, 2) = 0.0;
  tet4coords(0, 1, 0) = 1.0;
  tet4coords(0, 1, 1) = 0.0;
  tet4coords(0, 1, 2) = 0.0;
  tet4coords(0, 2, 0) = 0.0;
  tet4coords(0, 2, 1) = 1.0;
  tet4coords(0, 2, 2) = 0.0;
  tet4coords(0, 3, 0) = 0.0;
  tet4coords(0, 3, 1) = 0.0;
  tet4coords(0, 3, 2) = 1.0;

  points(0, 0, 0) = 0.2; //this point is inside the tet4
  points(0, 0, 1) = 0.2;
  points(0, 0, 2) = 0.2;

  points(0, 1, 0) = 1.2; //this point is NOT inside the tet4
  points(0, 1, 1) = 1.2;
  points(0, 1, 2) = 1.2;

  shards::CellTopology tet4 = stk::mesh::get_cell_topology(stk::topology::TET_4);
  ASSERT_TRUE(tet4.isValid());

  Intrepid2IntView inCell("inCell", numCells, numPoints);

  intrepid2ct::mapToReferenceFrame(refPoints, points, tet4coords, tet4);
  intrepid2ct::checkPointwiseInclusion(inCell, refPoints, tet4);

  EXPECT_EQ(1, inCell(0,0)) << "got wrong answer: point 0 should be in Tet4!";
  EXPECT_EQ(0, inCell(0,1)) << "got wrong answer: point 1 should not be in Tet4!";
}

#endif // STK_HAVE_INTREPID2

