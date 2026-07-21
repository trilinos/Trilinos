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
#include <stk_search_util/Intrepid2BasisFactory.hpp>
#include <stk_search_util/Intrepid2_HasParametricDistance.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Intrepid2_CellTools.hpp>
#include <array>

using I2ExecSpace = stk::search::BasisExecSpace;
using Intrepid2ct = Intrepid2::CellTools<typename I2ExecSpace::device_type>;
using Intrepid2rst = Intrepid2::RealSpaceTools<typename I2ExecSpace::device_type>;
using Intrepid2DoubleView = Kokkos::DynRankView<double,typename I2ExecSpace::memory_space>;
using Intrepid2IntView = Kokkos::DynRankView<int,typename I2ExecSpace::memory_space>;

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
      const bool hasReferenceCell = Intrepid2ct::hasReferenceCell(shardsCellTopo);
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
  constexpr unsigned spatialDim=3, numNodes=4, numCells=1, numPoints=2;
  Intrepid2DoubleView tet4coords("Tet4Cell", numCells, numNodes, spatialDim);
  Intrepid2DoubleView points("Point", numCells, numPoints, spatialDim);

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

  Intrepid2ct::checkPointwiseInclusion(inCell, points, tet4coords, tet4);

  EXPECT_EQ(1, inCell(0,0)) << "got wrong answer: point 0 should be in Tet4!";
  EXPECT_EQ(0, inCell(0,1)) << "got wrong answer: point 1 should not be in Tet4!";
}

TEST(SearchUtil, Intrepid2_mapToReferenceFrame_checkPointwiseInclusion_tet4)
{
  constexpr unsigned spatialDim=3, numNodes=4, numCells=1, numPoints=2;
  Intrepid2DoubleView tet4coords("Tet4Cell", numCells, numNodes, spatialDim);
  Intrepid2DoubleView points("Point", numCells, numPoints, spatialDim);
  Intrepid2DoubleView refPoints("Point", numCells, numPoints, spatialDim);

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

  points(0, 0, 0) = 0.25; //this point is inside the tet4
  points(0, 0, 1) = 0.25;
  points(0, 0, 2) = 0.25;

  points(0, 1, 0) = 1.25; //this point is NOT inside the tet4
  points(0, 1, 1) = 1.25;
  points(0, 1, 2) = 1.25;

  shards::CellTopology tet4 = stk::mesh::get_cell_topology(stk::topology::TET_4);
  ASSERT_TRUE(tet4.isValid());

  Intrepid2IntView inCell("inCell", numCells, numPoints);

  Intrepid2ct::mapToReferenceFrame(refPoints, points, tet4coords, tet4);
  Intrepid2ct::checkPointwiseInclusion(inCell, refPoints, tet4);

  EXPECT_EQ(1, inCell(0,0)) << "got wrong answer: point 0 should be in Tet4!";
  EXPECT_EQ(0, inCell(0,1)) << "got wrong answer: point 1 should not be in Tet4!";
}

template<typename ViewType>
void fill_test_tet10(ViewType& tet10coords)
{
  tet10coords(0, 0, 0) = 0.0; //3 base-triangle vertices
  tet10coords(0, 0, 1) = 0.0;
  tet10coords(0, 0, 2) = 0.0;
  tet10coords(0, 1, 0) = 1.0;
  tet10coords(0, 1, 1) = 0.0;
  tet10coords(0, 1, 2) = 0.0;
  tet10coords(0, 2, 0) = 0.0;
  tet10coords(0, 2, 1) = 1.0;
  tet10coords(0, 2, 2) = 0.0;

  tet10coords(0, 3, 0) = 0.0; // node at the top
  tet10coords(0, 3, 1) = 0.0;
  tet10coords(0, 3, 2) = 1.0;

  tet10coords(0, 4, 0) = 0.5; //3 base-triangle mid-edge nodes
  tet10coords(0, 4, 1) = 0.0;
  tet10coords(0, 4, 2) = 0.0;
  tet10coords(0, 5, 0) = 0.5;
  tet10coords(0, 5, 1) = 0.5;
  tet10coords(0, 5, 2) = 0.0;
  tet10coords(0, 6, 0) = 0.0;
  tet10coords(0, 6, 1) = 0.5;
  tet10coords(0, 6, 2) = 0.0;

  tet10coords(0, 7, 0) = 0.0; // mid-edge nodes of 'upright' edges
  tet10coords(0, 7, 1) = 0.0;
  tet10coords(0, 7, 2) = 0.5;
  tet10coords(0, 8, 0) = 0.5;
  tet10coords(0, 8, 1) = 0.0;
  tet10coords(0, 8, 2) = 0.5;
  tet10coords(0, 9, 0) = 0.0;
  tet10coords(0, 9, 1) = 0.5;
  tet10coords(0, 9, 2) = 0.5;
}

template<bool Intrepid2HasNewFunction, typename ViewType, typename BasisType>
void compute_parametric_coords(unsigned numParamCoords, unsigned numPoints,
                               ViewType& paramCoords,
                               ViewType& refCellCenter, ViewType& physCoordsView,
                               ViewType& elemNodesView, BasisType& basis)
{ 
  constexpr unsigned numCells=1;
  
  if constexpr (Intrepid2HasNewFunction) {
    Intrepid2rst::clone(paramCoords, refCellCenter);
    Intrepid2ct::mapToReferenceFrame(paramCoords, physCoordsView, elemNodesView, basis);
  }
  else {
    ViewType initGuess("init_guess", numCells, numPoints, numParamCoords);
    Intrepid2rst::clone(initGuess, refCellCenter);
    Intrepid2ct::mapToReferenceFrameInitGuess(paramCoords, initGuess, physCoordsView, elemNodesView, basis);
  }
}

void test_mapToReferenceFrameInitGuess_checkPointwiseInclusion_tet10(bool useComposite)
{
  constexpr unsigned spatialDim=3, numNodes=10, numCells=1, numPoints=2;
  Intrepid2DoubleView tet10coords("Tet4Cell", numCells, numNodes, spatialDim);
  Intrepid2DoubleView points("Point", numCells, numPoints, spatialDim);
  Intrepid2DoubleView refPoints("Point", numCells, numPoints, spatialDim);
  Intrepid2DoubleView refCellCenter("refCenter", spatialDim);

  stk::topology stkTopo = stk::topology::TET_10;
  shards::CellTopology tet10 = stk::mesh::get_cell_topology(stkTopo);
  Intrepid2ct::getReferenceCellCenter(refCellCenter, tet10);
  Intrepid2rst::clone(refPoints, refCellCenter); 

  fill_test_tet10(tet10coords);

  points(0, 0, 0) = 0.25; //this point is inside the tet10
  points(0, 0, 1) = 0.25;
  points(0, 0, 2) = 0.25;

  points(0, 1, 0) = 1.25; //this point is NOT inside the tet10
  points(0, 1, 1) = 1.25;
  points(0, 1, 2) = 1.25;

  ASSERT_TRUE(tet10.isValid());

  Teuchos::RCP<Intrepid2::Basis<typename I2ExecSpace::device_type,double,double>> basis = Teuchos::rcp(stk::search::lagrangeBasisFactory(stkTopo, useComposite));
  compute_parametric_coords<stk::search::intrepid2_has_param_dist>(spatialDim, numPoints,
        refPoints, refCellCenter, points, tet10coords, basis);

  Intrepid2IntView inCell("inCell", numCells, numPoints);
  Intrepid2ct::checkPointwiseInclusion(inCell, refPoints, tet10);

  EXPECT_EQ(1, inCell(0,0)) << "got wrong answer: point 0 should be in Tet10!";
  EXPECT_EQ(0, inCell(0,1)) << "got wrong answer: point 1 should not be in Tet10!";
}

TEST(SearchUtil, Intrepid2_mapToReferenceFrameInitGuess_checkPointwiseInclusion_tet10_composite)
{
  const bool useComposite = true;
  test_mapToReferenceFrameInitGuess_checkPointwiseInclusion_tet10(useComposite);
}

TEST(SearchUtil, Intrepid2_mapToReferenceFrameInitGuess_checkPointwiseInclusion_tet10_quadratic)
{
  const bool useComposite = false;
  test_mapToReferenceFrameInitGuess_checkPointwiseInclusion_tet10(useComposite);
}

#endif // STK_HAVE_INTREPID2

