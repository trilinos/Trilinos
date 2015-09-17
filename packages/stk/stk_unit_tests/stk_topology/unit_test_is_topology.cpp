// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_topology/topology.hpp>

TEST( stk_topology, isTri)
{
  using stk::topology;

  topology implementedTopologies[] =   // std::vector<topology> implementedTopologies(72)
 { stk::topology::INVALID_TOPOLOGY,
   stk::topology::BEGIN_TOPOLOGY, stk::topology::NODE,
   stk::topology::LINE_2, stk::topology::LINE_3,
   stk::topology::TRI_3, stk::topology::TRIANGLE_3, stk::topology::TRI_4, stk::topology::TRIANGLE_4, stk::topology::TRI_6, stk::topology::TRIANGLE_6,
   stk::topology::QUAD_4, stk::topology::QUADRILATERAL_4, stk::topology::QUAD_8, stk::topology::QUADRILATERAL_8, stk::topology::QUAD_9, stk::topology::QUADRILATERAL_9,
   stk::topology::PARTICLE,
   stk::topology::LINE_2_1D, stk::topology::LINE_3_1D,
   stk::topology::BEAM_2, stk::topology::BEAM_3,
   stk::topology::SHELL_LINE_2, stk::topology::SHELL_LINE_3,
   stk::topology::TRI_3_2D, stk::topology::TRIANGLE_3_2D, stk::topology::TRI_4_2D, stk::topology::TRIANGLE_4_2D, stk::topology::TRI_6_2D, stk::topology::TRIANGLE_6_2D,
   stk::topology::QUAD_4_2D, stk::topology::QUADRILATERAL_4_2D, stk::topology::QUAD_8_2D, stk::topology::QUADRILATERAL_8_2D, stk::topology::QUAD_9_2D, stk::topology::QUADRILATERAL_9_2D,
   stk::topology::SHELL_TRI_3, stk::topology::SHELL_TRIANGLE_3, stk::topology::SHELL_TRI_4, stk::topology::SHELL_TRIANGLE_4, stk::topology::SHELL_TRI_6, stk::topology::SHELL_TRIANGLE_6,
   stk::topology::SHELL_QUAD_4, stk::topology::SHELL_QUADRILATERAL_4, stk::topology::SHELL_QUAD_8, stk::topology::SHELL_QUADRILATERAL_8, stk::topology::SHELL_QUAD_9, stk::topology::SHELL_QUADRILATERAL_9,
   stk::topology::TET_4, stk::topology::TETRAHEDRON_4, stk::topology::TET_8, stk::topology:: TETRAHEDRON_8, stk::topology::TET_10, stk::topology::TETRAHEDRON_10, stk::topology::TET_11, stk::topology::TETRAHEDRON_11,
   stk::topology::PYRAMID_5, stk::topology::PYRAMID_13, stk::topology::PYRAMID_14,
   stk::topology::WEDGE_6, stk::topology::WEDGE_15, stk::topology::WEDGE_18,
   stk::topology::HEX_8, stk::topology:: HEXAHEDRON_8, stk::topology::HEX_20, stk::topology::HEXAHEDRON_20, stk::topology::HEX_27, stk::topology::HEXAHEDRON_27,
   stk::topology::END_TOPOLOGY, stk::topology::NUM_TOPOLOGIES,
   stk::topology::SUPERELEMENT_START,
   stk::topology::FORCE_TOPOLOGY_TO_UNSIGNED };

  for (uint i = 0; i<72;++i){
      bool amITri = stk::isTriangle(implementedTopologies[i]);
      if (i >23 && i<30 ){
        EXPECT_TRUE(amITri);
      }
      else{
        EXPECT_FALSE(amITri);
      }
  }

}


TEST( stk_topology, isQuad)
{
  using stk::topology;

  topology implementedTopologies[] =   // std::vector<topology> implementedTopologies(72)
 { stk::topology::INVALID_TOPOLOGY,
   stk::topology::BEGIN_TOPOLOGY, stk::topology::NODE,
   stk::topology::LINE_2, stk::topology::LINE_3,
   stk::topology::TRI_3, stk::topology::TRIANGLE_3, stk::topology::TRI_4, stk::topology::TRIANGLE_4, stk::topology::TRI_6, stk::topology::TRIANGLE_6,
   stk::topology::QUAD_4, stk::topology::QUADRILATERAL_4, stk::topology::QUAD_8, stk::topology::QUADRILATERAL_8, stk::topology::QUAD_9, stk::topology::QUADRILATERAL_9,
   stk::topology::PARTICLE,
   stk::topology::LINE_2_1D, stk::topology::LINE_3_1D,
   stk::topology::BEAM_2, stk::topology::BEAM_3,
   stk::topology::SHELL_LINE_2, stk::topology::SHELL_LINE_3,
   stk::topology::TRI_3_2D, stk::topology::TRIANGLE_3_2D, stk::topology::TRI_4_2D, stk::topology::TRIANGLE_4_2D, stk::topology::TRI_6_2D, stk::topology::TRIANGLE_6_2D,
   stk::topology::QUAD_4_2D, stk::topology::QUADRILATERAL_4_2D, stk::topology::QUAD_8_2D, stk::topology::QUADRILATERAL_8_2D, stk::topology::QUAD_9_2D, stk::topology::QUADRILATERAL_9_2D,
   stk::topology::SHELL_TRI_3, stk::topology::SHELL_TRIANGLE_3, stk::topology::SHELL_TRI_4, stk::topology::SHELL_TRIANGLE_4, stk::topology::SHELL_TRI_6, stk::topology::SHELL_TRIANGLE_6,
   stk::topology::SHELL_QUAD_4, stk::topology::SHELL_QUADRILATERAL_4, stk::topology::SHELL_QUAD_8, stk::topology::SHELL_QUADRILATERAL_8, stk::topology::SHELL_QUAD_9, stk::topology::SHELL_QUADRILATERAL_9,
   stk::topology::TET_4, stk::topology::TETRAHEDRON_4, stk::topology::TET_8, stk::topology:: TETRAHEDRON_8, stk::topology::TET_10, stk::topology::TETRAHEDRON_10, stk::topology::TET_11, stk::topology::TETRAHEDRON_11,
   stk::topology::PYRAMID_5, stk::topology::PYRAMID_13, stk::topology::PYRAMID_14,
   stk::topology::WEDGE_6, stk::topology::WEDGE_15, stk::topology::WEDGE_18,
   stk::topology::HEX_8, stk::topology:: HEXAHEDRON_8, stk::topology::HEX_20, stk::topology::HEXAHEDRON_20, stk::topology::HEX_27, stk::topology::HEXAHEDRON_27,
   stk::topology::END_TOPOLOGY, stk::topology::NUM_TOPOLOGIES,
   stk::topology::SUPERELEMENT_START,
   stk::topology::FORCE_TOPOLOGY_TO_UNSIGNED };

  for (uint i = 0; i<72;++i){
      bool amIQuad = stk::isQuadrilateral(implementedTopologies[i]);
      if (i >29 && i<36 ){
        EXPECT_TRUE(amIQuad);
      }
      else{
        EXPECT_FALSE(amIQuad);
      }
  }

}

TEST( stk_topology, isHex)
{
  using stk::topology;

  topology implementedTopologies[] =   // std::vector<topology> implementedTopologies(72)
 { stk::topology::INVALID_TOPOLOGY,
   stk::topology::BEGIN_TOPOLOGY, stk::topology::NODE,
   stk::topology::LINE_2, stk::topology::LINE_3,
   stk::topology::TRI_3, stk::topology::TRIANGLE_3, stk::topology::TRI_4, stk::topology::TRIANGLE_4, stk::topology::TRI_6, stk::topology::TRIANGLE_6,
   stk::topology::QUAD_4, stk::topology::QUADRILATERAL_4, stk::topology::QUAD_8, stk::topology::QUADRILATERAL_8, stk::topology::QUAD_9, stk::topology::QUADRILATERAL_9,
   stk::topology::PARTICLE,
   stk::topology::LINE_2_1D, stk::topology::LINE_3_1D,
   stk::topology::BEAM_2, stk::topology::BEAM_3,
   stk::topology::SHELL_LINE_2, stk::topology::SHELL_LINE_3,
   stk::topology::TRI_3_2D, stk::topology::TRIANGLE_3_2D, stk::topology::TRI_4_2D, stk::topology::TRIANGLE_4_2D, stk::topology::TRI_6_2D, stk::topology::TRIANGLE_6_2D,
   stk::topology::QUAD_4_2D, stk::topology::QUADRILATERAL_4_2D, stk::topology::QUAD_8_2D, stk::topology::QUADRILATERAL_8_2D, stk::topology::QUAD_9_2D, stk::topology::QUADRILATERAL_9_2D,
   stk::topology::SHELL_TRI_3, stk::topology::SHELL_TRIANGLE_3, stk::topology::SHELL_TRI_4, stk::topology::SHELL_TRIANGLE_4, stk::topology::SHELL_TRI_6, stk::topology::SHELL_TRIANGLE_6,
   stk::topology::SHELL_QUAD_4, stk::topology::SHELL_QUADRILATERAL_4, stk::topology::SHELL_QUAD_8, stk::topology::SHELL_QUADRILATERAL_8, stk::topology::SHELL_QUAD_9, stk::topology::SHELL_QUADRILATERAL_9,
   stk::topology::TET_4, stk::topology::TETRAHEDRON_4, stk::topology::TET_8, stk::topology:: TETRAHEDRON_8, stk::topology::TET_10, stk::topology::TETRAHEDRON_10, stk::topology::TET_11, stk::topology::TETRAHEDRON_11,
   stk::topology::PYRAMID_5, stk::topology::PYRAMID_13, stk::topology::PYRAMID_14,
   stk::topology::WEDGE_6, stk::topology::WEDGE_15, stk::topology::WEDGE_18,
   stk::topology::HEX_8, stk::topology:: HEXAHEDRON_8, stk::topology::HEX_20, stk::topology::HEXAHEDRON_20, stk::topology::HEX_27, stk::topology::HEXAHEDRON_27,
   stk::topology::END_TOPOLOGY, stk::topology::NUM_TOPOLOGIES,
   stk::topology::SUPERELEMENT_START,
   stk::topology::FORCE_TOPOLOGY_TO_UNSIGNED };

  for (uint i = 0; i<72;++i){
      bool amIHex = stk::isHexahedron(implementedTopologies[i]);
      if ((i >61 && i<68) || i == 69 ){                                                //not sure why one of the upper topologies is equivalent to HEXAHEDRON_27
        EXPECT_TRUE(amIHex);
      }
      else{
        EXPECT_FALSE(amIHex);
      }
  }

}

TEST( stk_topology, isTet)
{
  using stk::topology;

  topology implementedTopologies[] =   // std::vector<topology> implementedTopologies(72)
 { stk::topology::INVALID_TOPOLOGY,
   stk::topology::BEGIN_TOPOLOGY, stk::topology::NODE,
   stk::topology::LINE_2, stk::topology::LINE_3,
   stk::topology::TRI_3, stk::topology::TRIANGLE_3, stk::topology::TRI_4, stk::topology::TRIANGLE_4, stk::topology::TRI_6, stk::topology::TRIANGLE_6,
   stk::topology::QUAD_4, stk::topology::QUADRILATERAL_4, stk::topology::QUAD_8, stk::topology::QUADRILATERAL_8, stk::topology::QUAD_9, stk::topology::QUADRILATERAL_9,
   stk::topology::PARTICLE,
   stk::topology::LINE_2_1D, stk::topology::LINE_3_1D,
   stk::topology::BEAM_2, stk::topology::BEAM_3,
   stk::topology::SHELL_LINE_2, stk::topology::SHELL_LINE_3,
   stk::topology::TRI_3_2D, stk::topology::TRIANGLE_3_2D, stk::topology::TRI_4_2D, stk::topology::TRIANGLE_4_2D, stk::topology::TRI_6_2D, stk::topology::TRIANGLE_6_2D,
   stk::topology::QUAD_4_2D, stk::topology::QUADRILATERAL_4_2D, stk::topology::QUAD_8_2D, stk::topology::QUADRILATERAL_8_2D, stk::topology::QUAD_9_2D, stk::topology::QUADRILATERAL_9_2D,
   stk::topology::SHELL_TRI_3, stk::topology::SHELL_TRIANGLE_3, stk::topology::SHELL_TRI_4, stk::topology::SHELL_TRIANGLE_4, stk::topology::SHELL_TRI_6, stk::topology::SHELL_TRIANGLE_6,
   stk::topology::SHELL_QUAD_4, stk::topology::SHELL_QUADRILATERAL_4, stk::topology::SHELL_QUAD_8, stk::topology::SHELL_QUADRILATERAL_8, stk::topology::SHELL_QUAD_9, stk::topology::SHELL_QUADRILATERAL_9,
   stk::topology::TET_4, stk::topology::TETRAHEDRON_4, stk::topology::TET_8, stk::topology:: TETRAHEDRON_8, stk::topology::TET_10, stk::topology::TETRAHEDRON_10, stk::topology::TET_11, stk::topology::TETRAHEDRON_11,
   stk::topology::PYRAMID_5, stk::topology::PYRAMID_13, stk::topology::PYRAMID_14,
   stk::topology::WEDGE_6, stk::topology::WEDGE_15, stk::topology::WEDGE_18,
   stk::topology::HEX_8, stk::topology:: HEXAHEDRON_8, stk::topology::HEX_20, stk::topology::HEXAHEDRON_20, stk::topology::HEX_27, stk::topology::HEXAHEDRON_27,
   stk::topology::END_TOPOLOGY, stk::topology::NUM_TOPOLOGIES,
   stk::topology::SUPERELEMENT_START,
   stk::topology::FORCE_TOPOLOGY_TO_UNSIGNED };

  for (uint i = 0; i<72;++i){
      bool amITet = stk::isTetrahedron(implementedTopologies[i]);
      if (i >47 && i<56){
        EXPECT_TRUE(amITet);
      }
      else{
        EXPECT_FALSE(amITet);
      }
  }

}

