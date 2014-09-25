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

#include <stddef.h>                     // for size_t, NULL
#include <algorithm>                    // for count
#include <stdexcept>                    // for runtime_error, logic_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include "Shards_BasicTopologies.hpp"   // for getCellTopologyData, etc
#include "Shards_CellTopology.hpp"      // for operator==, CellTopology
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc


using stk::mesh::MetaData;

//----------------------------------------------------------------------------

TEST ( UnitTestMetaData, create )
{
  stk::mesh::MetaData fem_meta;
  EXPECT_TRUE ( true );
  EXPECT_FALSE( fem_meta.is_initialized() );
  EXPECT_EQ( fem_meta.spatial_dimension(), 0u );
  EXPECT_EQ( fem_meta.side_rank(), static_cast<unsigned>(stk::topology::INVALID_RANK) );

  // Verify throws/etc for FEM calls prior to initialization:
  stk::mesh::CellTopology invalid_cell_topology( NULL );
  stk::mesh::Part & universal_part = fem_meta.universal_part();
  ASSERT_THROW( fem_meta.register_cell_topology( invalid_cell_topology, stk::topology::INVALID_RANK ), std::logic_error );
  ASSERT_THROW( fem_meta.get_cell_topology_root_part( invalid_cell_topology), std::logic_error );
  ASSERT_THROW( fem_meta.get_cell_topology( universal_part), std::logic_error );
  ASSERT_THROW( stk::mesh::set_cell_topology( universal_part, invalid_cell_topology), std::logic_error );
}

TEST( UnitTestMetaData, initialize )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  EXPECT_TRUE( fem_meta.is_initialized() );
  EXPECT_EQ( fem_meta.spatial_dimension(), spatial_dimension );
}

TEST( UnitTestMetaData, initialize_only_once )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  ASSERT_THROW( fem_meta.initialize(2), std::runtime_error );
}

TEST( UnitTestMetaData, entity_ranks_1 )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 1;
  fem_meta.initialize(spatial_dimension);
  EXPECT_EQ( fem_meta.side_rank(), static_cast<unsigned>(stk::topology::NODE_RANK) );
}

TEST( UnitTestMetaData, entity_ranks_2 )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);
  EXPECT_EQ( fem_meta.side_rank(), static_cast<unsigned>(stk::topology::EDGE_RANK) );
}

TEST( UnitTestMetaData, entity_ranks_3 )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  EXPECT_EQ( fem_meta.side_rank(), static_cast<unsigned>(stk::topology::FACE_RANK) );
}

TEST( UnitTestMetaData, get_cell_topology_trivial )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.get_cell_topology_root_part(hex_top);

  EXPECT_TRUE( stk::mesh::is_auto_declared_part(hex_part) );
  EXPECT_EQ( hex_part.primary_entity_rank(), spatial_dimension );
  stk::mesh::CellTopology topology = fem_meta.get_cell_topology(hex_part);
  EXPECT_EQ( (topology == hex_top), true );
}

TEST( UnitTestMetaData, get_cell_topology_simple )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.get_cell_topology_root_part(hex_top);
  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  EXPECT_TRUE( stk::mesh::is_auto_declared_part(hex_part) );
  EXPECT_TRUE( !stk::mesh::is_auto_declared_part(A) );
  fem_meta.declare_part_subset( hex_part, A );
  stk::mesh::CellTopology topology = fem_meta.get_cell_topology(A);
  ASSERT_EQ( (topology == hex_top), true );
}

TEST( UnitTestMetaData, get_cell_topology_invalid )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);
  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  ASSERT_THROW( fem_meta.get_cell_topology_root_part(hex_top), std::runtime_error );
}

TEST( UnitTestMetaData, cell_topology_subsetting )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  const stk::mesh::EntityRank element_rank = static_cast<stk::mesh::EntityRank>(spatial_dimension);
  fem_meta.initialize(spatial_dimension);
  stk::mesh::Part & element_part = fem_meta.declare_part("element part", element_rank );
  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::set_cell_topology( element_part, hex_top );

  stk::mesh::Part & hex_part = fem_meta.get_cell_topology_root_part(hex_top);

  const stk::mesh::PartVector & element_part_supersets = element_part.supersets();
  EXPECT_EQ(
      std::count(element_part_supersets.begin(),element_part_supersets.end(),&hex_part), 1
      );
}

// 02/16/11:  Cell Topology Induced Membership
//
// Invariants:
// 1.  Root cell topology parts cannot be subsets of parts with cell topologies
// 2.  Incompatible cell topologies are prohibited.  I.e. parts with
//     different cell topologies of the same rank (as the part) cannot be subsets
//     of each other.
//
// Decision tree for declare_part_subset(superset,subset):
// Q:  Does the superset part have a topology?
// A:  No
//     [ Okay, go forward ] (Test 1)
// A:  Yes
//     Q:  Is the subset part a root cell topology part or are any of the subset's parts subsets a root cell topology part?
//     A:  Yes
//         [ Throw ] (Test 2 a,b)
//     A:  No
//         Q:  How many cell topologies of the same rank as the superset are defined for the subset part and the subset's subsets parts?
//         A:  0
//             [ Okay, go forward ] (Test 3 a,b,c)
//         A:  1+
//             Q:  Are all the cell topologies the same?
//             A:  Yes
//                 [ Okay, go forward ] (Test 4 a,b)
//             A:  No
//                 [ Throw ] (Test 5 a,b,c)
//
// Tests:
// The following parts are necessary for the tests [spatial_dimension = 3]:
// Part A, rank 3, no topology
// Part B, rank 3, no topology
// Part C, rank 2, no topology
// Part D, rank 2, no topolooy
// Part E, rank 3, no topology
// Part F, rank 3, no topolooy
// Part H, rank 3, Hex<8> topology, HR > H (HR = Hex<8> root cell topology part)
// Part Q, rank 2, Quad<4> topology, QR > Q (QR = Quad<4> root cell topology part)
// Part W, rank 3, Wedge<6> topology, WR > W (WR = Wedge<6> root cell topology part)
// 1:   Part A > HR -> Okay
// 2a:  HR > QR -> Throw
// 2b:  HR > A, B > QR, A > B -> Throw
// 3a:  Subset has no cell topology and subset's subsets have no cell topology
//      HR > A -> Okay
// 3b:  Different rank cell topology on subset
//      QR > D, HR > A, A > D -> Okay
// 3c:  Different rank cell topology on subset's subset
//      QR > C, B > C, HR > A, A > B -> Okay
// 4a:  Subset has same cell topology
//      HR > A, HR > B, A > B -> Okay
// 4b:  Subset's subsets have same cell topology
//      HR > C, B > C, HR > A, A > B -> Okay
// 5a:  Subset has different cell topology
//      WR > A, HR > B, A > B -> Throw
// 5b:  Subset's subsets have different cell topology
//      WR > E, B > E, HR > A, A > B -> Throw
// 5c:  Multiple different cell topologies in subset's subsets
//      HR > F, WR > E, B > F, B > E, HR > A, A > B -> Throw

// 1:   Part A > HR -> Okay
TEST( UnitTestMetaData, cell_topology_test_1 )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );

  fem_meta.declare_part_subset(A, HR);
  const stk::mesh::PartVector & HR_supersets = HR.supersets();
  EXPECT_EQ(
      std::count(HR_supersets.begin(),HR_supersets.end(),&A), 1
      );
}

// 2a:  HR > QR -> Throw
TEST( UnitTestMetaData, cell_topology_test_2a )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);
  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::CellTopology QR_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());

  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);
  stk::mesh::Part & QR = fem_meta.get_cell_topology_root_part(QR_top);

  ASSERT_THROW( fem_meta.declare_part_subset(HR, QR), std::runtime_error );
}

// 2b:  HR > A, B > QR, A > B -> Throw
TEST( UnitTestMetaData, cell_topology_test_2b )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology QR_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk::mesh::Part & QR = fem_meta.get_cell_topology_root_part(QR_top);

  fem_meta.declare_part_subset( HR, A );
  fem_meta.declare_part_subset( B, QR );
  ASSERT_THROW( fem_meta.declare_part_subset(A, B), std::runtime_error );
}

// 3a:  Subset has no cell topology and subset's subsets have no cell topology
//      HR > A -> Okay
TEST( UnitTestMetaData, cell_topology_test_3a )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  fem_meta.declare_part_subset( HR, A );

  const stk::mesh::PartVector & A_supersets = A.supersets();
  EXPECT_EQ(
      std::count(A_supersets.begin(),A_supersets.end(),&HR), 1
      );
}

// 3b:  Different rank cell topology on subset
//      QR > D, HR > A, A > D -> Okay
TEST( UnitTestMetaData, cell_topology_test_3b )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & D = fem_meta.declare_part("Part D", fem_meta.side_rank() );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology QR_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk::mesh::Part & QR = fem_meta.get_cell_topology_root_part(QR_top);

  fem_meta.declare_part_subset( QR, D );
  fem_meta.declare_part_subset( HR, A );
  fem_meta.declare_part_subset( A, D );

  const stk::mesh::PartVector & D_supersets = D.supersets();
  EXPECT_EQ(
      std::count(D_supersets.begin(),D_supersets.end(),&A), 1
      );
}

// 3c:  Different rank cell topology on subset's subset
//      QR > C, B > C, HR > A, A > B -> Okay
TEST( UnitTestMetaData, cell_topology_test_3c )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & C = fem_meta.declare_part("Part C", fem_meta.side_rank() );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology QR_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk::mesh::Part & QR = fem_meta.get_cell_topology_root_part(QR_top);

  fem_meta.declare_part_subset( QR, C );
  fem_meta.declare_part_subset( B, C );
  fem_meta.declare_part_subset( HR, A );
  fem_meta.declare_part_subset( A, B );

  const stk::mesh::PartVector & B_supersets = B.supersets();
  EXPECT_EQ(
      std::count(B_supersets.begin(),B_supersets.end(),&A), 1
      );
}

// 4a:  Subset has same cell topology
//      HR > A, HR > B, A > B -> Okay
TEST( UnitTestMetaData, cell_topology_test_4a )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  fem_meta.declare_part_subset( HR, A );
  fem_meta.declare_part_subset( HR, B );
  fem_meta.declare_part_subset( A, B );

  const stk::mesh::PartVector & B_supersets = B.supersets();
  EXPECT_EQ(
      std::count(B_supersets.begin(),B_supersets.end(),&A), 1
      );
}

// 4b:  Subset's subsets have same cell topology
//      HR > C, B > C, HR > A, A > B -> Okay
TEST( UnitTestMetaData, cell_topology_test_4b )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & C = fem_meta.declare_part("Part C", fem_meta.side_rank() );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  fem_meta.declare_part_subset( HR, C );
  fem_meta.declare_part_subset( B, C );
  fem_meta.declare_part_subset( HR, A );
  fem_meta.declare_part_subset( A, B );

  const stk::mesh::PartVector & B_supersets = B.supersets();
  EXPECT_EQ(
      std::count(B_supersets.begin(),B_supersets.end(),&A), 1
      );
}

// 5a:  Subset has different cell topology
//      WR > A, HR > B, A > B -> Throw
TEST( UnitTestMetaData, cell_topology_test_5a )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology WR_top(shards::getCellTopologyData<shards::Wedge<6> >());
  stk::mesh::Part & WR = fem_meta.get_cell_topology_root_part(WR_top);

  fem_meta.declare_part_subset( WR, A );
  fem_meta.declare_part_subset( HR, B );
  ASSERT_THROW( fem_meta.declare_part_subset(A, B), std::runtime_error );
}

// 5b:  Subset's subsets have different cell topology
//      WR > E, B > E, HR > A, A > B -> Throw
TEST( UnitTestMetaData, cell_topology_test_5b )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & E = fem_meta.declare_part("Part E", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology WR_top(shards::getCellTopologyData<shards::Wedge<6> >());
  stk::mesh::Part & WR = fem_meta.get_cell_topology_root_part(WR_top);

  fem_meta.declare_part_subset( WR, E );
  stk::mesh::CellTopology top = fem_meta.get_cell_topology(E);
  ASSERT_TRUE( top.isValid() );
  fem_meta.declare_part_subset( B, E );
  fem_meta.declare_part_subset( HR, A );
  ASSERT_THROW( fem_meta.declare_part_subset(A, B), std::runtime_error );
}

// 5c:  Multiple different cell topologies in subset's subsets
//      HR > F, WR > E, B > F, B > E, HR > A, A > B -> Throw
TEST( UnitTestMetaData, cell_topology_test_5c )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 3;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::Part & A = fem_meta.declare_part("Part A", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & B = fem_meta.declare_part("Part B", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & F = fem_meta.declare_part("Part F", stk::topology::ELEMENT_RANK );
  stk::mesh::Part & E = fem_meta.declare_part("Part E", stk::topology::ELEMENT_RANK );

  stk::mesh::CellTopology HR_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & HR = fem_meta.get_cell_topology_root_part(HR_top);

  stk::mesh::CellTopology WR_top(shards::getCellTopologyData<shards::Wedge<6> >());
  stk::mesh::Part & WR = fem_meta.get_cell_topology_root_part(WR_top);

  fem_meta.declare_part_subset( HR, F );
  fem_meta.declare_part_subset( WR, E );
  fem_meta.declare_part_subset( B, F );
  fem_meta.declare_part_subset( B, E );
  fem_meta.declare_part_subset( HR, A );
  ASSERT_THROW( fem_meta.declare_part_subset(A, B), std::runtime_error );
}

TEST( MetaData, register_cell_topology_duplicate )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);
  const stk::mesh::EntityRank hex_rank = stk::topology::ELEMENT_RANK;

  fem_meta.register_cell_topology( shards::getCellTopologyData<shards::Hexahedron<8> >(), hex_rank );
  ASSERT_NO_THROW( fem_meta.register_cell_topology( shards::getCellTopologyData<shards::Hexahedron<8> >(), hex_rank ) );
}

TEST( MetaData, register_cell_topology_duplicate_with_different_ranks )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);
  const stk::mesh::EntityRank hex_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank bad_rank = stk::topology::EDGE_RANK;

  fem_meta.register_cell_topology( shards::getCellTopologyData<shards::Hexahedron<8> >(), hex_rank );
  ASSERT_THROW( fem_meta.register_cell_topology( shards::getCellTopologyData<shards::Hexahedron<8> >(), bad_rank ), std::runtime_error );
}

TEST( MetaData, register_cell_topology_duplicate_with_invalid_rank )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);
  const stk::mesh::EntityRank invalid_rank = static_cast<stk::mesh::EntityRank>(4);

  ASSERT_THROW( fem_meta.register_cell_topology( shards::getCellTopologyData<shards::Hexahedron<8> >(), invalid_rank ), std::logic_error );
}

TEST( MetaData, get_cell_topology_root_part_invalid )
{
  stk::mesh::MetaData fem_meta;
  const size_t spatial_dimension = 2;
  fem_meta.initialize(spatial_dimension);

  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());

  ASSERT_THROW( fem_meta.get_cell_topology_root_part( hex_top ), std::runtime_error );
}
