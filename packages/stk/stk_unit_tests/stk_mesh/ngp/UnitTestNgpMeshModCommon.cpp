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
#include <algorithm>
#include <iterator>
#include <map>
#include "Kokkos_Core.hpp"
#include "ngp/NgpUnitTestUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SkinMesh.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"

#include "ngp/UnitTestNgpMeshModificationUtils.hpp"

#ifdef STK_USE_DEVICE_MESH
namespace
{
template<typename ViewType>
void init_sorted(ViewType& view)
{
  auto devicePolicy = stk::ngp::DeviceRangePolicy(0,1);
  Kokkos::parallel_for("init", devicePolicy, KOKKOS_LAMBDA(const int& /*idx*/) {
    view(0) = 0;
    view(1) = 2;
    view(2) = 2;
    view(3) = 4;
    view(4) = 9;
  });
}

template<typename ViewType>
void init_unsorted(ViewType& view)
{
  auto devicePolicy = stk::ngp::DeviceRangePolicy(0,1);
  Kokkos::parallel_for("init", devicePolicy, KOKKOS_LAMBDA(const int& /*idx*/) {
    view(0) = 0;
    view(1) = 2;
    view(2) = 1;
    view(3) = 4;
    view(4) = 9;
  });
}

TEST(NgpMeshImpl, is_sorted)
{
  Kokkos::View<int*,stk::ngp::ExecSpace> emptyView("emptyView", 0);
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},emptyView));

  Kokkos::View<int*,stk::ngp::ExecSpace> sortedView("sortedView", 5);
  init_sorted(sortedView);

  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},sortedView));

  Kokkos::View<int*,stk::ngp::ExecSpace> unsortedView("unsortedView", 5);
  init_unsorted(unsortedView);

  EXPECT_FALSE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},unsortedView));
}

TEST(NgpMeshImpl, get_sorted_view)
{
  Kokkos::View<int*,stk::ngp::ExecSpace> unsortedView("unsortedView", 5);
  init_unsorted(unsortedView);

  Kokkos::View<int*,stk::ngp::ExecSpace> sortedView = stk::mesh::impl::get_sorted_view(unsortedView);
  EXPECT_TRUE(Kokkos::Experimental::is_sorted(stk::ngp::ExecSpace{},sortedView));
}

TEST(NgpMeshConstructionTest, prevent_update_during_mesh_mod)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(3);
  auto bulk = builder.create();

  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));

  bulk->modification_begin();
  EXPECT_ANY_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));
  bulk->modification_end();

  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));
}

TEST_F(NgpMeshMod, PartCorrectnessAfterDeviceMeshUpdate)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  create_node(*m_bulk, 1, {&part1});

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
  }

  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::NODE);
  create_node(*m_bulk, 2, {&part2});

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::NODE_RANK);
  }

  stk::mesh::Part& part3 = m_meta->declare_part("part3");

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part3, stk::topology::INVALID_RANK);
  }
}

TEST_F(NgpMeshMod, PartCorrectnessAfterDeviceMeshUpdate_NoEntitiesAdded)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);
  create_node(*m_bulk, 1, {&part1});

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
  }

  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::HEX_8);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::ELEM_RANK);
  }

  stk::mesh::Part& part3 = m_meta->declare_part("part3");

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::ELEM_RANK);
    check_part_is_on_device(ngpMesh, part3, stk::topology::INVALID_RANK);
  }
}

TEST_F(NgpMeshMod, PartCorrectnessAfterDeviceMeshUpdate_NoEntitiesEver)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  build_empty_mesh(1, 1);
  commit_meta_data();

  stk::mesh::Part& part1 = m_meta->declare_part_with_topology("part1", stk::topology::NODE);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
  }

  stk::mesh::Part& part2 = m_meta->declare_part_with_topology("part2", stk::topology::HEX_8);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::ELEM_RANK);
  }

  stk::mesh::Part& part3 = m_meta->declare_part("part3");

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    check_part_is_on_device(ngpMesh, part1, stk::topology::NODE_RANK);
    check_part_is_on_device(ngpMesh, part2, stk::topology::ELEM_RANK);
    check_part_is_on_device(ngpMesh, part3, stk::topology::INVALID_RANK);
  }
}

}  // namespace
#endif
