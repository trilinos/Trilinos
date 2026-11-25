// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <memory>

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

// The ids of a block and a sideset declared, in that order, in a new mesh of its own.
std::array<int64_t, 2> ids_of_io_parts_declared_in_a_new_mesh()
{
  std::shared_ptr<stk::mesh::MetaData> meta = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create_meta_data();
  AuxMetaData & auxMeta = AuxMetaData::create(*meta);
  const stk::mesh::Part & block = auxMeta.declare_io_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  const stk::mesh::Part & sideset = auxMeta.declare_io_part("surface_1", meta->side_rank());
  return {{block.id(), sideset.id()}};
}

} // namespace

TEST(AuxMetaData, ioPartsDeclaredInOneMesh_haveDistinctIds)
{
  const std::array<int64_t, 2> ids = ids_of_io_parts_declared_in_a_new_mesh();
  EXPECT_NE(stk::mesh::Part::INVALID_ID, ids[0]);
  EXPECT_NE(stk::mesh::Part::INVALID_ID, ids[1]);
  EXPECT_NE(ids[0], ids[1]);
}

// Each mesh numbers its parts itself, so the meshes a process makes one after another (as a study
// writing one file per frame does) give the same parts the same ids.
TEST(AuxMetaData, samePartsDeclaredInTwoMeshes_haveTheSameIds)
{
  const std::array<int64_t, 2> first = ids_of_io_parts_declared_in_a_new_mesh();
  const std::array<int64_t, 2> second = ids_of_io_parts_declared_in_a_new_mesh();
  EXPECT_EQ(first, second);
}

} // namespace krino
