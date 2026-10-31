// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AuxMetaData.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshInterface.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Ionit_Initializer.h>
#include <gtest/gtest.h>

#include <string>
#include <vector>

// MeshFromFile opens an IOSS database; the IOSS db types (including the
// "generated" pseudo-database used below) must be registered. The main krino
// app does this at startup, but the unit-test main does not, so register here.
namespace { Ioss::Init::Initializer initialize_ioss; }

namespace krino {

TEST(MeshFromFile, GeneratedMesh_PopulatesMetaAndBulk)
{
  // 4x4x4 hex block over the unit cube: 125 nodes, 64 elements, independent of
  // decomposition. The IOSS "generated" database requires at least as many
  // Z-intervals as processors, so 4 lets this run on up to 4 ranks.
  const std::string meshSpec = "generated:4x4x4|bbox:0,0,0,1,1,1";

  MeshFromFile meshFromFile(meshSpec, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  meshFromFile.populate_mesh();

  EXPECT_EQ(3u, meshFromFile.meta_data().spatial_dimension());
  // AuxMetaData is created by the MeshFromFile ctor; confirm it is attached.
  EXPECT_TRUE(AuxMetaData::get(meshFromFile.meta_data()).active_part().name().size() > 0);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(meshFromFile.bulk_data(), counts);
  EXPECT_EQ(125u, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(64u,  counts[stk::topology::ELEM_RANK]);
}

// Exercises MeshInterface polymorphism through a base reference.
TEST(MeshInterface, MeshFromFileThroughBasePointer)
{
  // 4 Z-intervals so this distributes across up to 4 ranks (IOSS requires
  // Z-intervals >= processor count).
  const std::string meshSpec = "generated:1x1x4|bbox:0,0,0,1,1,4";

  MeshFromFile meshFromFile(meshSpec, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  MeshInterface & meshInterface = meshFromFile;
  meshInterface.populate_mesh();

  EXPECT_EQ(&meshFromFile.bulk_data(), &meshInterface.bulk_data());
  EXPECT_EQ(&meshFromFile.meta_data(), &meshInterface.meta_data());

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(meshInterface.bulk_data(), counts);
  EXPECT_EQ(20u, counts[stk::topology::NODE_RANK]); // 2x2x5 node lattice
  EXPECT_EQ(4u,  counts[stk::topology::ELEM_RANK]);
}

} // namespace krino
