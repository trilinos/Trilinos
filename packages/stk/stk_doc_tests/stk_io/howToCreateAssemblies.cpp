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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <stddef.h>                     // for size_t, NULL
#include <unistd.h>                     // for unlink
#include <stk_util/stk_config.h>
#include <stk_io/IossBridge.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Part.hpp>
#include <string>                       // for string
#include "stk_util/parallel/Parallel.hpp"

TEST(Assemblies, createAssemblyWithElementBlocks)
{
  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::Part& block1Part = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block2Part = meta.declare_part_with_topology("block_2", stk::topology::TRI_3_2D);

  stk::io::put_io_part_attribute(block1Part);
  stk::io::put_io_part_attribute(block2Part);

  stk::mesh::Part& assemblyPart = meta.declare_part("myAssembly");

  EXPECT_FALSE(stk::io::is_part_assembly_io_part(assemblyPart));
  EXPECT_FALSE(stk::io::is_part_io_part(assemblyPart));

  stk::io::put_assembly_io_part_attribute(assemblyPart);
  meta.set_part_id(assemblyPart, 100);

  EXPECT_TRUE(stk::io::is_part_assembly_io_part(assemblyPart));
  EXPECT_TRUE(stk::io::is_part_io_part(assemblyPart));

  meta.declare_part_subset(assemblyPart, block1Part);
  meta.declare_part_subset(assemblyPart, block2Part);

  EXPECT_TRUE(assemblyPart.contains(block1Part));
  EXPECT_TRUE(assemblyPart.contains(block2Part));
}

TEST(Assemblies, createAssemblyWithElementBlocksAndSurfaces)
{
  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::Part& block1Part = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  stk::mesh::Part& block2Part = meta.declare_part_with_topology("block_2", stk::topology::WEDGE_6);
  stk::mesh::Part& block3Part = meta.declare_part_with_topology("block_3", stk::topology::TET_4);

  stk::mesh::Part& surface1Part = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = meta.declare_part_with_topology("surface_2", stk::topology::TRI_3);

  stk::io::put_io_part_attribute(block1Part);
  stk::io::put_io_part_attribute(block2Part);
  stk::io::put_io_part_attribute(block3Part);
  stk::io::put_io_part_attribute(surface1Part);
  stk::io::put_io_part_attribute(surface2Part);

  std::string parentAssemblyName("parentAssembly");
  std::string subAssembly1Name("subAssembly1");
  std::string subAssembly2Name("subAssembly2");

  stk::mesh::Part& parentAssemblyPart = meta.declare_part(parentAssemblyName);
  stk::mesh::Part& subAssembly1Part = meta.declare_part(subAssembly1Name);
  stk::mesh::Part& subAssembly2Part = meta.declare_part(subAssembly2Name);

  stk::io::put_assembly_io_part_attribute(parentAssemblyPart);
  stk::io::put_assembly_io_part_attribute(subAssembly1Part);
  stk::io::put_assembly_io_part_attribute(subAssembly2Part);

  meta.declare_part_subset(parentAssemblyPart, subAssembly1Part);
  meta.declare_part_subset(parentAssemblyPart, subAssembly2Part);

  meta.declare_part_subset(subAssembly1Part, block1Part);
  meta.declare_part_subset(subAssembly1Part, block2Part);
  meta.declare_part_subset(subAssembly1Part, block3Part);

  meta.declare_part_subset(subAssembly2Part, surface1Part);
  meta.declare_part_subset(subAssembly2Part, surface2Part);

  EXPECT_TRUE(stk::io::has_sub_assemblies(meta, parentAssemblyName));

  stk::mesh::PartVector allLeafParts = stk::io::get_unique_leaf_parts(meta, parentAssemblyName);
  EXPECT_EQ(5u, allLeafParts.size());

  stk::mesh::PartVector subAssembly1LeafParts = stk::io::get_unique_leaf_parts(meta, subAssembly1Name);
  EXPECT_EQ(3u, subAssembly1LeafParts.size());

  stk::mesh::PartVector subAssembly2LeafParts = stk::io::get_unique_leaf_parts(meta, subAssembly2Name);
  EXPECT_EQ(2u, subAssembly2LeafParts.size());

  for(const stk::mesh::Part* subAssemblyPart : subAssembly1LeafParts) {
    EXPECT_TRUE(stk::mesh::contains(allLeafParts, *subAssemblyPart));
  }

  for(const stk::mesh::Part* subAssemblyPart : subAssembly2LeafParts) {
    EXPECT_TRUE(stk::mesh::contains(allLeafParts, *subAssemblyPart));
  }
}

TEST(Assemblies, cannotCreateAssemblyWithMixedRanks)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(3);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Part& block1Part = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  stk::mesh::Part& surface1Part = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);

  stk::io::put_io_part_attribute(block1Part);
  stk::io::put_io_part_attribute(surface1Part);

  std::string assemblyName("assembly");

  stk::mesh::Part& assemblyPart = meta.declare_part(assemblyName);
  stk::io::put_assembly_io_part_attribute(assemblyPart);

  meta.declare_part_subset(assemblyPart, block1Part);
  meta.declare_part_subset(assemblyPart, surface1Part);

  stk::io::fill_mesh("generated:2x2x2|sideset:x", *bulk);

  EXPECT_THROW(stk::io::write_mesh("outputFile.g", *bulk), std::runtime_error);
}

