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
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <fstream>
#include "Assembly.hpp"

namespace stk
{
namespace io
{
namespace unit_test
{

TEST_F(Assembly, readWriteAssembly_simple_blocks)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& block1Part = create_io_part(partNames[0]);
  stk::mesh::Part& block2Part = create_io_part(partNames[1]);
  declare_subsets(assemblyPart, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);

  test_write_then_read_block_assemblies(1);
}

TEST_F(Assembly, readWriteAssembly_simple_surfaces)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"surface_1", "surface_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& surface1Part = create_io_part(partNames[0], 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = create_io_part(partNames[1], 2, stk::topology::QUAD_4);
  declare_subsets(assemblyPart, {&surface1Part, &surface2Part});
  stk::io::fill_mesh("generated:2x2x2|sideset:xX", get_bulk());

  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface1Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface2Part, get_bulk().buckets(get_meta().side_rank())));

  stk::mesh::PartVector inputSideBlocksToIgnore;
  inputSideBlocksToIgnore.push_back(get_meta().get_part("surface_1_quad4"));
  inputSideBlocksToIgnore.push_back(get_meta().get_part("surface_2_quad4"));
  test_for_null_parts(inputSideBlocksToIgnore);

  stk::mesh::PartVector outputPartsToIgnore;

  test_write_then_read_surface_assemblies(1, outputPartsToIgnore, inputSideBlocksToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_nodesets)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"nodelist_1", "nodelist_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& nodeset1Part = create_io_part(partNames[0], 1, stk::topology::NODE);
  stk::mesh::Part& nodeset2Part = create_io_part(partNames[1], 2, stk::topology::NODE);
  declare_subsets(assemblyPart, {&nodeset1Part, &nodeset2Part});
  stk::io::fill_mesh("generated:2x2x2|nodeset:xX", get_bulk());

  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset1Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset2Part, get_bulk().buckets(stk::topology::NODE_RANK)));

  stk::mesh::PartVector inputPartsToIgnore;
  stk::mesh::PartVector outputPartsToIgnore;

  test_write_then_read_nodeset_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_emptyblock)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& block1Part = create_io_part(partNames[0], 1);
  stk::mesh::Part& block2Part = create_io_part(partNames[1], 2);
  declare_subsets(assemblyPart, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());

  test_write_then_read_block_assemblies(1, stk::mesh::PartVector{&block2Part});
}

TEST_F(Assembly, readWriteAssembly_simple_emptysurface)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"surface_1", "surface_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& surface1Part = create_io_part(partNames[0], 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = create_io_part(partNames[1], 2, stk::topology::QUAD_4);
  declare_subsets(assemblyPart, {&surface1Part, &surface2Part});
  stk::io::fill_mesh("generated:2x2x2|sideset:x", get_bulk());

  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface1Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(surface2Part, get_bulk().buckets(get_meta().side_rank())));

  stk::mesh::PartVector inputSideBlocksToIgnore;
  inputSideBlocksToIgnore.push_back(get_meta().get_part("surface_1_quad4"));
  test_for_null_parts(inputSideBlocksToIgnore);

  stk::mesh::PartVector outputPartsToIgnore;

  test_write_then_read_surface_assemblies(1, outputPartsToIgnore, inputSideBlocksToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_emptynodeset)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"nodelist_1", "nodelist_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& nodeset1Part = create_io_part(partNames[0], 1, stk::topology::NODE);
  stk::mesh::Part& nodeset2Part = create_io_part(partNames[1], 2, stk::topology::NODE);
  declare_subsets(assemblyPart, {&nodeset1Part, &nodeset2Part});
  stk::io::fill_mesh("generated:2x2x2|nodeset:x", get_bulk());

  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset1Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(nodeset2Part, get_bulk().buckets(stk::topology::NODE_RANK)));

  stk::mesh::PartVector inputPartsToIgnore;
  stk::mesh::PartVector outputPartsToIgnore;

  test_write_then_read_nodeset_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_excludeBlock2)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& block1Part = create_io_part(partNames[0], 1);
  stk::mesh::Part& block2Part = create_io_part(partNames[1], 2);
  declare_subsets(assemblyPart, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);

  test_write_then_read_block_assemblies(1, {&block2Part});
}

TEST_F(Assembly, readWriteAssembly_simple_excludeSurface2)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"surface_1", "surface_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& surface1Part = create_io_part(partNames[0], 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = create_io_part(partNames[1], 2, stk::topology::QUAD_4);
  declare_subsets(assemblyPart, {&surface1Part, &surface2Part});
  stk::io::fill_mesh("generated:2x2x2|sideset:xX", get_bulk());

  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface1Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface2Part, get_bulk().buckets(get_meta().side_rank())));

  stk::mesh::PartVector outputPartsToIgnore;
  outputPartsToIgnore.push_back(get_meta().get_part("surface_2"));
  test_for_null_parts(outputPartsToIgnore);

  stk::mesh::PartVector inputPartsToIgnore;
  inputPartsToIgnore.push_back(get_meta().get_part("surface_1_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_2_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_2"));
  test_for_null_parts(inputPartsToIgnore);

  test_write_then_read_surface_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_excludeNodeset2)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"nodelist_1", "nodelist_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& nodeset1Part = create_io_part(partNames[0], 1, stk::topology::NODE);
  stk::mesh::Part& nodeset2Part = create_io_part(partNames[1], 2, stk::topology::NODE);
  declare_subsets(assemblyPart, {&nodeset1Part, &nodeset2Part});
  stk::io::fill_mesh("generated:2x2x2|nodeset:xX", get_bulk());

  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset1Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset2Part, get_bulk().buckets(stk::topology::NODE_RANK)));

  stk::mesh::PartVector inputPartsToIgnore;
  inputPartsToIgnore.push_back(get_meta().get_part("nodelist_2"));
  test_for_null_parts(inputPartsToIgnore);

  stk::mesh::PartVector outputPartsToIgnore;
  outputPartsToIgnore.push_back(get_meta().get_part("nodelist_2"));
  test_for_null_parts(outputPartsToIgnore);

  test_write_then_read_nodeset_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_excludeBothAssemblyBlocks)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames {"block_1", "block_2", "block_3"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& block1Part = create_io_part(partNames[0], 1);
  stk::mesh::Part& block2Part = create_io_part(partNames[1], 2);
  stk::mesh::Part& block3Part = create_io_part(partNames[2], 3);
  declare_subsets(assemblyPart, {&block2Part, &block3Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);
  move_element(2, block1Part, block3Part);

  test_write_then_read_block_assemblies(0, {&block2Part, &block3Part});
}

TEST_F(Assembly, readWriteAssembly_simple_excludeBothAssemblySurfaces)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"surface_1", "surface_2", "surface_3"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& surface1Part = create_io_part(partNames[0], 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = create_io_part(partNames[1], 2, stk::topology::QUAD_4);
  stk::mesh::Part& surface3Part = create_io_part(partNames[2], 3, stk::topology::QUAD_4);
  declare_subsets(assemblyPart, {&surface2Part, &surface3Part});
  stk::io::fill_mesh("generated:2x2x2|sideset:xyz", get_bulk());

  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface1Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface2Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface3Part, get_bulk().buckets(get_meta().side_rank())));

  stk::mesh::PartVector outputPartsToIgnore;
  outputPartsToIgnore.push_back(get_meta().get_part("surface_2"));
  outputPartsToIgnore.push_back(get_meta().get_part("surface_3"));
  test_for_null_parts(outputPartsToIgnore);

  stk::mesh::PartVector inputPartsToIgnore;
  inputPartsToIgnore.push_back(get_meta().get_part("surface_1_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_2_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_3_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_2"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_3"));
  test_for_null_parts(inputPartsToIgnore);

  test_write_then_read_surface_assemblies(0, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_excludeBothAssemblyNodesets)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"nodelist_1", "nodelist_2", "nodelist_3"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& nodeset1Part = create_io_part(partNames[0], 1, stk::topology::NODE);
  stk::mesh::Part& nodeset2Part = create_io_part(partNames[1], 2, stk::topology::NODE);
  stk::mesh::Part& nodeset3Part = create_io_part(partNames[2], 3, stk::topology::NODE);
  declare_subsets(assemblyPart, {&nodeset2Part, &nodeset3Part});
  stk::io::fill_mesh("generated:2x2x2|nodeset:xyz", get_bulk());

  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset1Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset2Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset3Part, get_bulk().buckets(stk::topology::NODE_RANK)));

  stk::mesh::PartVector inputPartsToIgnore;
  inputPartsToIgnore.push_back(get_meta().get_part("nodelist_2"));
  inputPartsToIgnore.push_back(get_meta().get_part("nodelist_3"));
  test_for_null_parts(inputPartsToIgnore);

  stk::mesh::PartVector outputPartsToIgnore;
  outputPartsToIgnore.push_back(get_meta().get_part("nodelist_2"));
  outputPartsToIgnore.push_back(get_meta().get_part("nodelist_3"));
  test_for_null_parts(outputPartsToIgnore);

  test_write_then_read_nodeset_assemblies(0, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_includeBothAssemblySurfaces)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"surface_1", "surface_2", "surface_3"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& surface1Part = create_io_part(partNames[0], 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2Part = create_io_part(partNames[1], 2, stk::topology::QUAD_4);
  stk::mesh::Part& surface3Part = create_io_part(partNames[2], 3, stk::topology::QUAD_4);
  declare_subsets(assemblyPart, {&surface2Part, &surface3Part});
  stk::io::fill_mesh("generated:2x2x2|sideset:xyz", get_bulk());

  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface1Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface2Part, get_bulk().buckets(get_meta().side_rank())));
  EXPECT_EQ(4u, stk::mesh::count_selected_entities(surface3Part, get_bulk().buckets(get_meta().side_rank())));

  stk::mesh::PartVector outputPartsToIgnore;
  stk::mesh::PartVector inputPartsToIgnore;
  inputPartsToIgnore.push_back(get_meta().get_part("surface_2_quad4"));
  inputPartsToIgnore.push_back(get_meta().get_part("surface_3_quad4"));
  test_for_null_parts(inputPartsToIgnore);

  test_write_then_read_surface_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_simple_includeBothAssemblyNodesets)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames{"nodelist_1", "nodelist_2", "nodelist_3"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 10);
  stk::mesh::Part& nodeset1Part = create_io_part(partNames[0], 1, stk::topology::NODE);
  stk::mesh::Part& nodeset2Part = create_io_part(partNames[1], 2, stk::topology::NODE);
  stk::mesh::Part& nodeset3Part = create_io_part(partNames[2], 3, stk::topology::NODE);
  declare_subsets(assemblyPart, {&nodeset2Part, &nodeset3Part});
  stk::io::fill_mesh("generated:2x2x2|nodeset:xyz", get_bulk());

  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset1Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset2Part, get_bulk().buckets(stk::topology::NODE_RANK)));
  EXPECT_EQ(9u, stk::mesh::count_selected_entities(nodeset3Part, get_bulk().buckets(stk::topology::NODE_RANK)));

  stk::mesh::PartVector inputPartsToIgnore;
  stk::mesh::PartVector outputPartsToIgnore;

  test_write_then_read_nodeset_assemblies(1, outputPartsToIgnore, inputPartsToIgnore);
}

TEST_F(Assembly, readWriteAssembly_multiple)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assembly1Name("assembly1");
  const std::string assembly2Name("assembly2");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assembly1Part = create_assembly(assembly1Name, 100);
  stk::mesh::Part& assembly2Part = create_assembly(assembly2Name, 200);
  stk::mesh::Part& block1Part = create_io_part(partNames[0]);
  stk::mesh::Part& block2Part = create_io_part(partNames[1]);
  declare_subsets(assembly1Part, {&block1Part, &block2Part});
  declare_subsets(assembly2Part, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);

  test_write_then_read_block_assemblies(2);
}

TEST_F(Assembly, readWriteAssembly_hierarchy)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");
  std::vector<std::string> leafPartNames {"block_1", "block_2", "block_3"};

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName, 1000,
                                                     {subAssembly1Name, subAssembly2Name},
                                                     {2000, 2001});
  stk::mesh::Part& block1Part = create_io_part(leafPartNames[0]);
  stk::mesh::Part& block2Part = create_io_part(leafPartNames[1]);
  stk::mesh::Part& block3Part = create_io_part(leafPartNames[2]);

  ASSERT_EQ(2u, subAssemblyParts.size());
  declare_subsets(*subAssemblyParts[0], {&block1Part, &block2Part});
  declare_subsets(*subAssemblyParts[1], {&block2Part, &block3Part});

  stk::io::fill_mesh("generated:2x2x2", get_bulk());

  move_element(1, block1Part, block2Part);
  move_element(2, block1Part, block3Part);

  test_write_then_read_block_assemblies(3);
}

TEST_F(Assembly, readWriteAssembly_deeperHierarchy)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  create_deep_assembly_hierarchy();

  stk::io::fill_mesh("generated:2x2x2", get_bulk());

  stk::mesh::Part& block1Part = *get_meta().get_part("block_1");
  stk::mesh::Part& block2Part = *get_meta().get_part("block_2");
  stk::mesh::Part& block3Part = *get_meta().get_part("block_3");

  move_element(1, block1Part, block2Part);
  move_element(2, block1Part, block3Part);

  test_write_then_read_block_assemblies(6);
}

TEST_F(AssemblyFilter, omitLeafAssembly)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string fileName = "assembly.g";
  unsigned nBlock = 5;
  create_multi_block_mesh_with_deep_assembly(nBlock, fileName);
  read_mesh_with_assembly_filter(fileName, {"assembly_9005"}, {});

  test_omit_leaf_assembly();
}

TEST_F(AssemblyFilter, omitBranchAssembly)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string fileName = "assembly.g";
  unsigned nBlock = 5;
  create_multi_block_mesh_with_deep_assembly(nBlock, fileName);
  read_mesh_with_assembly_filter(fileName, {"assembly_9003"}, {});

  test_omit_branch_assembly();
}

TEST_F(AssemblyFilter, omitAndIncludeWillThrow)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string fileName = "assembly.g";
  unsigned nBlock = 5;
  create_multi_block_mesh_with_deep_assembly(nBlock, fileName);

  EXPECT_THROW(read_mesh_with_assembly_filter(fileName, {"assembly_9001"}, {"assembly_9002"}), std::runtime_error);
}

TEST_F(AssemblyFilter, includeLeafAssembly)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string fileName = "assembly.g";
  unsigned nBlock = 5;
  create_multi_block_mesh_with_deep_assembly(nBlock, fileName);
  read_mesh_with_assembly_filter(fileName, {}, {"assembly_9002"});

  test_include_leaf_assembly();
}

TEST_F(AssemblyFilter, includeBranchAssembly)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string fileName = "assembly.g";
  unsigned nBlock = 5;
  create_multi_block_mesh_with_deep_assembly(nBlock, fileName);
  read_mesh_with_assembly_filter(fileName, {}, {"assembly_9003"});

  test_include_branch_assembly();
}

}  // namespace unit_test
}  // namespace io
}  // namespace stk
