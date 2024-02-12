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
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_topology/topology.hpp>
#include <stk_io/IossBridge.hpp>
#include "Assembly.hpp"

namespace stk
{
namespace io
{
namespace unit_test
{

TEST_F(Assembly, ioPartIsNotAssembly)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Part& part = create_io_part("myPart");

  EXPECT_FALSE(stk::io::is_part_assembly_io_part(part));
}

TEST_F(Assembly, createAssembly)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("myAssembly");
  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 100);

  test_assembly_part_attributes(assemblyPart);

  EXPECT_TRUE(stk::io::get_unique_leaf_parts(get_meta(), assemblyName).empty());
}

TEST_F(Assembly, getAssemblyNames_empty)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  EXPECT_TRUE(assemblyNames.empty());
}

TEST_F(Assembly, getAssemblyNames_singleAssembly)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string assemblyName("myAssembly");
  create_assembly(assemblyName, 100);

  create_io_part("myPart");

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(1u, assemblyNames.size());
  EXPECT_EQ(assemblyName, assemblyNames[0]);
}

TEST_F(Assembly, getAssemblyNames_multipleAssemblies_inOrderOfCreation)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string assemblyBBB("BBB_assembly");
  create_assembly(assemblyBBB, 100);

  std::string assemblyZZZ("ZZZ_assembly");
  create_assembly(assemblyZZZ, 200);

  std::string assemblyAAA("AAA_assembly");
  create_assembly(assemblyAAA, 300);

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(3u, assemblyNames.size());
  EXPECT_EQ(assemblyBBB, assemblyNames[0]);
  EXPECT_EQ(assemblyZZZ, assemblyNames[1]);
  EXPECT_EQ(assemblyAAA, assemblyNames[2]);
}

TEST_F(Assembly, createAssemblyWithSubAssemblies)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName, 100,
                                                     {subAssembly1Name, subAssembly2Name},
                                                     {200, 201});

  EXPECT_TRUE(stk::io::has_sub_assemblies(get_meta(), parentAssemblyName));

  test_assembly_part_attributes(*parentAssemblyPart);
  test_assembly_part_attributes(subAssemblyParts);
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_getAssemblyNames)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  create_assembly_hierarchy(parentAssemblyName, 100, {subAssembly1Name, subAssembly2Name}, {200, 201});

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(3u, assemblyNames.size());
  EXPECT_EQ(parentAssemblyName, assemblyNames[0]);
  EXPECT_EQ(subAssembly1Name, assemblyNames[1]);
  EXPECT_EQ(subAssembly2Name, assemblyNames[2]);
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_getSubAssemblyNames)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  create_assembly_hierarchy(parentAssemblyName, 100, {subAssembly1Name, subAssembly2Name}, {200, 201});

  std::vector<std::string> subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), parentAssemblyName);
  ASSERT_EQ(2u, subAssemblyNames.size());
  EXPECT_EQ(subAssembly1Name, subAssemblyNames[0]);
  EXPECT_EQ(subAssembly2Name, subAssemblyNames[1]);

  EXPECT_TRUE(stk::io::get_sub_assembly_names(get_meta(), subAssembly1Name).empty());
  EXPECT_TRUE(stk::io::get_sub_assembly_names(get_meta(), subAssembly2Name).empty());
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_leafParts)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName, 100,
                                                     {subAssembly1Name, subAssembly2Name}, {200, 201});
  stk::mesh::Part& block1Part = create_io_part("block_1");
  stk::mesh::Part& block2Part = create_io_part("block_2");
  stk::mesh::Part& block3Part = create_io_part("block_3");

  ASSERT_EQ(2u, subAssemblyParts.size());
  declare_subsets(*subAssemblyParts[0], {&block1Part, &block2Part});
  declare_subsets(*subAssemblyParts[1], {&block2Part, &block3Part});

  EXPECT_FALSE(stk::io::has_sub_assemblies(get_meta(), subAssembly1Name));
  EXPECT_FALSE(stk::io::has_sub_assemblies(get_meta(), subAssembly2Name));

  EXPECT_TRUE(subAssemblyParts[0]->contains(block1Part));
  EXPECT_TRUE(subAssemblyParts[0]->contains(block2Part));
  EXPECT_TRUE(subAssemblyParts[1]->contains(block2Part));
  EXPECT_TRUE(subAssemblyParts[1]->contains(block3Part));

  EXPECT_EQ(5u, parentAssemblyPart->subsets().size());

  stk::mesh::PartVector leafParts = stk::io::get_unique_leaf_parts(get_meta(), parentAssemblyName);

  ASSERT_EQ(3u, leafParts.size());
  EXPECT_EQ(block1Part.mesh_meta_data_ordinal(), leafParts[0]->mesh_meta_data_ordinal());
  EXPECT_EQ(block2Part.mesh_meta_data_ordinal(), leafParts[1]->mesh_meta_data_ordinal());
  EXPECT_EQ(block3Part.mesh_meta_data_ordinal(), leafParts[2]->mesh_meta_data_ordinal());
}

TEST_F(Assembly, deeperHierarchy_names)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  create_deep_assembly_hierarchy();

  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");
  std::string subSubAssemblyName("mySubSubAssembly");
  std::string subSubSubAssembly1Name("mySubSubSubAssembly1");
  std::string subSubSubAssembly2Name("mySubSubSubAssembly2");

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(6u, assemblyNames.size());
  EXPECT_EQ(parentAssemblyName,     assemblyNames[0]);
  EXPECT_EQ(subAssembly1Name,       assemblyNames[1]);
  EXPECT_EQ(subAssembly2Name,       assemblyNames[2]);
  EXPECT_EQ(subSubAssemblyName,     assemblyNames[3]);
  EXPECT_EQ(subSubSubAssembly1Name, assemblyNames[4]);
  EXPECT_EQ(subSubSubAssembly2Name, assemblyNames[5]);

  test_sub_assembly_names(parentAssemblyName, {subAssembly1Name, subAssembly2Name});
  test_sub_assembly_names(subAssembly1Name, {subSubAssemblyName});
  test_sub_assembly_names(subAssembly2Name, {});
  test_sub_assembly_names(subSubAssemblyName, {subSubSubAssembly1Name, subSubSubAssembly2Name});
  test_sub_assembly_names(subSubSubAssembly1Name, {});
  test_sub_assembly_names(subSubSubAssembly2Name, {});
}

TEST_F(Assembly, deeperHierarchy_parts)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  create_deep_assembly_hierarchy();

  stk::mesh::PartVector leafParts = stk::io::get_unique_leaf_parts(get_meta(), "myParentAssembly");
  ASSERT_EQ(3u, leafParts.size());
  EXPECT_EQ("block_1", leafParts[0]->name());
  EXPECT_EQ("block_2", leafParts[1]->name());
  EXPECT_EQ("block_3", leafParts[2]->name());
}

}  // namespace unit_test
}  // namespace io
}  // namespace stk
