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

class Assembly : public ::ngp_testing::Test {
protected:
  Assembly()
  : m_meta(2)
  {}

  stk::mesh::Part& create_assembly(const std::string& assemblyName)
  {
    stk::mesh::Part& assemblyPart = get_meta().declare_part(assemblyName);
    stk::io::put_assembly_io_part_attribute(assemblyPart);
    return assemblyPart;
  }

  stk::mesh::Part& create_io_part(const std::string& partName)
  {
    stk::mesh::Part& part = get_meta().declare_part(partName, stk::topology::ELEM_RANK);
    stk::io::put_io_part_attribute(part);
    return part;
  }

  void declare_subsets(stk::mesh::Part& parentPart, const stk::mesh::PartVector& subsetParts)
  {
    for(stk::mesh::Part* subsetPart : subsetParts) {
      get_meta().declare_part_subset(parentPart, *subsetPart);
    }
  }

  std::pair<stk::mesh::Part*, stk::mesh::PartVector>
  create_assembly_hierarchy(const std::string& parentAssemblyName,
                            const std::vector<std::string>& subAssemblyNames)
  {
    stk::mesh::Part& parentAssemblyPart = create_assembly(parentAssemblyName);

    stk::mesh::PartVector subAssemblyParts;
    for(const std::string& subAssemblyName : subAssemblyNames) {
      subAssemblyParts.push_back(&create_assembly(subAssemblyName));
    }

    declare_subsets(parentAssemblyPart, subAssemblyParts);

    return std::make_pair(&parentAssemblyPart, subAssemblyParts);
  }

  void test_assembly_part_attributes(const stk::mesh::Part& part)
  {
    EXPECT_TRUE(stk::io::is_part_assembly_io_part(part));
    EXPECT_TRUE(stk::io::is_part_io_part(part));
  }

  void test_assembly_part_attributes(const stk::mesh::PartVector& parts)
  {
    for(const stk::mesh::Part* part : parts) {
      test_assembly_part_attributes(*part);
    }
  }

  stk::mesh::MetaData& get_meta() { return m_meta; }

  stk::mesh::MetaData m_meta;
};

TEST_F(Assembly, ioPartIsNotAssembly)
{
  stk::mesh::Part& part = create_io_part("myPart");

  EXPECT_FALSE(stk::io::is_part_assembly_io_part(part));
}

TEST_F(Assembly, createAssembly)
{
  stk::mesh::Part& assemblyPart = create_assembly("myAssembly");

  test_assembly_part_attributes({&assemblyPart});
}

TEST_F(Assembly, getAssemblyNames_empty)
{
  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  EXPECT_TRUE(assemblyNames.empty());
}

TEST_F(Assembly, getAssemblyNames_singleAssembly)
{
  std::string assemblyName("myAssembly");
  create_assembly(assemblyName);

  create_io_part("myPart");

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(1u, assemblyNames.size());
  EXPECT_EQ(assemblyName, assemblyNames[0]);

  EXPECT_TRUE(stk::io::get_unique_leaf_parts(get_meta(), assemblyName).empty());
}

TEST_F(Assembly, getAssemblyNames_multipleAssemblies)
{
  std::string assemblyBBB("BBB_assembly");
  create_assembly(assemblyBBB);

  std::string assemblyZZZ("ZZZ_assembly");
  create_assembly(assemblyZZZ);

  std::string assemblyAAA("AAA_assembly");
  create_assembly(assemblyAAA);

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(3u, assemblyNames.size());
  EXPECT_EQ(assemblyBBB, assemblyNames[0]);
  EXPECT_EQ(assemblyZZZ, assemblyNames[1]);
  EXPECT_EQ(assemblyAAA, assemblyNames[2]);

  EXPECT_TRUE(stk::io::get_unique_leaf_parts(get_meta(), assemblyAAA).empty());
  EXPECT_TRUE(stk::io::get_unique_leaf_parts(get_meta(), assemblyBBB).empty());
  EXPECT_TRUE(stk::io::get_unique_leaf_parts(get_meta(), assemblyZZZ).empty());
}

TEST_F(Assembly, createAssemblyWithSubAssemblies)
{
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName,
                                                     {subAssembly1Name, subAssembly2Name});

  EXPECT_TRUE(stk::io::has_sub_assemblies(get_meta(), parentAssemblyName));

  test_assembly_part_attributes(*parentAssemblyPart);
  test_assembly_part_attributes(subAssemblyParts);
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_getAssemblyNames)
{
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  create_assembly_hierarchy(parentAssemblyName, {subAssembly1Name, subAssembly2Name});

  std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
  ASSERT_EQ(3u, assemblyNames.size());
  EXPECT_EQ(parentAssemblyName, assemblyNames[0]);
  EXPECT_EQ(subAssembly1Name, assemblyNames[1]);
  EXPECT_EQ(subAssembly2Name, assemblyNames[2]);
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_getSubAssemblyNames)
{
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  create_assembly_hierarchy(parentAssemblyName, {subAssembly1Name, subAssembly2Name});

  std::vector<std::string> subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), parentAssemblyName);
  ASSERT_EQ(2u, subAssemblyNames.size());
  EXPECT_EQ(subAssembly1Name, subAssemblyNames[0]);
  EXPECT_EQ(subAssembly2Name, subAssemblyNames[1]);

  EXPECT_TRUE(stk::io::get_sub_assembly_names(get_meta(), subAssembly1Name).empty());
  EXPECT_TRUE(stk::io::get_sub_assembly_names(get_meta(), subAssembly2Name).empty());
}

TEST_F(Assembly, createAssemblyWithSubAssemblies_leafParts)
{
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName,
                                                     {subAssembly1Name, subAssembly2Name});
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

TEST_F(Assembly, deeperHierarchy)
{
  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName,
                                                     {subAssembly1Name, subAssembly2Name});

  std::string subSubAssemblyName("mySubSubAssembly");
  std::string subSubSubAssembly1Name("mySubSubSubAssembly1");
  std::string subSubSubAssembly2Name("mySubSubSubAssembly2");

  stk::mesh::Part* subSubAssemblyPart;
  stk::mesh::PartVector subSubSubAssemblyParts;
  std::tie(subSubAssemblyPart, subSubSubAssemblyParts) = create_assembly_hierarchy(subSubAssemblyName,
                                                     {subSubSubAssembly1Name, subSubSubAssembly2Name});

  ASSERT_EQ(2u, subAssemblyParts.size());
  declare_subsets(*subAssemblyParts[0], {subSubAssemblyPart});
  declare_subsets(*subAssemblyParts[1], {subSubAssemblyPart});

  stk::mesh::Part& block1Part = create_io_part("block_1");
  stk::mesh::Part& block2Part = create_io_part("block_2");

  ASSERT_EQ(2u, subSubSubAssemblyParts.size());
  declare_subsets(*subSubSubAssemblyParts[0], {&block1Part, &block2Part});
  declare_subsets(*subSubSubAssemblyParts[1], {&block1Part, &block2Part});

  const bool bottomSubset1IsAddedToTopSuperset = parentAssemblyPart->contains(block1Part);
  const bool bottomSubset2IsAddedToTopSuperset = parentAssemblyPart->contains(block2Part);
  EXPECT_TRUE(bottomSubset1IsAddedToTopSuperset);
  EXPECT_TRUE(bottomSubset2IsAddedToTopSuperset);

  stk::mesh::PartVector leafParts = stk::io::get_unique_leaf_parts(get_meta(), parentAssemblyName);
  ASSERT_EQ(2u, leafParts.size());
  EXPECT_EQ(block1Part.mesh_meta_data_ordinal(), leafParts[0]->mesh_meta_data_ordinal());
  EXPECT_EQ(block2Part.mesh_meta_data_ordinal(), leafParts[1]->mesh_meta_data_ordinal());
}

