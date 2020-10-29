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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <fstream>

class Assembly : public stk::unit_test_util::MeshFixture
{
protected:
  stk::mesh::Part& create_assembly(const std::string& assemblyName)
  {
    stk::mesh::Part& assemblyPart = get_meta().declare_part(assemblyName);
    stk::io::put_assembly_io_part_attribute(assemblyPart);
    return assemblyPart;
  }

  stk::mesh::Part& create_io_part(const std::string& partName)
  {
    stk::mesh::Part& part = get_meta().declare_part_with_topology(partName, stk::topology::HEX_8);
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

  void create_deep_assembly_hierarchy()
  {
    stk::mesh::Part& block1Part = create_io_part("block_1");
    stk::mesh::Part& block2Part = create_io_part("block_2");
    stk::mesh::Part& block3Part = create_io_part("block_3");
  
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
    declare_subsets(*subAssemblyParts[1], {&block3Part});
  
    ASSERT_EQ(2u, subSubSubAssemblyParts.size());
    declare_subsets(*subSubSubAssemblyParts[0], {&block1Part, &block2Part});
    declare_subsets(*subSubSubAssemblyParts[1], {&block1Part, &block2Part});

    const bool leafPart1IsAddedToTopSuperset = parentAssemblyPart->contains(block1Part);
    const bool leafPart2IsAddedToTopSuperset = parentAssemblyPart->contains(block2Part);
    const bool leafPart3IsAddedToTopSuperset = parentAssemblyPart->contains(block3Part);
    EXPECT_TRUE(leafPart1IsAddedToTopSuperset);
    EXPECT_TRUE(leafPart2IsAddedToTopSuperset);
    EXPECT_TRUE(leafPart3IsAddedToTopSuperset);
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

  void move_element(const stk::mesh::EntityId elemId,
                    stk::mesh::Part& sourcePart,
                    stk::mesh::Part& destPart)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    get_bulk().batch_change_entity_parts({elem}, {&destPart}, {&sourcePart});
  }

  void test_sub_assembly_names(const std::string& assemblyName,
                               const std::vector<std::string>& expectedSubAssemblyNames)
  {
    std::vector<std::string> subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), assemblyName);
    ASSERT_EQ(expectedSubAssemblyNames.size(), subAssemblyNames.size());
    for(size_t i=0; i<expectedSubAssemblyNames.size(); ++i) {
       EXPECT_EQ(expectedSubAssemblyNames[i], subAssemblyNames[i]);
    }
  }

  void compare_assemblies(const stk::mesh::MetaData& meta1,
                          const stk::mesh::MetaData& meta2)
  {
    std::vector<std::string> assemblyNames1 = stk::io::get_assembly_names(meta1);
    std::vector<std::string> assemblyNames2 = stk::io::get_assembly_names(meta2);

    std::sort(assemblyNames1.begin(), assemblyNames1.end());
    std::sort(assemblyNames2.begin(), assemblyNames2.end());
    EXPECT_EQ(assemblyNames1, assemblyNames2);

    for(size_t i=0; i<assemblyNames1.size(); ++i) {
      const stk::mesh::Part* assemblyPart1 = meta1.get_part(assemblyNames1[i]);
      const stk::mesh::Part* assemblyPart2 = meta2.get_part(assemblyNames2[i]);
      EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart1));
      EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart2));

      EXPECT_EQ(stk::io::get_sub_assembly_names(meta1, assemblyNames1[i]),
                stk::io::get_sub_assembly_names(meta2, assemblyNames2[i]));

      stk::mesh::PartVector leafParts1 = stk::io::get_unique_leaf_parts(meta1, assemblyNames1[i]);
      stk::mesh::PartVector leafParts2 = stk::io::get_unique_leaf_parts(meta2, assemblyNames2[i]);
      ASSERT_EQ(leafParts1.size(), leafParts2.size());
      for(size_t j=0; j<leafParts1.size(); ++j) {
        EXPECT_TRUE(stk::mesh::find(leafParts2, leafParts1[j]->name()) != nullptr);
      }
    }
  }

  void test_write_then_read_assemblies(size_t expectedNumAssemblies)
  {
    const std::string fileName("meshWithAssemblies.e");
    stk::io::write_mesh(fileName, get_bulk());

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh(fileName, bulk);

    EXPECT_EQ(expectedNumAssemblies, stk::io::get_assembly_names(meta).size());
    compare_assemblies(get_meta(), meta);

    unlink(fileName.c_str());
  }
};

