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
#ifndef STK_UNIT_TESTS_STK_IO_Assembly_hpp
#define STK_UNIT_TESTS_STK_IO_Assembly_hpp

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "IOMeshFixture.hpp"
#include <algorithm>
#include <fstream>

namespace stk
{
namespace io
{
namespace unit_test
{
class Assembly : public IOMeshFixture
{
protected:
  stk::mesh::Part& create_assembly(const std::string& assemblyName, int id)
  {
    stk::mesh::Part& assemblyPart = get_meta().declare_part(assemblyName);
    stk::io::put_assembly_io_part_attribute(assemblyPart);
    get_meta().set_part_id(assemblyPart, id);
    return assemblyPart;
  }

  void declare_subsets(stk::mesh::Part& parentPart, const stk::mesh::PartVector& subsetParts)
  {
    for(stk::mesh::Part* subsetPart : subsetParts) {
      get_meta().declare_part_subset(parentPart, *subsetPart);
    }
  }

  std::pair<stk::mesh::Part*, stk::mesh::PartVector>
  create_assembly_hierarchy(const std::string& parentAssemblyName, int parentId,
                            const std::vector<std::string>& subAssemblyNames,
                            const std::vector<int>& subAssemblyIds)
  {
    stk::mesh::Part& parentAssemblyPart = create_assembly(parentAssemblyName, parentId);

    stk::mesh::PartVector subAssemblyParts;
    unsigned counter = 0;
    for(const std::string& subAssemblyName : subAssemblyNames) {
      subAssemblyParts.push_back(&create_assembly(subAssemblyName, subAssemblyIds[counter++]));
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
    std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName, 100,
                                                       {subAssembly1Name, subAssembly2Name},
                                                       {200, 201});
  
    std::string subSubAssemblyName("mySubSubAssembly");
    std::string subSubSubAssembly1Name("mySubSubSubAssembly1");
    std::string subSubSubAssembly2Name("mySubSubSubAssembly2");
  
    stk::mesh::Part* subSubAssemblyPart;
    stk::mesh::PartVector subSubSubAssemblyParts;
    std::tie(subSubAssemblyPart, subSubSubAssemblyParts) = create_assembly_hierarchy(subSubAssemblyName, 101,
                                                       {subSubSubAssembly1Name, subSubSubAssembly2Name},
                                                       {300, 301});
  
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

  void test_sub_assembly_names(const std::string& assemblyName,
                               const std::vector<std::string>& expectedSubAssemblyNames)
  {
    std::vector<std::string> subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), assemblyName);
    ASSERT_EQ(expectedSubAssemblyNames.size(), subAssemblyNames.size());
    for(size_t i=0; i<expectedSubAssemblyNames.size(); ++i) {
       EXPECT_EQ(expectedSubAssemblyNames[i], subAssemblyNames[i]);
    }
  }

  bool find_by_name(const stk::mesh::PartVector& parts, const std::string& partName)
  {
    auto nameMatches = [&](const stk::mesh::Part* part){return partName==part->name();};
    return parts.end() != std::find_if(parts.begin(), parts.end(), nameMatches);
  }

  void compare_assemblies(const stk::mesh::MetaData& meta1,
      const stk::mesh::MetaData& meta2,
      const stk::mesh::PartVector& excludedParts = stk::mesh::PartVector())
  {
    std::vector<std::string> assemblyNames1 = stk::io::get_assembly_names(meta1);
    std::vector<std::string> assemblyNames2 = stk::io::get_assembly_names(meta2);

    std::sort(assemblyNames1.begin(), assemblyNames1.end());
    std::sort(assemblyNames2.begin(), assemblyNames2.end());
    EXPECT_EQ(assemblyNames1, assemblyNames2);

    const stk::mesh::BulkData& bulk1 = meta1.mesh_bulk_data();
    const stk::mesh::BulkData& bulk2 = meta2.mesh_bulk_data();

    for(size_t i=0; i<assemblyNames1.size(); ++i) {
      const stk::mesh::Part* assemblyPart1 = meta1.get_part(assemblyNames1[i]);
      const stk::mesh::Part* assemblyPart2 = meta2.get_part(assemblyNames2[i]);
      EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart1));
      EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart2));
      EXPECT_EQ(assemblyPart1->id(), assemblyPart2->id());

      EXPECT_EQ(stk::io::get_sub_assembly_names(meta1, assemblyNames1[i]),
                stk::io::get_sub_assembly_names(meta2, assemblyNames2[i]));

      stk::mesh::PartVector leafParts1 = stk::io::get_unique_leaf_parts(meta1, assemblyNames1[i]);
      stk::mesh::PartVector leafParts2 = stk::io::get_unique_leaf_parts(meta2, assemblyNames2[i]);
      ASSERT_GE(leafParts1.size(), leafParts2.size());
      for(size_t j=0; j<leafParts1.size(); ++j) {
        if (excludedParts.empty() || !find_by_name(excludedParts, leafParts1[j]->name())) {
          EXPECT_TRUE(stk::mesh::find(leafParts2, leafParts1[j]->name()) != nullptr)
              << "Could not find " << leafParts1[j]->name();

          stk::mesh::EntityRank rank1 = leafParts1[j]->primary_entity_rank();
          stk::mesh::EntityRank rank2 = leafParts2[j]->primary_entity_rank();
          EXPECT_EQ(rank1, rank2);

          unsigned numFaces1 = stk::mesh::count_selected_entities(*leafParts1[j], bulk1.buckets(rank1));
          unsigned numFaces2 = stk::mesh::count_selected_entities(*leafParts2[j], bulk2.buckets(rank2));
          EXPECT_EQ(numFaces1, numFaces2);
        }
      }
    }
  }

  void test_write_then_read_block_assemblies(
      size_t expectedNumAssemblies, const stk::mesh::PartVector& partsToExclude = stk::mesh::PartVector())
  {
    const std::string fileName("meshWithBlockAssemblies.e");
    stk::mesh::Selector meshSubsetSelector = create_block_subset_selector(partsToExclude);
    stk::io::write_mesh_subset(fileName, get_bulk(), meshSubsetSelector);

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh(fileName, bulk);

    EXPECT_EQ(expectedNumAssemblies, stk::io::get_assembly_names(meta).size());
    if (expectedNumAssemblies > 0) {
      compare_assemblies(get_meta(), meta, partsToExclude);
    }

    unlink(fileName.c_str());
  }

  void test_write_then_read_non_element_assemblies(size_t expectedNumAssemblies,
      const stk::mesh::PartVector& outputPartsToExclude,
      const stk::mesh::PartVector& inputPartsToExclude,
      const std::string& fileName,
      stk::mesh::EntityRank rank)
  {
    stk::mesh::Selector meshSubsetSelector = !stk::mesh::Selector(stk::mesh::selectUnion(outputPartsToExclude));
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());
    size_t outputFileIndex = stkIo.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
    stkIo.set_output_selector(outputFileIndex, rank, meshSubsetSelector);
    stkIo.write_output_mesh(outputFileIndex);

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh(fileName, bulk);

    EXPECT_EQ(expectedNumAssemblies, stk::io::get_assembly_names(meta).size());
    if (expectedNumAssemblies > 0) {
      compare_assemblies(get_meta(), meta, inputPartsToExclude);
    }

    unlink(fileName.c_str());
  }

  void test_write_then_read_surface_assemblies(size_t expectedNumAssemblies,
      const stk::mesh::PartVector& outputPartsToExclude = stk::mesh::PartVector(),
      const stk::mesh::PartVector& inputPartsToExclude = stk::mesh::PartVector())
  {
    const std::string fileName("meshWithSurfaceAssemblies.e");
    test_write_then_read_non_element_assemblies(
        expectedNumAssemblies, outputPartsToExclude, inputPartsToExclude, fileName, get_meta().side_rank());
  }

  void test_write_then_read_nodeset_assemblies(size_t expectedNumAssemblies,
      const stk::mesh::PartVector& outputPartsToExclude = stk::mesh::PartVector(),
      const stk::mesh::PartVector& inputPartsToExclude = stk::mesh::PartVector())
  {
    const std::string fileName("meshWithNodesetAssemblies.e");
    test_write_then_read_non_element_assemblies(
        expectedNumAssemblies, outputPartsToExclude, inputPartsToExclude, fileName, stk::topology::NODE_RANK);
  }

  void test_for_null_parts(const stk::mesh::PartVector& parts)
  {
    for (const stk::mesh::Part* part : parts) {
      EXPECT_TRUE(part != nullptr);
    }
  }
};

}  // namespace unit_test
}  // namespace io
}  // namespace stk

#endif

