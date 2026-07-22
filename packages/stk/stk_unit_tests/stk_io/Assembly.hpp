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
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/AssemblyUtils.hpp>
#include <stk_util/environment/Env.hpp>               // for parallel_size
#include "IOMeshFixture.hpp"
#include <algorithm>
#include <fstream>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include <Ioss_DBUsage.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_IOFactory.h>
#include <Ioss_DatabaseIO.h>

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

namespace stk
{
namespace io
{
namespace unit_test
{

class Assembly : public IOMeshFixture
{
protected:
  stk::mesh::Part& create_assembly(const std::string& assemblyName, int id,
                                   stk::mesh::EntityRank rank = stk::topology::INVALID_RANK)
  {
    stk::mesh::Part& assemblyPart = rank==stk::topology::INVALID_RANK ?
           get_meta().declare_part(assemblyName) : get_meta().declare_part(assemblyName, rank);
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
    EXPECT_FALSE(stk::mesh::is_element_block(part));
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

    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(fileName, *bulk);

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

    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(fileName, *bulk);

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

class AssemblyFilter : public IOMeshFixture
{
protected:

  std::string get_multi_block_mesh_desc(unsigned numBlocks)
  {
    std::ostringstream oss;
    unsigned proc = 0;
    for(unsigned i = 0; i < numBlocks; ++i) {
      unsigned elemId = i + 1;
      unsigned firstNodeId = i * 4 + 1;
      oss << proc << "," << elemId << ",HEX_8,";
      for(unsigned node = firstNodeId; node < firstNodeId + 8; ++node) {
        oss << node << ",";
      }
      unsigned blockId = i + 1;
      oss << "block_" << blockId;

      if(i < numBlocks - 1) {
        oss << "\n";
      }
    }

    return oss.str();
  }

  std::vector<double> get_multi_block_coordinates(unsigned numBlocks)
  {
    std::vector<double> planeCoords = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
    std::vector<double> coordinates;

    coordinates.insert(coordinates.end(), planeCoords.begin(), planeCoords.end());

    for(unsigned i = 1; i <= numBlocks; ++i) {
      for(unsigned point = 0; point < 4; ++point) {
        planeCoords[3 * point + 2] += 1;
      }

      coordinates.insert(coordinates.end(), planeCoords.begin(), planeCoords.end());
    }

    return coordinates;
  }

  void create_multi_block_mesh(const unsigned nBlock, stk::mesh::BulkData& bulk)
  {
    std::string meshDesc = stk::unit_test_util::get_full_text_mesh_desc(get_multi_block_mesh_desc(nBlock),
                                                                        get_multi_block_coordinates(nBlock));
    meshDesc = "textmesh:" + meshDesc;
    stk::io::fill_mesh(meshDesc, bulk);
  }

  stk::mesh::Part& create_assembly_part(stk::mesh::MetaData& meta,
                                        stk::unit_test_util::AssemblyManager& assemblyManager,
                                        stk::unit_test_util::AssemblyDescription& assembly)
  {
    const std::string& assemblyName = assemblyManager.get_assembly_name(assembly.id);
    stk::mesh::Part& assemblyPart = meta.declare_part(assemblyName);
    stk::io::put_assembly_io_part_attribute(assemblyPart);
    meta.set_part_id(assemblyPart, assembly.id);

    return assemblyPart;
  }

  void create_assembly_child_parts(stk::mesh::MetaData& meta, stk::mesh::Part& assemblyPart,
                                   stk::unit_test_util::AssemblyManager& assemblyManager,
                                   stk::unit_test_util::AssemblyDescription& assembly)
  {
    for(const std::string& partName : assembly.partNames) {
      stk::mesh::Part* part = meta.get_part(partName);
      EXPECT_TRUE(nullptr != part);
      meta.declare_part_subset(assemblyPart, *part);
    }

    for(unsigned subAssemblyId : assembly.subAssemblies) {
      const std::string& partName = assemblyManager.get_assembly_name(subAssemblyId);

      stk::mesh::Part* part = meta.get_part(partName);
      EXPECT_TRUE(nullptr != part);
      meta.declare_part_subset(assemblyPart, *part);
    }
  }

  void add_deep_assembly_to_mesh(stk::mesh::MetaData& meta)
  {
    stk::unit_test_util::AssemblyManager assemblyManager;
    unsigned parentAssembly = stk::unit_test_util::create_deep_assembly(assemblyManager);

    std::vector<unsigned> traversalList = assemblyManager.get_assembly_reverse_traversal_list(parentAssembly);

    for(unsigned leafId : traversalList) {
      stk::unit_test_util::AssemblyDescription& assembly = assemblyManager.get_assembly(leafId);
      stk::mesh::Part& assemblyRootPart = create_assembly_part(meta, assemblyManager, assembly);
      create_assembly_child_parts(meta, assemblyRootPart, assemblyManager, assembly);
    }
  }

  void create_multi_block_mesh_with_deep_assembly(const unsigned nBlock, const std::string& fileName)
  {
    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> outputBulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD);
    stk::mesh::MetaData& outputMeta = outputBulk->mesh_meta_data();

    create_multi_block_mesh(nBlock, *outputBulk);
    add_deep_assembly_to_mesh(outputMeta);

    stk::io::write_mesh(fileName, *outputBulk);
  }

  void read_mesh_with_assembly_filter(const std::string& fileName,
                                      const std::vector<std::string>& exclusions,
                                      const std::vector<std::string>& inclusions)
  {
    stk::io::StkMeshIoBroker broker;

    Ioss::PropertyManager properties;
    Ioss::DatabaseIO * db = Ioss::IOFactory::create("exodus", fileName, Ioss::READ_MODEL,
                                                    sierra::Env::parallel_comm(), properties);

    try {
      db->set_assembly_omissions(exclusions, inclusions);
    } catch (...) {
      delete db;
      unlink(fileName.c_str());
      throw;
    }

    Ioss::Region *region = nullptr;

    try {
      region = new Ioss::Region(db, "region");
    } catch (...) {
      unlink(fileName.c_str());
      throw;
    }

    std::shared_ptr<Ioss::Region> io = std::shared_ptr<Ioss::Region>(region);

    size_t db_index = broker.add_mesh_database(io);

    broker.set_active_mesh(db_index);
    broker.set_bulk_data(get_bulk());
    broker.create_input_mesh();
    broker.add_all_mesh_fields_as_input_fields();
    broker.populate_bulk_data();

    unlink(fileName.c_str());
  }

  void test_omit_leaf_assembly()
  {
    std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
    std::sort(assemblyNames.begin(), assemblyNames.end(), std::less<std::string>());

    EXPECT_EQ(5u, assemblyNames.size());
    EXPECT_EQ("assembly_9000", assemblyNames[0]);
    EXPECT_EQ("assembly_9001", assemblyNames[1]);
    EXPECT_EQ("assembly_9002", assemblyNames[2]);
    EXPECT_EQ("assembly_9003", assemblyNames[3]);
    EXPECT_EQ("assembly_9004", assemblyNames[4]);

    stk::mesh::Part* assembly9000 = get_meta().get_part(assemblyNames[0]);
    stk::mesh::Part* assembly9001 = get_meta().get_part(assemblyNames[1]);
    stk::mesh::Part* assembly9002 = get_meta().get_part(assemblyNames[2]);
    stk::mesh::Part* assembly9003 = get_meta().get_part(assemblyNames[3]);
    stk::mesh::Part* assembly9004 = get_meta().get_part(assemblyNames[4]);

    stk::mesh::Part* block_1 = get_meta().get_part("block_1");
    stk::mesh::Part* block_2 = get_meta().get_part("block_2");
    stk::mesh::Part* block_3 = get_meta().get_part("block_3");
    stk::mesh::Part* block_4 = get_meta().get_part("block_4");
    stk::mesh::Part* block_5 = get_meta().get_part("block_5");

    EXPECT_TRUE(nullptr != block_1);
    EXPECT_TRUE(nullptr != block_2);
    EXPECT_TRUE(nullptr != block_3);
    EXPECT_TRUE(nullptr == block_4);
    EXPECT_TRUE(nullptr != block_5);

    std::vector<std::string> subAssemblyNames;

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9000");
    std::sort(subAssemblyNames.begin(), subAssemblyNames.end(), std::less<std::string>());
    EXPECT_EQ(2u, subAssemblyNames.size());
    EXPECT_EQ("assembly_9001", subAssemblyNames[0]);
    EXPECT_EQ("assembly_9002", subAssemblyNames[1]);
    EXPECT_TRUE(assembly9000->contains(*block_1));
    EXPECT_TRUE(assembly9000->contains(*block_2));
    EXPECT_TRUE(assembly9000->contains(*block_3));
    EXPECT_FALSE(assembly9000->contains(*block_5));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9001");
    EXPECT_EQ(1u, subAssemblyNames.size());
    EXPECT_EQ("assembly_9003", subAssemblyNames[0]);
    EXPECT_TRUE(assembly9001->contains(*block_1));
    EXPECT_TRUE(assembly9001->contains(*block_2));
    EXPECT_FALSE(assembly9001->contains(*block_5));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9002");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9002->contains(*block_3));
    EXPECT_FALSE(assembly9002->contains(*block_5));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9003");
    EXPECT_EQ(1u, subAssemblyNames.size());
    EXPECT_EQ("assembly_9004", subAssemblyNames[0]);
    EXPECT_TRUE(assembly9003->contains(*block_1));
    EXPECT_TRUE(assembly9003->contains(*block_2));
    EXPECT_FALSE(assembly9003->contains(*block_5));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9004");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9004->contains(*block_1));
    EXPECT_TRUE(assembly9004->contains(*block_2));
    EXPECT_FALSE(assembly9004->contains(*block_5));
  }

  void test_omit_branch_assembly()
  {
    std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
    std::sort(assemblyNames.begin(), assemblyNames.end(), std::less<std::string>());

    EXPECT_EQ(3u, assemblyNames.size());
    EXPECT_EQ("assembly_9000", assemblyNames[0]);
    EXPECT_EQ("assembly_9001", assemblyNames[1]);
    EXPECT_EQ("assembly_9002", assemblyNames[2]);

    stk::mesh::Part* assembly9000 = get_meta().get_part(assemblyNames[0]);
    stk::mesh::Part* assembly9001 = get_meta().get_part(assemblyNames[1]);
    stk::mesh::Part* assembly9002 = get_meta().get_part(assemblyNames[2]);

    stk::mesh::Part* block_1 = get_meta().get_part("block_1");
    stk::mesh::Part* block_2 = get_meta().get_part("block_2");
    stk::mesh::Part* block_3 = get_meta().get_part("block_3");
    stk::mesh::Part* block_4 = get_meta().get_part("block_4");
    stk::mesh::Part* block_5 = get_meta().get_part("block_5");

    EXPECT_TRUE(nullptr == block_1);
    EXPECT_TRUE(nullptr == block_2);
    EXPECT_TRUE(nullptr != block_3);
    EXPECT_TRUE(nullptr == block_4);
    EXPECT_TRUE(nullptr != block_5);

    std::vector<std::string> subAssemblyNames;

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9000");
    std::sort(subAssemblyNames.begin(), subAssemblyNames.end(), std::less<std::string>());
    EXPECT_EQ(2u, subAssemblyNames.size());
    EXPECT_EQ("assembly_9001", subAssemblyNames[0]);
    EXPECT_EQ("assembly_9002", subAssemblyNames[1]);
    EXPECT_TRUE(assembly9000->contains(*block_3));
    EXPECT_FALSE(assembly9000->contains(*block_5));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9001");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9001->subsets().empty());

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9002");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9002->contains(*block_3));
    EXPECT_FALSE(assembly9002->contains(*block_5));
  }

  void test_include_leaf_assembly()
  {
    std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());

    EXPECT_EQ(1u, assemblyNames.size());
    EXPECT_EQ("assembly_9002", assemblyNames[0]);

    stk::mesh::Part* assembly9002 = get_meta().get_part(assemblyNames[0]);

    stk::mesh::Part* block_1 = get_meta().get_part("block_1");
    stk::mesh::Part* block_2 = get_meta().get_part("block_2");
    stk::mesh::Part* block_3 = get_meta().get_part("block_3");
    stk::mesh::Part* block_4 = get_meta().get_part("block_4");
    stk::mesh::Part* block_5 = get_meta().get_part("block_5");

    EXPECT_TRUE(nullptr == block_1);
    EXPECT_TRUE(nullptr == block_2);
    EXPECT_TRUE(nullptr != block_3);
    EXPECT_TRUE(nullptr == block_4);
    EXPECT_TRUE(nullptr == block_5);

    std::vector<std::string> subAssemblyNames;

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9002");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9002->contains(*block_3));
  }

  void test_include_branch_assembly()
  {
    std::vector<std::string> assemblyNames = stk::io::get_assembly_names(get_meta());
    std::sort(assemblyNames.begin(), assemblyNames.end(), std::less<std::string>());

    EXPECT_EQ(3u, assemblyNames.size());
    EXPECT_EQ("assembly_9003", assemblyNames[0]);
    EXPECT_EQ("assembly_9004", assemblyNames[1]);
    EXPECT_EQ("assembly_9005", assemblyNames[2]);

    stk::mesh::Part* assembly9003 = get_meta().get_part(assemblyNames[0]);
    stk::mesh::Part* assembly9004 = get_meta().get_part(assemblyNames[1]);
    stk::mesh::Part* assembly9005 = get_meta().get_part(assemblyNames[2]);

    stk::mesh::Part* block_1 = get_meta().get_part("block_1");
    stk::mesh::Part* block_2 = get_meta().get_part("block_2");
    stk::mesh::Part* block_3 = get_meta().get_part("block_3");
    stk::mesh::Part* block_4 = get_meta().get_part("block_4");
    stk::mesh::Part* block_5 = get_meta().get_part("block_5");

    EXPECT_TRUE(nullptr != block_1);
    EXPECT_TRUE(nullptr != block_2);
    EXPECT_TRUE(nullptr == block_3);
    EXPECT_TRUE(nullptr != block_4);
    EXPECT_TRUE(nullptr == block_5);

    std::vector<std::string> subAssemblyNames;

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9003");
    std::sort(subAssemblyNames.begin(), subAssemblyNames.end(), std::less<std::string>());
    EXPECT_EQ(2u, subAssemblyNames.size());
    EXPECT_EQ("assembly_9004", subAssemblyNames[0]);
    EXPECT_EQ("assembly_9005", subAssemblyNames[1]);
    EXPECT_TRUE(assembly9003->contains(*block_1));
    EXPECT_TRUE(assembly9003->contains(*block_2));
    EXPECT_TRUE(assembly9003->contains(*block_4));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9004");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9004->contains(*block_1));
    EXPECT_TRUE(assembly9004->contains(*block_2));

    subAssemblyNames = stk::io::get_sub_assembly_names(get_meta(), "assembly_9005");
    EXPECT_TRUE(subAssemblyNames.empty());
    EXPECT_TRUE(assembly9005->contains(*block_4));
  }
};

}  // namespace unit_test
}  // namespace io
}  // namespace stk

#endif

