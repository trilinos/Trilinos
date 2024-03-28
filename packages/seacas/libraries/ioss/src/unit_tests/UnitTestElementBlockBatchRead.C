// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>

#include <unistd.h> // for unlink

#include "Ionit_Initializer.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_Field.h"      // for Field, etc
#include "Ioss_IOFactory.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"

#include "Ioss_ElementBlock.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Region.h"

#include "Ioss_Utils.h"

namespace {

  std::string get_randomized_many_block_mesh_desc(unsigned numBlocks)
  {
    std::ostringstream    oss;
    std::vector<unsigned> elementIds(numBlocks);
    std::iota(elementIds.begin(), elementIds.end(), 1);

    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(elementIds.begin(), elementIds.end(), g);

    unsigned proc = 0;
    for (unsigned i = 0; i < numBlocks; ++i) {
      unsigned elemId      = elementIds[i];
      unsigned firstNodeId = i * 4 + 1;
      oss << proc << "," << elemId << ",HEX_8,";
      for (unsigned node = firstNodeId; node < firstNodeId + 8; ++node) {
        oss << node << ",";
      }
      unsigned blockId = i + 1;
      oss << "block_" << blockId;

      if (i < numBlocks - 1) {
        oss << "\n";
      }
    }

    oss << "|coordinates:";

    std::vector<double> planeCoords = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};

    for (double coord : planeCoords) {
      oss << coord << ",";
    }

    for (unsigned i = 1; i <= numBlocks; ++i) {
      for (unsigned point = 0; point < 4; ++point) {
        planeCoords[3 * point + 2] += 1;
      }

      for (double coord : planeCoords) {
        oss << coord << ",";
      }
    }

    return oss.str();
  }

  void define_model(const Ioss::Region &i_region, Ioss::Region &o_region)
  {
    Ioss::DatabaseIO *o_database = o_region.get_database();

    o_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    Ioss::NodeBlock *i_nb        = i_region.get_node_blocks()[0];
    int64_t          spatial_dim = 3;
    int64_t          num_nodes   = i_nb->entity_count();
    Ioss::NodeBlock *o_nb = new Ioss::NodeBlock(o_database, "nodeblock_1", num_nodes, spatial_dim);
    o_region.add(o_nb);

    for (Ioss::ElementBlock *i_eb : i_region.get_element_blocks()) {
      Ioss::ElementBlock *o_eb = new Ioss::ElementBlock(
          o_database, i_eb->name(), i_eb->topology()->name(), i_eb->entity_count());
      o_eb->property_add(i_eb->get_property("id"));
      o_region.add(o_eb);
    }

    o_region.end_mode(Ioss::STATE_DEFINE_MODEL);
  }

  void write_model(const Ioss::Region &i_region, Ioss::Region &o_region)
  {
    Ioss::NodeBlock *i_nb = i_region.get_node_blocks()[0];
    Ioss::NodeBlock *o_nb = o_region.get_node_blocks()[0];

    o_region.begin_mode(Ioss::STATE_MODEL);
    std::vector<double> coordinates;
    std::vector<int>    node_ids;
    i_nb->get_field_data("ids", node_ids);
    i_nb->get_field_data("mesh_model_coordinates", coordinates);

    o_nb->put_field_data("ids", node_ids);
    o_nb->put_field_data("mesh_model_coordinates", coordinates);

    for (Ioss::ElementBlock *i_eb : i_region.get_element_blocks()) {
      Ioss::ElementBlock *o_eb = o_region.get_element_block(i_eb->name());
      std::vector<int>    elem_ids;
      std::vector<int>    connectivity;

      i_eb->get_field_data("ids", elem_ids);
      i_eb->get_field_data("connectivity", connectivity);

      o_eb->put_field_data("ids", elem_ids);
      o_eb->put_field_data("connectivity", connectivity);
    }

    o_region.end_mode(Ioss::STATE_MODEL);
  }

  void define_transient(const Ioss::Region &i_region, Ioss::Region &o_region,
                        const std::string &elemFieldName)
  {
    o_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t      num_elem = o_eb->get_property("entity_count").get_int();
      std::string storage  = "scalar";

      Ioss::Field field(elemFieldName, Ioss::Field::REAL, storage, 1, Ioss::Field::Field::TRANSIENT,
                        num_elem);
      o_eb->field_add(field);
    }
    o_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  void write_transient(Ioss::Region &o_region, const std::string &elemFieldName)
  {
    o_region.begin_mode(Ioss::STATE_TRANSIENT);
    int step = o_region.add_state(0.0);
    o_region.begin_state(step);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t num_elem = o_eb->get_property("entity_count").get_int();

      std::vector<double> field_data(num_elem);
      std::vector<int>    elem_ids;

      o_eb->get_field_data("ids", elem_ids);
      for (size_t i = 0; i < elem_ids.size(); i++) {
        field_data[i] = (double)elem_ids[i];
      }

      o_eb->put_field_data(elemFieldName, field_data);
    }

    o_region.end_state(step);
    o_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  void generate_randomized_many_block_mesh_from_textmesh(int                numBlocks,
                                                         const std::string &elemFieldName,
                                                         const std::string &outFile)
  {
    Ioss::Init::Initializer io;

    Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());

    if (util.parallel_rank() == 0) {
      std::string meshDesc = get_randomized_many_block_mesh_desc(numBlocks);

      Ioss::PropertyManager propertyManager;

      Ioss::DatabaseIO *i_database =
          Ioss::IOFactory::create("textmesh", meshDesc, Ioss::READ_MODEL,
                                  Ioss::ParallelUtils::comm_self(), propertyManager);
      Ioss::Region i_region(i_database, "input_model");
      EXPECT_TRUE(i_database != nullptr);
      EXPECT_TRUE(i_database->ok(true));

      Ioss::DatabaseIO *o_database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_self(), propertyManager);
      Ioss::Region o_region(o_database, "output_model");
      EXPECT_TRUE(o_database != nullptr);
      EXPECT_TRUE(o_database->ok(true));

      define_model(i_region, o_region);
      write_model(i_region, o_region);

      define_transient(i_region, o_region, elemFieldName);
      write_transient(o_region, elemFieldName);

      o_database->finalize_database();
      o_database->flush_database();
      o_database->closeDatabase();
    }
  }

  template <typename INT>
  std::pair<double, double>
  do_connectivity_timing_impl(const Ioss::Region                &region,
                              const Ioss::ElementBlockContainer &elemBlocks)
  {
    std::vector<INT> blockBatchConnectivity;

    double              startTime = Ioss::Utils::timer();
    std::vector<size_t> dataOffset =
        region.get_entity_field_data("connectivity", elemBlocks, blockBatchConnectivity);
    double elapsedBlockBatchConnectivityTime = Ioss::Utils::timer() - startTime;

    double elapsedConnectivityTime = 0.0;
    for (unsigned i = 0; i < elemBlocks.size(); i++) {
      const Ioss::ElementBlock *entity = elemBlocks[i];
      int64_t                   iblk   = int64_t(i);

      std::vector<INT> connectivity;
      startTime        = Ioss::Utils::timer();
      int64_t numToGet = entity->get_field_data("connectivity", connectivity);
      elapsedConnectivityTime += Ioss::Utils::timer() - startTime;

      unsigned numComponents = entity->topology()->number_nodes();
      int64_t  numEntities   = entity->entity_count();
      EXPECT_EQ(numToGet, numEntities);

      int64_t blockBatchNumToGet = dataOffset[iblk + 1] - dataOffset[iblk];
      EXPECT_EQ(blockBatchNumToGet, numEntities * numComponents);

      for (unsigned eIndex = 0; eIndex < numEntities; eIndex++) {
        for (unsigned comp = 0; comp < numComponents; comp++) {
          size_t connIndex           = eIndex * numComponents + comp;
          size_t blockBatchConnIndex = dataOffset[iblk] + connIndex;

          INT connValue           = connectivity[connIndex];
          INT blockBatchConnValue = blockBatchConnectivity[blockBatchConnIndex];

          EXPECT_EQ(connValue, blockBatchConnValue);
        }
      }
    }

    return std::make_pair(elapsedConnectivityTime, elapsedBlockBatchConnectivityTime);
  }

  std::pair<double, double> do_connectivity_timing(const Ioss::Region                &region,
                                                   const Ioss::ElementBlockContainer &elemBlocks)
  {
    auto result = std::pair<double, double>(0.0, 0.0);

    bool is64Bit = (region.get_database()->int_byte_size_api() == 8);
    if (is64Bit) {
      result = do_connectivity_timing_impl<int64_t>(region, elemBlocks);
    }
    else {
      result = do_connectivity_timing_impl<int>(region, elemBlocks);
    }

    return result;
  }

  std::pair<double, double>
  do_field_timing_and_verification(const Ioss::Region                &region,
                                   const Ioss::ElementBlockContainer &elemBlocks,
                                   const std::string                 &fieldName)
  {
    std::vector<double> blockFieldData;

    double              startBlockFieldTime = Ioss::Utils::timer();
    std::vector<size_t> dataOffset =
        region.get_entity_field_data(fieldName, elemBlocks, blockFieldData);
    double elapsedBlockFieldTime = Ioss::Utils::timer() - startBlockFieldTime;

    double elapsedFieldTime = 0.0;

    std::vector<double> fieldData;

    for (unsigned i = 0; i < elemBlocks.size(); i++) {
      const Ioss::ElementBlock *entity = elemBlocks[i];

      int64_t iblk = int64_t(i);

      if (entity->field_exists(fieldName)) {
        double  startFieldTime = Ioss::Utils::timer();
        int64_t numToGet       = entity->get_field_data(fieldName, fieldData);
        elapsedFieldTime += Ioss::Utils::timer() - startFieldTime;

        Ioss::Field iossField        = entity->get_field(fieldName);
        unsigned    numComponents    = iossField.raw_storage()->component_count();
        int64_t     numEntities      = entity->entity_count();
        int64_t     expectedNumToGet = numEntities * numComponents;
        EXPECT_EQ(numToGet, expectedNumToGet);

        int64_t blockNumToGet = dataOffset[iblk + 1] - dataOffset[iblk];
        EXPECT_EQ(blockNumToGet, expectedNumToGet);

        for (unsigned eIndex = 0; eIndex < numEntities; eIndex++) {
          for (unsigned comp = 0; comp < numComponents; comp++) {
            size_t fieldIndex      = eIndex * numComponents + comp;
            size_t blockFieldIndex = dataOffset[iblk] + fieldIndex;

            double fieldValue      = fieldData[fieldIndex];
            double blockFieldValue = blockFieldData[blockFieldIndex];

            EXPECT_NEAR(fieldValue, blockFieldValue, 1.0e-5);
          }
        }
      }
    }

    return std::make_pair(elapsedFieldTime, elapsedBlockFieldTime);
  }

  Ioss::ElementBlockContainer get_all_blocks(const Ioss::Region &region)
  {
    return region.get_element_blocks();
  }

  Ioss::ElementBlockContainer get_alternate_blocks(const Ioss::Region &region)
  {
    const Ioss::ElementBlockContainer &elemBlocks = region.get_element_blocks();

    Ioss::ElementBlockContainer alternateElemBlocks;
    alternateElemBlocks.reserve(elemBlocks.size() / 2);

    for (unsigned i = 0; i < elemBlocks.size(); i++) {
      if (i % 2 == 0) {
        alternateElemBlocks.push_back(elemBlocks[i]);
      }
    }

    return alternateElemBlocks;
  }

  using ElemBlockFunc = std::function<Ioss::ElementBlockContainer(const Ioss::Region &region)>;

  void run_block_batch_test(int numBlocks, ElemBlockFunc func)
  {
    std::ostringstream oss;
    oss << "randomizedManyBlocks_";
    oss << numBlocks << ".g";

    std::string outFile(oss.str());
    std::string elemFieldName = "elem_field";
    generate_randomized_many_block_mesh_from_textmesh(numBlocks, elemFieldName, outFile);

    {
      Ioss_MPI_Comm       comm = Ioss::ParallelUtils::comm_world();
      Ioss::ParallelUtils util(comm);
      util.barrier();

      std::string           decompMethod = "RCB";
      Ioss::Property        decompProp("DECOMPOSITION_METHOD", decompMethod);
      Ioss::PropertyManager propertyManager;
      propertyManager.add(decompProp);

      Ioss::DatabaseIO *database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::READ_MODEL, comm, propertyManager);
      Ioss::Region region(database, "input_model");
      EXPECT_TRUE(database != nullptr);
      EXPECT_TRUE(database->ok(true));

      region.property_add(decompProp);
      region.begin_state(1);

      Ioss::ElementBlockContainer elemBlocks = func(region);

      double elapsedConnectivityTime;
      double elapsedBlockBatchConnectivityTime;
      std::tie(elapsedConnectivityTime, elapsedBlockBatchConnectivityTime) =
          do_connectivity_timing(region, elemBlocks);

      double maxBlockBatchConnectivityDuration =
          util.global_minmax(elapsedBlockBatchConnectivityTime, Ioss::ParallelUtils::DO_MAX);
      double minBlockBatchConnectivityDuration =
          util.global_minmax(elapsedBlockBatchConnectivityTime, Ioss::ParallelUtils::DO_MIN);

      double maxConnectivityDuration =
          util.global_minmax(elapsedConnectivityTime, Ioss::ParallelUtils::DO_MAX);
      double minConnectivityDuration =
          util.global_minmax(elapsedConnectivityTime, Ioss::ParallelUtils::DO_MIN);

      if (util.parallel_rank() == 0) {
        std::cout << std::endl;
        std::cout << "MAX       connectivity read time = " << std::setw(8) << std::fixed
                  << maxConnectivityDuration
                  << " : MIN       connectivity read time = " << std::setw(8) << std::fixed
                  << minConnectivityDuration << std::endl;
        std::cout << "MAX BATCH connectivity read time = " << std::setw(8) << std::fixed
                  << maxBlockBatchConnectivityDuration
                  << " : MIN BATCH connectivity read time = " << std::setw(8) << std::fixed
                  << minBlockBatchConnectivityDuration << std::endl;
        std::cout << std::endl;
      }

      double elapsedFieldTime;
      double elapsedBlockBatchFieldTime;

      std::tie(elapsedFieldTime, elapsedBlockBatchFieldTime) =
          do_field_timing_and_verification(region, elemBlocks, elemFieldName);

      double maxBlockBatchFieldDuration =
          util.global_minmax(elapsedBlockBatchFieldTime, Ioss::ParallelUtils::DO_MAX);
      double minBlockBatchFieldDuration =
          util.global_minmax(elapsedBlockBatchFieldTime, Ioss::ParallelUtils::DO_MIN);

      double maxFieldDuration = util.global_minmax(elapsedFieldTime, Ioss::ParallelUtils::DO_MAX);
      double minFieldDuration = util.global_minmax(elapsedFieldTime, Ioss::ParallelUtils::DO_MIN);

      if (util.parallel_rank() == 0) {
        std::cout << std::endl;
        std::cout << "MAX       field read time = " << std::setw(8) << std::fixed
                  << maxFieldDuration << " : MIN       field read time = " << std::setw(8)
                  << std::fixed << minFieldDuration << std::endl;
        std::cout << "MAX BATCH field read time = " << std::setw(8) << std::fixed
                  << maxBlockBatchFieldDuration << " : MIN BATCH field read time = " << std::setw(8)
                  << std::fixed << minBlockBatchFieldDuration << std::endl;
        std::cout << std::endl;
      }
    }

    unlink(outFile.c_str());
  }

  TEST(TestReadBlockBatch, readAllBlocks)
  {
    int numBlocks = 1000;
    run_block_batch_test(numBlocks, get_all_blocks);
  }

  TEST(TestReadBlockBatch, readAlternateBlocks)
  {
    int numBlocks = 1000;
    run_block_batch_test(numBlocks, get_alternate_blocks);
  }

} // namespace
