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

#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for rand, srand, RAND_MAX

namespace {

class CustomBulkData : public stk::mesh::BulkData
{
protected:
  friend class CustomMeshBuilder;

  CustomBulkData(std::shared_ptr<stk::mesh::MetaData> metaData,
                 stk::ParallelMachine parallel,
                 enum AutomaticAuraOption autoAuraOption = AUTO_AURA,
                 std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager = std::unique_ptr<stk::mesh::FieldDataManager>(),
                 unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                 unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity(),
                 std::shared_ptr<stk::mesh::impl::AuraGhosting> auraGhosting = std::shared_ptr<stk::mesh::impl::AuraGhosting>(),
                 bool createUpwardConnectivity = true)
#ifdef SIERRA_MIGRATION
    : BulkData(metaData, parallel, autoAuraOption, false, std::move(fieldDataManager), initialBucketCapacity,
               maximumBucketCapacity, auraGhosting, createUpwardConnectivity)
#else
    : BulkData(metaData, parallel, autoAuraOption, std::move(fieldDataManager), initialBucketCapacity,
               maximumBucketCapacity, auraGhosting, createUpwardConnectivity)
#endif
  {
  }
};

class CustomMeshBuilder : public stk::mesh::MeshBuilder
{
public:
  CustomMeshBuilder() = default;
  CustomMeshBuilder(stk::ParallelMachine comm)
    : stk::mesh::MeshBuilder(comm)
  {}

  virtual ~CustomMeshBuilder() override = default;

  //using statement to avoid compile-warning about 'only partially overridden'
  using stk::mesh::MeshBuilder::create;

  virtual std::unique_ptr<stk::mesh::BulkData> create(std::shared_ptr<stk::mesh::MetaData> metaData) override
  {
    STK_ThrowRequireMsg(m_haveComm, "MeshBuilder must be given an MPI communicator before creating BulkData");

    return std::unique_ptr<CustomBulkData>(new CustomBulkData(metaData,
                                                              m_comm,
                                                              m_auraOption,
                                                              std::move(m_fieldDataManager),
                                                              m_initialBucketCapacity,
                                                              m_maximumBucketCapacity));
  }

  virtual std::shared_ptr<stk::mesh::MetaData> create_meta_data() override
  {
    std::shared_ptr<stk::mesh::MetaData> meta;
    if (m_spatialDimension > 0 || !m_entityRankNames.empty()) {
      meta = std::make_shared<stk::mesh::MetaData>(m_spatialDimension, m_entityRankNames);
    }
    else {
      meta = std::make_shared<stk::mesh::MetaData>();
    }

    meta->declare_part("CUSTOM_PART");

    return meta;
  }
};


TEST(StkMeshIoBroker, useExternalBulkData)
{
  stk::mesh::MeshBuilder meshBuilder(MPI_COMM_WORLD);
  std::shared_ptr<stk::mesh::BulkData> bulk = meshBuilder.set_spatial_dimension(3).create();

  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  stkMeshIoBroker.set_bulk_data(bulk);
  stkMeshIoBroker.add_mesh_database("generated:1x1x8", stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  EXPECT_EQ(stkMeshIoBroker.bulk_data_ptr().get(), bulk.get());
}

TEST(StkMeshIoBroker, useDefaultMeshBuilder)
{
  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  stkMeshIoBroker.add_mesh_database("generated:1x1x8", stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  ASSERT_TRUE(stkMeshIoBroker.meta_data_ptr());
  ASSERT_TRUE(stkMeshIoBroker.bulk_data_ptr());
}

TEST(StkMeshIoBroker, useCustomMeshBuilder)
{
  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  stkMeshIoBroker.set_mesh_builder(std::make_shared<CustomMeshBuilder>(MPI_COMM_WORLD));
  stkMeshIoBroker.add_mesh_database("generated:1x1x8", stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  ASSERT_TRUE(stkMeshIoBroker.meta_data_ptr());
  ASSERT_TRUE(stkMeshIoBroker.bulk_data_ptr());

  EXPECT_TRUE(stkMeshIoBroker.meta_data().get_part("CUSTOM_PART") != nullptr);
  EXPECT_TRUE(dynamic_cast<CustomBulkData*>(stkMeshIoBroker.bulk_data_ptr().get()) != nullptr);
}

TEST(StkMeshIoBroker, useCustomMeshBuilder_afterInternalAlreadyGenerated)
{
  {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database("generated:1x1x8", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    EXPECT_ANY_THROW(stkMeshIoBroker.set_mesh_builder(std::make_shared<CustomMeshBuilder>(MPI_COMM_WORLD)));
  }
  {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database("generated:1x1x8", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    EXPECT_ANY_THROW(stkMeshIoBroker.set_mesh_builder(std::make_shared<CustomMeshBuilder>(MPI_COMM_WORLD)));
  }
}

TEST(StkMeshIoBroker, useCustomMeshBuilder_afterSetExternal)
{
  stk::mesh::MeshBuilder meshBuilder(MPI_COMM_WORLD);
  std::shared_ptr<stk::mesh::BulkData> bulk = meshBuilder.set_spatial_dimension(3).create();

  stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
  stkMeshIoBroker.set_bulk_data(bulk);
  EXPECT_ANY_THROW(stkMeshIoBroker.set_mesh_builder(std::make_shared<CustomMeshBuilder>(MPI_COMM_WORLD)));
}

}
