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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <gtest/gtest.h>

void verify_comm(stk::ParallelMachine expectedComm, const stk::mesh::BulkData& bulk)
{
  int commCompareResult = 0;
  MPI_Comm_compare(expectedComm, bulk.parallel(), &commCompareResult);
  EXPECT_EQ(MPI_IDENT, commCompareResult);

  EXPECT_EQ(stk::parallel_machine_size(expectedComm), bulk.parallel_size());
  EXPECT_EQ(stk::parallel_machine_rank(expectedComm), bulk.parallel_rank());
}

TEST(MeshBuilder, create_meta_no_comm)
{
  std::shared_ptr<stk::mesh::MetaData> meta = stk::mesh::MeshBuilder().create_meta_data();
  EXPECT_TRUE(nullptr != meta);
}

TEST(MeshBuilder, create_bulkdata_no_comm_throws)
{
  EXPECT_ANY_THROW(stk::mesh::MeshBuilder().create());
}

TEST(MeshBuilder, construct_builder_then_set_comm)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::mesh::MeshBuilder builder;
  builder.set_communicator(comm);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_TRUE(nullptr != bulk);
  verify_comm(comm, *bulk);
}

TEST(MeshBuilder, create_simplest_comm_world)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();

  EXPECT_TRUE(nullptr != bulk);
  verify_comm(comm, *bulk);
}

TEST(MeshBuilder, create_simplest_comm_self)
{
  stk::ParallelMachine comm = MPI_COMM_SELF;
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();

  EXPECT_TRUE(nullptr != bulk);
  verify_comm(comm, *bulk);
}

TEST(MeshBuilder, bulkdata_and_metadata_outlive_builder)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();

  EXPECT_TRUE(nullptr != bulk);
  verify_comm(comm, *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::BulkData& metaBulk = meta.mesh_bulk_data();
  EXPECT_EQ(&metaBulk, bulk.get());
}

TEST(MeshBuilder, bulkdata_aura_default)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  EXPECT_TRUE(bulk->is_automatic_aura_on());
}

TEST(MeshBuilder, bulkdata_aura_on)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_TRUE(bulk->is_automatic_aura_on());
}

TEST(MeshBuilder, bulkdata_aura_off)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_FALSE(bulk->is_automatic_aura_on());
}

TEST(MeshBuilder, set_spatial_dimension)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned expectedSpatialDim = 2;
  builder.set_spatial_dimension(expectedSpatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(expectedSpatialDim, bulk->mesh_meta_data().spatial_dimension());
}

TEST(MeshBuilder, spatial_dimension_default_then_initialize)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  EXPECT_EQ(0u, bulk->mesh_meta_data().spatial_dimension());
  const unsigned expectedSpatialDim = 2;
  bulk->mesh_meta_data().initialize(expectedSpatialDim);

  EXPECT_EQ(expectedSpatialDim, bulk->mesh_meta_data().spatial_dimension());
}

TEST(MeshBuilder, set_entity_rank_names_without_spatial_dimension)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_entity_rank_names({"node","edge","face","elem","constraint"});
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_spatial_dimension_zero_and_empty_entity_rank_names)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned expectedSpatialDim = 0;
  builder.set_spatial_dimension(expectedSpatialDim);
  std::vector<std::string> expectedRankNames = {};
  builder.set_entity_rank_names(expectedRankNames);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(expectedSpatialDim, bulk->mesh_meta_data().spatial_dimension());
  EXPECT_EQ(expectedRankNames, bulk->mesh_meta_data().entity_rank_names());
}

TEST(MeshBuilder, set_spatial_dimension_and_entity_rank_names)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned expectedSpatialDim = 3;
  builder.set_spatial_dimension(expectedSpatialDim);
  std::vector<std::string> expectedRankNames = {"node","edge","face","elem","constraint"};
  builder.set_entity_rank_names(expectedRankNames);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(expectedSpatialDim, bulk->mesh_meta_data().spatial_dimension());
  EXPECT_EQ(expectedRankNames, bulk->mesh_meta_data().entity_rank_names());
}

TEST(MeshBuilder, bulkdata_add_fmwk_data_default)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  EXPECT_FALSE(bulk->add_fmwk_data());
}

TEST(MeshBuilder, bulkdata_add_fmwk_data)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_add_fmwk_data(true);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

#ifdef SIERRA_MIGRATION
  EXPECT_TRUE(bulk->add_fmwk_data());
#else
  EXPECT_FALSE(bulk->add_fmwk_data());
#endif
}

TEST(MeshBuilder, bulkdata_add_fmwk_data_false)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_add_fmwk_data(false);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_FALSE(bulk->add_fmwk_data());
}

TEST(MeshBuilder, default_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), stk::mesh::get_default_initial_bucket_capacity());
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), stk::mesh::get_default_maximum_bucket_capacity());
}

TEST(MeshBuilder, set_invalid_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_bucket_capacity(0);
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned setCapacity = 256;
  builder.set_bucket_capacity(setCapacity);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), setCapacity);
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), setCapacity);
}

TEST(MeshBuilder, set_invalid_initial_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_initial_bucket_capacity(0);
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_initial_bucket_capacity_bigger_than_maximum)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned initialCapacity = stk::mesh::get_default_maximum_bucket_capacity() * 2;
  builder.set_initial_bucket_capacity(initialCapacity);
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_initial_bucket_capacity_smaller_than_maximum)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned initialCapacity = stk::mesh::get_default_maximum_bucket_capacity() / 2;
  builder.set_initial_bucket_capacity(initialCapacity);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), initialCapacity);
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), stk::mesh::get_default_maximum_bucket_capacity());
}

TEST(MeshBuilder, set_invalid_maximum_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_maximum_bucket_capacity(0);
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_maximum_bucket_capacity_bigger_than_initial)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned maxCapacity = stk::mesh::get_default_initial_bucket_capacity() * 2;
  builder.set_maximum_bucket_capacity(maxCapacity);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), stk::mesh::get_default_initial_bucket_capacity());
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), maxCapacity);
}

TEST(MeshBuilder, set_maximum_bucket_capacity_smaller_than_initial)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned maxCapacity = stk::mesh::get_default_initial_bucket_capacity() / 2;
  builder.set_maximum_bucket_capacity(maxCapacity);
  EXPECT_ANY_THROW(builder.create());
}

TEST(MeshBuilder, set_identical_initial_and_maximum_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned initialCapacity = 256;
  const unsigned maximumCapacity = 256;
  builder.set_initial_bucket_capacity(initialCapacity);
  builder.set_maximum_bucket_capacity(maximumCapacity);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), initialCapacity);
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), maximumCapacity);
}

TEST(MeshBuilder, set_different_initial_and_maximum_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned initialCapacity = 32;
  const unsigned maximumCapacity = 256;
  builder.set_initial_bucket_capacity(initialCapacity);
  builder.set_maximum_bucket_capacity(maximumCapacity);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  EXPECT_EQ(bulk->get_initial_bucket_capacity(), initialCapacity);
  EXPECT_EQ(bulk->get_maximum_bucket_capacity(), maximumCapacity);
}

TEST(MeshBuilder, set_faulty_initial_and_maximum_bucket_capacity)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  const unsigned initialCapacity = 64;
  const unsigned maximumCapacity = 32;
  builder.set_initial_bucket_capacity(initialCapacity);
  builder.set_maximum_bucket_capacity(maximumCapacity);
  EXPECT_ANY_THROW(builder.create());
}
