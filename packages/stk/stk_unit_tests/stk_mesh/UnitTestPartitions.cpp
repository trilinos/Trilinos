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
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for copy, reverse
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, Bucket::iterator
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator&, etc
#include <stk_mesh/base/Types.hpp>      // for PartOrdinal, BucketVector, etc
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator, etc

#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository, etc
#include "stk_mesh/baseImpl/GlobalIdEntitySorter.hpp"  // for BucketRepository, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/SelectorFixture.hpp"  // for SelectorFixture

namespace {

struct ReverseSorter : public stk::mesh::EntitySorterBase
{
  virtual void sort(stk::mesh::BulkData &bulk, stk::mesh::EntityVector& entityVector) const
  {
    std::sort(entityVector.begin(), entityVector.end(),
              [&bulk](const stk::mesh::Entity lhs, const stk::mesh::Entity rhs)->bool
    { return bulk.entity_key(lhs) > bulk.entity_key(rhs); });
  }
};

using stk::mesh::fixtures::SelectorFixture;


// Borrow a lot from UnitTestSelector.  Bulk up the SelectorFixture to have parts
// with enough entities so that each partition (bucket family) comprises multiple
// buckets.

stk::mesh::EntityId
addEntitiesToFixture(SelectorFixture& fixture, stk::mesh::EntityId start_id, size_t num_to_add,
                     const std::vector<stk::mesh::Part*> &partMembership,
                     std::vector<stk::mesh::Entity> &collector)
{
  stk::mesh::EntityId ent_id = start_id;
  for (size_t i = 0; i < num_to_add; ++i)
  {
    stk::mesh::Entity ent = fixture.m_bulk_data.declare_node(ent_id, partMembership);
    collector.push_back(ent);
    ++ent_id;
  }
  return ent_id;
}

void initializeFiveEntityCollections(SelectorFixture& fixture,
                                     std::vector<stk::mesh::Entity> &ec1_arg,
                                     std::vector<stk::mesh::Entity> &ec2_arg,
                                     std::vector<stk::mesh::Entity> &ec3_arg,
                                     std::vector<stk::mesh::Entity> &ec4_arg,
                                     std::vector<stk::mesh::Entity> &ec5_arg
                                     )
{
  fixture.m_meta_data.commit();
  fixture.m_bulk_data.modification_begin();
  fixture.generate_mesh();

  const size_t bucket_size = 1000;     // Default value for BucketRepository constructor.
  const size_t lb_num_buckets_per_partition = 3;

  const size_t bf_size = bucket_size * lb_num_buckets_per_partition;
  stk::mesh::EntityId ent_id = 1001;   // Want to keep numerical alignment.
  std::vector<stk::mesh::Part*> partMembership;

  // Note that the loop variables start at 1 because SelectorFixture::generate_mesh() has
  // already created an Entity in each partition.

  // Entities in collection 1 are contained in PartA
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  ent_id = addEntitiesToFixture(fixture, ent_id, bf_size - 1, partMembership, ec1_arg);

  // Entities in collection 2 are in PartA and PartB
  ++ent_id;     // For numerical alignment.
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  partMembership.push_back( & fixture.m_partB );
  ent_id = addEntitiesToFixture(fixture, ent_id, bf_size - 1, partMembership, ec2_arg);

  // Entities in collection 3 are in PartB and PartC
  ++ent_id;
  partMembership.clear();
  partMembership.push_back( & fixture.m_partB );
  partMembership.push_back( & fixture.m_partC );
  ent_id = addEntitiesToFixture(fixture, ent_id, bf_size - 1, partMembership, ec3_arg);

  // Entities in collection 4 are in PartC
  ++ent_id;
  partMembership.clear();
  partMembership.push_back( & fixture.m_partC );
  ent_id = addEntitiesToFixture(fixture, ent_id, bf_size - 1, partMembership, ec4_arg);

  // Entities in collection 5 are not contained in any Part
  ++ent_id;
  partMembership.clear();
  ent_id = addEntitiesToFixture(fixture, ent_id, bf_size - 1, partMembership, ec5_arg);

  bool me_result = fixture.m_bulk_data.modification_end();
  ASSERT_TRUE(me_result);
}

void initializeFivePartitionsWithSixBucketsEach(SelectorFixture& fixture)
{
  std::vector<stk::mesh::Entity> ec1;
  std::vector<stk::mesh::Entity> ec2;
  std::vector<stk::mesh::Entity> ec3;
  std::vector<stk::mesh::Entity> ec4;
  std::vector<stk::mesh::Entity> ec5;

  initializeFiveEntityCollections(fixture, ec1, ec2, ec3, ec4, ec5);
}

void initialize_unsorted(SelectorFixture& fixture)
{
  initializeFivePartitionsWithSixBucketsEach(fixture);

  stk::mesh::impl::BucketRepository &bucket_repository = fixture.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    partition.sort(ReverseSorter());
  }
}

// Initialize field data in the fixture to correspond to the EntityIds of the entities.
void setFieldDataUsingEntityIDs(SelectorFixture& fix)
{
  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  stk::mesh::BulkData &mesh = fix.m_bulk_data;

  std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();

  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    for (stk::mesh::BucketVector::iterator bkt_i= partition.begin(); bkt_i != partition.end(); ++bkt_i)
    {
      stk::mesh::Bucket &bkt = **bkt_i;
      double *field_data = stk::mesh::field_data(*fix.m_fieldABC, bkt, 0);
      if (!field_data)
      {
        continue;
      }

      size_t bkt_size = bkt.size();
      for (size_t k = 0; k < bkt_size; ++k)
      {
        stk::mesh::Entity curr_ent = bkt[k];
        field_data[k] = mesh.identifier(curr_ent);
      }
    }
  }
}


bool areAllEntitiesSelected(const stk::mesh::BulkData& mesh, const stk::mesh::Selector selects,
                            const std::vector<stk::mesh::Entity> &cands)
{
  for (size_t i = 0; i < cands.size(); ++i)
  {
    if (!selects(mesh.bucket(cands[i])))
      return false;
  }
  return true;
}

bool check_bucket_ptrs(const stk::mesh::Bucket &bucket)
{
  if (bucket.size() == 0 )
  {
    std::cout << "Bucket has size zero!" << std::endl;
    return false;
  }

  stk::mesh::Bucket::iterator e_i = bucket.begin();
  stk::mesh::Bucket::iterator e_e = bucket.end();
  for(; e_i != e_e; ++e_i)
  {
    if (&bucket != bucket.mesh().bucket_ptr(*e_i))
      return false;
  }
  return true;
}

template <typename Data_T>
bool check_nonempty_strictly_ordered(Data_T data[], size_t sz, bool reject_0_lt_0 = true )
{
  if (sz == 0)
    return false;

  for (size_t i = 0; i < sz - 1; ++i)
  {
    if ((data[i] >= data[i + 1]) && reject_0_lt_0)
    {
      std::cout << "i = " << i << ": data[i] = " << data[i]
                   << ", data[i + 1] = " << data[i + 1] << std::endl;
      return false;
    }
  }
  return true;
}

void check_test_partition_invariant(SelectorFixture& fix,
                                    const stk::mesh::impl::Partition &partition)
{
  const std::vector<unsigned> &partition_key = partition.get_legacy_partition_id();
  for (stk::mesh::BucketVector::const_iterator bkt_i= partition.begin();
       bkt_i != partition.end(); ++bkt_i)
  {
    const stk::mesh::Bucket &bkt = **bkt_i;
    EXPECT_EQ(&partition, bkt.getPartition() );
    EXPECT_TRUE(check_bucket_ptrs(bkt));

    double *field_data = stk::mesh::field_data(*fix.m_fieldABC, bkt, 0);
    if (field_data)
    {
      EXPECT_TRUE(check_nonempty_strictly_ordered(field_data, bkt.size()));
    }
    const unsigned *bucket_key = bkt.key();
    for (size_t k = 0; k < partition_key.size() - 1; ++k)
    {
      EXPECT_EQ(partition_key[k], bucket_key[k]);
    }

    stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
    EXPECT_EQ(bucket_repository.get_bucket(bkt.entity_rank(), bkt.bucket_id()), &bkt);
  }
}

unsigned checkGetBucketAndCountNonEmpty(stk::mesh::impl::BucketRepository &bucket_repo,
                                        stk::mesh::EntityRank entity_rank, unsigned max_id_to_try)
{
  unsigned count = 0;
  for (unsigned i = 0; i < max_id_to_try; ++i)
  {
    stk::mesh::Bucket *b = bucket_repo.get_bucket(entity_rank, i);
    if (b)
    {
      ++count;
      EXPECT_EQ(i, b->bucket_id());
    }
  }
  return count;
}

void check_bucket_ids_testset_A(stk::mesh::impl::BucketRepository &bucket_repository)
{
  const unsigned num_node_buckets = bucket_repository.buckets(stk::topology::NODE_RANK).size();
  unsigned numNonEmptyNodeBucketsInRepository = checkGetBucketAndCountNonEmpty(bucket_repository, stk::topology::NODE_RANK, num_node_buckets);
  EXPECT_EQ(numNonEmptyNodeBucketsInRepository, num_node_buckets);
  const unsigned num_face_buckets = bucket_repository.buckets(stk::topology::FACE_RANK).size();
  unsigned numNonEmptyFaceBucketsInRepository = checkGetBucketAndCountNonEmpty(bucket_repository, stk::topology::FACE_RANK, num_face_buckets);
  EXPECT_EQ(numNonEmptyFaceBucketsInRepository, num_face_buckets);
}


/** \defgroup stk_mesh_partition_unit "stk::mesh::Partition Unit Testing"
 * \addtogroup stk_mesh_partition_unit
 * \{
 *
 * Selector unit testing environment. <br>
 * A special set of mesh parts and entities are set up in the
 * following configuration for the Selector unit testing.<br>
 * Parts:  PartA, PartB, PartC, PartD, PartU <br>
 * PartU = MetaData.universal_part() <br>
 * Entities:  Entity1, Entity2, Entity3, Entity4, Entity5 <br>
 *
 * PartA contains Entity1, Entity2 <br>
 * PartB contains Entity2, Entity3 <br>
 * PartC contains Entity3, Entity4 <br>
 * PartD contains no entities <br>
 * Entity5 is not contained in any Part <br>
 *
 * <PRE>
 * |----------|--|-------|--|----------|    |-------------|
 * |<--PartA---->|       |<--PartC---->|    |   PartD     |
 * |          |<---PartB--->|          |    |             |
 * |  1       |2 |       |3 |       4  | 5  |             |
 * |          |  |       |  |          |    |             |
 * |          |  |       |  |          |    |             |
 * |----------|--|-------|--|----------|    |-------------|
 * </PRE>
 *
 * Note:  The unit test names use the convention of "i" for
 * intersection, "u" for union, and "c" for complement.
 *
 * */


/// Verify we can construct the unit testing fixture.
TEST( UnitTestPartition, Partition_testInitialize )
{
  std::vector<stk::mesh::Entity> ec1;
  std::vector<stk::mesh::Entity> ec2;
  std::vector<stk::mesh::Entity> ec3;
  std::vector<stk::mesh::Entity> ec4;
  std::vector<stk::mesh::Entity> ec5;

  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }

  initializeFiveEntityCollections(fix, ec1, ec2, ec3, ec4, ec5);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  for (size_t i = 0; i < num_partitions; ++i)
  {
    const stk::mesh::impl::Partition &partition = *partitions[i];
    check_test_partition_invariant(fix, partition);
  }

  stk::mesh::Selector selector;

  selector = fix.m_partA & !fix.m_partB;
  EXPECT_TRUE(areAllEntitiesSelected(fix.m_bulk_data, selector, ec1));

  selector = fix.m_partA & fix.m_partB;
  EXPECT_TRUE(areAllEntitiesSelected(fix.m_bulk_data, selector, ec2));

  selector = fix.m_partB & fix.m_partC;
  EXPECT_TRUE(areAllEntitiesSelected(fix.m_bulk_data, selector, ec3));

  selector = (!fix.m_partB) & fix.m_partC;
  EXPECT_TRUE(areAllEntitiesSelected(fix.m_bulk_data, selector, ec4));

  selector = !(fix.m_partA | fix.m_partB | fix.m_partC | fix.m_partD);
  EXPECT_TRUE(areAllEntitiesSelected(fix.m_bulk_data, selector, ec5));

  check_bucket_ids_testset_A(bucket_repository);
}

/// Test Partition::sort()
TEST( UnitTestPartition, Partition_testSort)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initialize_unsorted(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    // check_test_partition_invariant(fix, partition);
    partition.sort(stk::mesh::impl::GlobalIdEntitySorter());
    check_test_partition_invariant(fix, partition);
  }

  check_bucket_ids_testset_A(bucket_repository);
}

/// Test Partition::remove(.)
TEST( UnitTestPartition, Partition_testRemove)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initializeFivePartitionsWithSixBucketsEach(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];

    size_t old_size = partition.size();
    size_t num_removed = 0;

    // Remove non-last entity in a bucket.
    stk::mesh::Bucket &bkt_0 = **partition.begin();
    stk::mesh::Entity entity = bkt_0[0];
    partition.remove(entity);
    //    fix.m_bulk_data.my_set_mesh_index(entity, 0, 0);
    ++num_removed;

    // Remove last entity in a bucket.
    stk::mesh::Bucket &bkt_1 = **partition.begin();
    stk::mesh::Entity e_last_in_1 = bkt_1[bkt_1.size() - 1];
    partition.remove(e_last_in_1);
    //    fix.m_bulk_data.my_set_mesh_index(e_last_in_1, 0, 0);
    ++num_removed;

    // Need to sort before checking whether the invariant holds.
    partition.default_sort_if_needed();

    EXPECT_EQ(old_size,  partition.size() + num_removed);
    check_test_partition_invariant(fix, partition);
  }

  check_bucket_ids_testset_A(bucket_repository);

  //
  // Now's a good time to exercise BucketRepository::sync_from_partitions().
  //
  bucket_repository.sync_from_partitions();
  partitions = bucket_repository.get_partitions(stk::topology::NODE_RANK);
  num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    check_test_partition_invariant(fix, partition);
  }

  check_bucket_ids_testset_A(bucket_repository);
}

/// Test Partition::add(.).
TEST( UnitTestPartition, Partition_testAdd)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initializeFivePartitionsWithSixBucketsEach(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_partitions = partitions.size();
  EXPECT_EQ(num_partitions, 5u);

  std::vector<stk::mesh::Entity> first_entities(num_partitions);
  std::vector<size_t> old_sizes(num_partitions);

  // Get the first entity in each partition.
  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    stk::mesh::Bucket &bkt = **partition.begin();
    first_entities[i] = bkt[0];
    old_sizes[i] = partition.size();
  }

  // Now remove them from those partitions.
  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::Entity entity = first_entities[i];
    partitions[i]->remove(entity);
    fix.m_bulk_data.my_set_mesh_index(entity, 0, 0);
  }

  const stk::mesh::BulkData& mesh = fix.m_bulk_data;

  // "Rotate" each former first entity to another partition.
  for (size_t i = 0; i < num_partitions; ++i)
  {
    size_t dst_partition_idx = (i + 1) % num_partitions;
    stk::mesh::impl::Partition &dst_partition = *partitions[dst_partition_idx];
    stk::mesh::Entity entity = first_entities[i];
    dst_partition.add(entity);
    stk::mesh::Bucket &bkt = mesh.bucket(entity);
    double *field_data = stk::mesh::field_data(*fix.m_fieldABC, bkt, 0);
    if (field_data)
    {
      field_data[mesh.bucket_ordinal(entity)] = mesh.identifier(entity);
    }
  }

  // Because the each first entity had a low identifier, it will again be
  // the first entity in its new partition.
  for (size_t i = 0; i < num_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *partitions[i];
    partition.default_sort_if_needed();
    EXPECT_EQ(old_sizes[i], partition.size());
    check_test_partition_invariant(fix, partition);
  }

  check_bucket_ids_testset_A(bucket_repository);
}

/// Test Partition::move_to(..)
TEST( UnitTestPartition, Partition_testMoveTo)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initializeFivePartitionsWithSixBucketsEach(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> all_partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  const size_t num_all_partitions = all_partitions.size();
  EXPECT_EQ(num_all_partitions, 5u);

  const size_t num_data_partitions = num_all_partitions - 1;
  std::vector<stk::mesh::impl::Partition *> data_partitions(num_data_partitions);
  std::copy(all_partitions.begin() + 1, all_partitions.end(), &data_partitions[0]);

  std::vector<stk::mesh::Entity> first_entities(num_data_partitions);
  std::vector<size_t> old_sizes(num_data_partitions);

  // Get the first entity in each partition.
  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *data_partitions[i];
    stk::mesh::Bucket &bkt = **partition.begin();
    first_entities[i] = bkt[0];
    old_sizes[i] = partition.size();
  }

  // "Rotate" each first entity to another partition.
  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    size_t dst_partition_idx = (i + 1) % num_data_partitions;
    stk::mesh::impl::Partition &dst_partition = *data_partitions[dst_partition_idx];
    data_partitions[i]->move_to(first_entities[i], dst_partition);
  }

  // Because the each first entity had a low identifier, it will again be
  // the first entity in its new partition.
  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *data_partitions[i];
    partition.default_sort_if_needed();

    EXPECT_EQ(old_sizes[i], partition.size());
    check_test_partition_invariant(fix, partition);
  }

  check_bucket_ids_testset_A(bucket_repository);
}

// Test the OrdinalVector version of get_or_create_partition(..).
TEST( UnitTestPartition, Partition_testGetOrCreateOV)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initializeFivePartitionsWithSixBucketsEach(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  stk::mesh::OrdinalVector scratch;
  std::vector<stk::mesh::PartOrdinal> parts;
  parts.push_back(fix.m_meta_data.universal_part().mesh_meta_data_ordinal());
  parts.push_back(fix.m_meta_data.locally_owned_part().mesh_meta_data_ordinal());
  parts.push_back(fix.m_partA.mesh_meta_data_ordinal() );
  stk::mesh::impl::Partition *partitionA =
      bucket_repository.get_or_create_partition(stk::topology::NODE_RANK, parts);
  ASSERT_TRUE(0 != partitionA);
  size_t numEntitiesPerPartition = 3000;
  size_t bucketCapacity = stk::mesh::get_default_bucket_capacity();
  size_t expectedNumBucketsPerPartition = (numEntitiesPerPartition + (bucketCapacity - 1u)) / bucketCapacity;
  EXPECT_EQ(expectedNumBucketsPerPartition, partitionA->num_buckets());

  parts.push_back(fix.m_partC.mesh_meta_data_ordinal());
  stk::mesh::impl::Partition *partitionAC =
      bucket_repository.get_or_create_partition(stk::topology::NODE_RANK, parts);
  ASSERT_TRUE(0 != partitionAC);
  EXPECT_EQ(0u, partitionAC->num_buckets());

  stk::mesh::impl::Partition *partitionAC_again =
      bucket_repository.get_or_create_partition(stk::topology::NODE_RANK, parts);
  ASSERT_TRUE(partitionAC == partitionAC_again);

  check_bucket_ids_testset_A(bucket_repository);
}

/** \} */

} // namespace

/// Test Partition::move_to(..) more rigorously
TEST( UnitTestPartition, Partition_testMoveToBetter)
{
  SelectorFixture fix;

  if (fix.m_bulk_data.parallel_size() > 1)
  {
    return;
  }
  initializeFivePartitionsWithSixBucketsEach(fix);
  setFieldDataUsingEntityIDs(fix);

  stk::mesh::impl::BucketRepository &bucket_repository = fix.m_bulk_data.my_get_bucket_repository();
  bucket_repository.sync_from_partitions();

  std::vector<stk::mesh::impl::Partition *> all_partitions =
      bucket_repository.get_partitions(stk::topology::NODE_RANK);
  size_t num_all_partitions = all_partitions.size();
  EXPECT_EQ(num_all_partitions, 5u);

  const size_t num_data_partitions = num_all_partitions - 1;
  std::vector<stk::mesh::impl::Partition *> data_partitions(num_data_partitions);
  std::copy(all_partitions.begin() + 1, all_partitions.end(), &data_partitions[0]);

  std::vector<std::vector<stk::mesh::Entity> > entities_to_move(num_data_partitions);
  std::vector<size_t> old_sizes(num_data_partitions);

  // Get the entities from the first bucket in each partition.
  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *data_partitions[i];
    stk::mesh::Bucket &bkt = **partition.begin();

    size_t bkt_sz = bkt.size();
    const size_t default_bucket_capacity = stk::mesh::get_default_bucket_capacity();
    EXPECT_EQ(bkt_sz, default_bucket_capacity);
    for (size_t j = 0; j < bkt_sz; ++j)
    {
      entities_to_move[i].push_back(bkt[j]);
    }
    old_sizes[i] = partition.size();
    EXPECT_EQ(old_sizes[i], 3000u);
  }

  // "Rotate" the entities to another partition.
  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    size_t dst_partition_idx = (i + 1) % num_data_partitions;
    stk::mesh::impl::Partition &dst_partition = *data_partitions[dst_partition_idx];

    size_t num_to_move = entities_to_move[i].size();
    // std::cout << "Moving entities from " << *data_partitions[i] << " to " << dst_partition
    //           << "\n starting with " << entities_to_move[i][0] << " and ending with "
    //           << entities_to_move[i][num_to_move - 1] << std::endl;
    for (size_t j = 0; j < num_to_move; ++j)
    {
      data_partitions[i]->move_to(entities_to_move[i][j], dst_partition);
    }
  }

  for (size_t i = 0; i < num_data_partitions; ++i)
  {
    stk::mesh::impl::Partition &partition = *data_partitions[i];
    partition.default_sort_if_needed();

    EXPECT_EQ(partition.size(), old_sizes[i]);

    // std::cout << "Check " << partition << std::endl;
    // std::cout << "source partition was " << src_partition << std::endl;

    EXPECT_EQ(old_sizes[i], partition.size());
    check_test_partition_invariant(fix, partition);
  }
}
