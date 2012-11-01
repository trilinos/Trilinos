

#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

// Borrow a lot from UnitTestSelector.  Bulk up the SelectorFixture to have parts
// with enough entities so that each partition (bucket family) comprises multiple
// buckets.

namespace {

using stk::mesh::fixtures::SelectorFixture ;

void initialize(SelectorFixture& fixture,
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

  const size_t bucket_size = 1000;  // Default value for BucketRepository constructor.
  const size_t lb_num_buckets_per_partition = 3;
  stk::mesh::EntityRank ent_type = 0; // rank

  const size_t bf_size = bucket_size * lb_num_buckets_per_partition;
  stk::mesh::EntityId ent_id = 1000;
  std::vector<stk::mesh::Part*> partMembership;

  // Entities in collection 1 are contained in PartA
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec1_arg.push_back(ent);
      ++ent_id;
  }

  // Entity2 is contained in PartA and PartB
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  partMembership.push_back( & fixture.m_partB );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec2_arg.push_back(ent);
      ++ent_id;
  }

  // Entity3 is contained in PartB and PartC
  partMembership.clear();
  partMembership.push_back( & fixture.m_partB );
  partMembership.push_back( & fixture.m_partC );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec3_arg.push_back(ent);
      ++ent_id;
  }

  // Entity4 is contained in PartC
  partMembership.clear();
  partMembership.push_back( & fixture.m_partC );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec4_arg.push_back(ent);
      ++ent_id;
  }

  // Entity5 is not contained in any Part
  partMembership.clear();
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec5_arg.push_back(ent);
      ++ent_id;
  }

  STKUNIT_ASSERT(fixture.m_bulk_data.modification_end());
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

/** \brief Verify we can construct the unit testing fixture.
 *
 * */
STKUNIT_UNIT_TEST( UnitTestPartition, UTBF_testInitialize )
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

    initialize(fix, ec1, ec2, ec3, ec4, ec5);

    stk::mesh::Selector selector;
  //std::cout << "Selector = " << selector << std::endl;

    {
        selector = fix.m_partA & !fix.m_partB;
        for (size_t i = 0; i < ec1.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec1[i]->bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = fix.m_partA & fix.m_partB;
        for (size_t i = 0; i < ec2.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec2[i]->bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = fix.m_partB & fix.m_partC;
        for (size_t i = 0; i < ec3.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec3[i]->bucket() ;
            bool result = selector(bucket);
          STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = !fix.m_partB & fix.m_partC;
        for (size_t i = 0; i < ec4.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec4[i]->bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
      selector = !(fix.m_partA | fix.m_partB | fix.m_partC | fix.m_partD);
      for (size_t i = 0; i < ec5.size(); ++i)
      {
          const stk::mesh::Bucket & bucket = ec5[i]->bucket() ;
          bool result = selector(bucket);
          STKUNIT_EXPECT_TRUE(result);
      }
    }
}


STKUNIT_UNIT_TEST( UnitTestPartition, UTBF_testForReal)
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
    initialize(fix, ec1, ec2, ec3, ec4, ec5);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.update_partitions();

    // Correct means:
    //   - A partition for each bucket
    //   - Each partition contains the right number of buckets (didn't lose any buckets!)

    stk::mesh::Selector selector_A_only = fix.m_partA & !fix.m_partB;
    stk::mesh::Selector selector_A_and_B = fix.m_partA & fix.m_partB;
    stk::mesh::Selector selector_B_and_C = fix.m_partB & fix.m_partC;
    stk::mesh::Selector selector_C_only = !fix.m_partB & fix.m_partC;
    stk::mesh::Selector selector_no_letter = !(fix.m_partA | fix.m_partB | fix.m_partC | fix.m_partD);

    std::vector<stk::mesh::Selector> selectors;
    selectors.push_back(selector_A_only);
    selectors.push_back(selector_A_and_B);
    selectors.push_back(selector_B_and_C);
    selectors.push_back(selector_C_only);
    selectors.push_back(selector_no_letter);

    std::vector<bool> found_it(5, false);

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        // size_t count = 0;
        stk::mesh::impl::Partition &partition = *partitions[i];
        const std::vector<unsigned> &partition_key = partition.get_legacy_partition_id();
        for (std::vector<stk::mesh::Bucket *>::iterator bkt= partition.begin(); bkt != partition.end(); ++bkt)
        {
            STKUNIT_EXPECT_EQ(partitions[i], (*bkt)->getPartition() );
            const unsigned *bucket_key = (*bkt)->key();
            for (size_t k = 0; k < partition_key.size(); ++k)
            {
                STKUNIT_EXPECT_EQ(partition_key[k], bucket_key[k]);
            }
        }
    }
}



/** \} */


} // namespace
