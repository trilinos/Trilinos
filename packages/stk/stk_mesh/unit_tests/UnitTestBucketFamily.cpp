

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

// Borrow a lot from UnitTestBucketFamily.  Bulk up the SelectorFixture to have parts
// with enough entities so that each partition (bucket family) comprises multiple
// buckets.

namespace {

using stk::mesh::fixtures::SelectorFixture ;

void initialize(SelectorFixture& fixture,
                std::vector<stk::mesh::Entity *> &ec1_arg,
                std::vector<stk::mesh::Entity *> &ec2_arg,
                std::vector<stk::mesh::Entity *> &ec3_arg,
                std::vector<stk::mesh::Entity *> &ec4_arg,
                std::vector<stk::mesh::Entity *> &ec5_arg
                )
{
  fixture.m_meta_data.commit();
  fixture.m_bulk_data.modification_begin();
  fixture.generate_mesh();

  const size_t bucket_size = 1000;  // Default value for BucketRepository constructor.
  const size_t lb_num_buckets_per_family = 3;
  stk::mesh::EntityRank ent_type = 0; // rank

  const size_t bf_size = bucket_size * lb_num_buckets_per_family;
  stk::mesh::EntityId ent_id = 1000;
  std::vector<stk::mesh::Part*> partMembership;

  // Entities in collection 1 are contained in PartA
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity *ent = & fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec1_arg.push_back(ent);
      ++ent_id;
  }

  // Entity2 is contained in PartA and PartB
  partMembership.clear();
  partMembership.push_back( & fixture.m_partA );
  partMembership.push_back( & fixture.m_partB );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity *ent = & fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec2_arg.push_back(ent);
      ++ent_id;
  }

  // Entity3 is contained in PartB and PartC
  partMembership.clear();
  partMembership.push_back( & fixture.m_partB );
  partMembership.push_back( & fixture.m_partC );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity *ent = & fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec3_arg.push_back(ent);
      ++ent_id;
  }

  // Entity4 is contained in PartC
  partMembership.clear();
  partMembership.push_back( & fixture.m_partC );
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity *ent = & fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec4_arg.push_back(ent);
      ++ent_id;
  }

  // Entity5 is not contained in any Part
  partMembership.clear();
  for (size_t i = 1; i < bf_size; ++i)
  {
      stk::mesh::Entity *ent = & fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
      ec5_arg.push_back(ent);
      ++ent_id;
  }

  STKUNIT_ASSERT(fixture.m_bulk_data.modification_end());
}

/** \defgroup stk_mesh_bucket_family_unit "stk::mesh::BucketFamily Unit Testing"
  * \addtogroup stk_mesh_bucket_family_unit
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

/** \brief Verify we can construct the selector unit testing fixture.
 *
 * */
STKUNIT_UNIT_TEST( UnitTestBucketFamily, one_FatSelectorFixture )
{
    std::vector<stk::mesh::Entity *> ec1;
    std::vector<stk::mesh::Entity *> ec2;
    std::vector<stk::mesh::Entity *> ec3;
    std::vector<stk::mesh::Entity *> ec4;
    std::vector<stk::mesh::Entity *> ec5;

    SelectorFixture fix;
    initialize(fix, ec1, ec2, ec3, ec4, ec5);

    STKUNIT_EXPECT_TRUE(true);
}



/** \brief Test containment directly.
 *
 * Verify PartA contains Entity1 and Entity2, and does not contain
 * Entity3, Entity4, or Entity5.
 * */
STKUNIT_UNIT_TEST( UnitTestBucketFamily, UTBF_A_12345 )
{
    std::vector<stk::mesh::Entity *> ec1;
    std::vector<stk::mesh::Entity *> ec2;
    std::vector<stk::mesh::Entity *> ec3;
    std::vector<stk::mesh::Entity *> ec4;
    std::vector<stk::mesh::Entity *> ec5;

    SelectorFixture fix;
    initialize(fix, ec1, ec2, ec3, ec4, ec5);

  stk::mesh::Selector selector( fix.m_partA );
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket() ;
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}

STKUNIT_UNIT_TEST( UnitTestBucketFamily, UTBF_B_12345 )
{
    std::vector<stk::mesh::Entity *> ec1;
    std::vector<stk::mesh::Entity *> ec2;
    std::vector<stk::mesh::Entity *> ec3;
    std::vector<stk::mesh::Entity *> ec4;
    std::vector<stk::mesh::Entity *> ec5;

    SelectorFixture fix;
    initialize(fix, ec1, ec2, ec3, ec4, ec5);

  stk::mesh::Selector selector( fix.m_partA );
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket() ;
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}



/** \} */


} // namespace
