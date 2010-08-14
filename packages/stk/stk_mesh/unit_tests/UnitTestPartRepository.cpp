/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/parallel/Parallel.hpp>



namespace {

STKUNIT_UNIT_TEST( UnitTestPartRepository, construct )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * universal_part = partRepo.universal_part();
  STKUNIT_EXPECT_TRUE( universal_part != NULL );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, universal_in_subset )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * universal_part = partRepo.universal_part();
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  STKUNIT_ASSERT_THROW(
    partRepo.declare_subset(*part_A,*universal_part),
    std::runtime_error
    );
}


STKUNIT_UNIT_TEST( UnitTestPartRepository, subset_equal_superset )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  STKUNIT_ASSERT_THROW(
    partRepo.declare_subset(*part_A,*part_A),
    std::runtime_error
    );
}


STKUNIT_UNIT_TEST( UnitTestPartRepository, universal_in_intersection )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * universal_part = partRepo.universal_part();
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  std::vector<stk::mesh::Part *> intersection;
  intersection.push_back(universal_part);
  STKUNIT_ASSERT_THROW(
    partRepo.declare_part(intersection),
    std::runtime_error
    );
  intersection.push_back(part_A);
  STKUNIT_ASSERT_THROW(
    partRepo.declare_part(intersection),
    std::runtime_error
    );
}


STKUNIT_UNIT_TEST( UnitTestPartRepository, circular_subset )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  stk::mesh::Part * part_B = partRepo.declare_part("B",0);
  stk::mesh::Part * part_C = partRepo.declare_part("C",0);
  partRepo.declare_subset(*part_A,*part_B);
  partRepo.declare_subset(*part_B,*part_C);
  STKUNIT_ASSERT_THROW(
    partRepo.declare_subset(*part_C,*part_A),
    std::runtime_error
    );
}


STKUNIT_UNIT_TEST( UnitTestPartRepository, inconsistent_rank_subset )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  stk::mesh::Part * part_B = partRepo.declare_part("B",1);
  // lower rank cannot contain higher rank:
  STKUNIT_ASSERT_THROW(
    partRepo.declare_subset(*part_A,*part_B),
    std::runtime_error
    );
  // higher rank can contain lower rank:
  STKUNIT_ASSERT_NO_THROW(
    partRepo.declare_subset(*part_B,*part_A)
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, two_part_repositories )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo_1(&meta);
  stk::mesh::Part * part_1_A = partRepo_1.declare_part("A",0);
  stk::mesh::Part * part_1_B = partRepo_1.declare_part("B",0);
  stk::mesh::impl::PartRepository partRepo_2(&meta);
  stk::mesh::Part * part_2_A = partRepo_2.declare_part("A",0);
  // subset/superset parts must come from same part repository
  STKUNIT_ASSERT_THROW(
    partRepo_1.declare_subset(*part_1_A,*part_2_A),
    std::runtime_error
    );

  // intersection contains parts from another part repository
  {
    std::vector<stk::mesh::Part *> intersection;
    intersection.push_back(part_1_A);
    intersection.push_back(part_1_B);
    STKUNIT_ASSERT_THROW(
        partRepo_2.declare_part( intersection ),
        std::runtime_error
        );
  }
  // intersection contains parts from multiple part repositories
  {
    std::vector<stk::mesh::Part *> intersection;
    intersection.push_back(part_1_A);
    intersection.push_back(part_2_A);
    STKUNIT_ASSERT_THROW(
        partRepo_1.declare_part( intersection ),
        std::runtime_error
        );
  }
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, invalid_relation )
{
  std::vector<std::string> entity_rank_names;
  entity_rank_names.push_back("node");
  entity_rank_names.push_back("element");
  stk::mesh::MetaData meta(entity_rank_names);
  stk::mesh::impl::PartRepository partRepo(&meta);
  stk::mesh::Part * part_A = partRepo.declare_part("A",0);
  stk::mesh::Part * part_B = partRepo.declare_part("B",0);
  stk::mesh::PartRelation relation;
  relation.m_root = part_B;
  relation.m_target = part_A;
  STKUNIT_ASSERT_THROW(
      partRepo.declare_part_relation( *part_A, relation, *part_B ),
      std::runtime_error
      );
}
} // namespace 
