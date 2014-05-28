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

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::impl::PartRepository;

class UnitTestPartRepository
{
public:
  UnitTestPartRepository();
  ~UnitTestPartRepository() {}

   const int spatial_dimension;
   MetaData meta;
   stk::mesh::impl::PartRepository partRepo;
   stk::mesh::impl::PartRepository partRepo_1;
   stk::mesh::impl::PartRepository partRepo_2;

   stk::mesh::Part * universal_part;
   stk::mesh::Part * part_A;
   stk::mesh::Part * part_B;
   stk::mesh::Part * part_C;
   stk::mesh::Part * part_D;
   stk::mesh::Part * part_1_A;
   stk::mesh::Part * part_1_B;
   stk::mesh::Part * part_2_A;
   std::vector<stk::mesh::Part *> intersection;
   stk::mesh::PartRelation relation;
   const CellTopologyData * singleton;
};

UnitTestPartRepository::UnitTestPartRepository()
  : spatial_dimension(3)
  , meta( stk::mesh::fem::entity_rank_names(spatial_dimension) )
  , partRepo( &meta )
  , partRepo_1( &meta )
  , partRepo_2( &meta )
  , universal_part(      partRepo.universal_part()    )
  , part_A   (           partRepo.declare_part("A",0) )
  , part_B   (           partRepo.declare_part("B",0) )
  , part_C   (           partRepo.declare_part("C",0) )
  , part_D   (           partRepo.declare_part("D",1) )
  , part_1_A (           partRepo_1.declare_part("A",0) )
  , part_1_B (           partRepo_1.declare_part("B",0) )
  , part_2_A (           partRepo_2.declare_part("A",0) )
  , intersection ( )
  , relation ( )
  , singleton ( NULL )
{
 meta.commit();
}

namespace {

STKUNIT_UNIT_TEST( UnitTestPartRepository, construct )
{
  UnitTestPartRepository upr;
  STKUNIT_EXPECT_TRUE( upr.universal_part != NULL );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, universal_in_subset )
{
  UnitTestPartRepository upr;
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_subset(*upr.part_A,*upr.universal_part),
    std::runtime_error
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, subset_equal_superset )
{
  UnitTestPartRepository upr;
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_subset(*upr.part_A,*upr.part_A),
    std::runtime_error
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, universal_in_intersection )
{
  UnitTestPartRepository upr;
  upr.intersection.push_back(upr.universal_part);
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_part(upr.intersection),
    std::runtime_error
    );
  upr.intersection.push_back(upr.part_A);
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_part(upr.intersection),
    std::runtime_error
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, circular_subset )
{
  UnitTestPartRepository upr;
  upr.partRepo.declare_subset(*upr.part_A,*upr.part_B);
  upr.partRepo.declare_subset(*upr.part_B,*upr.part_C);
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_subset(*upr.part_C,*upr.part_A),
    std::runtime_error
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, inconsistent_rank_subset )
{
  UnitTestPartRepository upr;
  // lower rank cannot contain higher rank:
  STKUNIT_ASSERT_THROW(
    upr.partRepo.declare_subset(*upr.part_A,*upr.part_D),
    std::runtime_error
    );
  // higher rank can contain lower rank:
  STKUNIT_ASSERT_NO_THROW(
    upr.partRepo.declare_subset(*upr.part_D,*upr.part_A)
    );
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, two_part_repositories )
{
  UnitTestPartRepository upr;
  // subset/superset parts must come from same part repository
  STKUNIT_ASSERT_THROW(
    upr.partRepo_1.declare_subset(*upr.part_1_A,*upr.part_2_A),
    std::runtime_error
    );

  // intersection contains parts from another part repository
  {
    upr.intersection.push_back(upr.part_1_A);
    upr.intersection.push_back(upr.part_1_B);
    STKUNIT_ASSERT_THROW(
        upr.partRepo_2.declare_part( upr.intersection ),
        std::runtime_error
        );
  }
  // intersection contains parts from multiple part repositories
  {
    upr.intersection.push_back(upr.part_1_A);
    upr.intersection.push_back(upr.part_2_A);
    STKUNIT_ASSERT_THROW(
        upr.partRepo_1.declare_part( upr.intersection ),
        std::runtime_error
        );
  }
}

STKUNIT_UNIT_TEST( UnitTestPartRepository, invalid_relation )
{
  UnitTestPartRepository upr;
  upr.relation.m_root = upr.part_B;
  upr.relation.m_target = upr.part_A;
  STKUNIT_ASSERT_THROW(
      upr.partRepo.declare_part_relation( *upr.part_A, upr.relation, *upr.part_B ),
      std::runtime_error
      );
}

//Test covers declare_attribute_no_delete in PartRepository.hpp and PartImpl.hpp
STKUNIT_UNIT_TEST( UnitTestPartRepository, declare_attribute_no_delete )
{
  UnitTestPartRepository upr;
  upr.partRepo.declare_attribute_no_delete(*upr.part_A, upr.singleton);

}

}
