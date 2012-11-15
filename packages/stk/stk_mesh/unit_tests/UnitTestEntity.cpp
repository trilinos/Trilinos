/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>
#include <map>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/baseImpl/EntityRepository.hpp>

using stk::ParallelMachine;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::EntityKey;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::impl::PartRepository;
using stk::mesh::impl::EntityRepository;

namespace {

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestEntity,testEntityKey)
{
  EntityKey key_bad_zero = EntityKey();
  EntityKey key_good_0_1 = EntityKey( 0 , 1 );
  EntityKey key_good_1_1 = EntityKey( 1 , 1 );
  EntityKey key_good_2_10 = EntityKey( 2 , 10);
  EntityKey key_order_1_12 = EntityKey( 1 , 12 );
  EntityKey key_order_2_10 = EntityKey( 2 , 10 );

  STKUNIT_ASSERT( ! entity_key_valid(  key_bad_zero ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_0_1 ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_1_1 ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_2_10 ) );

  STKUNIT_ASSERT( 0  == entity_rank( key_good_0_1));
  STKUNIT_ASSERT( 1  == entity_rank( key_good_1_1) );
  STKUNIT_ASSERT( 2  == entity_rank( key_good_2_10) );
  STKUNIT_ASSERT( 1  == entity_id( key_good_0_1) );
  STKUNIT_ASSERT( 1  == entity_id( key_good_1_1) );
  STKUNIT_ASSERT( 10 == entity_id( key_good_2_10) );

  STKUNIT_ASSERT(  key_order_1_12 <  key_order_2_10);
  STKUNIT_ASSERT( !( key_order_1_12 >  key_order_2_10));

#ifndef NDEBUG
  STKUNIT_ASSERT_THROW( EntityKey( ~0u , 1 ) , std::logic_error );
  STKUNIT_ASSERT_THROW( EntityKey( 0 , ~stk::mesh::EntityKey::raw_key_type(0) ) , std::logic_error );
#endif // NDEBUG
}

//----------------------------------------------------------------------
}//namespace <anonymous>

