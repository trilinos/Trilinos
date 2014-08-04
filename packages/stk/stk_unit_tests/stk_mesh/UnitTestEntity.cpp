/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/EntityKey.hpp>  // for EntityKey
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stdexcept>

namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { namespace impl { class EntityRepository; } } }
namespace stk { namespace mesh { namespace impl { class PartRepository; } } }
namespace stk { namespace mesh { struct Entity; } }





using stk::ParallelMachine;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityKey;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::impl::PartRepository;
using stk::mesh::impl::EntityRepository;

namespace {

//----------------------------------------------------------------------

TEST(UnitTestEntity,testEntityKey)
{
  EntityKey key_bad_zero = EntityKey();
  EntityKey key_good_0_1 = EntityKey( stk::topology::NODE_RANK , 1 );
  EntityKey key_good_1_1 = EntityKey( stk::topology::EDGE_RANK , 1 );
  EntityKey key_good_2_10 = EntityKey( stk::topology::FACE_RANK , 10);
  EntityKey key_order_1_12 = EntityKey( stk::topology::EDGE_RANK , 12 );
  EntityKey key_order_2_10 = EntityKey( stk::topology::FACE_RANK , 10 );

  ASSERT_TRUE( ! key_bad_zero.is_valid() );
  ASSERT_TRUE(   key_good_0_1.is_valid() );
  ASSERT_TRUE(   key_good_1_1.is_valid() );
  ASSERT_TRUE(   key_good_2_10.is_valid() );

  ASSERT_TRUE( stk::topology::NODE_RANK  == key_good_0_1.rank());
  ASSERT_TRUE( stk::topology::EDGE_RANK  == key_good_1_1.rank() );
  ASSERT_TRUE( stk::topology::FACE_RANK  == key_good_2_10.rank() );
  ASSERT_TRUE( 1  == key_good_0_1.id() );
  ASSERT_TRUE( 1  == key_good_1_1.id() );
  ASSERT_TRUE( 10 == key_good_2_10.id() );

  ASSERT_TRUE(  key_order_1_12 <  key_order_2_10);
  ASSERT_TRUE( !( key_order_1_12 >  key_order_2_10));

#ifndef NDEBUG
  ASSERT_THROW( EntityKey( stk::topology::INVALID_RANK, 1 ) , std::logic_error );
  ASSERT_THROW( EntityKey( stk::topology::NODE_RANK , ~0ull ) , std::logic_error );
#endif // NDEBUG
}

//----------------------------------------------------------------------
}//namespace <anonymous>

