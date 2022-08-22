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

#include <gtest/gtest.h>                // for ASSERT_TRUE, AssertHelper, etc
#include <stk_mesh/base/EntityKey.hpp>  // for EntityKey
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stdexcept>
#include <unordered_set>

namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
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

namespace {

//----------------------------------------------------------------------

TEST(UnitTestEntity, testHashEntity)
{
  std::unordered_set<stk::mesh::Entity, std::hash<stk::mesh::Entity>> set_of_entities;
  stk::mesh::Entity entity1, entity2, entity3, entity4;
  entity1.set_local_offset(1);
  entity2.set_local_offset(2);
  entity3.set_local_offset(3);
  entity4.set_local_offset(4);
  set_of_entities.insert(entity2);
  set_of_entities.insert(entity3);
  set_of_entities.insert(entity1);
  set_of_entities.insert(entity2);//insert entity2 redundantly
  EXPECT_EQ(3u, set_of_entities.size());
  EXPECT_FALSE(set_of_entities.find(entity1) == set_of_entities.end());
  EXPECT_TRUE(set_of_entities.find(entity4) == set_of_entities.end());//entity4 not in set
}

TEST(UnitTestEntity, testHashEntityKey)
{
  std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>> set_of_keys;
  stk::mesh::EntityKey
      key1(stk::topology::NODE_RANK,1),
      key2(stk::topology::EDGE_RANK,1),
      key3(stk::topology::FACE_RANK,1),
      key4(stk::topology::ELEM_RANK,1);
  set_of_keys.insert(key2);
  set_of_keys.insert(key3);
  set_of_keys.insert(key1);
  set_of_keys.insert(key2);//insert key2 redundantly
  EXPECT_EQ(3u, set_of_keys.size());
  EXPECT_FALSE(set_of_keys.find(key1) == set_of_keys.end());
  EXPECT_TRUE(set_of_keys.find(key4) == set_of_keys.end());//key4 not in set
}

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

