// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_MESH_ENTITYKEY_HPP
#define STK_MESH_ENTITYKEY_HPP

#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for uint64_t
#include <boost/static_assert.hpp>      // for BOOST_STATIC_ASSERT
#include <boost/type_traits/is_same.hpp>  // for is_same
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityRank
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowAssertMsg
#include "boost/functional/hash/hash.hpp"  // for hash_value
#include "boost/mpl/bool.hpp"           // for bool_<>::value





namespace stk {
namespace mesh {

struct EntityKey
{
  BOOST_STATIC_ASSERT(( boost::is_same<EntityId, uint64_t>::value ));

  enum entity_key_t {
      RANK_SHIFT = 56ULL
    , MIN_ID = 0ULL
    , MAX_ID = (1ULL << RANK_SHIFT) - 1ULL
    , ID_MASK = MAX_ID
    , INVALID = ~0ULL
  };

  static bool is_valid_id( EntityId id )
  {
    return id > MIN_ID && id <=  EntityKey::MAX_ID;
  }

  EntityKey()
    : m_value(INVALID)
  {}

  EntityKey( entity_key_t value )
    : m_value(value)
  {}

  EntityKey( EntityRank arg_rank, EntityId arg_id )
    : m_value( static_cast<entity_key_t>( static_cast<uint64_t>(arg_rank) << RANK_SHIFT | arg_id) )
  {
    ThrowAssertMsg( arg_rank <= static_cast<EntityRank>(255), "Error: given an out of range entity rank " << arg_rank);
    ThrowAssertMsg( arg_id <= MAX_ID, "Error: given an out of range entity id " << arg_id);
  }

  EntityId   id() const   { return m_value & ID_MASK; }
  EntityRank rank() const { return static_cast<EntityRank>(m_value >> RANK_SHIFT); }

  bool is_valid() const { return m_value != INVALID; }

  operator entity_key_t() const { return m_value; }

  entity_key_t m_value;
};

std::ostream & operator << ( std::ostream & out, EntityKey  key);

inline size_t hash_value(EntityKey k)
{
  return boost::hash_value(static_cast<size_t>(k.m_value));
}

}} // namespace stk::mesh

#endif /* STK_MESH_ENTITYKEY_HPP */
