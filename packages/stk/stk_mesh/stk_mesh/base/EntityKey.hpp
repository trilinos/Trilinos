/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_ENTITYKEY_HPP
#define STK_MESH_ENTITYKEY_HPP

#include <iosfwd>

#include <stdint.h>

#include <stk_mesh/base/Types.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <stk_util/environment/ReportHandler.hpp>

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
  EntityRank rank() const { return m_value >> RANK_SHIFT; }

  bool is_valid() const { return m_value != INVALID; }

  operator entity_key_t() const { return m_value; }

  entity_key_t m_value;
};

std::ostream & operator << ( std::ostream & out, EntityKey  key);

}} // namespace stk::mesh

#endif /* STK_MESH_ENTITYKEY_HPP */
