/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_baseImpl_FieldRestriction_hpp
#define stk_mesh_baseImpl_FieldRestriction_hpp

#include <vector>
#include <Shards_Array.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/util/SimpleArrayOps.hpp>

namespace stk {
namespace mesh {

struct FieldRestriction {
  typedef shards::array_traits::int_t size_type ;

  EntityKey key ;
  size_type stride[ MaximumFieldDimension ];

  FieldRestriction() : key() {
    Copy<MaximumFieldDimension>( stride , size_type(0) );
  }

  FieldRestriction( const FieldRestriction & rhs ): key( rhs.key ) {
    Copy< MaximumFieldDimension >( stride , rhs.stride );
  }

  FieldRestriction & operator = ( const FieldRestriction & rhs ) {
    key = rhs.key ;
    Copy< MaximumFieldDimension >( stride , rhs.stride );
    return *this ;
  }

  FieldRestriction( EntityRank rank , EntityId ordinal)
    : key( EntityKey(rank, ordinal))
  {
    Copy< MaximumFieldDimension >( stride , size_type(0) );
  }

  EntityRank type()    const { return entity_rank( key ); }
  EntityId ordinal() const { return entity_id( key ); }
};

typedef std::vector<FieldRestriction> FieldRestrictionVector;

} // namespace mesh
} // namespace stk



#endif // stk_mesh_baseImpl_FieldRestriction_hpp
