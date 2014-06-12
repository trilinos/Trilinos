#ifndef stk_mesh_fem_CellTopology_hpp
#define stk_mesh_fem_CellTopology_hpp

#ifdef HAVE_SHARDS_DEBUG
#define STK_MESH_FEM_CHECK_REQUIRE( S )  S
#else
#define STK_MESH_FEM_CHECK_REQUIRE( S ) /* empty */
#endif

#include <Shards_CellTopologyTraits.hpp>
#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>

namespace stk_classic {
namespace mesh {
namespace fem {

typedef shards::CellTopology CellTopology;



template< typename id_type >
int findPermutation( const CellTopology top ,
                     const id_type * const expected_node ,
                     const id_type * const actual_node )
{
  return shards::findPermutation( *top.getCellTopologyData() , expected_node , actual_node );
}

/** \} */

} // namespace fem
} // namespace mesh
} // namespace stk_classic

#undef STK_MESH_FEM_CHECK_REQUIRE

#endif // stk_mesh_fem_CellTopology_hpp

