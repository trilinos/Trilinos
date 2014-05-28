#include <stk_mesh/fem/Stencils.hpp>

#include <stk_util/environment/ReportHandler.hpp>

namespace stk {
namespace mesh {
namespace fem {

int
element_node_stencil_2d(
  EntityRank            from_type ,
  EntityRank            to_type ,
  unsigned              identifier )
{
  static const size_t spatial_dimension = 2;

  int ordinal = -1 ;

  if ( spatial_dimension == from_type && FEMMetaData::NODE_RANK == to_type ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}


int
element_node_stencil_3d(
  EntityRank            from_type ,
  EntityRank            to_type ,
  unsigned              identifier )
{
  static const size_t spatial_dimension = 3;

  int ordinal = -1 ;

  if ( spatial_dimension == from_type && FEMMetaData::NODE_RANK == to_type ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}


relation_stencil_ptr
get_element_node_stencil(
  size_t                spatial_dimension) 
{
  ThrowRequire(spatial_dimension == 2 || spatial_dimension == 3);
  
  if (spatial_dimension == 3)
    return & element_node_stencil_3d;
  else // if (spatial_dimension == 2)
    return & element_node_stencil_2d;
}

} // namespace fem
} // namespace mesh
} // namespace stk
