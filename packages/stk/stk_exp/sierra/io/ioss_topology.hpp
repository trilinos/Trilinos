#ifndef STK_SIERRA_MESH_IO_IOSS_TOPOLOGY_HPP
#define STK_SIERRA_MESH_IO_IOSS_TOPOLOGY_HPP

#include <Shards_CellTopology.hpp>
#include <Ioss_SubSystem.h>

namespace sierra {
namespace mesh {
namespace io {

inline
const CellTopologyData *map_ioss_to_topology( const std::string &element_type ,
          const int node_count)
{
  /// \todo REFACTOR Is there a good type to return for an "unknown"
  /// topology type other than NULL?

  const CellTopologyData* celltopo = NULL ;

  const char *etype = element_type.c_str();
  if ( 0 == strncasecmp( "circle" , etype , 6 ) ) {
    celltopo = shards::getCellTopologyData< shards::Particle >();
  }
  else if ( 0 == strncasecmp( "sphere" , etype , 6) ) {
    celltopo = shards::getCellTopologyData< shards::Particle >();
  }
  // bar, beam, truss, rod...
  else if ( 0 == strncasecmp( "bar" , etype , 3 ) ) {
    if ( node_count == 2 ) {
      celltopo = shards::getCellTopologyData< shards::Beam<2> >();
    }
    else if ( node_count == 3 ) {
      celltopo = shards::getCellTopologyData< shards::Beam<3> >();
    }
  }
  else if ( 0 == strncasecmp( "shellline2d" , etype , 11 ) ) {
    if ( node_count == 2) {
      celltopo = shards::getCellTopologyData< shards::ShellLine<2> >();
    }
    else if ( node_count == 3) {
      celltopo = shards::getCellTopologyData< shards::ShellLine<3> >();
    }
  } else if ( 0 == strncasecmp( "shell" , etype , 5 ) ) {
    // shell4, shell8, shell9
    if ( node_count == 4 ) {
      celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<4> >();
    }
    else if ( node_count == 8 ) {
      celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<8> >();
    }
    else if ( node_count == 9 ) {
      celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<9> >();
    }
  }
  else if ( 0 == strncasecmp( "quad" , etype , 3 ) ) {
    // The 2D types would be quad4, quad8, and quad9.
    // The 3D types would be quad faces of a hex... quadface4,
    // quadface8, quadface9.
    if ( node_count == 4 ) {
      celltopo = shards::getCellTopologyData< shards::Quadrilateral<4> >();
    }
    else if ( node_count == 8 ) {
      celltopo = shards::getCellTopologyData< shards::Quadrilateral<8> >();
    }
    else if ( node_count == 9 ) {
      celltopo = shards::getCellTopologyData< shards::Quadrilateral<9> >();
    }
  }
  else if ( 0 == strncasecmp( "trishell" , etype , 8 ) ) {
    if ( node_count == 3 ) {
      celltopo = shards::getCellTopologyData< shards::ShellTriangle<3> >();
    }
    else if ( node_count == 6 ) {
      celltopo = shards::getCellTopologyData< shards::ShellTriangle<6> >();
    }
  }

  else if (0 == strncasecmp("triface", etype, 7) ||
      0 == strncasecmp("tri",     etype, 3)) {
    if ( node_count == 3 ) {
      celltopo = shards::getCellTopologyData< shards::Triangle<3> >();
    }
    else if ( node_count == 4 ) {
      celltopo = shards::getCellTopologyData< shards::Triangle<4> >();
    }
    else if ( node_count == 6 ) {
      celltopo = shards::getCellTopologyData< shards::Triangle<6> >();
    }
  }

  else if ( 0 == strncasecmp( "pyramid" , etype , 7 ) ) {
    if ( node_count == 5 ) {
      celltopo = shards::getCellTopologyData< shards::Pyramid<5> >();
    }
    else if ( node_count == 13 ) {
      celltopo = shards::getCellTopologyData< shards::Pyramid<13> >();
    }
    else if ( node_count == 14 ) {
      celltopo = shards::getCellTopologyData< shards::Pyramid<14> >();
    }
  }

  /// \todo REFACTOR Need to handle 8-node tet...
  else if ( 0 == strncasecmp( "tetra" , etype , 5 ) ) {
    if ( node_count == 4 ) {
      celltopo = shards::getCellTopologyData< shards::Tetrahedron<4> >();
    }
    else if ( node_count ==  8 ) {
      celltopo = shards::getCellTopologyData< shards::Tetrahedron<8> >();
    }
    else if ( node_count == 10 ) {
      celltopo = shards::getCellTopologyData< shards::Tetrahedron<10> >();
    }
  }
  else if ( 0 == strncasecmp( "wedge" , etype , 5 ) ) {
    if ( node_count == 6 ) {
      celltopo = shards::getCellTopologyData< shards::Wedge<6> >();
    }
    else if ( node_count == 15 ) {
      celltopo = shards::getCellTopologyData< shards::Wedge<15> >();
    }
    else if ( node_count == 18 ) {
      celltopo = shards::getCellTopologyData< shards::Wedge<18> >();
    }
  }
  else if ( 0 == strncasecmp( "hex" , etype , 3 ) ) {
    if ( node_count == 8 ) {
      celltopo = shards::getCellTopologyData< shards::Hexahedron<8> >();
    }
    else if ( node_count == 20 ) {
      celltopo = shards::getCellTopologyData< shards::Hexahedron<20> >();
    }
    else if ( node_count == 27 ) {
      celltopo = shards::getCellTopologyData< shards::Hexahedron<27> >();
    }
  }

  else if (0 == strncasecmp("edge", etype, 4)) {
    if ( node_count == 2) {
      // edge2, edge2d2, edge3d2
      celltopo = shards::getCellTopologyData< shards::Line<2> >();
    }
    else if ( node_count == 3) {
      // edge3, edge2d3, edge3d3
      celltopo = shards::getCellTopologyData< shards::Line<3> >();
    }
  }

  else if (0 == strncasecmp("node", etype, 4)) {
    celltopo = shards::getCellTopologyData< shards::Node >();
  }

  if ( NULL == celltopo ) {
    std::ostringstream oss;
    oss << "ERROR, unsupported topology name = '" << element_type
      << "' , node_count = " << node_count;
    throw std::runtime_error(oss.str());
  }

  return celltopo;
}

/// \todo QUESTION Should this function be at the application level,
/// or provided by stk_io? In either case, applications should have
/// capabilty to register new mappings.
// ========================================================================
inline
const CellTopologyData *map_topology_ioss_to_cell(const Ioss::ElementTopology *topology)
{
  /// \todo REFACTOR Consider either using or integrating the
  /// Trilinos CellTopology package into or with the
  /// Ioss::ElementTopology classes. That would then totally
  /// eliminate the need for these fragile mapping functions.
  /// However, it would still need to be extensible via application
  /// registration of new type mappings.

  std::string name         = topology->name();
  int io_nodes_per_element = topology->number_nodes();

  const CellTopologyData *cell_topology = map_ioss_to_topology(name, io_nodes_per_element);

  return cell_topology;
}

}//namespace io
}//namespace mesh
}//namespace sierra

#endif // STK_SIERRA_MESH_IO_IOSS_TOPOLOGY_HPP

