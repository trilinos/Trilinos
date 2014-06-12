#ifndef STK_SAMBA_SAMBA_IO_IOSS_TOPOLOGY_MAP_HPP
#define STK_SAMBA_SAMBA_IO_IOSS_TOPOLOGY_MAP_HPP

#include <Ioss_SubSystem.h>
#include <samba/entity_topology.hpp>

namespace samba {
namespace io {


inline
std::string map_topology_to_ioss(entity_topology topology)
{
  std::string ioss_string("unknown");
  switch(topology())
  {
  case entity_topology::tetrahedron_4_type::value: 
    ioss_string = "tetra4"; 
    break;
  case entity_topology::hexahedron_8_type::value: 
    ioss_string = "hex8";
    break;
  case entity_topology::triangle_3_type::value:
    ioss_string = "triface";
    break;
  case entity_topology::quadrilateral_4_type::value:
    ioss_string = "quadface4";
    break;
  default: 
    break;
  }
  return ioss_string;
}


/// \todo Look for a data structure that enables this to be done
/// in sublinear time.
inline
entity_topology map_ioss_to_topology( const std::string &element_type ,
                                      const int node_count)
{
  /// \todo REFACTOR Is there a good type to return for an "unknown"
  /// topology type other than entity_topology::invalid()?

  entity_topology celltopo = entity_topology::invalid();

  const char *etype = element_type.c_str();
  if ( 0 == strncasecmp( "circle" , etype , 6 ) ) {
    celltopo = entity_topology::particle();
  }
  else if ( 0 == strncasecmp( "sphere" , etype , 6) ) {
    celltopo = entity_topology::particle();
  }
  // bar, beam, truss, rod...
  else if ( 0 == strncasecmp( "bar" , etype , 3 ) ) {
    if ( node_count == 2 ) {
      celltopo = entity_topology::beam_2();
    }
    else if ( node_count == 3 ) {
      celltopo = entity_topology::beam_3();
    }
  }
  else if ( 0 == strncasecmp( "shellline2d" , etype , 11 ) ) {
    if ( node_count == 2) {
      celltopo = entity_topology::shell_line_2();
    }
    else if ( node_count == 3) {
      celltopo = entity_topology::shell_line_3();
    }
  } else if ( 0 == strncasecmp( "shell" , etype , 5 ) ) {
    // shell4, shell8, shell9
    if ( node_count == 4 ) {
      celltopo =  entity_topology::shell_quadrilateral_4();
    }
    else if ( node_count == 8 ) {
      celltopo =  entity_topology::shell_quadrilateral_8();
    }
    else if ( node_count == 9 ) {
      celltopo =  entity_topology::shell_quadrilateral_9();
    }
  }
  else if ( 0 == strncasecmp( "quad" , etype , 3 ) ) {
    // The 2D types would be quad4, quad8, and quad9.
    // The 3D types would be quad faces of a hex... quadface4,
    // quadface8, quadface9.
    if ( node_count == 4 ) {
      celltopo =  entity_topology::quadrilateral_4();
    }
    else if ( node_count == 8 ) {
      celltopo =  entity_topology::quadrilateral_8();
    }
    else if ( node_count == 9 ) {
      celltopo =  entity_topology::quadrilateral_9();
    }
  }
  else if ( 0 == strncasecmp( "trishell" , etype , 8 ) ) {
    if ( node_count == 3 ) {
      celltopo =  entity_topology::shell_triangle_3();
    }
    else if ( node_count == 6 ) {
      celltopo =  entity_topology::shell_triangle_6();
    }
  }

  else if (0 == strncasecmp("triface", etype, 7) ||
      0 == strncasecmp("tri",     etype, 3)) {
    if ( node_count == 3 ) {
      celltopo =  entity_topology::triangle_3();
    }
    else if ( node_count == 4 ) {
      celltopo =  entity_topology::triangle_4();
    }
    else if ( node_count == 6 ) {
      celltopo =  entity_topology::triangle_6();
    }
  }

  else if ( 0 == strncasecmp( "pyramid" , etype , 7 ) ) {
    if ( node_count == 5 ) {
      celltopo =  entity_topology::pyramid_5();
    }
    else if ( node_count == 13 ) {
      celltopo =  entity_topology::pyramid_13();
    }
    else if ( node_count == 14 ) {
      celltopo =  entity_topology::pyramid_14();
    }
  }

  else if ( 0 == strncasecmp( "tetra" , etype , 5 ) ) {
    if ( node_count == 4 ) {
      celltopo =  entity_topology::tetrahedron_4();
    }
    else if ( node_count ==  8 ) {
      celltopo =  entity_topology::tetrahedron_8();
    }
    else if ( node_count == 10 ) {
      celltopo =  entity_topology::tetrahedron_10();
    }
  }
  else if ( 0 == strncasecmp( "wedge" , etype , 5 ) ) {
    if ( node_count == 6 ) {
      celltopo =  entity_topology::wedge_6();
    }
    else if ( node_count == 15 ) {
      celltopo =  entity_topology::wedge_15();
    }
    else if ( node_count == 18 ) {
      celltopo =  entity_topology::wedge_18();
    }
  }
  else if ( 0 == strncasecmp( "hex" , etype , 3 ) ) {
    if ( node_count == 8 ) {
      celltopo =  entity_topology::hexahedron_8();
    }
    else if ( node_count == 20 ) {
      celltopo =  entity_topology::hexahedron_20();
    }
    else if ( node_count == 27 ) {
      celltopo =  entity_topology::hexahedron_27();
    }
  }

  else if (0 == strncasecmp("edge", etype, 4)) {
    if ( node_count == 2) {
      // edge2, edge2d2, edge3d2
      celltopo =  entity_topology::line_2();
    }
    else if ( node_count == 3) {
      // edge3, edge2d3, edge3d3
      celltopo =  entity_topology::line_3();
    }
  }

  else if (0 == strncasecmp("node", etype, 4)) {
    celltopo = entity_topology::node();
  }

 if ( entity_topology::invalid() == celltopo ) {
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
entity_topology map_topology_ioss_to_cell(const Ioss::ElementTopology *topology)
{
  /// \todo REFACTOR Consider either using or integrating the
  /// Trilinos CellTopology package into or with the
  /// Ioss::ElementTopology classes. That would then totally
  /// eliminate the need for these fragile mapping functions.
  /// However, it would still need to be extensible via application
  /// registration of new type mappings.

  std::string name         = topology->name();
  int io_nodes_per_element = topology->number_nodes();

  return map_ioss_to_topology(name, io_nodes_per_element);
}

} // namespace io
} // namespace samba


#endif
