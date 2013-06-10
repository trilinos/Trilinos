#ifndef STK_SIERRA_MESH_IO_ARRAY_MESH_IOSS_TOPOLOGY_HPP
#define STK_SIERRA_MESH_IO_ARRAY_MESH_IOSS_TOPOLOGY_HPP

#include <Ioss_SubSystem.h>

namespace sierra {
namespace mesh {
namespace io {

inline
std::string map_array_mesh_side_topology_to_ioss(stk::topology elem_topology, int node_count)
{
	std::string ioss_string("unknown");

	switch(elem_topology)
	{
    case stk::topology::TET_4: ioss_string = "triface"; break;
	  case stk::topology::HEX_8: ioss_string = "quadface4"; break;
	  default: break;
	}

	return ioss_string;
}

inline
std::string map_array_mesh_topology_to_ioss(stk::topology elem_topology)
{
	std::string ioss_string("unknown");

	switch(elem_topology)
	{
    case stk::topology::NODE:  ioss_string = "node"; break;
    case stk::topology::TET_4: ioss_string = "tetra"; break;
	  case stk::topology::HEX_8: ioss_string = "hex"; break;
	  default: break;
	}

	return ioss_string;
}

inline
stk::topology map_ioss_topology_to_array_mesh( const std::string &element_type ,
                                     int node_count)
{
//  int topology = sierra::mesh::InvalidTopology::value;
  stk::topology topology = stk::topology::INVALID_TOPOLOGY;

  const char *etype = element_type.c_str();
  if ( 0 == strncasecmp( "circle" , etype , 6 ) ) {
  }
  else if ( 0 == strncasecmp( "sphere" , etype , 6) ) {
  }
  // bar, beam, truss, rod...
  else if ( 0 == strncasecmp( "bar" , etype , 3 ) ) {
    if ( node_count == 2 ) {
    }
    else if ( node_count == 3 ) {
    }
  }
  else if ( 0 == strncasecmp( "shellline2d" , etype , 11 ) ) {
    if ( node_count == 2) {
    }
    else if ( node_count == 3) {
    }
  } else if ( 0 == strncasecmp( "shell" , etype , 5 ) ) {
    // shell4, shell8, shell9
    if ( node_count == 4 ) {
    }
    else if ( node_count == 8 ) {
    }
    else if ( node_count == 9 ) {
    }
  }
  else if ( 0 == strncasecmp( "quad" , etype , 3 ) ) {
    // The 2D types would be quad4, quad8, and quad9.
    // The 3D types would be quad faces of a hex... quadface4,
    // quadface8, quadface9.
    if ( node_count == 4 ) {
    }
    else if ( node_count == 8 ) {
    }
    else if ( node_count == 9 ) {
    }
  }
  else if ( 0 == strncasecmp( "trishell" , etype , 8 ) ) {
    if ( node_count == 3 ) {
    }
    else if ( node_count == 6 ) {
    }
  }

  else if (0 == strncasecmp("triface", etype, 7) ||
      0 == strncasecmp("tri",     etype, 3)) {
    if ( node_count == 3 ) {
    }
    else if ( node_count == 4 ) {
    }
    else if ( node_count == 6 ) {
    }
  }

  else if ( 0 == strncasecmp( "pyramid" , etype , 7 ) ) {
    if ( node_count == 5 ) {
    }
    else if ( node_count == 13 ) {
    }
    else if ( node_count == 14 ) {
    }
  }

  else if ( 0 == strncasecmp( "tetra" , etype , 5 ) ) {
    if ( node_count == 4 ) {
//      topology = sierra::mesh::Tet4::value;
      topology = stk::topology::TET_4;
    }
    else if ( node_count ==  8 ) {
    }
    else if ( node_count == 10 ) {
    }
  }
  else if ( 0 == strncasecmp( "wedge" , etype , 5 ) ) {
    if ( node_count == 6 ) {
    }
    else if ( node_count == 15 ) {
    }
    else if ( node_count == 18 ) {
    }
  }
  else if ( 0 == strncasecmp( "hex" , etype , 3 ) ) {
    if ( node_count == 8 ) {
//      topology = sierra::mesh::Hex8::value;
      topology = stk::topology::HEX_8;
    }
    else if ( node_count == 20 ) {
    }
    else if ( node_count == 27 ) {
    }
  }

  else if (0 == strncasecmp("edge", etype, 4)) {
    if ( node_count == 2) {
      // edge2, edge2d2, edge3d2
    }
    else if ( node_count == 3) {
      // edge3, edge2d3, edge3d3
    }
  }

  else if (0 == strncasecmp("node", etype, 4)) {
//    topology = sierra::mesh::Node::value;
    topology = stk::topology::NODE;
  }

  if ( topology == stk::topology::INVALID_TOPOLOGY ) {
    std::ostringstream oss;
    oss << "ERROR, unsupported topology name = '" << element_type
      << "' , node_count = " << node_count;
    throw std::runtime_error(oss.str());
  }

  std::cout << topology << std::endl;

  return topology;
}

}//namespace io
}//namespace mesh
}//namespace sierra

#endif // STK_SIERRA_MESH_IO_ARRAY_MESH_IOSS_TOPOLOGY_HPP

