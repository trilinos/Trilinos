#include <stk_topology/topology.hpp>

namespace stk {

const char * topology::rank_names[] =
{
    "NODE_RANK"
  , "EDGE_RANK"
  , "FACE_RANK"
  , "ELEMENT_RANK"
  , "INVALID_RANK"
};

const char * topology::topology_names[] =
{
    "NODE"
  , "LINE_2"
  , "LINE_3"
  , "TRIANGLE_3"
  , "TRIANGLE_4"
  , "TRIANGLE_6"
  , "QUADRILATERAL_4"
  , "QUADRILATERAL_8"
  , "QUADRILATERAL_9"
  , "PARTICLE"
  , "LINE_2_1D"
  , "LINE_3_1D"
  , "BEAM_2"
  , "BEAM_3"
  , "SHELL_LINE_2"
  , "SHELL_LINE_3"
  , "TRIANGLE_3_2D"
  , "TRIANGLE_4_2D"
  , "TRIANGLE_6_2D"
  , "QUADRILATERAL_4_2D"
  , "QUADRILATERAL_8_2D"
  , "QUADRILATERAL_9_2D"
  , "SHELL_TRIANGLE_3"
  , "SHELL_TRIANGLE_4"
  , "SHELL_TRIANGLE_6"
  , "SHELL_QUADRILATERAL_4"
  , "SHELL_QUADRILATERAL_8"
  , "SHELL_QUADRILATERAL_9"
  , "TETRAHEDRON_4"
  , "TETRAHEDRON_8"
  , "TETRAHEDRON_10"
  , "TETRAHEDRON_11"
  , "PYRAMID_5"
  , "PYRAMID_13"
  , "PYRAMID_14"
  , "WEDGE_6"
  , "WEDGE_15"
  , "WEDGE_18"
  , "HEXAHEDRON_8"
  , "HEXAHEDRON_20"
  , "HEXAHEDRON_27"
  , "INVALID_TOPOLOGY"
};

} //namespace stk


