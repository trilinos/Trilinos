#ifndef STKTOPOLOGY_TYPES_TCC
#define STKTOPOLOGY_TYPES_TCC

#include <stk_topology/topology_type.tcc>

namespace stk {

struct topology::types {
  typedef topology_type<NODE>          node;
  typedef topology_type<LINE_2>        line_2;
  typedef topology_type<LINE_3>        line_3;
  typedef topology_type<TRI_3>         tri_3;
  typedef topology_type<TRI_4>         tri_4;
  typedef topology_type<TRI_6>         tri_6;
  typedef topology_type<QUAD_4>        quad_4;
  typedef topology_type<QUAD_8>        quad_8;
  typedef topology_type<QUAD_9>        quad_9;
  typedef topology_type<PARTICLE>      particle;
  typedef topology_type<LINE_2_1D>     line_2_1d;
  typedef topology_type<LINE_3_1D>     line_3_1d;
  typedef topology_type<BEAM_2>        beam_2;
  typedef topology_type<BEAM_3>        beam_3;
  typedef topology_type<SHELL_LINE_2>  shell_line_2;
  typedef topology_type<SHELL_LINE_3>  shell_line_3;
  typedef topology_type<TRI_3_2D>      tri_3_2d;
  typedef topology_type<TRI_4_2D>      tri_4_2d;
  typedef topology_type<TRI_6_2D>      tri_6_2d;
  typedef topology_type<QUAD_4_2D>     quad_4_2d;
  typedef topology_type<QUAD_8_2D>     quad_8_2d;
  typedef topology_type<QUAD_9_2D>     quad_9_2d;
  typedef topology_type<SHELL_TRI_3>   shell_tri_3;
  typedef topology_type<SHELL_TRI_4>   shell_tri_4;
  typedef topology_type<SHELL_TRI_6>   shell_tri_6;
  typedef topology_type<SHELL_QUAD_4>  shell_quad_4;
  typedef topology_type<SHELL_QUAD_8>  shell_quad_8;
  typedef topology_type<SHELL_QUAD_9>  shell_quad_9;
  typedef topology_type<TET_4>         tet_4;
  typedef topology_type<TET_8>         tet_8;
  typedef topology_type<TET_10>        tet_10;
  typedef topology_type<TET_11>        tet_11;
  typedef topology_type<PYRAMID_5>     pyramid_5;
  typedef topology_type<PYRAMID_13>    pyramid_13;
  typedef topology_type<PYRAMID_14>    pyramid_14;
  typedef topology_type<WEDGE_6>       wedge_6;
  typedef topology_type<WEDGE_15>      wedge_15;
  typedef topology_type<WEDGE_18>      wedge_18;
  typedef topology_type<HEX_8>         hex_8;
  typedef topology_type<HEX_20>        hex_20;
  typedef topology_type<HEX_27>        hex_27;
};

} //namespace stk

#endif //STKTOPOLOGY_TYPES_TCC


