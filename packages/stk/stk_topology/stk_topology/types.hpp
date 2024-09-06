// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
#ifndef STKTOPOLOGY_TYPES_TCC
#define STKTOPOLOGY_TYPES_TCC

#include "stk_topology/topology_type.hpp"

namespace stk {

struct topology::types {
  typedef topology_type<NODE>                          node;
  typedef topology_type<LINE_2>                        line_2;
  typedef topology_type<LINE_3>                        line_3;
  typedef topology_type<TRI_3>                         tri_3;
  typedef topology_type<TRI_4>                         tri_4;
  typedef topology_type<TRI_6>                         tri_6;
  typedef topology_type<QUAD_4>                        quad_4;
  typedef topology_type<QUAD_6>                        quad_6;
  typedef topology_type<QUAD_8>                        quad_8;
  typedef topology_type<QUAD_9>                        quad_9;
  typedef topology_type<PARTICLE>                      particle;
  typedef topology_type<LINE_2_1D>                     line_2_1d;
  typedef topology_type<LINE_3_1D>                     line_3_1d;
  typedef topology_type<BEAM_2>                        beam_2;
  typedef topology_type<BEAM_3>                        beam_3;
  typedef topology_type<SHELL_LINE_2>                  shell_line_2;
  typedef topology_type<SHELL_LINE_3>                  shell_line_3;
  typedef topology_type<SHELL_SIDE_BEAM_2>             shell_side_beam_2;
  typedef topology_type<SHELL_SIDE_BEAM_3>             shell_side_beam_3;
  typedef topology_type<TRI_3_2D>                      tri_3_2d;
  typedef topology_type<TRI_4_2D>                      tri_4_2d;
  typedef topology_type<TRI_6_2D>                      tri_6_2d;
  typedef topology_type<QUAD_4_2D>                     quad_4_2d;
  typedef topology_type<QUAD_8_2D>                     quad_8_2d;
  typedef topology_type<QUAD_9_2D>                     quad_9_2d;
  typedef topology_type<SHELL_TRI_3>                   shell_tri_3;
  typedef topology_type<SHELL_TRI_4>                   shell_tri_4;
  typedef topology_type<SHELL_TRI_6>                   shell_tri_6;
  typedef topology_type<SHELL_QUAD_4>                  shell_quad_4;
  typedef topology_type<SHELL_QUAD_8>                  shell_quad_8;
  typedef topology_type<SHELL_QUAD_9>                  shell_quad_9;
  typedef topology_type<SHELL_QUAD_4_ALL_FACE_SIDES>   shell_quad_4_all_face_sides;
  typedef topology_type<SHELL_QUAD_8_ALL_FACE_SIDES>   shell_quad_8_all_face_sides;
  typedef topology_type<SHELL_QUAD_9_ALL_FACE_SIDES>   shell_quad_9_all_face_sides;
  typedef topology_type<TET_4>                         tet_4;
  typedef topology_type<TET_8>                         tet_8;
  typedef topology_type<TET_10>                        tet_10;
  typedef topology_type<TET_11>                        tet_11;
  typedef topology_type<PYRAMID_5>                     pyramid_5;
  typedef topology_type<PYRAMID_13>                    pyramid_13;
  typedef topology_type<PYRAMID_14>                    pyramid_14;
  typedef topology_type<WEDGE_6>                       wedge_6;
  typedef topology_type<WEDGE_15>                      wedge_15;
  typedef topology_type<WEDGE_18>                      wedge_18;
  typedef topology_type<HEX_8>                         hex_8;
  typedef topology_type<HEX_20>                        hex_20;
  typedef topology_type<HEX_27>                        hex_27;
};

} //namespace stk

#endif //STKTOPOLOGY_TYPES_TCC


