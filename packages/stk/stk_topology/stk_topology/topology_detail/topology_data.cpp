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

#include "stk_topology/topology_detail/topology_data.hpp"

namespace stk { namespace topology_detail {

// Declaration of static class member arrays does not count as a definition,
// so they must be initialized outside the class.  However, if they are also
// constexpr, then they must be initialized at declaration time.  Defining
// them as empty here works around this conflict, preserving their initialization
// at declaration time.
//
// In C++17, these external definitions are no longer necessary, meaning that
// this entire file can be deleted.

constexpr bool topology_data<topology::INVALID_TOPOLOGY  >::spatial_dimension_vector[];
constexpr bool topology_data<topology::NODE              >::spatial_dimension_vector[];
constexpr bool topology_data<topology::PARTICLE          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_2            >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_2_1D         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_3_1D         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_LINE_2      >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_LINE_3      >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_SIDE_BEAM_2 >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_SIDE_BEAM_3 >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SPRING_2          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SPRING_3          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_3             >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_3_2D          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_4_2D          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_6_2D          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_4            >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_4_2D         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_8_2D         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_9_2D         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TET_4             >::spatial_dimension_vector[];
constexpr bool topology_data<topology::PYRAMID_5         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::WEDGE_6           >::spatial_dimension_vector[];
constexpr bool topology_data<topology::HEX_8             >::spatial_dimension_vector[];

constexpr topology::topology_t topology_data<topology::INVALID_TOPOLOGY  >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::PARTICLE          >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::LINE_2            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::BEAM_2            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::BEAM_3            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_LINE_2      >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_LINE_3      >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_SIDE_BEAM_2 >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_SIDE_BEAM_3 >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SPRING_2          >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::SPRING_3          >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::TRI_3             >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::TRI_6             >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::QUAD_4            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::QUAD_6            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::QUAD_8            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::TET_4             >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::TET_10            >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::PYRAMID_5         >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::PYRAMID_13        >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_6           >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_12          >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_15          >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::HEX_8             >::edge_topology_vector[];
constexpr topology::topology_t topology_data<topology::HEX_20            >::edge_topology_vector[];

constexpr topology::topology_t topology_data<topology::INVALID_TOPOLOGY             >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::NODE                         >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::LINE_2                       >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::TRI_3                        >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::QUAD_4                       >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::PARTICLE                     >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_3                  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_4                  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_6                  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_3_ALL_FACE_SIDES   >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_4_ALL_FACE_SIDES   >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_6_ALL_FACE_SIDES   >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_4                 >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_8                 >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_9                 >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_4_ALL_FACE_SIDES  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_8_ALL_FACE_SIDES  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_9_ALL_FACE_SIDES  >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::TET_4                        >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::TET_8                        >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::TET_10                       >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::PYRAMID_5                    >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::PYRAMID_13                   >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::PYRAMID_14                   >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_6                      >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_12                     >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_15                     >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::WEDGE_18                     >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::HEX_8                        >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::HEX_20                       >::face_topology_vector[];
constexpr topology::topology_t topology_data<topology::HEX_27                       >::face_topology_vector[];

constexpr topology::topology_t topology_data<topology::SHELL_TRI_3                  >::shell_side_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_4                  >::shell_side_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_TRI_6                  >::shell_side_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_4                 >::shell_side_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_8                 >::shell_side_topology_vector[];
constexpr topology::topology_t topology_data<topology::SHELL_QUAD_9                 >::shell_side_topology_vector[];

constexpr uint8_t topology_data<topology::INVALID_TOPOLOGY             >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::NODE                         >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::LINE_2                       >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TRI_3                        >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::QUAD_4                       >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PARTICLE                     >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_3                  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_4                  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_6                  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_3_ALL_FACE_SIDES   >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_4_ALL_FACE_SIDES   >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_TRI_6_ALL_FACE_SIDES   >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_4                 >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_8                 >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_9                 >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_4_ALL_FACE_SIDES  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_8_ALL_FACE_SIDES  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_9_ALL_FACE_SIDES  >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TET_4                        >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TET_8                        >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TET_10                       >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PYRAMID_5                    >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PYRAMID_13                   >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PYRAMID_14                   >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_6                      >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_12                     >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_15                     >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_18                     >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::HEX_8                        >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::HEX_20                       >::face_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::HEX_27                       >::face_node_ordinals_offsets[];

constexpr uint8_t topology_data<topology::INVALID_TOPOLOGY             >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::NODE                         >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::LINE_2                       >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_3                        >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_4                       >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PARTICLE                     >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_3                  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_4                  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_6                  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_3_ALL_FACE_SIDES   >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_4_ALL_FACE_SIDES   >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_TRI_6_ALL_FACE_SIDES   >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_4                 >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_8                 >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_9                 >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_4_ALL_FACE_SIDES  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_8_ALL_FACE_SIDES  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_QUAD_9_ALL_FACE_SIDES  >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_4                        >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_8                        >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_10                       >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_5                    >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_13                   >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_14                   >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_6                      >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_12                     >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_15                     >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_18                     >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_8                        >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_20                       >::face_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_27                       >::face_node_ordinals_vector[];

constexpr uint8_t topology_data<topology::INVALID_TOPOLOGY  >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::NODE              >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::LINE_2            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TRI_3             >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TRI_6             >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::QUAD_4            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::QUAD_6            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::QUAD_8            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PARTICLE          >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::BEAM_2            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::BEAM_3            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_LINE_2      >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_LINE_3      >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_SIDE_BEAM_2 >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::SHELL_SIDE_BEAM_3 >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TET_4             >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::TET_10            >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PYRAMID_5         >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::PYRAMID_13        >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_6           >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_12          >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::WEDGE_15          >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::HEX_8             >::edge_node_ordinals_offsets[];
constexpr uint8_t topology_data<topology::HEX_20            >::edge_node_ordinals_offsets[];

constexpr uint8_t topology_data<topology::INVALID_TOPOLOGY  >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::NODE              >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::LINE_2            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_3             >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_6             >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_4            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_6            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_8            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PARTICLE          >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::BEAM_2            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::BEAM_3            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_LINE_2      >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_LINE_3      >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_SIDE_BEAM_2 >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::SHELL_SIDE_BEAM_3 >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_4             >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_10            >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_5         >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_13        >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_6           >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_12          >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_15          >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_8             >::edge_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_20            >::edge_node_ordinals_vector[];

constexpr uint8_t topology_data<topology::INVALID_TOPOLOGY>::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::NODE            >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::LINE_2          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::LINE_3          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_3           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_4           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TRI_6           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_4          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_6          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_8          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::QUAD_9          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PARTICLE        >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_4           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_8           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_10          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::TET_11          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_5       >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_13      >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::PYRAMID_14      >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_6         >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_12        >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_15        >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::WEDGE_18        >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_8           >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_20          >::permutation_node_ordinals_vector[];
constexpr uint8_t topology_data<topology::HEX_27          >::permutation_node_ordinals_vector[];

}} // namespace stk::topology_detail
