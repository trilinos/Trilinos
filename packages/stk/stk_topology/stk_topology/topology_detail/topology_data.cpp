// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

//#include <stk_topology/topology_detail/topology_data.hpp>
#include <stk_topology/topology.hpp>

namespace stk { namespace topology_detail {

constexpr bool topology_data<topology::INVALID_TOPOLOGY>::spatial_dimension_vector[];
constexpr bool topology_data<topology::NODE            >::spatial_dimension_vector[];
constexpr bool topology_data<topology::PARTICLE        >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_2          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_2_1D       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::LINE_3_1D       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_LINE_2    >::spatial_dimension_vector[];
constexpr bool topology_data<topology::SHELL_LINE_3    >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_3           >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_3_2D        >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_4_2D        >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TRI_6_2D        >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_4          >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_4_2D       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_8_2D       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::QUAD_9_2D       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::TET_4           >::spatial_dimension_vector[];
constexpr bool topology_data<topology::PYRAMID_5       >::spatial_dimension_vector[];
constexpr bool topology_data<topology::WEDGE_6         >::spatial_dimension_vector[];
constexpr bool topology_data<topology::HEX_8           >::spatial_dimension_vector[];

}} // namespace stk::topology_detail
