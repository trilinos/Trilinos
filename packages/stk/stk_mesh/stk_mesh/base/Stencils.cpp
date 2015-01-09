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


#include <stk_mesh/base/Stencils.hpp>
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequire
#include "stk_mesh/base/Types.hpp"      // for EntityRank, etc

namespace stk {
namespace mesh {

int
element_node_stencil_2d(
  EntityRank            from_type ,
  EntityRank            to_type ,
  unsigned              identifier )
{
  static const size_t spatial_dimension = 2;

  int ordinal = -1 ;

  if ( spatial_dimension == from_type && stk::topology::NODE_RANK == to_type ) {
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

  if ( spatial_dimension == from_type && stk::topology::NODE_RANK == to_type ) {
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

} // namespace mesh
} // namespace stk
