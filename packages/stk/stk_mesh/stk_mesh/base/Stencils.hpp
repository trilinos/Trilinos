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

#ifndef stk_mesh_Stencils_hpp
#define stk_mesh_Stencils_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Types.hpp>      // for EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc



namespace stk {
namespace mesh {

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after April 2021 and delete this header after trilinos/percept is updated

STK_DEPRECATED relation_stencil_ptr get_element_node_stencil(size_t spatial_dimension);

template<class TopologyTraits, EntityRank element_rank >
STK_DEPRECATED int element_node_stencil( EntityRank , EntityRank , unsigned );


template<class TopologyTraits, EntityRank element_rank >
STK_DEPRECATED int element_node_stencil( EntityRank from_type , EntityRank to_type , unsigned identifier )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( element_rank == from_type &&
       stk::topology::NODE_RANK == to_type &&
       identifier < number_node ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

#endif // STK_HIDE_DEPRECATED_CODE

} // namespace mesh
} // namespace stk

#endif //  stk_mesh_Stencils_hpp
