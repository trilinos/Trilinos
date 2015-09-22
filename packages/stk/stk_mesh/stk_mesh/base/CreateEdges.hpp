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

#ifndef stk_mesh_CreateEdges_hpp
#define stk_mesh_CreateEdges_hpp

#include <stk_mesh/base/Types.hpp>      // for EntityVector, etc
#include <unordered_map>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk {
  namespace mesh {
    class BulkData;
    class Selector;
    class Part;

    /** Create all the edges in the mesh and attach them to
     * existing elements and defined faces.
     *
     * This is a parallel collective function (it should be called on all
     * processors at the same time
     *
     */
    void create_edges(  BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges = 0 );

    void create_edges( BulkData & mesh );

    namespace impl {
      typedef std::unordered_map<EntityVector,Entity,stk::mesh::impl::HashValueForEntityVector> edge_map_type;
      void connect_faces_to_edges(BulkData & mesh,
                                  const Selector & element_selector,
                                  edge_map_type edge_map);
    }
  }
}
#endif

