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

#ifndef STK_MESH_BASE_CREATE_EDGES_ON_FACES_H
#define STK_MESH_BASE_CREATE_EDGES_ON_FACES_H

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"


namespace stk {
namespace mesh {
namespace impl {

class EdgesOnFacesCreator
{
  public:
    EdgesOnFacesCreator(BulkData& bulkData, const Selector& surfaceSelector, Part* edgePart=nullptr) :
      m_bulk(bulkData),
      m_surfaceSelector(surfaceSelector),
      m_surfaceAndOwned(m_surfaceSelector & bulkData.mesh_meta_data().locally_owned_part()),
      m_edgePart(edgePart)
    {}

    Part* create_edges();

  public:

    void create_surface_edges();

    void attach_edges_to_non_owned_faces();

    [[nodiscard]] stk::topology get_edge_topology() const;

    stk::mesh::Part* create_edge_part_if_needed() const;

    void check_edge_part_topology();

    int compute_upper_bound_on_num_edges() const;

    void check_face_topology(stk::topology faceTopo);

    [[nodiscard]] stk::mesh::Entity get_common_edge(stk::mesh::Entity entity1, stk::mesh::Entity entity2) const;

    void sort_nodes_for_global_consistency(std::vector<Entity>& edgeNodes) const;

    BulkData& m_bulk;
    Selector m_surfaceSelector;
    Selector m_surfaceAndOwned;
    Part* m_edgePart;

};


}
}
}

#endif