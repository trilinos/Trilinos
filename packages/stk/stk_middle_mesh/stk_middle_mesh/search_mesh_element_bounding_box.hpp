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

#ifndef  STK_MIDDLE_MESH_PARALLEL_SEARCH_HPP
#define  STK_MIDDLE_MESH_PARALLEL_SEARCH_HPP

#include "search_mesh_element_bounding_box_base.hpp"
#include "mesh.hpp"
#include "field.hpp"
#include <stk_search/CoarseSearch.hpp>
#include "stk_search/Box.hpp"

#include <limits>
#include <map>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {


class SearchMeshElementBoundingBox : public SearchMeshElementBoundingBoxBase {
 public:

  SearchMeshElementBoundingBox(std::shared_ptr<mesh::Mesh> inputMesh, MPI_Comm unionComm)
    : m_mesh(inputMesh)
    , m_unionComm(unionComm)
    {}

  void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const override;

  std::shared_ptr<mesh::Mesh> get_mesh() const override { return m_mesh; }

 private:
  using Point         = stk::search::Point<double>;

  void fill_bounding_box(mesh::MeshEntityPtr element, stk::search::Point<double>& minCorner,
                         stk::search::Point<double>& maxCorner) const;

  std::shared_ptr<mesh::Mesh> m_mesh;
  MPI_Comm m_unionComm;
};




}
}
}
}

#endif
