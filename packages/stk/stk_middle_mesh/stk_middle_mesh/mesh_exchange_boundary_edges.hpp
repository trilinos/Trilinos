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

#include <unordered_map>
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/parallel_exchange.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "stk_middle_mesh/entity_sorted_by_owner.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "stk_middle_mesh/field.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct BoundaryVertexInfo
{
  RemoteSharedEntity owner;
  utils::Point pts;
};

struct BoundaryEdgeInfo
{
  RemoteSharedEntity vertex1Owner, vertex2Owner;
};

struct VertexUpdateInfo
{
  int localId;
  utils::Point pts;
};


class MeshExchangeBoundaryEdges
{
 public:
  template <typename T>
  using Exchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<T>;

  using EdgeExchange           = Exchanger<BoundaryEdgeInfo>;
  using VertexExchange         = Exchanger<BoundaryVertexInfo>;
  using EdgeCountExchange      = Exchanger<int>;
  using VertexPointExchanger   = Exchanger<utils::Point>;
  using Bool                   = int_least8_t;

  explicit MeshExchangeBoundaryEdges(std::shared_ptr<Mesh> inputMesh, MPI_Comm unionComm, int rootRank)
	 : m_inputMesh(inputMesh),
     m_boundaryMesh(nullptr),
     m_unionComm(unionComm),
     m_root(rootRank),
     m_myrank(utils::impl::comm_rank(m_unionComm)),
     m_entityCountExchanger(m_unionComm),
     m_sortedVerticesByOwner(utils::impl::comm_size(m_unionComm))
	{
    if (m_myrank == m_root)
      m_boundaryMesh = make_empty_mesh(MPI_COMM_SELF);

    if (m_inputMesh)
      set_boundary_vertex_field();
  }

	std::shared_ptr<Mesh> get_boundary_edge_mesh();

  void update_remote_vertices();

 private:
  void update_shared_vertex_coordinates();

  void apply_snapped_boundaries(VertexPointExchanger& exchanger);

  void fill_updated_remote_vertices_info(VertexPointExchanger& exchanger);

  void fill_boundary_vertices();

  void fill_local_boundary_edges();

  void gather_boundary_entity_counts_to_root();

  void gather_boundary_edges_and_vertices_to_root(VertexExchange& boundaryVertexExchanger,
                                                  EdgeExchange& boundaryEdgeExchanger);

  void fill_boundary_vertex_info(MeshEntityPtr vertex, BoundaryVertexInfo& boundaryInfo);

  void send_boundary_edges_to_root(EdgeExchange& boundaryEdgeExchanger);

  void send_boundary_vertices_to_root(VertexExchange& boundaryVertexExchanger);

  void recv_boundary_edges_and_vertices_in_root(VertexExchange& boundaryVertexExchanger,
                                                EdgeExchange& boundaryEdgeExchanger);

  void gather_boundary_edge_and_vertex_count_info();

  void add_boundary_edges_on_root(VertexExchange& boundaryVertexExchanger, EdgeExchange& boundaryEdgeExchanger);

  void set_boundary_vertex_field();

  bool is_local_boundary_vertex(MeshEntityPtr vertex);

  bool is_boundary_vertex(MeshEntityPtr vertex);

  bool is_boundary_edge(MeshEntityPtr edge);

  std::shared_ptr<Mesh> m_inputMesh;
  FieldPtr<Bool> m_inputMeshIsBoundaryNode;
  std::shared_ptr<Mesh> m_boundaryMesh;
	MPI_Comm m_unionComm;
  int m_root;
  int m_myrank;
  EdgeCountExchange m_entityCountExchanger;
	std::vector<MeshEntityPtr> m_localBoundaryEdges; 
	std::vector<MeshEntityPtr> m_localBoundaryVertices;
  std::vector<MeshEntityPtr> m_nonOwnedBoundaryVertices;
  EntitySortedByOwner m_sortedVerticesByOwner;
};

}
}
}
}
