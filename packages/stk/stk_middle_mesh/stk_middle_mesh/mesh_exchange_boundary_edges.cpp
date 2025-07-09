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

#include "stk_middle_mesh/mesh_exchange_boundary_edges.hpp"
#include "field_shared_reduction.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshExchangeBoundaryEdges::update_shared_vertex_coordinates()
{
  stk::DataExchangeKnownPatternNonBlockingBuffer<VertexUpdateInfo> exchanger(m_inputMesh->get_comm());
  int commSize = utils::impl::comm_size(m_inputMesh->get_comm());

  std::vector<int> recvCounts(commSize, 0);
  for (auto& nonOwnedVert : m_nonOwnedBoundaryVertices)
  {
    int ownerRank = get_owner(m_inputMesh, nonOwnedVert);
    recvCounts[ownerRank]++;
  }

  for (int i=0; i < commSize; ++i)
  {
    if (recvCounts[i] > 0)
    {
      exchanger.get_recv_buf(i).resize(recvCounts[i]);
    }
  }

  for (auto& ownedVert : m_localBoundaryVertices)
    for (int i=0; i < ownedVert->count_remote_shared_entities(); ++i)
    {
      RemoteSharedEntity remote = ownedVert->get_remote_shared_entity(i);
      exchanger.get_send_buf(remote.remoteRank).push_back(VertexUpdateInfo{remote.remoteId, ownedVert->get_point_orig(0)});
    }

  exchanger.start_nonblocking();

  auto f = [&](int /*rank*/, const std::vector<VertexUpdateInfo>& buf)
  {
    for (const VertexUpdateInfo& info : buf)
    {
      m_inputMesh->get_vertices()[info.localId]->set_point_orig(0, info.pts);
    }
  };

  exchanger.complete_receives(f);
}

void MeshExchangeBoundaryEdges::apply_snapped_boundaries(VertexPointExchanger& exchanger)
{
  auto f = [&](int /*rank*/, const std::vector<utils::Point>& updatedVertexInfo)
  {
    for (size_t i=0; i < m_localBoundaryVertices.size(); ++i)
    {
      MeshEntityPtr vertex = m_localBoundaryVertices[i];
      vertex->set_point_orig(0, updatedVertexInfo[i]);
    }
  };

  exchanger.complete_receives(f);

  if (m_inputMesh)
    update_shared_vertex_coordinates();
}


void MeshExchangeBoundaryEdges::fill_updated_remote_vertices_info(VertexPointExchanger& exchanger)
{

  int nprocs = utils::impl::comm_size(m_unionComm);
  for (int rank=0; rank < nprocs; ++rank)
  {
    auto& sendBuf = exchanger.get_send_buf(rank);
    sendBuf.reserve(m_sortedVerticesByOwner.get_values(rank).size());
    for (auto& vertex : m_sortedVerticesByOwner.get_values(rank))
      sendBuf.push_back(vertex->get_point_orig(0));
  }
}

void MeshExchangeBoundaryEdges::update_remote_vertices()
{
  VertexPointExchanger exchanger(m_unionComm);

  if (m_inputMesh)
  {
    exchanger.get_recv_buf(m_root).resize(m_localBoundaryVertices.size());
  }

  if (m_myrank == m_root) {
    fill_updated_remote_vertices_info(exchanger);
  }

  exchanger.start_nonblocking();

  if (m_inputMesh)
    apply_snapped_boundaries(exchanger);
}

void MeshExchangeBoundaryEdges::add_boundary_edges_on_root(VertexExchange& boundaryVertexExchanger,
                                                           EdgeExchange& boundaryEdgeExchanger)
{
  
  auto unpackVerts = [&](int /*rank*/, const std::vector<BoundaryVertexInfo>& recvVertexBuffer)
  {
    for (unsigned j = 0; j < recvVertexBuffer.size(); ++j) {
      auto newVertex = m_boundaryMesh->create_vertex(recvVertexBuffer[j].pts);
      m_sortedVerticesByOwner.insert(recvVertexBuffer[j].owner, newVertex);
    }
  };

  auto unpackEdges = [&](int /*rank*/, const std::vector<BoundaryEdgeInfo>& recvEdgeBuffer)
  {
    for (unsigned j = 0; j < recvEdgeBuffer.size(); ++j) {
      BoundaryEdgeInfo edgeInfo = recvEdgeBuffer[j];
      m_boundaryMesh->create_edge(m_sortedVerticesByOwner.get_value(edgeInfo.vertex1Owner),
                                m_sortedVerticesByOwner.get_value(edgeInfo.vertex2Owner));
    }
  };

  boundaryVertexExchanger.complete_receives(unpackVerts);
  boundaryEdgeExchanger.complete_receives(unpackEdges);
}


void MeshExchangeBoundaryEdges::send_boundary_edges_to_root(EdgeExchange& boundaryEdgeExchanger)
{
  auto& sendEdgeBuffer = boundaryEdgeExchanger.get_send_buf(m_root);

  for (auto edge : m_localBoundaryEdges) {
    assert(edge->count_down() == 2); 

    BoundaryEdgeInfo edgeInfo;
    edgeInfo.vertex1Owner = get_owner_remote(m_inputMesh, edge->get_down(0));
    edgeInfo.vertex2Owner = get_owner_remote(m_inputMesh, edge->get_down(1));

    edgeInfo.vertex1Owner.remoteRank = utils::impl::get_rank_on_other_comm(m_inputMesh->get_comm(), m_unionComm, 
                                                                           edgeInfo.vertex1Owner.remoteRank);
    edgeInfo.vertex2Owner.remoteRank = utils::impl::get_rank_on_other_comm(m_inputMesh->get_comm(), m_unionComm, 
                                                                           edgeInfo.vertex2Owner.remoteRank);

    sendEdgeBuffer.push_back(edgeInfo);
  }
}


void MeshExchangeBoundaryEdges::send_boundary_vertices_to_root(VertexExchange& boundaryVertexExchanger)
{
  auto& sendVertexBuffer = boundaryVertexExchanger.get_send_buf(m_root);

  for (auto vertex : m_localBoundaryVertices) {
    assert(check_is_entity_owner(m_inputMesh, vertex));
    assert(vertex->count_remote_shared_entities() <= 1);

    BoundaryVertexInfo vertexInfo;
    vertexInfo.owner = get_owner_remote(m_inputMesh, vertex);
    vertexInfo.owner.remoteRank = utils::impl::get_rank_on_other_comm(m_inputMesh->get_comm(), m_unionComm, vertexInfo.owner.remoteRank);
    vertexInfo.pts = vertex->get_point_orig(0);

    sendVertexBuffer.push_back(vertexInfo);
  }
}

void MeshExchangeBoundaryEdges::recv_boundary_edges_and_vertices_in_root(VertexExchange& boundaryVertexExchanger,
                                                                         EdgeExchange& boundaryEdgeExchanger)
{
  const int nprocs = utils::impl::comm_size(m_unionComm);

  for (int i = 0; i < nprocs; ++i) {
    auto& entityCountPerProc = m_entityCountExchanger.get_recv_buf(i);
    int boundaryEdgeCount = entityCountPerProc[0];
    int boundaryVertexCount = entityCountPerProc[1];
    boundaryEdgeExchanger.get_recv_buf(i).resize(boundaryEdgeCount);
    boundaryVertexExchanger.get_recv_buf(i).resize(boundaryVertexCount);
  }
}

void MeshExchangeBoundaryEdges::gather_boundary_edge_and_vertex_count_info()
{
  const int nprocs = utils::impl::comm_size(m_unionComm);

  if (m_myrank == m_root) {
    for (int i = 0; i < nprocs; ++i) {
      m_entityCountExchanger.get_recv_buf(i).resize(2);
    }
  }

  auto& entityCountVector = m_entityCountExchanger.get_send_buf(m_root);
  entityCountVector.push_back( (m_inputMesh != nullptr) ? m_localBoundaryEdges.size() : 0 );
  entityCountVector.push_back( (m_inputMesh != nullptr) ? m_localBoundaryVertices.size() : 0 );
}

void MeshExchangeBoundaryEdges::gather_boundary_edges_and_vertices_to_root(VertexExchange& boundaryVertexExchanger, 
                                                                           EdgeExchange& boundaryEdgeExchanger)
{
  if (m_myrank == m_root)
    recv_boundary_edges_and_vertices_in_root(boundaryVertexExchanger, boundaryEdgeExchanger);

  send_boundary_vertices_to_root(boundaryVertexExchanger);
  send_boundary_edges_to_root(boundaryEdgeExchanger);

  boundaryEdgeExchanger.start_nonblocking();
  boundaryVertexExchanger.start_nonblocking();
}

void MeshExchangeBoundaryEdges::gather_boundary_entity_counts_to_root()
{
  gather_boundary_edge_and_vertex_count_info();

  m_entityCountExchanger.start_nonblocking();

  auto f = [&](int /*rank*/, const std::vector<int>& /*buf*/) {};
  m_entityCountExchanger.complete_receives(f);
}

void MeshExchangeBoundaryEdges::fill_local_boundary_edges()
{
  if (m_inputMesh != nullptr)
    for (auto& edge : m_inputMesh->get_edges())
      if (edge && is_boundary_edge(edge)) {
        m_localBoundaryEdges.push_back(edge);
      }
}

void MeshExchangeBoundaryEdges::fill_boundary_vertices()
{
  if (m_inputMesh != nullptr)
    for (auto& vertex : m_inputMesh->get_vertices()) {
       if (vertex && is_boundary_vertex(vertex))
       {
          if (check_is_entity_owner(m_inputMesh, vertex))
            m_localBoundaryVertices.push_back(vertex);
          else
            m_nonOwnedBoundaryVertices.push_back(vertex);
       }
    }
}


void MeshExchangeBoundaryEdges::set_boundary_vertex_field()
{
  m_inputMeshIsBoundaryNode = create_field<Bool>(m_inputMesh, FieldShape(1, 0, 0), 1, false);
  auto& field = *m_inputMeshIsBoundaryNode;
  for (auto& vertex : m_inputMesh->get_vertices())
    if (vertex)
      field(vertex, 0, 0) = is_local_boundary_vertex(vertex);

  
  ReductionOpLOR<Bool> op;
  FieldSharedReduction<Bool> reducer(m_inputMeshIsBoundaryNode, op);
  reducer.reduce();

}


bool MeshExchangeBoundaryEdges::is_local_boundary_vertex(MeshEntityPtr vertex)
{
  for (int i=0; i < vertex->count_up(); ++i)
  {
    if (is_boundary_edge(vertex->get_up(i))) {
      return true;
    }
  }

  return false;
}

bool MeshExchangeBoundaryEdges::is_boundary_vertex(MeshEntityPtr vertex)
{
  return (*m_inputMeshIsBoundaryNode)(vertex, 0, 0);
}

bool MeshExchangeBoundaryEdges::is_boundary_edge(MeshEntityPtr edge)
{
  return (edge && edge->count_up() == 1 && edge->count_remote_shared_entities() == 0);
}

std::shared_ptr<Mesh> MeshExchangeBoundaryEdges::get_boundary_edge_mesh()
{
  VertexExchange boundaryVertexExchanger(m_unionComm);
  EdgeExchange boundaryEdgeExchanger(m_unionComm);
  fill_boundary_vertices();
  fill_local_boundary_edges();

  if (utils::impl::comm_size(m_unionComm) > 1) {
    gather_boundary_entity_counts_to_root();
    gather_boundary_edges_and_vertices_to_root(boundaryVertexExchanger, boundaryEdgeExchanger);
  }

  if (m_myrank == m_root) {
    for (int i=0; i > utils::impl::comm_size(m_unionComm); ++i)
    {
      boundaryVertexExchanger.get_recv_buf(i).resize(m_entityCountExchanger.get_recv_buf(i)[0]);
      boundaryEdgeExchanger.get_recv_buf(i).resize(m_entityCountExchanger.get_recv_buf(i)[1]);
    }
    add_boundary_edges_on_root(boundaryVertexExchanger, boundaryEdgeExchanger);
  }

  return m_boundaryMesh;
}

}
}
}
}
