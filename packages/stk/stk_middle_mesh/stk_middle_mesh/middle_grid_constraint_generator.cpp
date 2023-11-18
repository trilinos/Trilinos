#include "middle_grid_constraint_generator.hpp"
#include "mesh.hpp"
#include "mesh_entity.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"
#include "variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void MiddleGridConstraintGenerator::generate()
{
  create_mesh1_vertices();
  create_mesh2_interior_vertices();
  create_mesh1_edges();
  split_edges();
  create_internal_edges();
}

void MiddleGridConstraintGenerator::create_mesh1_vertices()
{
  if (m_output)
    std::cout << "creating mesh1 vertices" << std::endl;

  auto& verts1ToFakeVerts   = *(m_relationalData->verts1ToFakeVerts);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);

  std::vector<mesh::MeshEntityPtr> els1;
  ExchangerKnown exchanger(m_mesh1->get_comm());

  for (auto& vert1 : m_mesh1->get_vertices())
    if (vert1)
    {
      FakeVert fv = verts1ToFakeVerts(vert1, 0, 0);
      auto vertIn = m_meshIn->create_vertex(fv.pt);

      if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
        std::cout << "created vert_in id " << vertIn->get_id() << " from fakevert id " << fv.id << std::endl;

      mesh::MeshEntityPtr el1 = vert1->get_up(0)->get_up(0);
      int localId             = predicates::impl::get_entity_id(el1, vert1);
      vertsInClassOnMesh1(vertIn, 0, 0) = m_pointClassifier->create_vert_record(el1, localId);
      fakeVertsToVertsIn[fv.id] = vertIn;

      int nels = get_upward(vert1, 2, els1);
      for (int i = 0; i < nels; ++i)
      {
        mesh1ElsToVertsIn.insert(els1[i], 0, vertIn);
        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "adding vert_in id " << vertIn->get_id() << " to el1 with id " << els1[i]->get_id() << std::endl;
      }

      for (int i=0; i < vert1->count_remote_shared_entities(); ++i)
      {
        const mesh::RemoteSharedEntity& remote = vert1->get_remote_shared_entity(i);
        exchanger.get_send_buf(remote.remoteRank).push_back(remote.remoteId);
        exchanger.get_send_buf(remote.remoteRank).push_back(vertIn->get_id());

        exchanger.get_recv_buf(remote.remoteRank).push_back(-1);
        exchanger.get_recv_buf(remote.remoteRank).push_back(-1);        
      }
    }

  set_mesh1_vert_shared_entities(exchanger);
}

void MiddleGridConstraintGenerator::set_mesh1_vert_shared_entities(ExchangerKnown& exchanger)
{
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);
  auto& fakeVertsToVertsIn = m_relationalData->fakeVertsToVertsIn;

  exchanger.start_nonblocking();

  auto unpacker = [&](int rank, const std::vector<int>& buf)
  {
    assert(buf.size() % 2 == 0);
    for (size_t i=0; i < buf.size(); i += 2)
    {
      int vert1Id = buf[i];
      int vertInRemoteId = buf[i+1];

      mesh::MeshEntityPtr vert1 = m_mesh1->get_vertices()[vert1Id];
      FakeVert fv = verts1ToFakeVerts(vert1, 0, 0);
      mesh::MeshEntityPtr vertIn = fakeVertsToVertsIn[fv.id];

      vertIn->add_remote_shared_entity(mesh::RemoteSharedEntity(rank, vertInRemoteId));
    }
  };

  exchanger.complete_receives(unpacker);
}

void MiddleGridConstraintGenerator::create_mesh2_interior_vertices()
{
  if (m_output)
    std::cout << "creating mesh2 interior vertices" << std::endl;
  auto& verts2ClassOnMesh1  = *(m_relationalData->verts2ClassOnMesh1);
  auto& verts2ToFakeVerts   = *(m_relationalData->verts2ToFakeVerts);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);

  for (auto& vert2 : m_mesh2->get_vertices())
    if (vert2)
    {
      predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert2, 0, 0);
      if (record.type == predicates::impl::PointClassification::Interior)
      {
        FakeVert fv = verts2ToFakeVerts(vert2, 0, 0);
        auto vertIn = m_meshIn->create_vertex(fv.pt);

        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "created vert_in id " << vertIn->get_id() << " from fakevert id " << fv.id << std::endl;

        vertsInClassOnMesh1(vertIn, 0, 0) = record;         
        fakeVertsToVertsIn[fv.id]         = vertIn;
        mesh1ElsToVertsIn.insert(record.el, 0, vertIn);
        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "adding vert_in id " << vertIn->get_id() << " to el1 with id " << record.el->get_id()
                    << std::endl;
      }
    }
}

void MiddleGridConstraintGenerator::create_mesh1_edges()
{
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);
  ExchangerKnown exchanger(m_mesh1->get_comm());
  for (auto& edge : m_mesh1->get_edges())
    if (edge)
    {
      FakeVert fv1             = verts1ToFakeVerts(edge->get_down(0), 0, 0);
      FakeVert fv2             = verts1ToFakeVerts(edge->get_down(1), 0, 0);
      mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
      mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];
      assert(!stk::middle_mesh::mesh::get_common_edge(v1In, v2In));

      mesh::MeshEntityPtr edgeIn = m_meshIn->create_edge(v1In, v2In);

      for (int i=0; i < edge->count_remote_shared_entities(); ++i)
      {
        const mesh::RemoteSharedEntity& edge1Remote = edge->get_remote_shared_entity(i);

        auto& sendBuf = exchanger.get_send_buf(edge1Remote.remoteRank);
        auto& recvBuf = exchanger.get_recv_buf(edge1Remote.remoteRank);
        mesh::RemoteSharedEntity v1InRemote = mesh::get_remote_shared_entity(v1In, edge1Remote.remoteRank);
        mesh::RemoteSharedEntity v2InRemote = mesh::get_remote_shared_entity(v2In, edge1Remote.remoteRank);

        sendBuf.push_back(v1InRemote.remoteId);  recvBuf.push_back(-1);
        sendBuf.push_back(v2InRemote.remoteId);  recvBuf.push_back(-1);
        sendBuf.push_back(edgeIn->get_id());     recvBuf.push_back(-1);
      }
    }

  set_mesh1_edge_shared_entities(exchanger);
}

void MiddleGridConstraintGenerator::set_mesh1_edge_shared_entities(ExchangerKnown& exchanger)
{
  exchanger.start_nonblocking();

  auto unpacker = [&](int rank, const std::vector<int>& buf)
  {
    assert(buf.size() % 3 == 0);
    for (size_t i=0; i < buf.size(); i += 3)
    {
      int vert1Id       = buf[i];
      int vert2Id       = buf[i+1];
      int edgeRemoteId  = buf[i+2];

      mesh::MeshEntityPtr vert1 = m_meshIn->get_vertices()[vert1Id];
      mesh::MeshEntityPtr vert2 = m_meshIn->get_vertices()[vert2Id];
      mesh::MeshEntityPtr edge  = mesh::get_common_edge(vert1, vert2);
      edge->add_remote_shared_entity(mesh::RemoteSharedEntity(rank, edgeRemoteId));
    }
  };

  exchanger.complete_receives(unpacker);
}

void MiddleGridConstraintGenerator::split_edges()
{
  if (m_output)
    std::cout << "splitting edges" << std::endl;
  auto& mesh1EdgesToSplit   = *(m_relationalData->mesh1EdgesToSplit);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);
  ExchangerKnown exchanger(m_mesh1->get_comm());

  auto sharedEntityInfoPtr = mesh::create_variable_size_field<int>(m_mesh1, mesh::FieldShape(0, 2, 0));

  for (auto& edge1 : m_mesh1->get_edges())
  {
    if (edge1 && mesh1EdgesToSplit.get_num_comp(edge1, 0) > 0)
    {
      sort_edge_splits(edge1);
      bool edge1HasRemote = edge1->count_remote_shared_entities() > 0;

      mesh::MeshEntityPtr el1 = edge1->get_up(0);
      int localId             = predicates::impl::get_entity_id(el1, edge1);

      auto currEdgeIn = get_mesh_in_edge_from_mesh1_edge(edge1);
      double xiStart  = 0;
      mesh::MeshEntityPtr newEdges[2];

      for (int i = 0; i < mesh1EdgesToSplit.get_num_comp(edge1, 0); ++i)
      {
        double xi   = mesh1EdgesToSplit(edge1, 0, i).xi;
        FakeVert fv = mesh1EdgesToSplit(edge1, 0, i).vert;

        // map xi from edge1 to curr_edge_in
        double xiNew = (xi - xiStart) / (1 - xiStart);
        auto newVert = m_meshIn->split_edge_broken(currEdgeIn, xiNew, newEdges);

        if (m_output && (m_vertIds.count(newVert->get_id()) > 0))
          std::cout << "created vert_in id " << newVert->get_id() << " from fakevert id " << fv.id << std::endl;


        fakeVertsToVertsIn[fv.id] = newVert;

        double edgeXiOnReferenceEl = el1->get_down_orientation(localId) == mesh::EntityOrientation::Standard ? xi : 1 - xi;
        vertsInClassOnMesh1(newVert, 0, 0) = m_pointClassifier->create_edge_record(el1, localId, edgeXiOnReferenceEl);
        for (int j = 0; j < edge1->count_up(); ++j)
        {
          mesh1ElsToVertsIn.insert(edge1->get_up(j), 0, newVert);
          if (m_output && m_vertIds.count(newVert->get_id()) > 0)
            std::cout << "adding vert_in id " << newVert->get_id() << " to el1 with id " << edge1->get_up(j)->get_id()
                      << std::endl;
        }

        if (edge1HasRemote)
        {
          sharedEntityInfoPtr->insert(edge1, 0, newVert->get_id());
          sharedEntityInfoPtr->insert(edge1, 1, newEdges[0]->get_id());
        }

        // update for next iteration
        xiStart    = xi;
        currEdgeIn = newEdges[1];
      }

      if (edge1HasRemote)
      {
        sharedEntityInfoPtr->insert(edge1, 1, currEdgeIn->get_id());
        pack_edge_split_shared_info(exchanger, sharedEntityInfoPtr, edge1);
      }
    }
  }

  set_edge_split_shared_entities(exchanger, sharedEntityInfoPtr);
}

void MiddleGridConstraintGenerator::pack_edge_split_shared_info(ExchangerKnown& exchanger, mesh::VariableSizeFieldPtr<int> sharedEntityInfoPtr, mesh::MeshEntityPtr edge1)
{
  auto& sharedEntityInfo = *sharedEntityInfoPtr;
  for (int i=0; i < edge1->count_remote_shared_entities(); ++i)
  {
    const mesh::RemoteSharedEntity& remote = edge1->get_remote_shared_entity(i);

    auto& sendBuf = exchanger.get_send_buf(remote.remoteRank);
    auto& recvBuf = exchanger.get_recv_buf(remote.remoteRank);

    sendBuf.push_back(remote.remoteId); recvBuf.push_back(-1);
    for (int dim=0; dim < 2; ++dim)
    {
      for (int entityId : sharedEntityInfo(edge1, dim))
      {
        sendBuf.push_back(entityId);
        recvBuf.push_back(-1);
      }
    }
  }
}

void MiddleGridConstraintGenerator::set_edge_split_shared_entities(ExchangerKnown& exchanger, mesh::VariableSizeFieldPtr<int> sharedEntityInfoPtr)
{
  auto& sharedEntityInfo = *sharedEntityInfoPtr;
  exchanger.start_nonblocking();

  auto unpacker = [&](int rank, const std::vector<int>& buf)
  {
    size_t idx = 0;
    while (idx < buf.size())
    {
      int edge1Id = buf[idx++];
      mesh::MeshEntityPtr edge1 = m_mesh1->get_edges()[edge1Id];
      for (int dim=0; dim < 2; ++dim)
      {
        for (int i=0; i < sharedEntityInfo.get_num_comp(edge1, dim); ++i)
        {
          assert(idx < buf.size());
          int remoteEntityId = buf[idx++];
          int localEntityId = sharedEntityInfo(edge1, dim, i);

          mesh::MeshEntityPtr localVert = m_meshIn->get_mesh_entities(dim)[localEntityId];
          localVert->add_remote_shared_entity(mesh::RemoteSharedEntity(rank, remoteEntityId));
        }
      }
    }
  };

  exchanger.complete_receives(unpacker);
}

void MiddleGridConstraintGenerator::sort_edge_splits(mesh::MeshEntityPtr edge1)
{
  auto& mesh1EdgesToSplit = *(m_relationalData->mesh1EdgesToSplit);
  auto cmp = [](const EdgeSplitRecord& split1, const EdgeSplitRecord& split2) { return split1.xi < split2.xi; };

  int nsplits                  = mesh1EdgesToSplit.get_num_comp(edge1, 0);
  EdgeSplitRecord* splitsStart = &(mesh1EdgesToSplit(edge1, 0, 0));
  EdgeSplitRecord* splitsEnd   = &(mesh1EdgesToSplit(edge1, 0, nsplits - 1)) + 1;
  std::sort(splitsStart, splitsEnd, cmp);
}

mesh::MeshEntityPtr MiddleGridConstraintGenerator::get_mesh_in_edge_from_mesh1_edge(mesh::MeshEntityPtr edge1)
{
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);

  FakeVert fv1             = verts1ToFakeVerts(edge1->get_down(0), 0, 0);
  FakeVert fv2             = verts1ToFakeVerts(edge1->get_down(1), 0, 0);
  mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
  mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];

  return get_common_edge(v1In, v2In);
}

void MiddleGridConstraintGenerator::create_internal_edges()
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);

  for (auto& edge2 : m_mesh2->get_edges())
    if (edge2)
    {
      assert(edges2ToFakeVertsIn.get_num_comp(edge2, 0) >= 2);
      // sortVertsOnEdge2(edge2);  //TODO: make this an assert?

      for (int i = 0; i < edges2ToFakeVertsIn.get_num_comp(edge2, 0) - 1; ++i)
      {
        FakeVert fv1             = edges2ToFakeVertsIn(edge2, 0, i).vert;
        FakeVert fv2             = edges2ToFakeVertsIn(edge2, 0, i + 1).vert;
        mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
        mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];

        if (m_output && (v1In->get_id() == 1916 || v2In->get_id() == 1916))
        {
          std::cout << "found possible edge involving vert 1916" << std::endl;
          std::cout << "vert ids = " << v1In->get_id() << ", " << v2In->get_id() << std::endl;

          for (int j = 0; j < edges2ToFakeVertsIn.get_num_comp(edge2, 0); ++j)
            std::cout << "fake vert " << j << " on specified edge has id " << edges2ToFakeVertsIn(edge2, 0, j).vert.id
                      << ", xi = " << edges2ToFakeVertsIn(edge2, 0, j).xi << std::endl;
        }

        if (v1In && v2In && !get_common_edge(v1In, v2In))
        {
          m_meshIn->create_edge(v1In, v2In);
        }
      }
    }
}

void MiddleGridConstraintGenerator::sort_verts_on_edge2(mesh::MeshEntityPtr edge2)
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto cmp                  = [](const VertOnEdge& v1, const VertOnEdge& v2) { return v1.xi < v2.xi; };

  int nsplits            = edges2ToFakeVertsIn.get_num_comp(edge2, 0);
  VertOnEdge* vertsStart = &(edges2ToFakeVertsIn(edge2, 0, 0));
  VertOnEdge* vertsEnd   = &(edges2ToFakeVertsIn(edge2, 0, nsplits - 1)) + 1;
  std::sort(vertsStart, vertsEnd, cmp);
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
