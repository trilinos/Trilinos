#include "stk_middle_mesh/mesh_relational_data_scatter.hpp"
#include "mesh_entity.hpp"
#include "predicates/intersection_common.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

namespace {
using predicates::impl::PointRecord;
using predicates::impl::PointRecordForTriangle;
using predicates::impl::PointClassification;
}


MeshRelationalDataScatter::MeshRelationalDataScatter(
                            const MeshRelationalDataScatterInput& scatterInput,     
                            std::shared_ptr<MeshRelationalData> meshRelationalData,
                            std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifierInput,
                            std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifierOutput,                 
                            MPI_Comm unionComm) :
  m_mesh1(scatterInput.mesh1),
  m_mesh2(scatterInput.mesh2),
  m_mesh1ScatteredToMesh2(scatterInput.mesh1ScatteredToMesh2),
  m_mesh2ScatteredToMesh1(scatterInput.mesh2ScatteredToMesh1),
  m_mesh1EntityOrigins(scatterInput.mesh1EntityOrigins),
  m_mesh1EntityDestinations(scatterInput.mesh1EntityDestinations),
  m_mesh2EntityOrigins(scatterInput.mesh2EntityOrigins),
  m_mesh2EntityDestinations(scatterInput.mesh2EntityDestinations),
  m_meshRelationDataInput(meshRelationalData),
  m_pointClassifierInput(pointClassifierInput),
  m_pointClassifierOutput(pointClassifierOutput),
  m_quadToTrianglesOrigin(pointClassifierInput ? &(pointClassifierInput->get_quad_to_triangles()) : nullptr),
  m_quadToTrianglesDest(pointClassifierOutput ? &(pointClassifierOutput->get_quad_to_triangles()) : nullptr),
  m_unionComm(unionComm)
{
  if (m_mesh1)
  {
    // TODO: a flexible CSR structure might be better for this
    m_senderFakeVertsToLocal.resize(utils::impl::comm_size(m_unionComm));

    m_meshIn = mesh::make_empty_mesh(m_mesh1->get_comm());
    m_meshRelationDataOutput = std::make_shared<MeshRelationalData>(m_mesh1, m_mesh2ScatteredToMesh1, m_meshIn);
  }
  check_fields_are_on_correct_meshes();
}

std::shared_ptr<MeshRelationalData> MeshRelationalDataScatter::scatter()
{
  {
    Exchanger mesh1VertExchanger(m_unionComm);
    pack_mesh1_vert_fields(mesh1VertExchanger);
    size_recv_buffer_mesh1_vert_fields(mesh1VertExchanger);
    unpack_mesh1_vert_fields(mesh1VertExchanger);
  }

  {
    ExchangerUnknown mesh2VertExchanger(m_unionComm);
    pack_mesh2_vert_fields(mesh2VertExchanger);
    unpack_mesh2_vert_fields(mesh2VertExchanger);
  }

  {
    ExchangerUnknown mesh1EdgeExchanger(m_unionComm);
    pack_mesh1_edge_fields(mesh1EdgeExchanger);
    unpack_mesh1_edge_fields(mesh1EdgeExchanger);
  }

  {
    ExchangerUnknown mesh2EdgeExchanger(m_unionComm);
    pack_mesh2_edge_fields(mesh2EdgeExchanger);
    unpack_mesh2_edge_fields(mesh2EdgeExchanger);
  }

  if (m_meshRelationDataOutput)
  {
    m_meshRelationDataOutput->fakeVertsToVertsIn.resize(m_fakeVertGenerator.get_num_verts());
  }

  return m_meshRelationDataOutput;
}


void MeshRelationalDataScatter::pack_mesh1_vert_fields(Exchanger& exchanger)
{
  if (m_mesh1ScatteredToMesh2)
  {
    for (int phase=0; phase < 2; ++phase)
    {
      auto& verts1ToFakeVerts = *(m_meshRelationDataInput->verts1ToFakeVerts);
      auto& dests             = *m_mesh1EntityOrigins;
      for (auto& vert : m_mesh1ScatteredToMesh2->get_vertices())
        if (vert)
        {
          FakeVert& fakeVert = verts1ToFakeVerts(vert, 0, 0);
          for (mesh::RemoteSharedEntity remote : dests(vert, 0))
          {
            exchanger.get_send_buf(remote.remoteRank).pack(remote.remoteId);
            exchanger.get_send_buf(remote.remoteRank).pack(fakeVert.id);
            exchanger.get_send_buf(remote.remoteRank).pack(fakeVert.pt);
          }              
        }

      if (phase == 0)
        exchanger.allocate_send_buffers();
    }
  }
}

void MeshRelationalDataScatter::size_recv_buffer_mesh1_vert_fields(Exchanger& exchanger)
{
  if (m_mesh1)
  {
    auto& senders = *m_mesh1EntityDestinations;
    for (auto& vert : m_mesh1->get_vertices())
      if (vert)
      {
        for (mesh::RemoteSharedEntity sender : senders(vert, 0))
        {
          int id = -1;
          utils::Point pt;
          exchanger.get_recv_buf(sender.remoteRank).pack(id);
          exchanger.get_recv_buf(sender.remoteRank).pack(id);
          exchanger.get_recv_buf(sender.remoteRank).pack(pt);
        }
      }

    for (int i=0; i < utils::impl::comm_size(m_unionComm); ++i)
      exchanger.set_recv_buffer_size(i, exchanger.get_recv_buf(i).size());
  }
    
  exchanger.allocate_recv_buffers();

}

void MeshRelationalDataScatter::unpack_mesh1_vert_fields(Exchanger& exchanger)
{
  exchanger.start_nonblocking();

  // don't unpack data as it comes in: this would give the 
  // fakeVertIds a non-deterministic order
  auto f = [](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  exchanger.complete_receives(f);

  if (m_mesh1)
  {
    auto& verts1ToFakeVerts = *(m_meshRelationDataOutput->verts1ToFakeVerts);

    using Bool = int_least8_t;
    auto seenVertsPtr = mesh::create_field<Bool>(m_mesh1, FieldShape(1, 0, 0), 1, false);
    auto& seenVerts = *seenVertsPtr;
    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
    {
      auto& buf = exchanger.get_recv_buf(rank);
      while (buf.remaining() > 0)
      {
        int localVertId, senderFakeVertId;
        utils::Point pt;
        buf.unpack(localVertId);
        buf.unpack(senderFakeVertId);
        buf.unpack(pt);

        mesh::MeshEntityPtr vert = m_mesh1->get_vertices()[localVertId];

        FakeVert fv;
        if (seenVerts(vert, 0, 0))
        {
          fv = verts1ToFakeVerts(vert, 0, 0);
        } else
        {
          fv = m_fakeVertGenerator.get_vert(pt);
          verts1ToFakeVerts(vert, 0, 0) = fv;
          seenVerts(vert, 0, 0) = true;
        }

        m_senderFakeVertsToLocal[rank].insert(std::make_pair(senderFakeVertId, fv));
      }
    }
  }
}


void MeshRelationalDataScatter::pack_mesh2_vert_fields(ExchangerUnknown& exchanger)
{
  if (m_mesh2)
  {
    auto& verts2ToFakeVerts = *(m_meshRelationDataInput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(m_meshRelationDataInput->verts2ClassOnMesh1);
    auto& dests = *m_mesh2EntityDestinations;
    for (int phase=0; phase < 2; ++phase)
    {
      for (auto& vert : m_mesh2->get_vertices())
        if (vert)
        {
          for (mesh::RemoteSharedEntity dest : dests(vert, 0))
          {
            auto& buf = exchanger.get_send_buf(dest.remoteRank);
            FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
            const PointRecord& pointRecord = verts2ClassOnMesh1(vert, 0, 0);

            buf.pack(dest.remoteId);
            buf.pack(fv.id);
            buf.pack(fv.pt);
            pack(buf, dest.remoteRank, pointRecord);
          }
        }

      if (phase == 0)
        exchanger.allocate_send_buffers();
    }
  }
}


void MeshRelationalDataScatter::pack(stk::CommBuffer& buf, 
                                     int destRank, const PointRecord& record)
{

  assert(record.type != predicates::impl::PointClassification::Exterior);
  buf.pack(record.type);
  if (record.type == predicates::impl::PointClassification::Vert)
  {
    mesh::MeshEntityPtr vertLocal = predicates::impl::get_entity(record);
    int vertIdOnDest = get_id_on_dest(vertLocal, destRank);
    buf.pack(vertIdOnDest);
  } else if (record.type == predicates::impl::PointClassification::Edge)
  {
    mesh::MeshEntityPtr edgeLocal = predicates::impl::get_entity(record);
    int edgeIdOnDest = get_id_on_dest(edgeLocal, destRank); 
    buf.pack(edgeIdOnDest);
    if (edgeIdOnDest == -1)
    {
      return;
    }    
    buf.pack(m_pointClassifierInput->get_edge_xi(record));
  } else if (record.type == predicates::impl::PointClassification::Interior)
  {
    int elIdOnDest = get_id_on_dest(record.el, destRank);

    buf.pack(elIdOnDest);
    if (elIdOnDest == -1)
    {
      return;
    }       
    pack(buf, record.m_r1);
    pack(buf, record.m_r2);
  }
}


int MeshRelationalDataScatter::get_id_on_dest(mesh::MeshEntityPtr entity1, int destRank)
{
  auto& mesh1EntityOrigins = *m_mesh1EntityOrigins;
  for (mesh::RemoteSharedEntity remote : mesh1EntityOrigins(entity1, 0))
  {
    if (remote.remoteRank == destRank)
    {
      return remote.remoteId;
    }
  }

  return -1;
}

void MeshRelationalDataScatter::pack(stk::CommBuffer& buf, 
                                     const PointRecordForTriangle& record)
{
  int elId;
  if (record.el == m_quadToTrianglesOrigin->el1)
  {
    elId = 0;
  } else if (record.el == m_quadToTrianglesOrigin->el2)
  {
    elId = 1;
  } else if (record.el == nullptr)
  {
    elId = 2;
  } else
  {
    elId = 3;
  }

  buf.pack(record.type);
  buf.pack(record.id);
  buf.pack(elId);
  buf.pack(record.m_ptXi);
}



PointRecord MeshRelationalDataScatter::unpack_point_record(stk::CommBuffer& buf)
{
  predicates::impl::PointClassification type;  
  int entityId;

  buf.unpack(type);
  buf.unpack(entityId);

  // elId can be -1 in the partial overlap case: if a mesh2 element spans two
  // (or more) mesh1 elements that are on different procs, the mesh2 element
  // will be set to both procs, but some of its verts will be classified on
  // elements not on one of the destination procs
  // Example: 
  // |-------/--------|
  // |       /        |
  // |  P0   /    P1  | 
  // |-------/--------|
  // If the slash marks are the inter-processor boundary of mesh1,
  // the top right vertex will be classified on a P1 element, which
  // is not present on P0.
  if (entityId == -1)
  {
    return PointRecord();
  }

  if (type == predicates::impl::PointClassification::Vert)
  {
    mesh::MeshEntityPtr vert = m_mesh1->get_vertices()[entityId];
    mesh::MeshEntityPtr el = vert->get_up(0)->get_up(0);
    int vertIdOnEl = mesh::get_local_id(el, vert);

    return m_pointClassifierOutput->create_vert_record(el, vertIdOnEl);
  } else if (type == predicates::impl::PointClassification::Edge)
  {
    double edgeXi;
    buf.unpack(edgeXi);

    mesh::MeshEntityPtr edge = m_mesh1->get_edges()[entityId];
    mesh::MeshEntityPtr el   = edge->get_up(0);
    int edgeIdOnEl = mesh::get_local_id(el, edge);
    double edgeXiOnReferenceEl = el->get_down_orientation(edgeIdOnEl) == mesh::EntityOrientation::Standard ? edgeXi : 1 - edgeXi;

    return m_pointClassifierOutput->create_edge_record(el, edgeIdOnEl, edgeXiOnReferenceEl);
  } else if (type == predicates::impl::PointClassification::Interior)
  {
    mesh::MeshEntityPtr el = m_mesh1->get_elements()[entityId];

    PointRecordForTriangle r1 = unpack_point_record_for_triangle(buf, el);
    PointRecordForTriangle r2 = unpack_point_record_for_triangle(buf, el);

    return PointRecord(type, 0, el, r1, r2);
  } else
  {
    throw std::runtime_error("type cannot be Exterior");
  }
}

PointRecordForTriangle MeshRelationalDataScatter::unpack_point_record_for_triangle(
                         stk::CommBuffer& buf, mesh::MeshEntityPtr parentEl)
{
  predicates::impl::PointClassification type;  
  int entityId, elId;
  utils::Point ptXi;

  buf.unpack(type);
  buf.unpack(entityId);
  buf.unpack(elId);
  buf.unpack(ptXi);

  mesh::MeshEntityPtr el;
  switch (elId)
  {
    case 0 : { el = m_quadToTrianglesDest->el1; break; }
    case 1 : { el = m_quadToTrianglesDest->el2; break; }
    case 2 : { el = nullptr; break; }
    case 3 : { el = parentEl; break; }
    default:
    {
      throw std::runtime_error("invalid value for elId");
    }
  }

  return PointRecordForTriangle(type, entityId, el, ptXi);
}

void MeshRelationalDataScatter::unpack_mesh2_vert_fields(ExchangerUnknown& exchanger)
{
  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();
  // Don't unpack data as it arrives because that would lead to non-deterministic
  // assignment of FakeVert ids
  auto f = [](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  exchanger.complete_receives(f);

  if (m_mesh2ScatteredToMesh1)
  {
    auto& verts2ToFakeVerts = *(m_meshRelationDataOutput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(m_meshRelationDataOutput->verts2ClassOnMesh1);
    verts2ToFakeVerts.set(FakeVert{-1, {0, 0, 0}});

    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
    {
      stk::CommBuffer& buf = exchanger.get_recv_buf(rank);
      while (buf.remaining() > 0)
      {
        int vertLocalId, senderFakeVertId;
        utils::Point pt;
        buf.unpack(vertLocalId);
        buf.unpack(senderFakeVertId);
        buf.unpack(pt);
        predicates::impl::PointRecord record = unpack_point_record(buf);

        mesh::MeshEntityPtr vert = m_mesh2ScatteredToMesh1->get_vertices()[vertLocalId];
        FakeVert fv = verts2ToFakeVerts(vert, 0, 0);

        if (fv.id == -1)
        {
          bool fakeVertExists = m_senderFakeVertsToLocal[rank].count(senderFakeVertId);
          if (fakeVertExists)
          {
            fv = m_senderFakeVertsToLocal[rank][senderFakeVertId];
          } else{
            fv = m_fakeVertGenerator.get_vert(pt);
          }
        }

        m_senderFakeVertsToLocal[rank].insert(std::make_pair(senderFakeVertId, fv));
        verts2ToFakeVerts(vert, 0, 0)  = fv;
        verts2ClassOnMesh1(vert, 0, 0) = record;
      }
    }
  }
}

void MeshRelationalDataScatter::pack_mesh1_edge_fields(ExchangerUnknown& exchanger)
{
  if (m_mesh1ScatteredToMesh2)
  {
    auto& mesh1Dests = *m_mesh1EntityOrigins;
    std::vector<int> ownedSplits;
    for (int phase=0; phase < 2; ++phase)
    {
      for (auto& edge : m_mesh1ScatteredToMesh2->get_edges())
      {
        if (edge)
        {
          get_owned_splits(edge, ownedSplits);
          if (ownedSplits.size() > 0)
          {
            for (mesh::RemoteSharedEntity dest : mesh1Dests(edge, 0))
            {
              stk::CommBuffer& buf = exchanger.get_send_buf(dest.remoteRank);
              pack_edge_splits(buf, dest, edge, ownedSplits);
            }
          }              
        }
      }

      if (phase == 0)
        exchanger.allocate_send_buffers();
    }
  }  
}

void MeshRelationalDataScatter::get_owned_splits(mesh::MeshEntityPtr edge, std::vector<int>& ownedSplits)
{
  int myrank = utils::impl::comm_rank(m_mesh1ScatteredToMesh2->get_comm());
  auto& mesh1EdgesToSplit = *(m_meshRelationDataInput->mesh1EdgesToSplit);

  ownedSplits.resize(0);
  int nsplits = mesh1EdgesToSplit.get_num_comp(edge, 0);
  for (int i=0; i < nsplits; ++i)
  {
    mesh::MeshEntityPtr mesh1Entity = mesh1EdgesToSplit(edge, 0, i).otherMeshEntity;
    if (mesh::get_owner(m_mesh1ScatteredToMesh2, mesh1Entity) == myrank)
    {
      ownedSplits.push_back(i);
    }
  }
}

void MeshRelationalDataScatter::pack_edge_splits(stk::CommBuffer& buf, const mesh::RemoteSharedEntity& dest,
                                                 mesh::MeshEntityPtr edge, const std::vector<int> ownedSplits)
{
  auto& mesh1EdgesToSplit = *(m_meshRelationDataInput->mesh1EdgesToSplit);
  auto& mesh2Dests = *m_mesh2EntityDestinations;

  buf.pack(dest.remoteId);
  buf.pack(int(ownedSplits.size()));
  for (int i : ownedSplits)
  {
    EdgeSplitRecord& edgeSplit = mesh1EdgesToSplit(edge, 0, i);

    int otherEntityIdOnDest = -1;
    int otherEntityDim = mesh::get_type_dimension(edgeSplit.otherMeshEntity->get_type());
    for (auto& mesh2Dest : mesh2Dests(edgeSplit.otherMeshEntity, 0))
      if (mesh2Dest.remoteRank == dest.remoteRank)
      {
        otherEntityIdOnDest = mesh2Dest.remoteId;
        break;                          
      }
      
    buf.pack(edgeSplit.vert.id);
    buf.pack(otherEntityDim);
    buf.pack(otherEntityIdOnDest);
    buf.pack(edgeSplit.vert.pt);
    buf.pack(edgeSplit.xi);
  }
}

void MeshRelationalDataScatter::unpack_mesh1_edge_fields(ExchangerUnknown& exchanger)
{
  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();
  auto f = [](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  exchanger.complete_receives(f);

  if (m_mesh1)
  {
    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
    {
      stk::CommBuffer& buf = exchanger.get_recv_buf(rank);
      while (buf.remaining() > 0)
      {
        unpack_edge_split(buf, rank);
      }
    }
  }
}

void MeshRelationalDataScatter::unpack_edge_split(stk::CommBuffer& buf, int senderRank)
{      
  int edgeLocalId, nsplits;
  buf.unpack(edgeLocalId);
  buf.unpack(nsplits);

  mesh::MeshEntityPtr edge = m_mesh1->get_edges()[edgeLocalId];
  auto& mesh1EdgesToSplit = *(m_meshRelationDataOutput->mesh1EdgesToSplit);
  for (int i=0; i < nsplits; ++i)
  {
    int senderFakeVertId, otherEntityDim, otherEntityId;
    utils::Point pt;
    double xi;

    buf.unpack(senderFakeVertId);
    buf.unpack(otherEntityDim);
    buf.unpack(otherEntityId);
    buf.unpack(pt);
    buf.unpack(xi);

    bool vertExists = m_senderFakeVertsToLocal[senderRank].count(senderFakeVertId) > 0;
    FakeVert fv;
    if (vertExists)
    {
      fv = m_senderFakeVertsToLocal[senderRank][senderFakeVertId];
    } else
    {
      fv = m_fakeVertGenerator.get_vert(pt);
      m_senderFakeVertsToLocal[senderRank].insert(std::make_pair(senderFakeVertId, fv));
    }

    mesh::MeshEntityPtr otherEntity = m_mesh2ScatteredToMesh1->get_mesh_entities(otherEntityDim)[otherEntityId];
    mesh1EdgesToSplit.insert(edge, 0, EdgeSplitRecord{fv, xi, otherEntity});
  }      
}

void MeshRelationalDataScatter::pack_mesh2_edge_fields(ExchangerUnknown& exchanger)
{
  if (m_mesh2)
  {
    auto& edges2ToFakeVertsIn = *(m_meshRelationDataInput->edges2ToFakeVertsIn);
    auto& dests = *m_mesh2EntityDestinations;
    int myrank = utils::impl::comm_rank(m_mesh2->get_comm());
    for (int phase=0; phase < 2; ++phase)
    {
      for (auto& edge : m_mesh2->get_edges())
        if (edge)
        {
          if (mesh::get_owner(m_mesh2, edge) != myrank)
            continue;

          for (mesh::RemoteSharedEntity dest : dests(edge, 0))
          {
            stk::CommBuffer& buf = exchanger.get_send_buf(dest.remoteRank);
            int nverts = edges2ToFakeVertsIn.get_num_comp(edge, 0);

            buf.pack(dest.remoteId);
            buf.pack(nverts);
            for (VertOnEdge& vertOnEdge : edges2ToFakeVertsIn(edge, 0))
            {
              buf.pack(vertOnEdge);
            }
          }
        }

      if (phase == 0)
        exchanger.allocate_send_buffers();
    }
  }
}

void MeshRelationalDataScatter::unpack_mesh2_edge_fields(ExchangerUnknown& exchanger)
{
  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();
  auto f = [](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  exchanger.complete_receives(f);

  if (m_mesh2ScatteredToMesh1)
  {
    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
    {
      stk::CommBuffer& buf = exchanger.get_recv_buf(rank);

      while (buf.remaining() > 0)
      {
        unpack_vert_on_edge(buf, rank);
      }
    }
  }
}

void MeshRelationalDataScatter::unpack_vert_on_edge(stk::CommBuffer& buf, int senderRank)
{ 
  int edgeLocalId, nverts;
  buf.unpack(edgeLocalId);
  buf.unpack(nverts);

  mesh::MeshEntityPtr edge = m_mesh2ScatteredToMesh1->get_edges()[edgeLocalId];
  auto& edges2ToFakeVertsIn = *(m_meshRelationDataOutput->edges2ToFakeVertsIn);
  for (int i=0; i < nverts; ++i)
  {
    VertOnEdge vertOnEdge;
    buf.unpack(vertOnEdge);

    bool vertExists = m_senderFakeVertsToLocal[senderRank].count(vertOnEdge.vert.id) > 0;
    FakeVert fv;
    if (vertExists)
    {
      fv = m_senderFakeVertsToLocal[senderRank][vertOnEdge.vert.id];
    } else 
    {
      fv = m_fakeVertGenerator.get_vert(vertOnEdge.vert.pt);
      m_senderFakeVertsToLocal[senderRank].insert(std::make_pair(vertOnEdge.vert.id, fv));
    }

    vertOnEdge.vert = fv;
    edges2ToFakeVertsIn.insert(edge, 0, vertOnEdge);
  }      
}


void MeshRelationalDataScatter::check_fields_are_on_correct_meshes()
{
  if (m_mesh1)
  {
    STK_ThrowAssertMsg(m_quadToTrianglesDest, 
                        "quadToTrianglesDest must be provided on all procs where mesh1 is present");

    STK_ThrowAssertMsg(m_mesh1EntityDestinations->get_mesh() == m_mesh1,
                        "mesh1EntityDestinations must be on mesh1");
  }
  if (m_mesh1ScatteredToMesh2)
  {        
    STK_ThrowRequireMsg(m_meshRelationDataInput->verts1ToFakeVerts->get_mesh() == m_mesh1ScatteredToMesh2,
                        "verts1ToFakeVerts must be defined on mesh1ScatteredToMesh2"); 

    STK_ThrowRequireMsg(m_meshRelationDataInput->mesh1EdgesToSplit->get_mesh() == m_mesh1ScatteredToMesh2,
                        "mesh1EdgesToSplit must be defined on mesh1ScatteredToMesh2");  

    STK_ThrowRequireMsg(m_mesh1EntityOrigins->get_mesh() == m_mesh1ScatteredToMesh2,
                        "mesh1EntityOrigins must be defined on mesh1ScatteredToMesh2");
  }

  if (m_mesh2)
  {
    STK_ThrowRequireMsg(m_mesh1ScatteredToMesh2,
                        "mesh1ScatteredToMesh2 must be present on all procs where mesh2 is present");

    STK_ThrowAssertMsg(m_quadToTrianglesOrigin, 
                        "quadToTrianglesOrigin must be provided on all procs where mesh2 is present");

    STK_ThrowRequireMsg(m_meshRelationDataInput->verts2ToFakeVerts->get_mesh() == m_mesh2,
                        "verts2ToFakeVerts must be defined on mesh2");
            
    STK_ThrowRequireMsg(m_meshRelationDataInput->verts2ClassOnMesh1->get_mesh() == m_mesh2,
                        "verts2ClassOnMesh1 must be defined on mesh2");

    STK_ThrowRequireMsg(m_meshRelationDataInput->edges2ToFakeVertsIn->get_mesh() == m_mesh2,
                        "edges2ToFakeVertsIn must be defined on mesh2");

    STK_ThrowRequireMsg(m_mesh2EntityDestinations->get_mesh() == m_mesh2,
                        "mesh2EntityDestinations must be defined on mesh2");                               
  } 

  if (m_mesh2ScatteredToMesh1)
  {
    STK_ThrowRequireMsg(m_mesh2EntityOrigins->get_mesh() == m_mesh2ScatteredToMesh1,
                        "mesh2EntityDestinations must be defined on mesh2ScatteredToMesh1");          
  }
}



}
}
}
}