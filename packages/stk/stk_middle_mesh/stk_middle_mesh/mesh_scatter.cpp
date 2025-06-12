#include "mesh_scatter.hpp"
#include "mesh_entity.hpp"
#include "variable_size_field.hpp"
#include "destination_field_synchronizer.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "create_sharing_from_verts.hpp"



namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

MeshScatter::MeshScatter(std::shared_ptr<impl::MeshScatterSpec> scatterSpec,
            std::shared_ptr<Mesh> inputMesh, MPI_Comm scatteredMeshComm,
            bool computeEntityCorrespondence) :
  m_unionComm(scatterSpec->get_comm()),
  m_scatterSpec(scatterSpec),
  m_inputMesh(inputMesh),
  m_computeEntityCorrespondence(computeEntityCorrespondence)
{
  if (scatteredMeshComm != MPI_COMM_NULL)
  {
    m_outputMesh = make_empty_mesh(scatteredMeshComm);
    m_elementOrigins = create_field<RemoteSharedEntity>(m_outputMesh, FieldShape(0, 0, 1), 1);
  }
  check_mesh_comm_is_subset_of_union_comm();

  m_vertsBySrcMeshOwner = std::make_shared<EntitySortedByOwner>(utils::impl::comm_size(m_unionComm));
  m_edgesBySrcMeshOwner = std::make_shared<EntitySortedByOwner>(utils::impl::comm_size(m_unionComm));

  if (m_computeEntityCorrespondence)
  {
    if (inputMesh)
      m_entityDestinations = create_variable_size_field<RemoteSharedEntity>(inputMesh, FieldShape(1, 1, 1));

    if (m_outputMesh)
      m_entityOrigins = create_variable_size_field<RemoteSharedEntity>(m_outputMesh, FieldShape(1, 1, 1));
  }
}

std::shared_ptr<mesh::Mesh> MeshScatter::scatter()
{
  VariableSizeFieldPtr<int> destRanksOnUnionComm;
  if (m_inputMesh)
  {
    DestinationFieldSynchronizer gatherer(m_inputMesh, m_scatterSpec);
    destRanksOnUnionComm = gatherer.synchronize();
  }

  //TODO: possible performance improvement: change the order:
  //         1. pack verts, 
  //         2. start sending edges
  //         3. pack edges
  //         4. start sending edges
  //         5. unpack verts
  //         6. unpack edges
  send_verts(destRanksOnUnionComm);
  send_edges(destRanksOnUnionComm);
  send_elements();  


  if (m_outputMesh)
  {
    for (int dim=1; dim <= 2; ++dim)
    {
      CreateSharingFromVert edgeSharingCreator(m_outputMesh, dim);
      edgeSharingCreator.create_sharing_from_verts();
    }
  }
        

  return m_outputMesh;
}

mesh::FieldPtr<mesh::RemoteSharedEntity> MeshScatter::get_element_origins() const
{
  return m_elementOrigins;
}

VariableSizeFieldPtr<RemoteSharedEntity> MeshScatter::get_entity_origins() const
{
  assert(m_computeEntityCorrespondence);
  return m_entityOrigins;
}

VariableSizeFieldPtr<RemoteSharedEntity> MeshScatter::get_entity_destinations() const
{
  assert(m_computeEntityCorrespondence);
  return m_entityDestinations;
}

void MeshScatter::check_mesh_comm_is_subset_of_union_comm()
{
  if (m_outputMesh)
  {
    int meshCommSize = utils::impl::comm_size(m_outputMesh->get_comm());
    std::vector<int> meshCommRanks(meshCommSize), unionCommRanks(meshCommSize);
    for (int rank=0; rank < meshCommSize; ++rank)
      meshCommRanks[rank] = rank;

    MPI_Group meshGroup, unionGroup;
    MPI_Comm_group(m_outputMesh->get_comm(), &meshGroup);
    MPI_Comm_group(m_unionComm, &unionGroup);

    MPI_Group_translate_ranks(meshGroup, meshCommRanks.size(), meshCommRanks.data(),
                              unionGroup, unionCommRanks.data());

    for (int i=0; i < meshCommSize; ++i)
      if (unionCommRanks[i] == MPI_UNDEFINED)
        throw std::runtime_error("outputMeshComm is not a subset of unionComm");
  }
}


void MeshScatter::send_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr)
{
  EntityExchanger entityExchanger(m_unionComm);
  EntityIdExchanger entityCorrespondenceExchanger(m_unionComm);
  std::shared_ptr<EntityIdExchanger> sharingExchanger;
  if (m_outputMesh)
    sharingExchanger = std::make_shared<EntityIdExchanger>(m_outputMesh->get_comm());

  if (m_inputMesh)
  {
    pack_verts(destRanksOnUnionCommPtr, entityExchanger, entityCorrespondenceExchanger);
  }

  unpack_verts_and_pack_sharing(entityExchanger, sharingExchanger, entityCorrespondenceExchanger);

  if (m_outputMesh)
  {
    unpack_vert_sharing(sharingExchanger);
  }

  if (m_computeEntityCorrespondence)
  {
    unpack_returned_entity_ids(destRanksOnUnionCommPtr, 0, entityCorrespondenceExchanger);
  }
}

void MeshScatter::pack_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger, EntityIdExchanger& entityCorrespondenceExchanger)
{
  assert(m_inputMesh);

  std::vector<int> vertDestRanksOnUnionComm;
  auto& destRanksOnUnionComm = *destRanksOnUnionCommPtr;
  for (int phase=0; phase < 2; ++phase)
  {
    for (auto& vert : m_inputMesh->get_vertices())
      if (vert)
      {

        RemoteSharedEntity owner = get_owner_remote(m_inputMesh, vert);
        vertDestRanksOnUnionComm.assign(destRanksOnUnionComm(vert, 0).begin(), 
                                        destRanksOnUnionComm(vert, 0).end());
        for (const auto& destRankOnUnionComm : vertDestRanksOnUnionComm)
        {
          // need to pack my id (so the receiver can setup their fields)
          auto& buf = entityExchanger.get_send_buf(destRankOnUnionComm);
          buf.pack(vert->get_point_orig(0));
          buf.pack(vert->get_id());
          buf.pack(owner);
          buf.pack(vertDestRanksOnUnionComm);

          if (phase == 1 && m_computeEntityCorrespondence)
            entityCorrespondenceExchanger.get_recv_buf(destRankOnUnionComm).push_back(-1);
        }
      }

    if (phase == 0)
      entityExchanger.allocate_send_buffers();
  }      
}

void MeshScatter::unpack_verts_and_pack_sharing(EntityExchanger& entityExchanger, 
                                std::shared_ptr<EntityIdExchanger> sharingExchanger,
                                EntityIdExchanger& entityCorrespondenceExchanger)
{
  entityExchanger.start_nonblocking();
  entityExchanger.post_nonblocking_receives();

  // don't unpack data as it comes in because that would make the vertex ordering
  // nondeterministic      
  auto unpackVerts = [&](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  entityExchanger.complete_receives(unpackVerts);

  for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
  {
    auto& buf = entityExchanger.get_recv_buf(rank);
    if (buf.remaining() > 0)
      unpack_vert_buffer(rank, buf, sharingExchanger, entityCorrespondenceExchanger);
  }
}

void MeshScatter::unpack_vert_buffer(int rank, stk::CommBuffer& buf, 
                     std::shared_ptr<EntityIdExchanger> sharingExchanger,
                     EntityIdExchanger& entityCorrespondenceExchanger)
{
  int myRankOnOutputComm = utils::impl::comm_rank(m_outputMesh->get_comm());
  std::vector<int> destRanksOnUnionComm, destRanksOnOutputComm;
  while (buf.remaining() > 0)
  {
    utils::Point pt;
    int idOnSender;
    RemoteSharedEntity ownerOnSrcMesh;
    destRanksOnUnionComm.clear();
    buf.unpack(pt);
    buf.unpack(idOnSender);
    buf.unpack(ownerOnSrcMesh);
    buf.unpack(destRanksOnUnionComm);

    MeshEntityPtr v = m_vertsBySrcMeshOwner->get_value(ownerOnSrcMesh);
    if (!v)
    {
      v = m_outputMesh->create_vertex(pt);
      m_vertsBySrcMeshOwner->insert(ownerOnSrcMesh, v);

      translate_union_comm_ranks_to_output_comm(destRanksOnUnionComm, destRanksOnOutputComm);

      for (auto& destRank : destRanksOnOutputComm)
        if (destRank != myRankOnOutputComm)
        {
          sharingExchanger->get_send_buf(destRank).push_back(ownerOnSrcMesh.remoteRank);
          sharingExchanger->get_send_buf(destRank).push_back(ownerOnSrcMesh.remoteId);
          sharingExchanger->get_send_buf(destRank).push_back(v->get_id());
          for (int i=0; i < 3; ++i)
          {
            sharingExchanger->get_recv_buf(destRank).push_back(-1);
          }
        }
    }

    if (m_computeEntityCorrespondence)
    {
      m_entityOrigins->insert(v, 0, RemoteSharedEntity(rank, idOnSender));
      entityCorrespondenceExchanger.get_send_buf(rank).push_back(v->get_id());
    }
  }
}

void MeshScatter::unpack_vert_sharing(std::shared_ptr<EntityIdExchanger> sharingExchanger)
{
  auto unpackShared = [&](int rank, const std::vector<int>& buf)
  {
    unpack_vert_sharing_buffer(rank, buf);
  };
  sharingExchanger->start_nonblocking();
  sharingExchanger->complete_receives(unpackShared);
}

void MeshScatter::unpack_returned_entity_ids(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, int dim,
                                             EntityIdExchanger& entityCorrespondenceExchanger)
{
  auto& destRanksOnUnionComm = *destRanksOnUnionCommPtr;

  entityCorrespondenceExchanger.start_nonblocking();

  auto unpacker = [&](int rank, const std::vector<int>& ids)
  {
    assert(m_inputMesh);
    int idIdx = 0;
    std::vector<int> vertDestRanksOnUnionComm;
    for (auto& entity : m_inputMesh->get_mesh_entities(dim))
      if (entity)
      {
        vertDestRanksOnUnionComm.assign(destRanksOnUnionComm(entity, 0).begin(), 
                                        destRanksOnUnionComm(entity, 0).end());

        bool foundDest = std::find(vertDestRanksOnUnionComm.begin(), vertDestRanksOnUnionComm.end(), rank) != vertDestRanksOnUnionComm.end();
        if (foundDest)
        {
          m_entityDestinations->insert(entity, 0, RemoteSharedEntity(rank, ids[idIdx++]));
        }
      }
  };

  entityCorrespondenceExchanger.complete_receives(unpacker);
}

void MeshScatter::translate_union_comm_ranks_to_output_comm(const std::vector<int>& unionCommRanks, std::vector<int>& outputCommRanks)
{
  assert(m_outputMesh);

  MPI_Group unionGroup, outputGroup;
  MPI_Comm_group(m_unionComm, &unionGroup);
  MPI_Comm_group(m_outputMesh->get_comm(), &outputGroup);

  outputCommRanks.resize(unionCommRanks.size());
  MPI_Group_translate_ranks(unionGroup, unionCommRanks.size(), unionCommRanks.data(), outputGroup, outputCommRanks.data());
}

void MeshScatter::unpack_vert_sharing_buffer(int rank, const std::vector<int>& buf)
{
  assert(buf.size() % 3 == 0);
  for (size_t i=0; i < buf.size(); i += 3)
  {
    RemoteSharedEntity ownerOnSrcMesh{buf[i], buf[i+1]};
    int remoteId = buf[i+2];

    MeshEntityPtr v = m_vertsBySrcMeshOwner->get_value(ownerOnSrcMesh);
    v->add_remote_shared_entity(RemoteSharedEntity{rank, remoteId});
  }
}

void MeshScatter::send_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr)
{
  EntityExchanger entityExchanger(m_unionComm);
  EntityIdExchanger entityCorrespondenceExchanger(m_unionComm);
  pack_edges(destRanksOnUnionCommPtr, entityExchanger, entityCorrespondenceExchanger);
  unpack_edges(entityExchanger, entityCorrespondenceExchanger);

  if (m_computeEntityCorrespondence)
  {
    unpack_returned_entity_ids(destRanksOnUnionCommPtr, 1, entityCorrespondenceExchanger);
  }
}


void MeshScatter::pack_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger,
                             EntityIdExchanger& entityCorrespondenceExchanger)
{
  auto& destRanksOnUnionComm = *destRanksOnUnionCommPtr;
  if (m_inputMesh)
  {
    for (int phase=0; phase < 2; ++phase)
    {
      for (auto& edge : m_inputMesh->get_edges())
        if (edge)
        {
          RemoteSharedEntity edgeOwner = get_owner_remote(m_inputMesh, edge);
          RemoteSharedEntity v1Owner   = get_owner_remote(m_inputMesh, edge->get_down(0));
          RemoteSharedEntity v2Owner   = get_owner_remote(m_inputMesh, edge->get_down(1));

          for (int destRank : destRanksOnUnionComm(edge, 0))
          {
            auto& buf = entityExchanger.get_send_buf(destRank);
            buf.pack(edge->get_id());
            buf.pack(edgeOwner);
            buf.pack(v1Owner);
            buf.pack(v2Owner);

            if (phase == 1 && m_computeEntityCorrespondence)
            {
              entityCorrespondenceExchanger.get_recv_buf(destRank).push_back(-1);
            }              
          }
        }

      if (phase == 0)
        entityExchanger.allocate_send_buffers();
    }
  }      
}

void MeshScatter::unpack_edges(EntityExchanger& entityExchanger, EntityIdExchanger& entityCorrespondenceExchanger)
{
  // dont process the buffers as they arrive because that would make
  // the edge ids non-deterministic
  auto f = [&](int /*rank*/, stk::CommBuffer& /*buf*/) {};
  entityExchanger.start_nonblocking();
  entityExchanger.post_nonblocking_receives();
  entityExchanger.complete_receives(f);

    for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
    {
      unpack_edge_buffer(rank, entityExchanger.get_recv_buf(rank), entityCorrespondenceExchanger);
    }
}

void MeshScatter::unpack_edge_buffer(int rank, stk::CommBuffer& buf, EntityIdExchanger& entityCorrespondenceExchanger)
{
  while (buf.remaining() > 0)
  {
    int idOnSender;
    RemoteSharedEntity edgeOwnerOnSrcMesh, v1OwnerOnSrcMesh, v2OwnerOnSrcMesh;
    buf.unpack(idOnSender);
    buf.unpack(edgeOwnerOnSrcMesh);
    buf.unpack(v1OwnerOnSrcMesh);
    buf.unpack(v2OwnerOnSrcMesh);

    MeshEntityPtr v1 = m_vertsBySrcMeshOwner->get_value(v1OwnerOnSrcMesh);
    MeshEntityPtr v2 = m_vertsBySrcMeshOwner->get_value(v2OwnerOnSrcMesh);
    MeshEntityPtr edge = m_edgesBySrcMeshOwner->get_value(edgeOwnerOnSrcMesh);
    if (!edge)
    {
      edge = m_outputMesh->create_edge(v1, v2);
      m_edgesBySrcMeshOwner->insert(edgeOwnerOnSrcMesh, edge);    
    }

    if (m_computeEntityCorrespondence)
    {
      m_entityOrigins->insert(edge, 0, RemoteSharedEntity(rank, idOnSender));
      entityCorrespondenceExchanger.get_send_buf(rank).push_back(edge->get_id());
    }
  }
}


void MeshScatter::send_elements()
{
  EntityExchanger exchanger(m_unionComm);
  EntityIdExchanger entityCorrespondenceExchanger(m_unionComm);
  pack_elements(exchanger, entityCorrespondenceExchanger);
  unpack_elements(exchanger, entityCorrespondenceExchanger);

  if (m_computeEntityCorrespondence)
  {
    unpack_returned_element_ids(entityCorrespondenceExchanger);
  }
}

void MeshScatter::pack_elements(EntityExchanger& entityExchanger, EntityIdExchanger& entityCorrespondenceExchanger)
{
  std::vector<int> destRanksOnUnionComm;
  if (m_inputMesh)
  {
    for (int phase=0; phase < 2; ++phase)
    {
      for (auto& el : m_inputMesh->get_elements())
        if (el)
        {
          std::array<RemoteSharedEntity, MAX_DOWN> owners;
          std::fill(owners.begin(), owners.end(), RemoteSharedEntity{-1, -1});
          for (int i=0; i < el->count_down(); ++i)
            owners[i] = get_owner_remote(m_inputMesh, el->get_down(i));

          EntityOrientation edge1Orient = el->get_down_orientation(0);

          m_scatterSpec->get_destinations(el, destRanksOnUnionComm);
          for (auto& destRank : destRanksOnUnionComm)
          {
            auto& buf = entityExchanger.get_send_buf(destRank);
            buf.pack(el->get_id());
            for (int i=0; i < MAX_DOWN; ++i)
            {
              buf.pack(owners[i]);
            }
            buf.pack(static_cast<int>(edge1Orient));
            buf.pack(el->get_id());

            if (phase == 1 && m_computeEntityCorrespondence)
            {
              entityCorrespondenceExchanger.get_recv_buf(destRank).push_back(-1);
            }     
          }
        }

      if (phase == 0)
        entityExchanger.allocate_send_buffers();
    }
  }
}

void MeshScatter::unpack_elements(EntityExchanger& exchanger, EntityIdExchanger& entityCorrespondenceExchanger)
{
  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();
  auto f = [&](int /*rank*/, stk::CommBuffer& /*buf*/) {};

  exchanger.complete_receives(f);
  for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
  {
    unpack_element_buffer(rank, exchanger.get_recv_buf(rank), entityCorrespondenceExchanger);
  }      
}

void MeshScatter::unpack_element_buffer(int rank, stk::CommBuffer& buf, EntityIdExchanger& entityCorrespondenceExchanger)
{
  auto& elementOrigins = *m_elementOrigins;
  while(buf.remaining() > 0)
  {
    std::array<RemoteSharedEntity, MAX_DOWN> edgeOwnersOnSrcMesh;
    int idOnSender, edge1Orient, srcMeshElementId;
    buf.unpack(idOnSender);
    for (int i=0; i < MAX_DOWN; ++i)
    {
      buf.unpack(edgeOwnersOnSrcMesh[i]);
    }
    buf.unpack(edge1Orient);
    buf.unpack(srcMeshElementId);

    MeshEntityPtr edge1 = m_edgesBySrcMeshOwner->get_value(edgeOwnersOnSrcMesh[0]);
    MeshEntityPtr edge2 = m_edgesBySrcMeshOwner->get_value(edgeOwnersOnSrcMesh[1]);
    MeshEntityPtr edge3 = m_edgesBySrcMeshOwner->get_value(edgeOwnersOnSrcMesh[2]);
    MeshEntityPtr el;
    if (edgeOwnersOnSrcMesh[3].remoteRank == -1)
    {
      el = m_outputMesh->create_triangle(edge1, edge2, edge3, static_cast<EntityOrientation>(edge1Orient));
    } else
    {
      MeshEntityPtr edge4 = m_edgesBySrcMeshOwner->get_value(edgeOwnersOnSrcMesh[3]);
      el = m_outputMesh->create_quad(edge1, edge2, edge3, edge4, static_cast<EntityOrientation>(edge1Orient));
    }

    elementOrigins(el, 0, 0) = RemoteSharedEntity{rank, srcMeshElementId};

    if (m_computeEntityCorrespondence)
    {
      m_entityOrigins->insert(el, 0, RemoteSharedEntity(rank, idOnSender));
      entityCorrespondenceExchanger.get_send_buf(rank).push_back(el->get_id());
    }
  }
}

void MeshScatter::unpack_returned_element_ids(EntityIdExchanger& entityCorrespondenceExchanger)
{
  entityCorrespondenceExchanger.start_nonblocking();

  auto unpacker = [&](int rank, const std::vector<int>& ids)
  {
    std::vector<int> destRanksOnUnionComm;
    int idIdx = 0;
    for (auto el : m_inputMesh->get_elements())
      if (el)
      {
        destRanksOnUnionComm.clear();
        m_scatterSpec->get_destinations(el, destRanksOnUnionComm);

        auto it = std::find(destRanksOnUnionComm.begin(), destRanksOnUnionComm.end(), rank);
        if (it != destRanksOnUnionComm.end())
        {
          int remoteId = ids[idIdx++];
          m_entityDestinations->insert(el, 0, RemoteSharedEntity(rank, remoteId));
        }
      }
  };

  entityCorrespondenceExchanger.complete_receives(unpacker);
}


}
}
}
}
