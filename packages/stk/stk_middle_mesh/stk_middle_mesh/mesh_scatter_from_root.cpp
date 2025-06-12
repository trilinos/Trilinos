#include "mesh_scatter_from_root.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

MeshScatterFromRoot::MeshScatterFromRoot(MPI_Comm unionComm, std::shared_ptr<Mesh> inputMesh, MPI_Comm meshComm,
                                         FieldPtr<int> elementDestinationRanks)
  : m_unionComm(unionComm)
  , m_inputMesh(inputMesh)
  , m_elementDestinationRanks(elementDestinationRanks)
{
  if (inputMesh && utils::impl::comm_size(inputMesh->get_comm()) != 1)
    throw std::runtime_error("input_mesh must be serial");

  // TODO: extract MeshGatherToRoot::checkMeshCommIsSubsetOfInputComm()
  //       to a free function
  // checkComm2IsSubsetOfComm1(input_comm, mesh_comm);

  m_amIRootRank   = inputMesh != nullptr;
  m_amIOnMeshComm = meshComm != MPI_COMM_NULL;

  verify_only_one_root_debug_only(unionComm, m_amIRootRank);

  if (m_amIOnMeshComm)
    m_outputMesh = make_empty_mesh(meshComm);

  if (m_amIRootRank)
  {
    m_localIds.resize(utils::impl::comm_size(m_unionComm));
    m_entityDestinations = create_variable_size_field<RemoteSharedEntity>(m_inputMesh, FieldShape(1, 1, 1));
  }
}

std::shared_ptr<Mesh> MeshScatterFromRoot::scatter()
{
  send_entities<VertexSendData>();
  send_entities<EdgeSendData>();
  send_entities<ElementSendData>();

  return m_outputMesh;
}

VariableSizeFieldPtr<RemoteSharedEntity> MeshScatterFromRoot::get_entity_destinations()
{
  return m_entityDestinations;
}


void MeshScatterFromRoot::verify_only_one_root_debug_only([[maybe_unused]] MPI_Comm unionComm, [[maybe_unused]] bool amIRoot)
{
#ifndef NDEBUG
  int numLocalRoots  = amIRoot ? 1 : 0;
  int numGlobalRoots = 0;
  MPI_Allreduce(&numLocalRoots, &numGlobalRoots, 1, MPI_INT, MPI_SUM, unionComm);

  if (numGlobalRoots > 1)
    throw std::runtime_error("more than one process has the input mesh.  This is not supported");
#endif
}

template <typename T>
void MeshScatterFromRoot::send_entities()
{
  Exchanger exchanger(m_unionComm);
  if (m_amIRootRank)
  {
    for (int i = 0; i < 2; ++i)
    {
      reset_for_new_pack(T::DIMENSION);
      for (auto& entity : m_inputMesh->get_mesh_entities(T::DIMENSION))
        if (entity)
        {
          T sendData;
          get_send_data(entity, sendData);
          pack_data(exchanger, sendData);
        }

      if (i == 0)
        exchanger.allocate_send_buffers();
    }
  }

  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();

  // TODO: this results in non-deterministic entity order
  auto f = [&](int rank, stk::CommBuffer& buf) { unpack_data<T>(rank, buf); };
  exchanger.complete_receives(f);
  exchanger.complete_sends();
}

void MeshScatterFromRoot::reset_for_new_pack(int dimension)
{
  std::fill(m_localIds.begin(), m_localIds.end(), 0);
  m_entityDestinations->clear(dimension);
}

template <typename T>
void MeshScatterFromRoot::pack_data(Exchanger& exchanger, const T& sendData)
{
  for (int i = 0; i < int(sendData.ranks.size()); ++i)
  {
    typename T::InfoType info;
    get_send_data_for_destination(sendData, i, info);

    int destRankOnInputComm = sendData.ranks[i];
    auto& buf               = exchanger.get_send_buf(destRankOnInputComm);
    pack_buffer(buf, info);
  }
}

template <typename T>
void MeshScatterFromRoot::unpack_data(int /*rank*/, stk::CommBuffer& buf)
{
  typename T::InfoType info;
  while (buf.remaining() > 0)
  {
    unpack_buffer(buf, info);
    unpack_info(info);
  }
}

void MeshScatterFromRoot::get_send_data(MeshEntityPtr vert, VertexSendData& sendData)
{
  assert(vert->get_type() == MeshEntityType::Vertex);

  sendData.pt = vert->get_point_orig(0);

  auto& elementDestinationRank = *m_elementDestinationRanks;
  std::vector<MeshEntityPtr> els;
  int nels = get_upward(vert, 2, els);
  std::set<int> seenDestRanks;
  for (int i = 0; i < nels; ++i)
  {
    int rankOnInputComm = elementDestinationRank(els[i], 0, 0);
    if (seenDestRanks.count(rankOnInputComm) == 0)
    {
      seenDestRanks.insert(rankOnInputComm);
      int localId = m_localIds[rankOnInputComm]++;

      sendData.ranks.push_back(rankOnInputComm);
      sendData.localIds.push_back(localId);
      m_entityDestinations->insert(vert, 0, {rankOnInputComm, localId});
    }
  }
}

void MeshScatterFromRoot::get_send_data_for_destination(const VertexSendData& sendData, int rankIdx, VertexInfo& info)
{
  assert(sendData.ranks.size() >= 1);
  assert(rankIdx >= 0 && rankIdx < int(sendData.ranks.size()));

  info.pt = sendData.pt;
  info.otherRanks.resize(sendData.ranks.size() - 1);
  info.otherLocalIds.resize(sendData.localIds.size() - 1);

  int destRank = sendData.ranks[rankIdx];
  int idx      = 0;
  for (int i = 0; i < int(sendData.ranks.size()); ++i)
    if (sendData.ranks[i] != destRank)
    {
      info.otherRanks[idx]    = sendData.ranks[i];
      info.otherLocalIds[idx] = sendData.localIds[i];
      idx++;
    }
}

void MeshScatterFromRoot::pack_buffer(stk::CommBuffer& buf, const VertexInfo& info)
{
  buf.pack(info.pt);
  buf.pack(info.otherLocalIds);
  buf.pack(info.otherRanks);
}

void MeshScatterFromRoot::unpack_buffer(stk::CommBuffer& buf, VertexInfo& info)
{
  buf.unpack(info.pt);
  buf.unpack(info.otherLocalIds);
  buf.unpack(info.otherRanks);
}

void MeshScatterFromRoot::unpack_info(const VertexInfo& info)
{
  auto vert = m_outputMesh->create_vertex(info.pt);
  for (int i = 0; i < int(info.otherRanks.size()); ++i)
  {
    int otherRank = translate_input_comm_rank_to_mesh_comm_rank(info.otherRanks[i]);
    vert->add_remote_shared_entity({otherRank, info.otherLocalIds[i]});
  }
}

int MeshScatterFromRoot::translate_input_comm_rank_to_mesh_comm_rank(int unionCommRank)
{
  MPI_Group unionCommGroup, meshCommGroup;
  MPI_Comm_group(m_unionComm, &unionCommGroup);
  MPI_Comm_group(m_outputMesh->get_comm(), &meshCommGroup);

  int meshCommRank;
  MPI_Group_translate_ranks(unionCommGroup, 1, &unionCommRank, meshCommGroup, &meshCommRank);

  return meshCommRank;
}

void MeshScatterFromRoot::get_send_data(MeshEntityPtr edge, EdgeSendData& sendData)
{
  assert(edge->get_type() == MeshEntityType::Edge);
  MeshEntityPtr vert1 = edge->get_down(0);
  MeshEntityPtr vert2 = edge->get_down(1);

  auto& elementDestinationRank = *m_elementDestinationRanks;
  std::set<int> seenDestRanks;
  for (int i = 0; i < edge->count_up(); ++i)
  {
    MeshEntityPtr el        = edge->get_up(i);
    int destRankOnInputComm = elementDestinationRank(el, 0, 0);
    if (seenDestRanks.count(destRankOnInputComm) == 0)
    {
      seenDestRanks.insert(destRankOnInputComm);
      int edgeLocalId = m_localIds[destRankOnInputComm]++;

      sendData.vert1LocalIds.push_back(get_local_id_on_rank(vert1, destRankOnInputComm));
      sendData.vert2LocalIds.push_back(get_local_id_on_rank(vert2, destRankOnInputComm));
      sendData.edgeLocalIds.push_back(edgeLocalId);
      sendData.ranks.push_back(destRankOnInputComm);

      m_entityDestinations->insert(edge, 0, {destRankOnInputComm, edgeLocalId});
    }
  }
}

int MeshScatterFromRoot::get_local_id_on_rank(MeshEntityPtr entity, int rankOnInputComm)
{
  auto& entityDestinations = *m_entityDestinations;
  for (int i = 0; i < entityDestinations.get_num_comp(entity, 0); ++i)
    if (entityDestinations(entity, 0, i).remoteRank == rankOnInputComm)
      return entityDestinations(entity, 0, i).remoteId;

  throw std::runtime_error("unable to find local id on given rank");
}

void MeshScatterFromRoot::get_send_data_for_destination(const EdgeSendData& sendData, int rankIdx, EdgeInfo& info)
{
  assert(rankIdx >= 0 && rankIdx <= 1);
  info.vert1LocalId = sendData.vert1LocalIds[rankIdx];
  info.vert2LocalId = sendData.vert2LocalIds[rankIdx];

  if (sendData.ranks.size() > 1)
  {
    info.remoteId   = sendData.edgeLocalIds[1 - rankIdx];
    info.remoteRank = sendData.ranks[1 - rankIdx];
  }
}

void MeshScatterFromRoot::pack_buffer(stk::CommBuffer& buf, const EdgeInfo& info)
{
  buf.pack(info);
}

void MeshScatterFromRoot::unpack_buffer(stk::CommBuffer& buf, EdgeInfo& info)
{
  buf.unpack(info);
}

void MeshScatterFromRoot::unpack_info(const EdgeInfo& info)
{
  auto& verts = m_outputMesh->get_vertices();
  auto edge   = m_outputMesh->create_edge(verts[info.vert1LocalId], verts[info.vert2LocalId]);

  if (info.remoteRank != -1)
  {
    int rankOnMeshComm = translate_input_comm_rank_to_mesh_comm_rank(info.remoteRank);
    edge->add_remote_shared_entity({rankOnMeshComm, info.remoteId});
  }
}

void MeshScatterFromRoot::get_send_data(MeshEntityPtr el, ElementSendData& sendData)
{
  assert(get_type_dimension(el->get_type()) == 2);

  auto& elementDestinationRanks = *m_elementDestinationRanks;
  int destRankOnInputComm       = elementDestinationRanks(el, 0, 0);
  int elLocalId                 = m_localIds[destRankOnInputComm]++;

  std::array<int, 4> edgeLocalIds;
  for (int i = 0; i < el->count_down(); ++i)
    edgeLocalIds[i] = get_local_id_on_rank(el->get_down(i), destRankOnInputComm);

  sendData.ranks.resize(1);
  sendData.ranks[0] = destRankOnInputComm;
  if (el->get_type() == MeshEntityType::Triangle)
  {
    sendData.edge1LocalId     = edgeLocalIds[0];
    sendData.edge2LocalId     = edgeLocalIds[1];
    sendData.edge3LocalId     = edgeLocalIds[2];
    sendData.edge4LocalId     = -1;
    sendData.edge1Orientation = el->get_down_orientation(0);
  } else if (el->get_type() == MeshEntityType::Quad)
  {
    sendData.edge1LocalId     = edgeLocalIds[0];
    sendData.edge2LocalId     = edgeLocalIds[1];
    sendData.edge3LocalId     = edgeLocalIds[2];
    sendData.edge4LocalId     = edgeLocalIds[3];
    sendData.edge1Orientation = el->get_down_orientation(0);
  }

  m_entityDestinations->insert(el, 0, {destRankOnInputComm, elLocalId});
}

void MeshScatterFromRoot::get_send_data_for_destination(const ElementSendData& sendData, [[maybe_unused]] int rankIdx, ElementInfo& info)
{
  assert(rankIdx == 0);
  info.edge1LocalId     = sendData.edge1LocalId;
  info.edge2LocalId     = sendData.edge2LocalId;
  info.edge3LocalId     = sendData.edge3LocalId;
  info.edge4LocalId     = sendData.edge4LocalId;
  info.edge1Orientation = sendData.edge1Orientation;
}

void MeshScatterFromRoot::pack_buffer(stk::CommBuffer& buf, const ElementInfo& info)
{
  buf.pack(info);
}

void MeshScatterFromRoot::unpack_buffer(stk::CommBuffer& buf, ElementInfo& info)
{
  buf.unpack(info);
}

void MeshScatterFromRoot::unpack_info(const ElementInfo& info)
{
  auto& edges = m_outputMesh->get_edges();
  if (info.edge4LocalId == -1)
  {
    m_outputMesh->create_triangle(edges[info.edge1LocalId], edges[info.edge2LocalId], edges[info.edge3LocalId],
                                  info.edge1Orientation);
  } else
  {
    m_outputMesh->create_quad(edges[info.edge1LocalId], edges[info.edge2LocalId], edges[info.edge3LocalId],
                              edges[info.edge4LocalId], info.edge1Orientation);
  }
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
