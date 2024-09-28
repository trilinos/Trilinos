#include "mesh_gather_to_root.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

MeshGatherToRoot::MeshGatherToRoot(MPI_Comm inputComm, int rootRankOnInputComm, std::shared_ptr<Mesh> mesh)
  : m_inputComm(inputComm)
  , m_inputMesh(mesh)
{
  check_mesh_comm_is_subset_of_input_comm(inputComm);

  bool amIRootRank   = utils::impl::comm_rank(inputComm) == rootRankOnInputComm;
  bool amIOnMeshComm = mesh != nullptr;

  int color, key;
  if (amIRootRank)
  {
    color = 1;
    key   = 0;
  } else if (amIOnMeshComm)
  {
    color = 1;
    key   = utils::impl::comm_rank(mesh->get_comm()) + 1;
  } else
  {
    color = MPI_UNDEFINED;
    key   = 0;
  }

  MPI_Comm_split(inputComm, color, key, &m_comm);

  if (amIRootRank)
  {
    m_outputMesh     = make_empty_mesh(MPI_COMM_SELF);
    m_elementOrigins = create_field<RemoteSharedEntity>(m_outputMesh, FieldShape(0, 0, 1), 1);
  }
}

// performs the mesh gather.  On the root process (see constructor)
// returns the newly created mesh object.  Returns nullptr otherwise
std::shared_ptr<Mesh> MeshGatherToRoot::gather()
{
  if (m_comm != MPI_COMM_NULL)
  {
    send_sizes_to_root();
    send_entities<VertexInfo>(28);
    send_entities<EdgeInfo>(29);
    send_entities<ElementInfo>(30);
  }

  return m_outputMesh;
}

FieldPtr<RemoteSharedEntity> MeshGatherToRoot::get_element_origins()
{
  return m_elementOrigins;
}

void MeshGatherToRoot::check_mesh_comm_is_subset_of_input_comm(MPI_Comm inputComm)
{
  if (m_inputMesh)
  {
    MPI_Comm meshComm = m_inputMesh->get_comm();
    int meshCommSize  = utils::impl::comm_size(meshComm);
    std::vector<int> meshCommRanks(meshCommSize), meshCommRanksOnInputComm(meshCommSize);
    for (int i = 0; i < meshCommSize; ++i)
      meshCommRanks[i] = i;

    MPI_Group meshCommGroup, inputCommGroup;
    MPI_Comm_group(meshComm, &meshCommGroup);
    MPI_Comm_group(inputComm, &inputCommGroup);

    MPI_Group_translate_ranks(meshCommGroup, meshCommSize, meshCommRanks.data(), inputCommGroup,
                              meshCommRanksOnInputComm.data());

    for (int i = 0; i < meshCommSize; ++i)
      if (meshCommRanksOnInputComm[i] == MPI_UNDEFINED)
        throw std::runtime_error("mesh comm is not a subset of input_comm");
  }
}

void MeshGatherToRoot::send_sizes_to_root()
{
  const int root = 0, tag = 27;
  int nprocs = utils::impl::comm_size(m_comm);

  utils::impl::ParallelExchange<MeshEntityCounts> exchanger(m_comm, tag);
  if (m_inputMesh)
  {
    MeshEntityCounts counts;
    counts.rankOnMeshComm = utils::impl::comm_rank(m_inputMesh->get_comm());
    for (int dim = 0; dim < 3; ++dim)
    {
      counts.entityCounts[dim] = count_valid(m_inputMesh->get_mesh_entities(dim));
      counts.entityMaxIds[dim] = m_inputMesh->get_mesh_entities(dim).size();
    }

    exchanger.get_send_buffer(root).push_back(counts);

    exchanger.start_sends();
  }

  if (m_outputMesh)
  {
    for (int i = 1; i < nprocs; ++i)
      exchanger.set_recv_buffer_size(i, 1);

    if (m_inputMesh)
      exchanger.set_recv_buffer_size(root, 1);

    exchanger.start_recvs();
    exchanger.complete_recvs();
    m_entityCounts.resize(nprocs);
    m_splitCommToMeshCommRank.resize(nprocs, -1);
    for (int rank = 0; rank < nprocs; ++rank)
    {
      const auto& recvBuf = exchanger.get_recv_buffer(rank);
      if (recvBuf.size() > 0)
      {
        assert(recvBuf.size() == 1);
        m_entityCounts[rank]            = recvBuf[0];
        m_splitCommToMeshCommRank[rank] = recvBuf[0].rankOnMeshComm;
      }
    }
  }

  exchanger.complete_sends();
}

template <typename T>
void MeshGatherToRoot::send_entities(int tag)
{
  const int root = 0;
  int nprocs     = utils::impl::comm_size(m_comm);
  utils::impl::ParallelExchange<T> exchanger(m_comm, tag);
  if (m_inputMesh)
  {
    int myrank = utils::impl::comm_rank(m_inputMesh->get_comm());

    auto& buf = exchanger.get_send_buffer(root);
    for (auto& entity : m_inputMesh->get_mesh_entities(T::DIMENSION))
      if (entity)
      {
        T info;
        get_info(entity, myrank, info);
        buf.push_back(info);
      }

    exchanger.start_sends();
  }

  if (m_outputMesh)
  {
    for (int i = 0; i < nprocs; ++i)
      exchanger.set_recv_buffer_size(i, m_entityCounts[i].entityCounts[T::DIMENSION]);

    exchanger.start_recvs();
    exchanger.complete_recvs();
    unpack_info(exchanger);
  }

  exchanger.complete_sends();
}

void MeshGatherToRoot::get_info(MeshEntityPtr vert, int myrank, VertexInfo& vertInfo)
{
  vertInfo.ownerRank    = get_owner(m_inputMesh, vert);
  vertInfo.ownerLocalId = vert->get_id();
  vertInfo.coords       = vert->get_point_orig(0);

  if (vertInfo.ownerRank != myrank)
    vertInfo.ownerLocalId = get_remote_shared_entity(vert, vertInfo.ownerRank).remoteId;
}

void MeshGatherToRoot::unpack_info(utils::impl::ParallelExchange<VertexInfo>& exchanger)
{
  int nprocs = utils::impl::comm_size(m_comm);
  m_ownerLocalIdToOutputVertex.resize(nprocs);
  for (int splitCommRank = 0; splitCommRank < nprocs; ++splitCommRank)
  {
    int rankOnMeshComm = m_splitCommToMeshCommRank[splitCommRank];

    if (rankOnMeshComm == -1)
      continue;

    const auto& recvBuf              = exchanger.get_recv_buffer(splitCommRank);
    auto& ownerLocalIdToOutputVertex = m_ownerLocalIdToOutputVertex[rankOnMeshComm];
    ownerLocalIdToOutputVertex.resize(m_entityCounts[splitCommRank].entityMaxIds[0], nullptr);

    for (auto& vertInfo : recvBuf)
      if (vertInfo.ownerRank == rankOnMeshComm)
      {
        MeshEntityPtr vert                                = m_outputMesh->create_vertex(vertInfo.coords);
        ownerLocalIdToOutputVertex[vertInfo.ownerLocalId] = vert;
      }
  }
}

void MeshGatherToRoot::get_info(MeshEntityPtr edge, int myrank, EdgeInfo& edgeInfo)
{
  edgeInfo.edgeOwnerRank    = get_owner(m_inputMesh, edge);
  edgeInfo.edgeOwnerLocalId = edge->get_id();

  if (edgeInfo.edgeOwnerRank != myrank)
    edgeInfo.edgeOwnerLocalId = get_remote_shared_entity(edge, edgeInfo.edgeOwnerRank).remoteId;

  VertexInfo vert1Info, vert2Info;
  get_info(edge->get_down(0), myrank, vert1Info);
  get_info(edge->get_down(1), myrank, vert2Info);

  edgeInfo.vert1OwnerRank    = vert1Info.ownerRank;
  edgeInfo.vert1OwnerLocalId = vert1Info.ownerLocalId;

  edgeInfo.vert2OwnerRank    = vert2Info.ownerRank;
  edgeInfo.vert2OwnerLocalId = vert2Info.ownerLocalId;
}

void MeshGatherToRoot::unpack_info(utils::impl::ParallelExchange<EdgeInfo>& exchanger)
{
  int nprocs = utils::impl::comm_size(m_comm);
  m_ownerLocalIdToOutputEdge.resize(nprocs);
  for (int splitCommRank = 0; splitCommRank < nprocs; ++splitCommRank)
  {
    int rankOnMeshComm = m_splitCommToMeshCommRank[splitCommRank];
    if (rankOnMeshComm == -1)
      continue;

    const auto& recvBuf            = exchanger.get_recv_buffer(splitCommRank);
    auto& ownerLocalIdToOutputEdge = m_ownerLocalIdToOutputEdge[rankOnMeshComm];
    ownerLocalIdToOutputEdge.resize(m_entityCounts[splitCommRank].entityMaxIds[1], nullptr);
    for (auto& edgeInfo : recvBuf)
    {
      if (edgeInfo.edgeOwnerRank == rankOnMeshComm)
      {
        MeshEntityPtr vert1 = m_ownerLocalIdToOutputVertex[edgeInfo.vert1OwnerRank][edgeInfo.vert1OwnerLocalId];
        MeshEntityPtr vert2 = m_ownerLocalIdToOutputVertex[edgeInfo.vert2OwnerRank][edgeInfo.vert2OwnerLocalId];
        assert(vert1);
        assert(vert2);

        MeshEntityPtr edge                                  = m_outputMesh->create_edge(vert1, vert2);
        ownerLocalIdToOutputEdge[edgeInfo.edgeOwnerLocalId] = edge;
      }
    }
  }
}

void MeshGatherToRoot::get_info(MeshEntityPtr el, int myrank, ElementInfo& elInfo)
{
  for (int i = 0; i < el->count_down(); ++i)
  {
    EdgeInfo edgeInfo;
    get_info(el->get_down(i), myrank, edgeInfo);
    elInfo.edgeOwnerRanks[i] = edgeInfo.edgeOwnerRank;
    elInfo.edgeOwnerIds[i]   = edgeInfo.edgeOwnerLocalId;
  }

  elInfo.edge1Orient = el->get_down_orientation(0);
  elInfo.elementLocalId = el->get_id();
}

void MeshGatherToRoot::unpack_info(utils::impl::ParallelExchange<ElementInfo>& exchanger)
{
  auto& elementOrigins = *m_elementOrigins;
  int nprocs           = utils::impl::comm_size(m_comm);
  for (int splitCommRank = 0; splitCommRank < nprocs; ++splitCommRank)
  {
    int inputCommRank   = utils::impl::get_rank_on_other_comm(m_comm, m_inputComm, splitCommRank);
    const auto& recvBuf = exchanger.get_recv_buffer(splitCommRank);

    for (auto& elInfo : recvBuf)
    {
      std::array<MeshEntityPtr, 4> edges = {0};
      for (int i = 0; i < 4; ++i)
        if (elInfo.edgeOwnerRanks[i] != -1)
        {
          int ownerRank    = elInfo.edgeOwnerRanks[i];
          int ownerLocalId = elInfo.edgeOwnerIds[i];
          edges[i]         = m_ownerLocalIdToOutputEdge[ownerRank][ownerLocalId];
          assert(edges[i]);
        }

      MeshEntityPtr el;
      if (elInfo.edgeOwnerRanks[3] != -1)
      {
        el = m_outputMesh->create_quad(edges[0], edges[1], edges[2], edges[3], elInfo.edge1Orient);
      } else
      {
        el = m_outputMesh->create_triangle(edges[0], edges[1], edges[2], elInfo.edge1Orient);
      }

      elementOrigins(el, 0, 0) = RemoteSharedEntity{inputCommRank, elInfo.elementLocalId};
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
