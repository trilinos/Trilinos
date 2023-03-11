#include "gathered_mesh_coordinate_update.hpp"
#include "mesh_entity.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

GatheredMeshCoordinateUpdate::GatheredMeshCoordinateUpdate(MPI_Comm inputComm, std::shared_ptr<Mesh> serialMesh,
                                                           FieldPtr<RemoteSharedEntity> serialMeshElementOrigins,
                                                           std::shared_ptr<Mesh> parallelMesh)
  : m_inputComm(inputComm)
  , m_serialMesh(serialMesh)
  , m_serialMeshElementOrigins(serialMeshElementOrigins)
  , m_parallelMesh(parallelMesh)
  , m_amIRoot(!(serialMesh == nullptr))
{
  check_serial_mesh_comm_size();
  check_only_one_root();
  broadcast_root_process();
  check_element_count();
}

void GatheredMeshCoordinateUpdate::update()
{
  send_coordinates();
}

void GatheredMeshCoordinateUpdate::check_only_one_root()
{
  int rootVal       = m_amIRoot ? 1 : 0;
  int numberOfRoots = 0;
  MPI_Allreduce(&rootVal, &numberOfRoots, 1, MPI_INT, MPI_SUM, m_inputComm);

  if (numberOfRoots != 1)
    throw std::runtime_error("there can be exactly one root process");
}

void GatheredMeshCoordinateUpdate::check_serial_mesh_comm_size()
{
  if (m_serialMesh && utils::impl::comm_size(m_serialMesh->get_comm()) != 1)
    throw std::runtime_error("serial mesh communicator must have exactly one process");
}

void GatheredMeshCoordinateUpdate::broadcast_root_process()
{
  int rootRankInput  = m_serialMesh ? utils::impl::comm_rank(m_inputComm) : 0;
  int rootRankOutput = 0;
  MPI_Allreduce(&rootRankInput, &rootRankOutput, 1, MPI_INT, MPI_SUM, m_inputComm);

  m_rootRank = rootRankOutput;
}

void GatheredMeshCoordinateUpdate::check_element_count()
{
  int numLocalParallelMeshElements  = m_parallelMesh ? count_valid(m_parallelMesh->get_elements()) : 0;
  int numGlobalParallelMeshElements = 0;
  MPI_Reduce(&numLocalParallelMeshElements, &numGlobalParallelMeshElements, 1, MPI_INT, MPI_SUM, m_rootRank,
             m_inputComm);

  if (m_serialMesh && numGlobalParallelMeshElements != count_valid(m_serialMesh->get_elements()))
    throw std::runtime_error("number of elements on serial and parallel meshes must be equal");
}

void GatheredMeshCoordinateUpdate::send_coordinates()
{
  utils::impl::ParallelExchange<ElementCoordinateInfo> exchanger(m_inputComm, 75);

  if (m_amIRoot)
  {
    pack_send_buffers(exchanger);
  }

  if (m_parallelMesh)
    exchanger.set_recv_buffer_size(m_rootRank, count_valid(m_parallelMesh->get_elements()));

  exchanger.start_sends();
  exchanger.start_recvs();
  exchanger.complete_recvs();

  if (m_parallelMesh)
  {
    unpack_recv_buffers(exchanger);
  }

  exchanger.complete_sends();
}

void GatheredMeshCoordinateUpdate::pack_send_buffers(utils::impl::ParallelExchange<ElementCoordinateInfo>& exchanger)
{
  assert(m_amIRoot);

  auto& elOrigins = *m_serialMeshElementOrigins;
  std::array<MeshEntityPtr, MAX_DOWN> verts;
  std::array<utils::Point, 4> vertCoords;
  ElementCoordinateInfo info;
  for (auto& el : m_serialMesh->get_elements())
    if (el)
    {
      int nverts = get_downward(el, 0, verts.data());

      vertCoords[3] = utils::Point(0, 0, 0);
      for (int i = 0; i < nverts; ++i)
        vertCoords[i] = verts[i]->get_point_orig(0);

      RemoteSharedEntity remoteEntity = elOrigins(el, 0, 0);
      info.localId                    = remoteEntity.remoteId;
      info.vertCoords                 = vertCoords;

      exchanger.get_send_buffer(remoteEntity.remoteRank).push_back(info);
    }
}

void GatheredMeshCoordinateUpdate::unpack_recv_buffers(utils::impl::ParallelExchange<ElementCoordinateInfo>& exchanger)
{
  assert(m_parallelMesh);

  auto& buf = exchanger.get_recv_buffer(m_rootRank);
  std::array<MeshEntityPtr, MAX_DOWN> verts;
  for (const ElementCoordinateInfo& info : buf)
  {
    auto& el   = m_parallelMesh->get_elements()[info.localId];
    int nverts = get_downward(el, 0, verts.data());
    for (int i = 0; i < nverts; ++i)
      verts[i]->set_point_orig(0, info.vertCoords[i]);
  }
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
