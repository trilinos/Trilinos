#ifndef GATHERED_MESH_COORDINATE_UPDATE_H
#define GATHERED_MESH_COORDINATE_UPDATE_H

#include "field.hpp"
#include "mesh.hpp"
#include "parallel_exchange.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Given a serial mesh that was created by MeshGatherToRoot on the parallel_mesh, overwrites
// the coordinates of parallel_mesh with those from serial_mesh
class GatheredMeshCoordinateUpdate
{
  public:
    GatheredMeshCoordinateUpdate(MPI_Comm inputComm, std::shared_ptr<Mesh> serialMesh,
                                 FieldPtr<RemoteSharedEntity> serialMeshElementOrigins,
                                 std::shared_ptr<Mesh> parallelMesh);

    void update();

  private:
    void check_only_one_root();

    void check_serial_mesh_comm_size();

    void broadcast_root_process();

    void check_element_count();

    struct ElementCoordinateInfo
    {
        int localId;
        std::array<utils::Point, 4> vertCoords;
    };

    void send_coordinates();

    void pack_send_buffers(utils::impl::ParallelExchange<ElementCoordinateInfo>& exchanger);

    void unpack_recv_buffers(utils::impl::ParallelExchange<ElementCoordinateInfo>& exchanger);

    MPI_Comm m_inputComm;
    std::shared_ptr<Mesh> m_serialMesh;
    FieldPtr<RemoteSharedEntity> m_serialMeshElementOrigins;
    std::shared_ptr<Mesh> m_parallelMesh;
    bool m_amIRoot;
    int m_rootRank;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif