#ifndef STK_MIDDLE_MESH_EDGE_SHARING_FROM_VERT_H
#define STK_MIDDLE_MESH_EDGE_SHARING_FROM_VERT_H

#include "mesh.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Given a mesh with the vertex RemoteSharedEntity information set,
// compute the RemoteSharedEntity info for the given dimension entity.
class CreateSharingFromVert
{
  public:
    CreateSharingFromVert(std::shared_ptr<Mesh> mesh, int dim=1);

    void create_sharing_from_verts();

  private:
    using Exchanger = stk::DataExchangeUnknownPatternNonBlockingBuffer<int>;

    void pack_buffers();

    bool are_all_verts_shared(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                              int ndown);

    void get_entity_vert_remotes(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                                 int ndown,
                                 std::vector<std::vector<RemoteSharedEntity>>& remotesByRank);

    void unpack_buffer(int rank, const std::vector<int>& buf);

    MeshEntityPtr get_common_entity(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                                    int ndown);

    std::shared_ptr<Mesh> m_mesh;
    int m_dim;
    Exchanger m_exchanger;
};

}
}
}
}

#endif