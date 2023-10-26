#ifndef STK_MIDDLE_MESH_EDGE_SHARING_FROM_VERT_H
#define STK_MIDDLE_MESH_EDGE_SHARING_FROM_VERT_H

#include "mesh.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Given a mesh with the vertex RemoteSharedEntity information set,
// compute the edge RemoteSharedEntity.
class EdgeSharingCreatorFromVerts
{
  public:
    EdgeSharingCreatorFromVerts(std::shared_ptr<Mesh> mesh);

    void create_sharing_from_verts();

  private:
    using Exchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<int>;

    void pack_buffers();

    void get_edge_vert_remotes(MeshEntityPtr v1, MeshEntityPtr v2,
              std::vector<std::pair<RemoteSharedEntity, RemoteSharedEntity>>& remotes);

    void unpack_buffer(int rank, const std::vector<int>& buf);

    std::shared_ptr<Mesh> m_mesh;
    Exchanger m_exchanger;
};

}
}
}
}

#endif