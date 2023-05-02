#ifndef STK_MIDDLE_MESH_REMOTE_COORDINATE_UPDATOR
#define STK_MIDDLE_MESH_REMOTE_COORDINATE_UPDATOR

#include "mesh.hpp"
#include "utils.hpp"

#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// overwrites the coordinates of shared entities with
// the coordinates of the entity on the owning process
class RemoteCoordinateUpdator
{
  public:
    explicit RemoteCoordinateUpdator(std::shared_ptr<Mesh> mesh) :
      m_mesh(mesh),
      m_exchanger(mesh->get_comm())
    {}

    void start_update();

    void finish_update();
    
  private:
    struct IdAndCoords
    {
      int localId;
      utils::Point pt;
    };

    void pack_buffers();

    void unpack_buffer(int rank, const std::vector<IdAndCoords>& buf);

    std::shared_ptr<Mesh> m_mesh;
    stk::DataExchangeKnownPatternNonBlockingBuffer<IdAndCoords> m_exchanger;
};

}
}
}
}

#endif