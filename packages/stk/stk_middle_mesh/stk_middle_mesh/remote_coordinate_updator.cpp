#include "stk_middle_mesh/remote_coordinate_updator.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void RemoteCoordinateUpdator::start_update()
{
  m_exchanger.clear_recv_bufs();
  m_exchanger.clear_send_bufs();
  pack_buffers();
  m_exchanger.start_nonblocking();
}

void RemoteCoordinateUpdator::finish_update()
{
  auto f = [&](int rank, const std::vector<IdAndCoords>& buf)
  {
    unpack_buffer(rank, buf);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


void RemoteCoordinateUpdator::pack_buffers()
{
  int myRank = utils::impl::comm_rank(m_mesh->get_comm());
  for (auto& vert : m_mesh->get_vertices())
    if (vert)
    {
      RemoteSharedEntity owner = get_owner_remote(m_mesh, vert);
      if (owner.remoteRank == myRank)
      {
        for (int i=0; i < vert->count_remote_shared_entities(); ++i)
        {
          RemoteSharedEntity remote = vert->get_remote_shared_entity(i);
          m_exchanger.get_send_buf(remote.remoteRank).push_back(IdAndCoords{remote.remoteId, vert->get_point_orig(0)});
        }
      } else
      {
        m_exchanger.get_recv_buf(owner.remoteRank).push_back(IdAndCoords{-1, {0, 0, 0}});
      }
    }
}

void RemoteCoordinateUpdator::unpack_buffer(int /*rank*/, const std::vector<IdAndCoords>& buf)
{
  for (const IdAndCoords& data : buf)
  {
    MeshEntityPtr vert = m_mesh->get_vertices()[data.localId];
    vert->set_point_orig(0, data.pt);
  }
}

}
}
}
}