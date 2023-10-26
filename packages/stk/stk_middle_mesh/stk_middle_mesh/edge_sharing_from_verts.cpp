#include "edge_sharing_from_verts.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

EdgeSharingCreatorFromVerts::EdgeSharingCreatorFromVerts(std::shared_ptr<Mesh> mesh) :
  m_mesh(mesh),
  m_exchanger(mesh->get_comm())
{}

void EdgeSharingCreatorFromVerts::create_sharing_from_verts()
{
  pack_buffers();
  m_exchanger.start_nonblocking();

  auto f = [&](int rank, const std::vector<int>& buf)
  {
    unpack_buffer(rank, buf);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


void EdgeSharingCreatorFromVerts::pack_buffers()
{
  std::vector<std::pair<RemoteSharedEntity, RemoteSharedEntity>> remotes;
  for (auto& edge : m_mesh->get_edges())
    if (edge)
    {
      MeshEntityPtr v1 = edge->get_down(0);
      MeshEntityPtr v2 = edge->get_down(1);

      if (v1->count_remote_shared_entities() > 0 && v2->count_remote_shared_entities() > 0)
      {
        get_edge_vert_remotes(v1, v2, remotes);


        for (size_t i=0; i < remotes.size(); ++i)
        {
          auto& remotes_i = remotes[i];
          int destRank = remotes_i.first.remoteRank;
          m_exchanger.get_send_buf(destRank).push_back(remotes_i.first.remoteId);
          m_exchanger.get_send_buf(destRank).push_back(remotes_i.second.remoteId);
          m_exchanger.get_send_buf(destRank).push_back(edge->get_id());

          for (int j=0; j < 3; ++j)
            m_exchanger.get_recv_buf(destRank).push_back(-1);
        }
      }
    }
}

void EdgeSharingCreatorFromVerts::get_edge_vert_remotes(MeshEntityPtr v1, MeshEntityPtr v2,
                                    std::vector<std::pair<RemoteSharedEntity, RemoteSharedEntity>>& remotes)
{
  remotes.clear();
  for (int i=0; i < v1->count_remote_shared_entities(); ++i)
    for (int j=0; j < v2->count_remote_shared_entities(); ++j)
    {
      RemoteSharedEntity remote1 = v1->get_remote_shared_entity(i);
      RemoteSharedEntity remote2 = v2->get_remote_shared_entity(j);
      if (remote1.remoteRank == remote2.remoteRank)
      {
        remotes.push_back(std::make_pair(remote1, remote2));
      }
    }

  //throw std::runtime_error("failed to find remote rank for edge");
}

void EdgeSharingCreatorFromVerts::unpack_buffer(int rank, const std::vector<int>& buf)
{
  assert(buf.size() % 3 == 0);
  for (size_t i=0; i < buf.size(); i += 3)
  {
    MeshEntityPtr v1 = m_mesh->get_vertices()[buf[i]];
    MeshEntityPtr v2 = m_mesh->get_vertices()[buf[i+1]];
    int edgeRemoteId = buf[i+2];

    MeshEntityPtr edge = get_common_edge(v1, v2);
    edge->add_remote_shared_entity({rank, edgeRemoteId});
  }
}


}
}
}
}
