#include "create_sharing_from_verts.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

CreateSharingFromVert::CreateSharingFromVert(std::shared_ptr<Mesh> mesh, int dim) :
  m_mesh(mesh),
  m_dim(dim),
  m_exchanger(mesh->get_comm())
{}

void CreateSharingFromVert::create_sharing_from_verts()
{
  pack_buffers();
  m_exchanger.start_nonblocking();
  m_exchanger.post_nonblocking_receives();

  auto f = [&](int rank, const std::vector<int>& buf)
  {
    unpack_buffer(rank, buf);
  };

  m_exchanger.complete_receives(f);
  m_exchanger.complete_sends();
}


void CreateSharingFromVert::pack_buffers()
{
  std::vector<std::vector<RemoteSharedEntity>> remotesByRank;
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  for (auto& entity : m_mesh->get_mesh_entities(m_dim))
    if (entity)
    {
      int ndown = mesh::get_downward(entity, 0, verts.data());

      if (are_all_verts_shared(verts, ndown))
      {
        get_entity_vert_remotes(verts, ndown, remotesByRank);

        for (size_t i=0; i < remotesByRank.size(); ++i)
        {
          auto& remotes_i = remotesByRank[i];
          assert(remotes_i.size() == size_t(ndown));

          int destRank = remotes_i[0].remoteRank;
          m_exchanger.get_send_buf(destRank).push_back(remotes_i.size());
          for (auto& remote : remotes_i)
          {
            m_exchanger.get_send_buf(destRank).push_back(remote.remoteId);
          }

          m_exchanger.get_send_buf(destRank).push_back(entity->get_id());
        }
      }
    }
}

bool CreateSharingFromVert::are_all_verts_shared(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                                                       int ndown)
{
  for (int i=0; i < ndown; ++i)
  {
    if (verts[i]->count_remote_shared_entities() == 0)
    {
      return false;
    }
  }

  return true;
}

void CreateSharingFromVert::get_entity_vert_remotes(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                                                        int ndown,
                                                        std::vector<std::vector<RemoteSharedEntity>>& remotesByRank)
{
  std::vector<std::vector<RemoteSharedEntity>> allRemotes(ndown);
  auto compareRemoteEntityByRank = [](const RemoteSharedEntity& lhs, const RemoteSharedEntity& rhs)
  {
    return lhs.remoteRank < rhs.remoteRank;
  };  
  int minNumRemotes = std::numeric_limits<int>::max();
  int idxWithMinRemotes = -1;
  for (int i=0; i < ndown; ++i)
  {
    mesh::MeshEntityPtr vert = verts[i];
    if (vert->count_remote_shared_entities() < minNumRemotes)
    {
      minNumRemotes = vert->count_remote_shared_entities();
      idxWithMinRemotes = i;
    }

    idxWithMinRemotes = std::min(idxWithMinRemotes, vert->count_remote_shared_entities());
    for (int j=0; j < vert->count_remote_shared_entities(); ++j)
    {
      allRemotes[i].push_back(vert->get_remote_shared_entity(j));
    }

    std::sort(allRemotes[i].begin(), allRemotes[i].end(), compareRemoteEntityByRank);
  }

  remotesByRank.clear();
  std::vector<RemoteSharedEntity> remotesOnSameRank(ndown);

  for (int i=0; i < verts[idxWithMinRemotes]->count_remote_shared_entities(); ++i)
  {
    int rank = allRemotes[idxWithMinRemotes][i].remoteRank;
    bool foundOnAllVerts = true;
    for (int j=0; j < ndown; ++j)
    {
      if (j == idxWithMinRemotes)
      {
        remotesOnSameRank[j] = allRemotes[j][i];
      } else
      {
        auto it = std::lower_bound(allRemotes[j].begin(), allRemotes[j].end(), RemoteSharedEntity(rank, -1), compareRemoteEntityByRank);

        if (it != allRemotes[j].end() && it->remoteRank == rank)
        {
          remotesOnSameRank[j] = *it;
        } else
        {
          foundOnAllVerts = false;
          break;
        }
      }
    }

    if (foundOnAllVerts)
    {
      remotesByRank.push_back(remotesOnSameRank);
    }
  }
}

void CreateSharingFromVert::unpack_buffer(int rank, const std::vector<int>& buf)
{
  size_t idx = 0;
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  while (idx < buf.size())
  {
    int ndown = buf[idx++];
    for (int i=0; i < ndown; ++i)
    {
      verts[i] = m_mesh->get_vertices()[buf[idx++]];
    }
    int remoteEntityId = buf[idx++];

    mesh::MeshEntityPtr entity = get_common_entity(verts, ndown);
    if (entity)
    {
      entity->add_remote_shared_entity({rank, remoteEntityId});
    }
  }
}

mesh::MeshEntityPtr CreateSharingFromVert::get_common_entity(const std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN>& verts,
                                                                   int ndown)
{
  MeshEntityCompare comp;
  std::array<std::vector<mesh::MeshEntityPtr>, mesh::MAX_DOWN> vertEls;
  for (int i=0; i < ndown; ++i)
  {
    mesh::get_upward(verts[i], m_dim, vertEls[i]);
    std::sort(vertEls[i].begin(), vertEls[i].end(), comp);
  }

  std::vector<MeshEntityPtr> prevIntersection = vertEls[0];
  std::vector<MeshEntityPtr> output;
  for (int i=1; i < ndown; ++i)
  {
    output.clear();
    std::set_intersection(vertEls[i].begin(), vertEls[i].end(), 
                          prevIntersection.begin(), prevIntersection.end(),
                          std::back_inserter(output), comp);
    prevIntersection = output;
  }

  assert(prevIntersection.size() <= 1);

  return prevIntersection.size() == 1 ? prevIntersection[0] : nullptr;
}


}
}
}
}
