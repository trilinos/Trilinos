#include "mesh_layers.hpp"
#include "mesh_entity.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshLayers::initialize_que(std::vector<MeshEntityPtr>& roots, std::vector<MeshEntityPtr>& que)
{
  for (auto& v : roots)
  {
    que.push_back(v);
    mark_entity_seen(v);
  }
}

void MeshLayers::add_vert_to_output(MeshEntityPtr v, std::vector<MeshEntityPtr>& output)
{
  output.push_back(v);
  for (int i=0; i < v->count_remote_shared_entities(); ++i)
  {
    RemoteSharedEntity remote = v->get_remote_shared_entity(i);
    m_exchanger.get_send_buf(remote.remoteRank).push_back(RemoteQueueUpdate{remote.remoteId, AddToQue::Output});
  }
}


void MeshLayers::que_adjacent_verts(MeshEntityPtr v, std::vector<MeshEntityPtr>& que, std::set<int>* queIds)
{
  for (int j = 0; j < v->count_up(); ++j)
  {
    auto edge            = v->get_up(j);
    MeshEntityPtr vOther = edge->get_down(0) == v ? edge->get_down(1) : edge->get_down(0);
    if (!is_entity_seen(vOther))
    {
      que.push_back(vOther);
      mark_entity_seen(vOther);

      if (queIds)
        queIds->insert(vOther->get_id());

      for (int i=0; i < vOther->count_remote_shared_entities(); ++i)
      {
        RemoteSharedEntity remote = vOther->get_remote_shared_entity(i);
        m_exchanger.get_send_buf(remote.remoteRank).push_back(RemoteQueueUpdate{remote.remoteId, AddToQue::Next});
      }
    }
  }
}

void MeshLayers::process_remote_updates(std::vector<MeshEntityPtr>& queNext, std::vector<MeshEntityPtr>& output)
{
  //TODO: could also use a Field<Bool>, which is (maybe) more memory but faster access
  //TODO: for the output queue, we reconstruct the set once per layer, even though
  //      it grows monotonically
  std::set<int> queNextIds, outputIds;

  for (auto& v : output)
    outputIds.insert(v->get_id());

  for (auto& v : queNext)
    queNextIds.insert(v->get_id());    

  m_exchanger.execute();
  for (int rank=0; rank < utils::impl::comm_size(m_mesh->get_comm()); ++rank)
  {
    auto& recvBuf = m_exchanger.get_recv_buf(rank);

    for (const RemoteQueueUpdate& remoteUpdate : recvBuf)
    {
      MeshEntityPtr v = m_mesh->get_vertices()[remoteUpdate.localId];
      if (remoteUpdate.queue == AddToQue::Output && outputIds.count(v->get_id()) == 0)
      {
        output.push_back(v);
        que_adjacent_verts(v, queNext, &queNextIds);
        outputIds.insert(v->get_id());
      } else if (remoteUpdate.queue == AddToQue::Next && queNextIds.count(v->get_id()) == 0)
      {
        queNext.push_back(v);
        queNextIds.insert(v->get_id());
      }
      mark_entity_seen(v);
    }
  }

  m_exchanger.clear_recv_bufs();
  m_exchanger.clear_send_bufs();
}                                        

void MeshLayers::mark_entity_seen(MeshEntityPtr e)
{
  (*m_seenEntities)(e, 0, 0) = true;
}

bool MeshLayers::is_entity_seen(MeshEntityPtr e)
{
  return (*m_seenEntities)(e, 0, 0);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
