#include "active_vert_container.hpp"
#include "mesh_entity.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void ActiveVertContainer::setup_remote_patches(const std::vector<opt::impl::ActiveVertData>& nonOwnedActiveVerts)
{
  Exchanger exchanger(m_mesh->get_comm());
  for (int phase=0; phase < 2; ++phase)
  {
    for (const auto& patch : nonOwnedActiveVerts)
    {
      set_patch_destinations(patch, exchanger);
    }

    if (phase == 0)
    {
      exchanger.allocate_send_buffers();
    }
  }

  exchanger.execute();
  get_remote_patch_contributions(exchanger);
  collect_remote_vertices();
  setup_send_comm_lists();
  update_remote_coords();
  set_remote_coords_orig();
}

void ActiveVertContainer::set_patch_destinations(const opt::impl::ActiveVertData& patch, Exchanger& exchanger)
{
  RemoteActiveVertData remoteData;
  for (int i=0; i < patch.get_num_verts(); ++i)
  {
    MeshEntityPtr v = patch.get_local_verts()[i];
    RemoteSharedEntity owner = get_owner_remote(m_mesh, v);
    remoteData.vertIds.push_back(owner);
  }

  for (int i=0; i < patch.get_num_elements(); ++i)
  {
    std::array<int, 3> triVerts = patch.get_element_vert_ids(i);
    for (int j=0; j < 3; ++j)
    {
      remoteData.triVertIndices.push_back(triVerts[j]);
    }
  }

  pack(exchanger, get_owner(m_mesh, patch.get_current_vert()), remoteData);
}

void ActiveVertContainer::get_remote_patch_contributions(Exchanger& exchanger)
{
  std::map<int, int> vertexToPatchIdx;  //TODO: consider using a Field for this?
  compute_vertex_to_patch_map(vertexToPatchIdx);

  for (int rank=0; rank < utils::impl::comm_size(m_mesh->get_comm()); ++rank)
  {
    auto& buf = exchanger.get_recv_buf(rank);
    while (buf.remaining() > 0)
    {
      RemoteActiveVertData remoteData = unpack(exchanger, rank);
      int patchIdx = vertexToPatchIdx[remoteData.vertIds[0].remoteId];
      create_local_verts_used_by_remote_patches(patchIdx, remoteData);
    }

    buf.reset();
  }

  for (int rank=0; rank < utils::impl::comm_size(m_mesh->get_comm()); ++rank)
  {
    while (exchanger.get_recv_buf(rank).remaining() > 0)
    {
      RemoteActiveVertData remoteData = unpack(exchanger, rank);
      int patchIdx = vertexToPatchIdx[remoteData.vertIds[0].remoteId];
      merge_patches(patchIdx, rank, remoteData);
    }
  }
}

void ActiveVertContainer::compute_vertex_to_patch_map(std::map<int, int>& vertexToPatch)
{
  vertexToPatch.clear();
  for (size_t i=0; i < m_activeVerts.size(); ++i)
    vertexToPatch[m_activeVerts[i].get_current_vert()->get_id()] = i;
}


void ActiveVertContainer::pack(Exchanger& exchanger, int destRank, const RemoteActiveVertData& remoteData)
{
  auto& buf = exchanger.get_send_buf(destRank);
  buf.pack(remoteData.vertIds);
  buf.pack(remoteData.triVertIndices);
}

ActiveVertContainer::RemoteActiveVertData ActiveVertContainer::unpack(Exchanger& exchanger, int sendRank)
{
  RemoteActiveVertData remoteData;
  auto& buf = exchanger.get_recv_buf(sendRank);
  buf.unpack(remoteData.vertIds);
  buf.unpack(remoteData.triVertIndices);

  return remoteData;
} 

void ActiveVertContainer::create_local_verts_used_by_remote_patches(int patchIdx, RemoteActiveVertData& remoteData)
{
  int myrank = utils::impl::comm_rank(m_mesh->get_comm());
  opt::impl::ActiveVertData& patch = m_activeVerts[patchIdx];
  for (size_t i=0; i < remoteData.vertIds.size(); ++i)
  {
    RemoteSharedEntity owner = remoteData.vertIds[i];
    if (owner.remoteRank == myrank)
    {
      bool foundExistingVert = false;
      for (int j=0; j < patch.get_num_local_verts(); ++j)
      {
        if (patch.get_vert_owner(j) == owner)
        {
          foundExistingVert = true;
          break;
        }
      }

      if (!foundExistingVert)
      {
        patch.add_local_vert(m_mesh->get_vertices()[owner.remoteId]);
      }
    }
  }

}

void ActiveVertContainer::merge_patches(int patchIdx, int /*senderRank*/, RemoteActiveVertData& remoteData)
{    
  // add vertices
  opt::impl::ActiveVertData& patch = m_activeVerts[patchIdx];

  std::vector<int> uniqueVertIdxsOnDest;
  for (size_t i=0; i < remoteData.vertIds.size(); ++i)
  {
    RemoteSharedEntity owner = remoteData.vertIds[i];

    bool foundExistingVert = false;
    for (int j=0; j < patch.get_num_verts(); ++j)
      if (patch.get_vert_owner(j) == owner)
      {
        uniqueVertIdxsOnDest.push_back(j);
        foundExistingVert = true;
        break;
      }

    if (!foundExistingVert)
    {
      patch.add_remote_vert(owner);
      uniqueVertIdxsOnDest.push_back(patch.get_num_verts() - 1);
    }
  }

  // add triangles
  for (size_t i=0; i < remoteData.triVertIndices.size(); i += 3)
  {
    std::array<int, 3> vertIdxs = {remoteData.triVertIndices[i],
                                   remoteData.triVertIndices[i + 1],
                                   remoteData.triVertIndices[i + 2]};

    for (int j=0; j < 3; ++j)
      vertIdxs[j] = uniqueVertIdxsOnDest[vertIdxs[j]];

    patch.add_remote_element(vertIdxs);
  }
}


void ActiveVertContainer::collect_remote_vertices()
{
  int myrank = utils::impl::comm_rank(m_mesh->get_comm());
  std::vector<VertexUse> remoteVertsInPatches;
  for (size_t i=0; i < m_activeVerts.size(); ++i)
  {
    opt::impl::ActiveVertData& patch = m_activeVerts[i];
    for (int j=patch.get_num_local_verts(); j < patch.get_num_verts(); ++j)
    {
      RemoteSharedEntity owner = patch.get_vert_owner(j);
      if (owner.remoteRank != myrank)
      {
        VertexUse use{owner, static_cast<int>(i), j};
        remoteVertsInPatches.push_back(use);
      }
    }
  }

  auto sortByOwner = [](const VertexUse& lhs, const VertexUse& rhs)
  {
    return lhs.owner < rhs.owner;
  };

  std::sort(remoteVertsInPatches.begin(), remoteVertsInPatches.end(), sortByOwner);

  size_t idx = 0;
  while (idx < remoteVertsInPatches.size())
  {
    RemoteSharedEntity owner = remoteVertsInPatches[idx].owner;
    VertexUses& allVertUses = m_recvVertexLists[owner.remoteRank].emplace_back();
    allVertUses.ownerLocalId = owner.remoteId;
    while (idx < remoteVertsInPatches.size() && remoteVertsInPatches[idx].owner == owner)
    {
      VertexUse& vertUse = remoteVertsInPatches[idx];
      allVertUses.patchAndVertIdxs.push_back(std::make_pair(vertUse.patchIdx, vertUse.vertIdxInPatch));
      idx++;
    }
  }
}

void ActiveVertContainer::setup_send_comm_lists()
{
  stk::DataExchangeUnknownPatternNonBlockingBuffer<int> exchanger(m_mesh->get_comm());

  for (size_t i=0; i < m_recvVertexLists.size(); ++i)
    for (size_t j=0; j < m_recvVertexLists[i].size(); ++j)
    {
      exchanger.get_send_buf(i).push_back(m_recvVertexLists[i][j].ownerLocalId);
    }

  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();

  auto unpack = [&](int rank, const std::vector<int>& buf)
  {
    m_sendVertexLists[rank] = std::move(buf);
  };

  exchanger.complete_receives(unpack);
}

void ActiveVertContainer::start_coord_update()
{
  m_exchanger.clear_send_bufs();
  m_exchanger.clear_recv_bufs();

  for (size_t rank=0; rank < m_recvVertexLists.size(); ++rank)
  {
    if (m_recvVertexLists[rank].size() > 0)
      m_exchanger.get_recv_buf(rank).resize(m_recvVertexLists[rank].size());

    for (auto& localId : m_sendVertexLists[rank])
      m_exchanger.get_send_buf(rank).push_back(m_mesh->get_vertices()[localId]->get_point_orig(0));
  }

  m_exchanger.start_nonblocking();
  m_sharedCoordUpdator.start_update();
}

void ActiveVertContainer::finish_coord_update()
{
  auto unpack = [&](int rank, const std::vector<utils::Point>& buf)
  {
    assert(buf.size() == m_recvVertexLists[rank].size());

    for (size_t i=0; i < m_recvVertexLists[rank].size(); ++i)
    {
      auto& vertexUses = m_recvVertexLists[rank][i];
      for (auto& patchAndVertIdxs : vertexUses.patchAndVertIdxs)
      {
        int patchIdx = patchAndVertIdxs.first;
        int vertIdx = patchAndVertIdxs.second;
        opt::impl::ActiveVertData& patch = m_activeVerts[patchIdx];
        if (vertIdx < patch.get_num_local_verts())
        {
          patch.get_local_verts()[vertIdx]->set_point_orig(0, buf[i]);
        } else
        {
          patch.get_unique_verts()[vertIdx] = buf[i];
        }
      }
    } 
  };

  m_exchanger.complete_receives(unpack);
  m_sharedCoordUpdator.finish_update();
  m_exchanger.complete_sends();
}

void ActiveVertContainer::set_remote_coords_orig()
{
  for (auto& patch : m_activeVerts)
    patch.finish_initialization();
}

}
}
}
}