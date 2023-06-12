#include "active_vert_container.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void ActiveVertContainer::set_patch_destinations(opt::impl::ActiveVertData& patch, Exchanger& exchanger)
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
  std::map<int, int> vertexToPatchIdx;
  compute_vertex_to_patch_map(vertexToPatchIdx);
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

void ActiveVertContainer::merge_patches(int patchIdx, int senderRank, RemoteActiveVertData& remoteData)
{
  // add vertices
  opt::impl::ActiveVertData& patch = m_activeVerts[patchIdx];
  PatchRemoteInfo& remoteMapping   = m_recvPatchMapping[senderRank].emplace_back(patchIdx);

  for (size_t i=0; i < remoteData.vertIds.size(); ++i)
  {
    RemoteSharedEntity owner = remoteData.vertIds[i];

    bool foundExistingVert = false;
    for (int j=0; j < patch.get_num_verts(); ++j)
      if (patch.get_vert_owner(j) == owner)
      {
        remoteMapping.uniqueVertIdxsOnDest.push_back(j);
        foundExistingVert = true;
        break;
      }

    if (!foundExistingVert)
    {
      patch.add_remote_vert(owner);
      remoteMapping.uniqueVertIdxsOnDest.push_back(patch.get_num_verts() - 1);
    }
  }

  // add triangles
  for (size_t i=0; i < remoteData.triVertIndices.size(); i += 3)
  {
    std::array<int, 3> vertIdxs = {remoteData.triVertIndices[i],
                                   remoteData.triVertIndices[i + 1],
                                   remoteData.triVertIndices[i + 2]};

    for (int j=0; j < 3; ++j)
      vertIdxs[j] = remoteMapping.uniqueVertIdxsOnDest[vertIdxs[j]];

    patch.add_remote_element(vertIdxs);
  }
}

void ActiveVertContainer::start_coord_update()
{    
  m_coordExchanger.clear_recv_bufs();
  m_coordExchanger.clear_send_bufs();
  m_coordUpdator.start_update();

  for (size_t rank=0; rank < m_recvPatchMapping.size(); ++rank)
  {
    int numVerts = 0;
    for (auto& mapping : m_recvPatchMapping[rank])
      numVerts += mapping.uniqueVertIdxsOnDest.size();

    m_coordExchanger.get_recv_buf(rank).resize(numVerts);
  }

  for (opt::impl::ActiveVertData& patch : m_nonOwnedActiveVerts)
  {
    MeshEntityPtr v = patch.get_current_vert();
    int ownerRank = get_owner(m_mesh, v);
    for (int i=0; i < patch.get_num_verts(); ++i)
    {
      m_coordExchanger.get_send_buf(ownerRank).push_back(patch.get_unique_verts()[i]);
    }
  }

  m_coordExchanger.start_nonblocking();
}

void ActiveVertContainer::finish_coord_update()
{
  m_coordUpdator.finish_update();

  auto f = [&](int rank, const std::vector<utils::Point>& vertcoords)
  {
    unpack_buffer(rank, vertcoords);
  };

  m_coordExchanger.complete_receives(f);
  m_coordExchanger.complete_sends();
}

void ActiveVertContainer::unpack_buffer(int rank, const std::vector<utils::Point>& vertCoords)
{
  auto& mappings = m_recvPatchMapping[rank];

  int recvBufIdx = 0;
  for (PatchRemoteInfo& mapping : mappings)
  {
    opt::impl::ActiveVertData& patch = m_activeVerts[mapping.patchIdx];
    for (size_t i=0; i < mapping.uniqueVertIdxsOnDest.size(); ++i)
    {
      int patchVertIdx = mapping.uniqueVertIdxsOnDest[i];
      utils::Point pt = vertCoords[recvBufIdx++];

      if (patch.get_vert_owner(patchVertIdx).remoteRank == rank)
      {
        patch.get_unique_verts()[patchVertIdx] = pt;
      }
    }
  }
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