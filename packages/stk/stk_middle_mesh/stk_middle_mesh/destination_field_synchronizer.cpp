#include "destination_field_synchronizer.hpp"
#include "stk_util/util/SortAndUnique.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

DestinationFieldSynchronizer::DestinationFieldSynchronizer(std::shared_ptr<Mesh> mesh, std::shared_ptr<MeshScatterSpec> scatterSpec) :
  m_mesh(mesh),
  m_scatterSpec(scatterSpec)
{}

VariableSizeFieldPtr<int> DestinationFieldSynchronizer::synchronize()
{
  Exchanger vertExchanger(m_mesh->get_comm()), edgeExchanger(m_mesh->get_comm());
  auto fieldPtr = create_variable_size_field<int>(m_mesh, FieldShape(1, 1, 0));

  get_local_destinations_and_pack_buffers(vertExchanger, 0, nullptr);
  vertExchanger.allocate_send_buffers();
  get_local_destinations_and_pack_buffers(vertExchanger, 0, fieldPtr);
  vertExchanger.start_nonblocking();
  vertExchanger.post_nonblocking_receives();

  get_local_destinations_and_pack_buffers(edgeExchanger, 1, nullptr);
  edgeExchanger.allocate_send_buffers();
  get_local_destinations_and_pack_buffers(edgeExchanger, 1, fieldPtr); 
  edgeExchanger.start_nonblocking();
  edgeExchanger.post_nonblocking_receives();

  unpack_remote_element_destinations(vertExchanger, 0, fieldPtr);
  unpack_remote_element_destinations(edgeExchanger, 1, fieldPtr);
  sort_and_unique(fieldPtr);

  return fieldPtr;
}

void DestinationFieldSynchronizer::get_local_destinations_and_pack_buffers(Exchanger& exchanger,
                                                                       int dim, VariableSizeFieldPtr<int> fieldPtr)
{
  std::vector<int> destRanks;
  for (auto& entity : m_mesh->get_mesh_entities(dim))
    if (entity)
    {
      get_destinations_from_scatterspec(entity, destRanks);


      for (int i=0; i < entity->count_remote_shared_entities(); ++i)      
      {
        RemoteSharedEntity remote = entity->get_remote_shared_entity(i);
        exchanger.get_send_buf(remote.remoteRank).pack(remote.remoteId);
        exchanger.get_send_buf(remote.remoteRank).pack(destRanks);
      }

      if (fieldPtr)
      {
        for (auto& destRank : destRanks)
        {
          fieldPtr->insert(entity, 0, destRank);
        }
      }
    }
}

void DestinationFieldSynchronizer::get_destinations_from_scatterspec(MeshEntityPtr vert, std::vector<int>& destRanks)
{
  std::vector<int> destRanksForEl;
  destRanks.clear();
  std::vector<MeshEntityPtr> els;
  get_upward(vert, 2, els);
  for (auto& el : els)
  {
    m_scatterSpec->get_destinations(el, destRanksForEl);
    for (int destRank : destRanksForEl)
      destRanks.push_back(destRank);
  }
  stk::util::sort_and_unique(destRanks); 
}

void DestinationFieldSynchronizer::unpack_remote_element_destinations(Exchanger& exchanger, int dim, VariableSizeFieldPtr<int> fieldPtr)
{
  auto f = [&](int rank, stk::CommBuffer& buf)
  {
    unpack_buffer(rank, dim, buf, fieldPtr);
  };

  exchanger.complete_receives(f);
}

void DestinationFieldSynchronizer::unpack_buffer(int /*rank*/, int dim, stk::CommBuffer& buf, VariableSizeFieldPtr<int> fieldPtr)
{
  auto& field = *fieldPtr;
  std::vector<int> destRanks;
  while (buf.remaining() > 0)
  {
    int localId;
    destRanks.clear();
    buf.unpack(localId);
    buf.unpack(destRanks);

    MeshEntityPtr entity = m_mesh->get_mesh_entities(dim)[localId];
    for (auto& destRank : destRanks)
    {  
      field.insert(entity, 0, destRank);
    }
  }
}

void DestinationFieldSynchronizer::sort_and_unique(VariableSizeFieldPtr<int> fieldPtr)
{
  auto& field = *fieldPtr;
  for (int dim=0; dim < 2; ++dim)
    for (auto& entity : m_mesh->get_mesh_entities(dim))
      if (entity)
      {
        auto begin = field(entity, 0).begin();
        auto end   = field(entity, 0).end();
        std::sort(begin, end);
        auto newEnd = std::unique(begin, end);

        field.resize(entity, 0, newEnd - begin);     
      }
}


}
}
}
}
