#ifndef STK_MIDDLE_MESH_DESTINATION_FIELD_GATHERER_H
#define STK_MIDDLE_MESH_DESTINATION_FIELD_GATHERER_H

#include "mesh_scatter_spec.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "mesh.hpp"
#include "variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {


class DestinationFieldSynchronizer
{
  using Exchanger = stk::DataExchangeUnknownPatternNonBlockingCommBuffer;

  public:
    DestinationFieldSynchronizer(std::shared_ptr<Mesh> mesh, std::shared_ptr<MeshScatterSpec> scatterSpec);

    VariableSizeFieldPtr<int> synchronize();

  private:

    void get_local_destinations_and_pack_buffers(Exchanger& exchanger, int dim, VariableSizeFieldPtr<int> fieldPtr);

    void get_destinations_from_scatterspec(MeshEntityPtr vert, std::vector<int>& destRanks);

    void unpack_remote_element_destinations(Exchanger& exchanger, int dim, VariableSizeFieldPtr<int> fieldPtr);

    void unpack_buffer(int rank, int dim, stk::CommBuffer& buf, VariableSizeFieldPtr<int> fieldPtr);

    void sort_and_unique(VariableSizeFieldPtr<int> fieldPtr);

  private:
    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<MeshScatterSpec> m_scatterSpec;
};


}
}
}
}

#endif