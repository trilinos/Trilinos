#ifndef STK_MIDDLE_MESH_MESH_SCATTER
#define STK_MIDDLE_MESH_MESH_SCATTER

#include "mesh_scatter_spec.hpp"
#include "mesh.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "variable_size_field.hpp"
#include "field.hpp"
#include "edge_sharing_from_verts.hpp"
#include "destination_field_gatherer.hpp"
#include "entity_sorted_by_owner.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Class to be implemented
class MeshScatter
{
  public:
    MeshScatter(std::shared_ptr<impl::MeshScatterSpec> scatterSpec,
                std::shared_ptr<Mesh> inputMesh, MPI_Comm scatteredMeshComm);

    std::shared_ptr<mesh::Mesh> scatter();

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_element_origins();
    

  private:
    using EntityExchanger      = stk::DataExchangeUnknownPatternNonBlockingCommBuffer;
    using VertSharingExchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<int>;

    void check_mesh_comm_is_subset_of_union_comm();


    int compute_dest_comm_size();
    

    void send_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr);

    void pack_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger);

    void unpack_verts_and_pack_sharing(EntityExchanger& entityExchanger, 
                                   std::shared_ptr<VertSharingExchanger> sharingExchanger);

    void unpack_vert_buffer(int rank, stk::CommBuffer& buf, std::shared_ptr<VertSharingExchanger> sharingExchanger);

    void unpack_vert_sharing(std::shared_ptr<VertSharingExchanger> sharingExchanger);

    void translate_union_comm_ranks_to_output_comm(const std::vector<int>& unionCommRanks, std::vector<int>& outputCommRanks);

    void unpack_vert_sharing_buffer(int rank, const std::vector<int>& buf);

    void send_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr);


    void pack_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger);

    void unpack_edges(EntityExchanger& entityExchanger);

    void unpack_edge_buffer(int rank, stk::CommBuffer& buf);

    void send_elements();

    void pack_elements(EntityExchanger& entityExchanger);

    void unpack_elements(EntityExchanger& exchanger);

    void unpack_element_buffer(int rank, stk::CommBuffer& buf);

    MPI_Comm m_unionComm;
    std::shared_ptr<impl::MeshScatterSpec> m_scatterSpec;
    std::shared_ptr<Mesh> m_inputMesh;
    std::shared_ptr<Mesh> m_outputMesh;
    FieldPtr<RemoteSharedEntity> m_elementOrigins;
    std::shared_ptr<EntitySortedByOwner> m_vertsBySrcMeshOwner;
    std::shared_ptr<EntitySortedByOwner> m_edgesBySrcMeshOwner;
};




}
}
}
}

#endif