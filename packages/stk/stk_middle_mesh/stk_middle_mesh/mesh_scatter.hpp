#ifndef STK_MIDDLE_MESH_MESH_SCATTER
#define STK_MIDDLE_MESH_MESH_SCATTER

#include "mesh_entity.hpp"
#include "mesh_scatter_spec.hpp"
#include "mesh.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "variable_size_field.hpp"
#include "field.hpp"
#include "entity_sorted_by_owner.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Class to be implemented
class MeshScatter
{
  public:
    // computeEntityCorrespondence: if true, get_entity_origins and get_entity_destinations will
    //                         return fields describeing the corespondence between entities on the 
    //                         input mesh and the scattered mesh
    MeshScatter(std::shared_ptr<impl::MeshScatterSpec> scatterSpec,
                std::shared_ptr<Mesh> inputMesh, MPI_Comm scatteredMeshComm,
                bool computeEntityCorrespondence=false);

    std::shared_ptr<Mesh> scatter();

    FieldPtr<RemoteSharedEntity> get_element_origins() const;

    VariableSizeFieldPtr<RemoteSharedEntity> get_entity_origins() const;

    VariableSizeFieldPtr<RemoteSharedEntity> get_entity_destinations() const;
    

  private:
    using EntityExchanger   = stk::DataExchangeUnknownPatternNonBlockingCommBuffer;
    using EntityIdExchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<int>;

    void check_mesh_comm_is_subset_of_union_comm();    

    void send_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr);

    void pack_verts(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger, 
                    EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_verts_and_pack_sharing(EntityExchanger& entityExchanger, 
                                   std::shared_ptr<EntityIdExchanger> sharingExchanger,
                                   EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_vert_buffer(int rank, stk::CommBuffer& buf, std::shared_ptr<EntityIdExchanger> sharingExchanger,
                            EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_vert_sharing(std::shared_ptr<EntityIdExchanger> sharingExchanger);

    void unpack_returned_entity_ids(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, int dim,
                                    EntityIdExchanger& entityCorrespondenceExchanger);

    void translate_union_comm_ranks_to_output_comm(const std::vector<int>& unionCommRanks, std::vector<int>& outputCommRanks);

    void unpack_vert_sharing_buffer(int rank, const std::vector<int>& buf);

    void send_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr);

    void pack_edges(VariableSizeFieldPtr<int> destRanksOnUnionCommPtr, EntityExchanger& entityExchanger,
                    EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_edges(EntityExchanger& entityExchanger, EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_edge_buffer(int rank, stk::CommBuffer& buf, EntityIdExchanger& entityCorrespondenceExchanger);

    void send_elements();

    void pack_elements(EntityExchanger& entityExchanger, EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_elements(EntityExchanger& exchanger, EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_element_buffer(int rank, stk::CommBuffer& buf, EntityIdExchanger& entityCorrespondenceExchanger);

    void unpack_returned_element_ids(EntityIdExchanger& entityCorrespondenceExchanger);


    MPI_Comm m_unionComm;
    std::shared_ptr<impl::MeshScatterSpec> m_scatterSpec;
    std::shared_ptr<Mesh> m_inputMesh;
    bool m_computeEntityCorrespondence;
    std::shared_ptr<Mesh> m_outputMesh;
    FieldPtr<RemoteSharedEntity> m_elementOrigins;
    VariableSizeFieldPtr<RemoteSharedEntity> m_entityOrigins;
    VariableSizeFieldPtr<RemoteSharedEntity> m_entityDestinations;
    std::shared_ptr<EntitySortedByOwner> m_vertsBySrcMeshOwner;
    std::shared_ptr<EntitySortedByOwner> m_edgesBySrcMeshOwner;
};




}
}
}
}

#endif