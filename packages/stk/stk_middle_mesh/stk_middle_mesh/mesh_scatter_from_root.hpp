#ifndef MESH_SCATTER_FROM_ROOT
#define MESH_SCATTER_FROM_ROOT

#include "field.hpp"
#include "mesh.hpp"

#include "parallel_exchange.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {


class MeshScatterFromRoot
{
  public:
    // input_comm is a communicator that is at least the union of
    // input_mesh->get_comm() and mesh_comm.  It can include other,
    // unrelated processes as well.
    // input_mesh is a serial mesh
    // mesh_comm is the communicator the newly-created mesh will be
    // created on.  If the process the input_mesh is defined on
    // is not part of the mesh_comm, it *must* pass in MPI_COMM_NULL
    // for mesh_comm
    // element_destination_rank gives the rank on input_comm that
    // each element in mesh_in will be sent to
    MeshScatterFromRoot(MPI_Comm unionComm, std::shared_ptr<Mesh> inputMesh, MPI_Comm meshComm,
                        FieldPtr<int> elementDestinationRanks);

    std::shared_ptr<Mesh> scatter();

    VariableSizeFieldPtr<RemoteSharedEntity> get_entity_destinations();

  private:
    void verify_only_one_root_debug_only(MPI_Comm unionComm, bool amIRoot);

    using Exchanger = DataExchangeUnknownPatternNonBlockingCommBuffer;

    template <typename T>
    void send_entities();

    void reset_for_new_pack(int dimension);

    template <typename T>
    void pack_data(Exchanger& exchanger, const T& sendData);

    template <typename T>
    void unpack_data(int rank, CommBuffer& buf);

    struct VertexInfo
    {
        utils::Point pt;
        std::vector<int> otherRanks;
        std::vector<int> otherLocalIds;
    };

    struct VertexSendData
    {
        const static int DIMENSION = 0;
        using InfoType             = VertexInfo;
        utils::Point pt;
        std::vector<int> ranks;
        std::vector<int> localIds;
    };

    void get_send_data(MeshEntityPtr vert, VertexSendData& sendData);

    void get_send_data_for_destination(const VertexSendData& sendData, int rankIdx, VertexInfo& info);

    void pack_buffer(CommBuffer& buf, const VertexInfo& info);

    void unpack_buffer(CommBuffer& buf, VertexInfo& info);

    void unpack_info(const VertexInfo& info);

    int translate_input_comm_rank_to_mesh_comm_rank(int unionCommRank);

    struct EdgeInfo
    {
        int vert1LocalId;
        int vert2LocalId;
        int remoteId   = -1;
        int remoteRank = -1;
    };

    struct EdgeSendData
    {
        const static int DIMENSION = 1;
        using InfoType             = EdgeInfo;
        std::vector<int> vert1LocalIds;
        std::vector<int> vert2LocalIds;
        std::vector<int> edgeLocalIds;
        std::vector<int> ranks;
    };

    void get_send_data(MeshEntityPtr edge, EdgeSendData& sendData);

    int get_local_id_on_rank(MeshEntityPtr entity, int rankOnMeshComm);

    void get_send_data_for_destination(const EdgeSendData& sendData, int rankIdx, EdgeInfo& info);

    void pack_buffer(CommBuffer& buf, const EdgeInfo& info);

    void unpack_buffer(CommBuffer& buf, EdgeInfo& info);

    void unpack_info(const EdgeInfo& info);

    struct ElementInfo
    {
        int edge1LocalId                   = -1;
        int edge2LocalId                   = -1;
        int edge3LocalId                   = -1;
        int edge4LocalId                   = -1;
        EntityOrientation edge1Orientation = EntityOrientation::Standard;
    };

    struct ElementSendData
    {
        const static int DIMENSION = 2;
        using InfoType             = ElementInfo;
        std::vector<int> ranks;
        int edge1LocalId                   = -1;
        int edge2LocalId                   = -1;
        int edge3LocalId                   = -1;
        int edge4LocalId                   = -1;
        EntityOrientation edge1Orientation = EntityOrientation::Standard;
    };

    void get_send_data(MeshEntityPtr el, ElementSendData& sendData);

    void get_send_data_for_destination(const ElementSendData& sendData, int rankIdx, ElementInfo& info);

    void pack_buffer(CommBuffer& buf, const ElementInfo& info);

    void unpack_buffer(CommBuffer& buf, ElementInfo& info);

    void unpack_info(const ElementInfo& info);

    bool m_amIRootRank;
    bool m_amIOnMeshComm;
    MPI_Comm m_unionComm;
    std::shared_ptr<Mesh> m_inputMesh;
    std::shared_ptr<Mesh> m_outputMesh;
    FieldPtr<int> m_elementDestinationRanks;
    std::vector<int> m_localIds;
    VariableSizeFieldPtr<RemoteSharedEntity> m_entityDestinations;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
