#ifndef MESH_GATHER_TO_ROOT
#define MESH_GATHER_TO_ROOT

#include "field.hpp"
#include "mesh.hpp"
#include "parallel_exchange.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// Gathers a parallel mesh onto a single process, the specified rank of the specified
// communicator.  The specified communicator (input_comm) must contain all the processes
// in the mesh communicator.  The constructor for this class is collective on input_comm.
// The member functions are collective on the union of the specified rank of input_comm and
// the mesh communicator.  Processes on input_comm but not the mesh communicator should
// pass in a nullptr for the mesh.
class MeshGatherToRoot
{
  public:
    MeshGatherToRoot(MPI_Comm inputComm, int rootRankOnInputComm, std::shared_ptr<Mesh> mesh);

    std::shared_ptr<Mesh> gather();

    // the ranks are the ranks on the input_comm
    FieldPtr<RemoteSharedEntity> get_element_origins();

  private:
    void check_mesh_comm_is_subset_of_input_comm(MPI_Comm inputComm);

    struct MeshEntityCounts
    {
        std::array<int, 3> entityCounts = {0, 0, 0};
        std::array<int, 3> entityMaxIds = {0, 0, 0};
        int rankOnMeshComm              = -1;
    };

    void send_sizes_to_root();

    template <typename T>
    void send_entities(int tag);

    struct VertexInfo
    {
        const static int DIMENSION = 0;
        int ownerRank              = -1;
        int ownerLocalId           = -1;
        utils::Point coords;
    };

    void get_info(MeshEntityPtr vert, int myrank, VertexInfo& vertInfo);

    void unpack_info(utils::impl::ParallelExchange<VertexInfo>& exchanger);

    struct EdgeInfo
    {
        const static int DIMENSION = 1;
        int edgeOwnerRank          = -1;
        int edgeOwnerLocalId       = -1;
        int vert1OwnerRank         = -1;
        int vert1OwnerLocalId      = -1;
        int vert2OwnerRank         = -1;
        int vert2OwnerLocalId      = -1;
    };

    void get_info(MeshEntityPtr edge, int myrank, EdgeInfo& edgeInfo);

    void unpack_info(utils::impl::ParallelExchange<EdgeInfo>& exchanger);

    struct ElementInfo
    {
        const static int DIMENSION        = 2;
        std::array<int, 4> edgeOwnerRanks = {-1, -1, -1, -1};
        std::array<int, 4> edgeOwnerIds   = {-1, -1, -1, -1};
        EntityOrientation edge1Orient;
        int elementLocalId;
    };

    void get_info(MeshEntityPtr el, int myrank, ElementInfo& elInfo);

    void unpack_info(utils::impl::ParallelExchange<ElementInfo>& exchanger);

    MPI_Comm m_inputComm;
    MPI_Comm m_comm;
    std::shared_ptr<Mesh> m_inputMesh;
    std::shared_ptr<Mesh> m_outputMesh;
    std::vector<MeshEntityCounts> m_entityCounts;
    std::vector<int> m_splitCommToMeshCommRank;
    std::vector<std::vector<MeshEntityPtr>> m_ownerLocalIdToOutputVertex;
    std::vector<std::vector<MeshEntityPtr>> m_ownerLocalIdToOutputEdge;
    FieldPtr<RemoteSharedEntity> m_elementOrigins;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif