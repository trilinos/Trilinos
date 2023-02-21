#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHBUILDER_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHBUILDER_HPP_
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_TopologyData.hpp>

namespace krino {

class AuxMetaData;
class Phase_Support;

template<stk::topology::topology_t TOPO>
class StkMeshBuilder
{
public:
    static constexpr unsigned DIM = TopologyData<TOPO>::spatial_dimension();
    static constexpr unsigned NPE = TopologyData<TOPO>::num_nodes();

    StkMeshBuilder(stk::mesh::BulkData & mesh, const stk::ParallelMachine comm);

    void build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
        const std::vector<std::vector<std::array<unsigned,NPE>>> &elemConnPerProc,
        const unsigned blockId=1u);

    void build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
        const std::vector<std::array<unsigned, NPE>> &elementConn,
        const std::vector<unsigned> &elementBlockIDs,
        const std::vector<int> &specifiedElementProcOwners = {});

    void build_mesh_with_all_needed_block_ids(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
        const std::vector<std::array<unsigned, NPE>> &elementConn,
        const std::vector<unsigned> &elementBlockIDs,
        const std::vector<unsigned> &allBlocksIncludingThoseThatDontHaveElements,
        const std::vector<int> &specifiedElementProcOwners);

    void build_mesh_nodes_and_elements(
        const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
        const std::vector<std::array<unsigned, NPE>> &elementConn,
        const std::vector<unsigned> &elementBlockIDs,
        const std::vector<int> &specifiedElementProcOwners);

    const std::vector<stk::mesh::Entity> & get_owned_elements() const { return mOwnedElems; }
    const std::vector<stk::mesh::EntityId> & get_assigned_node_global_ids() const { return mAssignedGlobalNodeIdsforAllNodes; }
    const std::vector<stk::mesh::EntityId> & get_assigned_element_global_ids() const { return mAssignedGlobalElementIdsforAllElements; }
    std::vector<stk::mesh::EntityId> get_ids_of_elements_with_given_indices(const std::vector<unsigned> & elemIndices) const;

    bool check_boundary_sides() const;
    bool check_block_boundary_sides() const;

    AuxMetaData & get_aux_meta() { return mAuxMeta; }
    const AuxMetaData & get_aux_meta() const { return mAuxMeta; }

    Phase_Support & get_phase_support() { return mPhaseSupport; }
    const Phase_Support & get_phase_support() const { return mPhaseSupport; }

    void create_sideset_parts(const std::vector<unsigned> &sidesetIds);
    void add_sides_to_sidesets(const std::vector<stk::mesh::Entity> &sides, const std::vector<std::vector<unsigned>> &sidesetIdsPerSide);
    stk::mesh::Entity get_side_with_nodes(const std::vector<stk::mesh::Entity> &nodesOfSide) const;
    void create_block_parts(const std::vector<unsigned> &elementBlockIDs);
    stk::math::Vector3d get_node_coordinates(const stk::mesh::Entity node) const;

private:
    stk::mesh::BulkData & mMesh;
    AuxMetaData & mAuxMeta;
    Phase_Support & mPhaseSupport;
    const stk::ParallelMachine mComm;
    std::vector<stk::mesh::EntityId> mAssignedGlobalNodeIdsforAllNodes;
    std::vector<stk::mesh::EntityId> mAssignedGlobalElementIdsforAllElements;
    std::vector<stk::mesh::Entity> mOwnedElems;

    void declare_coordinates();

    void create_boundary_sides();
    void create_block_boundary_sides();

    void set_node_coordinates(const stk::mesh::Entity node, const stk::math::Vector3d &newLoc);
    stk::mesh::Entity create_node(const stk::math::Vector3d &loc, const std::vector<int> &sharingProcs, stk::mesh::EntityId nodeId);
    stk::mesh::Entity create_element(const std::vector<stk::mesh::Entity> &nodes, stk::mesh::EntityId elementId, unsigned blockId);

    std::vector<stk::mesh::Entity> create_parallel_nodes(const std::vector<stk::math::Vec<double,DIM>>& nodeLocs,
        const std::map<unsigned,std::vector<int>> &nodeIndicesWithSharingProcs,
        const std::vector<stk::mesh::EntityId> & assignedGlobalNodeIdsforAllNodes);

    std::vector<stk::mesh::Entity> create_parallel_elements(const std::vector<std::array<unsigned, NPE>> &elementConn,
        const std::vector<unsigned> &elementBlockIDs,
        const std::vector<int> &elementProcOwners,
        const std::vector<stk::mesh::Entity>& nodesWhichAreValidIfTheyExistOnProc,
        const std::vector<stk::mesh::EntityId> & assignedGlobalElementIdsforAllElements);

    std::map<unsigned,std::vector<int>> build_node_sharing_procs(const std::vector<std::array<unsigned, NPE>> &elementConn,
        const std::vector<int> &elementProcOwners) const;
    std::map<unsigned,std::vector<int>> build_node_sharing_procs_for_all_nodes_on_all_procs(const unsigned numNodes, const unsigned numProcs) const;
};

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHBUILDER_HPP_ */
