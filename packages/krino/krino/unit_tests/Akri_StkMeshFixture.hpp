#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_

#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <array>
#include <vector>
#include <Akri_StkMeshBuilder.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

template<stk::topology::topology_t TOPO>
class StkMeshFixture : public ::testing::Test
{
protected:
    static constexpr unsigned NUM_DIM = TopologyData<TOPO>::spatial_dimension();
    static constexpr unsigned NPE = TopologyData<TOPO>::num_nodes();
    static constexpr unsigned theBlockId = 1;
    const stk::ParallelMachine mComm = MPI_COMM_WORLD;
    const int mProc{stk::parallel_machine_rank(mComm)};
    std::vector< std::string > entityRankNames{ std::string("NODE"), std::string("EDGE"), std::string("FACE"), std::string("ELEMENT"), std::string("FAMILY_TREE") };
    std::unique_ptr<stk::mesh::BulkData> mMeshPtr{stk::mesh::MeshBuilder(mComm).set_spatial_dimension(NUM_DIM).set_entity_rank_names(entityRankNames).create()};
    stk::mesh::BulkData & mMesh{*mMeshPtr};
    StkMeshBuilder<TOPO> mBuilder{mMesh, mComm};

    AuxMetaData & get_aux_meta() { return mBuilder.get_aux_meta(); }
    const AuxMetaData & get_aux_meta() const { return mBuilder.get_aux_meta(); }
    const std::vector<stk::mesh::EntityId> & get_assigned_node_global_ids() const { return mBuilder.get_assigned_node_global_ids(); }
    stk::mesh::Entity get_assigned_node_for_index(const size_t nodeIndex) const { return mMesh.get_entity(stk::topology::NODE_RANK, get_assigned_node_global_ids()[nodeIndex]); }
    const std::vector<stk::mesh::EntityId> & get_assigned_element_global_ids() const { return mBuilder.get_assigned_element_global_ids(); }
    std::vector<stk::mesh::EntityId> get_ids_of_elements_with_given_indices(const std::vector<unsigned> & elemIndices) const { return mBuilder.get_ids_of_elements_with_given_indices(elemIndices); }
    const std::vector<stk::mesh::Entity> & get_owned_elements() const { return mBuilder.get_owned_elements(); }
    stk::math::Vector3d get_node_coordinates(const stk::mesh::Entity node) const { return mBuilder.get_node_coordinates(node); }

    template <typename MeshSpecType>
    void build_full_np1_mesh(const MeshSpecType &meshSpec)
    {
        build_mesh(meshSpec.mNodeLocs, {meshSpec.mAllTetConn});
    }

    void build_mesh(const std::vector<stk::math::Vec<double,NUM_DIM>> &nodeLocs,
                    const std::vector<std::vector<std::array<unsigned,NPE>>> &elemConnPerProc)
    {
      mMesh.mesh_meta_data().use_simple_fields();
      mBuilder.build_mesh(nodeLocs, elemConnPerProc, theBlockId);
    }

    void build_mesh(const std::vector<stk::math::Vec<double,NUM_DIM>> &nodeLocs,
                    const std::vector<std::array<unsigned, NPE>> &elementConn,
                    const std::vector<unsigned> &elementBlockIDs,
                    const std::vector<int> &specifiedElementProcOwners = {})
    {
      mMesh.mesh_meta_data().use_simple_fields();
      mBuilder.create_block_parts(elementBlockIDs);
      mBuilder.build_mesh(nodeLocs, elementConn, elementBlockIDs, specifiedElementProcOwners);
    }
};

typedef StkMeshFixture<stk::topology::TETRAHEDRON_4> StkMeshTetFixture;
typedef StkMeshFixture<stk::topology::TRIANGLE_3_2D> StkMeshTriFixture;

}


#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_ */
