#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_

#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <array>
#include <vector>
#include <Akri_StkMeshBuilder.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

template<stk::topology::topology_t TOPO>
class StkMeshAndBuilder
{
protected:
    static constexpr unsigned NUM_DIM = TopologyData<TOPO>::spatial_dimension();
    static constexpr unsigned NPE = TopologyData<TOPO>::num_nodes();
    static constexpr unsigned theBlockId = 1;
public:
    stk::mesh::BulkData & get_mesh() { return mMesh; }
    const stk::mesh::BulkData & get_mesh() const { return mMesh; }
    stk::ParallelMachine get_comm() const { return mComm; }
    int parallel_size() const { return stk::parallel_machine_size(mComm); }
    int parallel_rank() const { return stk::parallel_machine_rank(mComm); }
    AuxMetaData & get_aux_meta() { return mBuilder.get_aux_meta(); }
    const AuxMetaData & get_aux_meta() const { return mBuilder.get_aux_meta(); }
    const std::vector<stk::mesh::EntityId> & get_assigned_node_global_ids() const { return mBuilder.get_assigned_node_global_ids(); }
    stk::mesh::Entity get_assigned_node_for_index(const size_t nodeIndex) const { return mBuilder.get_assigned_node_for_index(nodeIndex); }
    const std::vector<stk::mesh::EntityId> & get_assigned_element_global_ids() const { return mBuilder.get_assigned_element_global_ids(); }
    std::vector<stk::mesh::EntityId> get_ids_of_elements_with_given_indices(const std::vector<unsigned> & elemIndices) const { return mBuilder.get_ids_of_elements_with_given_indices(elemIndices); }
    std::vector<stk::mesh::EntityId> get_ids_of_nodes_with_given_indices(const std::vector<unsigned> & nodeIndices) const { return mBuilder.get_ids_of_nodes_with_given_indices(nodeIndices); }
    const std::vector<stk::mesh::Entity> & get_owned_elements() const { return mBuilder.get_owned_elements(); }
    stk::math::Vector3d get_node_coordinates(const stk::mesh::Entity node) const { return mBuilder.get_node_coordinates(node); }
    const stk::mesh::FieldBase & get_coordinates_field() const { return mBuilder.get_coordinates_field(); }

    template <typename MeshSpecType>
    void build_full_np1_mesh(const MeshSpecType &meshSpec)
    {
        build_mesh(meshSpec.mNodeLocs, {meshSpec.mAllTetConn});
    }

    void build_mesh(const std::vector<stk::math::Vec<double,NUM_DIM>> &nodeLocs,
                    const std::vector<std::vector<std::array<unsigned,NPE>>> &elemConnPerProc)
    {
      mBuilder.build_mesh(nodeLocs, elemConnPerProc, theBlockId);
    }

    void build_mesh(const std::vector<stk::math::Vec<double,NUM_DIM>> &nodeLocs,
                    const std::vector<std::array<unsigned, NPE>> &elementConn,
                    const std::vector<unsigned> &elementBlockIDs,
                    const std::vector<int> &specifiedElementProcOwners = {})
    {
      mBuilder.build_mesh(nodeLocs, elementConn, elementBlockIDs, specifiedElementProcOwners);
    }

    void write_mesh(const std::string &fileName)
    {
      mBuilder.write_mesh(fileName);
    }

protected:
    const stk::ParallelMachine mComm = MPI_COMM_WORLD;
    const int mProc{stk::parallel_machine_rank(mComm)};
    std::unique_ptr<stk::mesh::BulkData> mMeshPtr{stk::mesh::MeshBuilder(mComm).set_spatial_dimension(NUM_DIM).create()};
    stk::mesh::BulkData & mMesh{*mMeshPtr};
    StkMeshBuilder<TOPO> mBuilder{mMesh, mComm};
};

template<stk::topology::topology_t TOPO>
class StkMeshFixture : public ::testing::Test, public StkMeshAndBuilder<TOPO>
{
public:
  void set_valid_proc_sizes_for_test(const std::vector<int> & procSizes) { mTestProcSizes = procSizes; }
  bool is_valid_proc_size_for_test() const { STK_ThrowRequireMsg(!mTestProcSizes.empty(), "Valid proc sizes not set for test."); return std::find(mTestProcSizes.begin(), mTestProcSizes.end(), this->parallel_size()) != mTestProcSizes.end(); }
protected:
    std::vector<int> mTestProcSizes;
};

typedef StkMeshFixture<stk::topology::BEAM_2> StkMeshBeamFixture;
typedef StkMeshFixture<stk::topology::TETRAHEDRON_4> StkMeshTetFixture;
typedef StkMeshFixture<stk::topology::TRIANGLE_3_2D> StkMeshTriFixture;
typedef StkMeshFixture<stk::topology::QUADRILATERAL_4_2D> StkMeshQuadFixture;
typedef StkMeshFixture<stk::topology::HEXAHEDRON_8> StkMeshHexFixture;
typedef StkMeshAndBuilder<stk::topology::TETRAHEDRON_4> StkTetMeshAndBuilder;
typedef StkMeshAndBuilder<stk::topology::TRIANGLE_3_2D> StkTriMeshAndBuilder;

}


#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_ */
