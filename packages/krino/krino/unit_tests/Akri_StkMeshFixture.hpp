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

template<int DIM>
class StkMeshFixture : public ::testing::Test
{
protected:
    static constexpr int NPE = DIM+1;
    static constexpr unsigned theBlockId = 1;
    const stk::ParallelMachine mComm = MPI_COMM_WORLD;
    const int mProc{stk::parallel_machine_rank(mComm)};
    std::unique_ptr<stk::mesh::BulkData> mMeshPtr{stk::mesh::MeshBuilder(mComm).set_spatial_dimension(DIM).create()};
    stk::mesh::BulkData & mMesh{*mMeshPtr};
    StkMeshBuilder<DIM> mBuilder{mMesh, mComm};

    const std::vector<stk::mesh::EntityId> & get_assigned_node_global_ids() const { return mBuilder.get_assigned_node_global_ids(); }
    stk::mesh::Entity get_assigned_node_for_index(const size_t nodeIndex) const { return mMesh.get_entity(stk::topology::NODE_RANK, get_assigned_node_global_ids()[nodeIndex]); }
    const std::vector<stk::mesh::Entity> & get_owned_elements() const { return mBuilder.get_owned_elements(); }

    template <typename MeshSpecType>
    void build_full_np1_mesh(const MeshSpecType &meshSpec)
    {
        build_mesh(meshSpec.mNodeLocs, {meshSpec.mAllTetConn});
    }

    void build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
                    const std::vector<std::vector<std::array<unsigned,NPE>>> &elemConnPerProc)
    {
      mMesh.mesh_meta_data().use_simple_fields();
      mBuilder.build_mesh(nodeLocs, elemConnPerProc, theBlockId);
    }

    void build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
                    const std::vector<std::array<unsigned, NPE>> &elementConn,
                    const std::vector<unsigned> &elementBlockIDs,
                    const std::vector<int> &specifiedElementProcOwners = {})
    {
      mMesh.mesh_meta_data().use_simple_fields();
      mBuilder.create_block_parts(elementBlockIDs);
      mBuilder.build_mesh(nodeLocs, elementConn, elementBlockIDs, specifiedElementProcOwners);
    }
};

typedef StkMeshFixture<3> StkMeshTetFixture;
typedef StkMeshFixture<2> StkMeshTriFixture;

}


#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_ */
