#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_

#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <algorithm>
#include <array>
#include <map>
#include <vector>
#include <Akri_StkMeshBuilder.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino
{

template<stk::topology::topology_t TOPO>
class StkMeshAndBuilder
{
protected:
    static constexpr unsigned NUM_DIM = TopologyData<TOPO>::spatial_dimension();
    static constexpr unsigned NPE = TopologyData<TOPO>::num_nodes();
    static constexpr unsigned theBlockId = 1;

    static stk::topology element_topology() { return TOPO; }
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
    stk::mesh::Entity get_assigned_element_for_index(const size_t elemIndex) const { return mBuilder.get_assigned_element_for_index(elemIndex); }
    // Whether the element or node with the given index in the mesh spec is present on this processor (owned, shared
    // or ghost), or owned by it.
    bool is_element_present_on_proc(const size_t elemIndex) const { return mMesh.is_valid(get_assigned_element_for_index(elemIndex)); }
    bool is_owned_element(const size_t elemIndex) const { return is_element_present_on_proc(elemIndex) && mMesh.bucket(get_assigned_element_for_index(elemIndex)).owned(); }
    bool is_node_present_on_proc(const size_t nodeIndex) const { return mMesh.is_valid(get_assigned_node_for_index(nodeIndex)); }
    bool is_owned_node(const size_t nodeIndex) const { return is_node_present_on_proc(nodeIndex) && mMesh.bucket(get_assigned_node_for_index(nodeIndex)).owned(); }
    const std::vector<stk::mesh::EntityId> & get_assigned_element_global_ids() const { return mBuilder.get_assigned_element_global_ids(); }
    std::vector<stk::mesh::EntityId> get_ids_of_elements_with_given_indices(const std::vector<unsigned> & elemIndices) const { return mBuilder.get_ids_of_elements_with_given_indices(elemIndices); }
    std::vector<stk::mesh::EntityId> get_ids_of_nodes_with_given_indices(const std::vector<unsigned> & nodeIndices) const { return mBuilder.get_ids_of_nodes_with_given_indices(nodeIndices); }
    const std::vector<stk::mesh::Entity> & get_owned_elements() const { return mBuilder.get_owned_elements(); }
    stk::math::Vector3d get_node_coordinates(const stk::mesh::Entity node) const { return mBuilder.get_node_coordinates(node); }
    void set_node_coordinates(const stk::mesh::Entity node, const stk::math::Vector3d &newLoc) { mBuilder.set_node_coordinates(node, newLoc); }
    const stk::mesh::FieldBase & get_coordinates_field() const { return mBuilder.get_coordinates_field(); }

    // Element owners by processor count; counts without an entry use the default decomposition.
    using ElementOwnersByProcCount = std::map<int, std::vector<int>>;

    // Builds the mesh of meshSpec (nodeLocs, allElementConn) in block 1, decomposed over the processors.  By default
    // element i is owned by processor i % np; give owners for specific processor counts in elementOwnersByProcCount
    // to choose the decomposition, e.g. {{2, {0,0,1,1}}, {4, {0,1,2,3}}}.  Without owners for np, nothing is built
    // when np exceeds the number of elements.
    template <typename MeshSpecType>
    void build_mesh_from_spec(const MeshSpecType &meshSpec, const ElementOwnersByProcCount & elementOwnersByProcCount = {})
    {
      build_mesh_from_spec_in_blocks(meshSpec, std::vector<unsigned>(meshSpec.allElementConn.size(), theBlockId), elementOwnersByProcCount);
    }

    // As build_mesh_from_spec, with element i in block elementBlockIDs[i].
    template <typename MeshSpecType>
    void build_mesh_from_spec_in_blocks(const MeshSpecType &meshSpec, const std::vector<unsigned> & elementBlockIDs,
        const ElementOwnersByProcCount & elementOwnersByProcCount = {})
    {
      const size_t numElements = meshSpec.allElementConn.size();
      STK_ThrowRequireMsg(elementBlockIDs.size() == numElements,
          "Expected " << numElements << " element block ids, got " << elementBlockIDs.size() << ".");
      mMaxProcSizeForDefaultDecomposition = numElements;
      std::vector<int> elementOwners;
      const auto iter = elementOwnersByProcCount.find(parallel_size());
      if (iter != elementOwnersByProcCount.end())
      {
        elementOwners = iter->second;
        STK_ThrowRequireMsg(elementOwners.size() == numElements,
            "Expected " << numElements << " element owners for " << parallel_size() << " processors, got " << elementOwners.size() << ".");
        for (int owner : elementOwners)
          STK_ThrowRequireMsg(owner >= 0 && owner < parallel_size(),
              "Element owner " << owner << " is not a valid processor for " << parallel_size() << " processors.");
      }
      else if ((size_t)parallel_size() <= numElements)
      {
        elementOwners = mBuilder.get_processor_distribution_for_num_elements(numElements);
      }
      else
      {
        return;
      }
      build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, elementOwners);
    }

    // Number of elements of the mesh built from a spec (0 if none): the most processors the default can use.
    size_t get_max_proc_size_for_default_decomposition() const { return mMaxProcSizeForDefaultDecomposition; }

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
    size_t mMaxProcSizeForDefaultDecomposition{0};
};

template<stk::topology::topology_t TOPO>
class StkMeshFixture : public ::testing::Test, public StkMeshAndBuilder<TOPO>
{
public:
  void set_valid_proc_sizes_for_test(const std::vector<int> & procSizes) { mTestProcSizes = procSizes; }

  // As in StkMeshAndBuilder, but these build nothing on processor counts excluded by set_valid_proc_sizes_for_test.
  template <typename MeshSpecType>
  void build_mesh_from_spec(const MeshSpecType &meshSpec, const typename StkMeshAndBuilder<TOPO>::ElementOwnersByProcCount & elementOwnersByProcCount = {})
  {
    if (is_excluded_proc_size()) return;
    StkMeshAndBuilder<TOPO>::build_mesh_from_spec(meshSpec, elementOwnersByProcCount);
  }

  template <typename MeshSpecType>
  void build_mesh_from_spec_in_blocks(const MeshSpecType &meshSpec, const std::vector<unsigned> & elementBlockIDs,
      const typename StkMeshAndBuilder<TOPO>::ElementOwnersByProcCount & elementOwnersByProcCount = {})
  {
    if (is_excluded_proc_size()) return;
    StkMeshAndBuilder<TOPO>::build_mesh_from_spec_in_blocks(meshSpec, elementBlockIDs, elementOwnersByProcCount);
  }

  // Whether the test runs on this number of processors: the sizes given to set_valid_proc_sizes_for_test if called,
  // otherwise 1 up to the number of elements of the mesh built by build_mesh_from_spec.
  bool is_valid_proc_size_for_test() const
  {
    if (mTestProcSizes.empty())
    {
      STK_ThrowRequireMsg(this->get_max_proc_size_for_default_decomposition() > 0, "Valid proc sizes not set for test.");
      return (size_t)this->parallel_size() <= this->get_max_proc_size_for_default_decomposition();
    }
    return std::find(mTestProcSizes.begin(), mTestProcSizes.end(), this->parallel_size()) != mTestProcSizes.end();
  }
protected:
    bool is_excluded_proc_size() const { return !mTestProcSizes.empty() && std::find(mTestProcSizes.begin(), mTestProcSizes.end(), this->parallel_size()) == mTestProcSizes.end(); }
    std::vector<int> mTestProcSizes;
};

typedef StkMeshFixture<stk::topology::BEAM_2> StkMeshBeamFixture;
typedef StkMeshFixture<stk::topology::TETRAHEDRON_4> StkMeshTetFixture;
typedef StkMeshFixture<stk::topology::TETRAHEDRON_10> StkMeshTet10Fixture;
typedef StkMeshFixture<stk::topology::TRIANGLE_3_2D> StkMeshTriFixture;
typedef StkMeshFixture<stk::topology::TRIANGLE_6_2D> StkMeshTri6Fixture;
typedef StkMeshFixture<stk::topology::QUADRILATERAL_4_2D> StkMeshQuadFixture;
typedef StkMeshFixture<stk::topology::HEXAHEDRON_8> StkMeshHexFixture;
typedef StkMeshAndBuilder<stk::topology::TETRAHEDRON_4> StkTetMeshAndBuilder;
typedef StkMeshAndBuilder<stk::topology::TRIANGLE_3_2D> StkTriMeshAndBuilder;

}


#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_STKMESHFIXTURE_HPP_ */
