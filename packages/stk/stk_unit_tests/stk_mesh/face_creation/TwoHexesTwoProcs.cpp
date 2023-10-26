#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

#include "FaceCreatorFixture.hpp"

namespace
{

class FaceCreatorUsingBulkDataFaceSharingTester : public FaceCreatorFixture
{
protected:

  ~FaceCreatorUsingBulkDataFaceSharingTester()
  {
    bulkData.reset();
    metaData.reset();
  }

  virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                             unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                             unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity())
  {
    STK_ThrowRequireMsg((initialBucketCapacity == stk::mesh::get_default_initial_bucket_capacity()) &&
                        (maximumBucketCapacity == stk::mesh::get_default_maximum_bucket_capacity()),
                        "allocate_bulk: BulkDataFaceSharingTester doesn't recognize non-default bucket capacity");

    metaData = stk::mesh::MeshBuilder(communicator)
                .set_spatial_dimension(m_spatialDim)
                .set_entity_rank_names(m_entityRankNames)
                .set_aura_option(auraOption).create_meta_data();

    set_bulk(std::make_shared<stk::unit_test_util::BulkDataFaceSharingTester>(get_meta(), get_comm(), auraOption));
  }

  void test_elem_elem_graph_for_face_connection_info()
  {
    stk::mesh::ElemElemGraph egraph(get_bulk());
    test_sides_using_graph(egraph);
  }

private:

  stk::mesh::EntityVector get_shared_sides()
  {
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(get_meta().globally_shared_part(), get_bulk().buckets(get_meta().side_rank()), sides);
    EXPECT_EQ(1u, sides.size());
    return sides;
  }

  void test_sides_using_graph(const stk::mesh::ElemElemGraph & egraph)
  {
    stk::mesh::EntityVector sides = get_shared_sides();
    for(size_t i = 0;i < sides.size();++i)
      test_graph_info_on_side(egraph, sides[i]);
  }

  void test_graph_info_on_side(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity side)
  {
    unsigned num_elements = get_bulk().num_elements(side);
    const stk::mesh::Entity* elements = get_bulk().begin_elements(side);
    for(unsigned j=0;j<num_elements;++j)
      check_if_remote_sides_have_consistent_parallel_info(egraph, elements[j], side);
  }

  void check_if_remote_sides_have_consistent_parallel_info(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity element, stk::mesh::Entity side)
  {
    if(get_bulk().bucket(element).owned())
    {
      int num_connected_elems = egraph.get_num_connected_elems(element);
      for(int k=0;k<num_connected_elems;++k)
        check_if_remote_side_has_consistent_parallel_info(egraph, element, k, side);
    }
  }

  void check_if_remote_side_has_consistent_parallel_info(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity element, int element_index, stk::mesh::Entity side)
  {
    if(!egraph.is_connected_elem_locally_owned(element, element_index))
    {
      stk::mesh::impl::IdViaSidePair sidePair = egraph.get_connected_remote_id_and_via_side(element, element_index);
      int remoteSideOrdinal = egraph.get_connected_elements_side(element, element_index);
      const stk::mesh::impl::ParallelInfo &p_info = egraph.get_const_parallel_edge_info(element, sidePair.side, sidePair.id, remoteSideOrdinal);
      test_parallel_info_for_side(p_info, element, side, sidePair.side);
    }
  }

  void test_parallel_info_for_side(const stk::mesh::impl::ParallelInfo &p_info, stk::mesh::Entity element, stk::mesh::Entity side, int side_ordinal)
  {
    test_parallel_info(p_info);
    stk::mesh::EntityVector side_nodes(get_bulk().begin_nodes(side),get_bulk().end_nodes(side));
    stk::mesh::EntityVector permuted_element_side_nodes = get_permuted_element_side_nodes(element, side_ordinal, p_info.get_proc_rank_of_neighbor(), p_info.m_permutation);
    EXPECT_TRUE(side_nodes == permuted_element_side_nodes);
  }

  stk::mesh::EntityVector get_permuted_element_side_nodes(stk::mesh::Entity element, int side_ordinal, int other_proc_id, unsigned permutation_other_proc)
  {
    stk::mesh::EntityVector element_nodes(get_bulk().begin_nodes(element),get_bulk().end_nodes(element));
    stk::topology side_topology = get_side_topology(element, side_ordinal);
    stk::mesh::EntityVector element_side_nodes = get_element_side_nodes(element, element_nodes, side_ordinal);
    return get_permuted_nodes(side_topology, element_side_nodes, other_proc_id, permutation_other_proc);
  }

  stk::mesh::EntityVector get_element_side_nodes(stk::mesh::Entity element, const stk::mesh::EntityVector& element_nodes, int side_ordinal)
  {
    stk::topology element_topology = get_bulk().bucket(element).topology();
    stk::mesh::EntityVector element_side_nodes(get_side_topology(element, side_ordinal).num_nodes());
    element_topology.side_nodes(element_nodes.data(), side_ordinal, element_side_nodes.data());
    return element_side_nodes;
  }

  stk::topology get_side_topology(stk::mesh::Entity element, int side_ordinal)
  {
    stk::topology element_topology = get_bulk().bucket(element).topology();
    return element_topology.side_topology(side_ordinal);
  }

  stk::mesh::EntityVector get_permuted_nodes(stk::topology side_topology, const stk::mesh::EntityVector &element_side_nodes, int other_proc_id, unsigned permutation_other_proc)
  {
    stk::mesh::EntityVector permuted_element_side_nodes(side_topology.num_nodes());
    unsigned permutation = get_bulk().parallel_rank()>other_proc_id ? permutation_other_proc : 0;
    side_topology.permutation_nodes(element_side_nodes.data(), permutation, permuted_element_side_nodes.data());
    return permuted_element_side_nodes;
  }

  void test_parallel_info(const stk::mesh::impl::ParallelInfo &p_info)
  {
    EXPECT_EQ(1-get_bulk().parallel_rank(), p_info.get_proc_rank_of_neighbor());
    EXPECT_EQ(4, p_info.m_permutation);
    EXPECT_EQ(stk::topology::HEXAHEDRON_8, p_info.m_remote_element_topology);
  }
};

TEST_F(FaceCreatorUsingBulkDataFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }
}

//  These two tests demonstrate that if faces are created with a permutation other than 0 then
//  the graph parallel info is incorrect because it assumes all faces are created with permutation
//  0 on the originating element/side.  If we can guarantee all faces are created through
//  declare_element_side, then this issue goes away.
TEST_F(FaceCreatorUsingBulkDataFaceSharingTester, testFaceDataUsingElemElemGraphWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    test_elem_elem_graph_for_face_connection_info();
  }
}

TEST_F(FaceCreatorUsingBulkDataFaceSharingTester, testFaceDataUsingElemElemGraphWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    test_elem_elem_graph_for_face_connection_info();
  }
}

}
