#include <set>
#include <vector>

#include <Akri_AllReduce.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture_Edges.hpp>
#include <Akri_Unit_RefinementFixture_Tet.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

// Refinement of explicit edge entities (not sides) of tet4 and tet10 meshes (LINE_2 and LINE_3
// edges): each element edge of an active element is attached to an edge with the same nodes, edges
// split by refinement have child edges with the parent edge's parts, and unrefinement removes the
// child edges.  Hex8 meshes with edges are tested in Akri_Unit_Refine_Hex.cpp.

namespace krino
{

using RightTet10SurroundedByFaceTets =
    HigherOrderSimplexMeshSpec<RightTetSurroundedByFaceTets, stk::topology::TETRAHEDRON_10>;

template <typename MESHSPEC>
class RefineEdgesOfTetsFixture : public WithExplicitEdges<RefinementFixture<MESHSPEC>>
{
public:
  using RefinementFixture<MESHSPEC>::mMesh;
  using RefinementFixture<MESHSPEC>::myRefinement;

  RefineEdgesOfTetsFixture()
  {
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    const stk::topology edgeTopology = stk::topology(MESHSPEC::TOPOLOGY).edge_topology(0);
    myEdgePart = &meta.declare_part_with_topology("edge_part", edgeTopology);
    myEdgeSetPart = &meta.declare_part("edge_set", stk::topology::EDGE_RANK);
    this->set_valid_proc_sizes_for_test({1, 2, 3, 4});
    this->build_mesh_from_spec_in_blocks(this->meshSpec, {1, 2, 2, 2, 2});
    if (!this->is_valid_proc_size_for_test()) return;
    this->create_edges(myEdgePart);
    put_edges_of_central_element_in_edge_set();
    myAllOriginalEdgeSetSegments =
        this->gather_segments(this->get_owned_edge_segments(*myEdgeSetPart));
  }

protected:
  // The edges of element 0, the tet surrounded by the others, also go into edge_set, so that a part
  // carried only by some edges is checked to be inherited by exactly their children.
  void put_edges_of_central_element_in_edge_set()
  {
    std::vector<stk::mesh::Entity> edges;
    if (this->is_owned_element(0))
    {
      const stk::mesh::Entity elem = this->get_assigned_element_for_index(0);
      edges.assign(mMesh.begin_edges(elem), mMesh.end_edges(elem));
    }
    std::vector<stk::mesh::Entity> ownedEdges;
    for (auto edge : edges)
      if (mMesh.bucket(edge).owned()) ownedEdges.push_back(edge);
    mMesh.batch_change_entity_parts(ownedEdges,
        std::vector<stk::mesh::PartVector>(ownedEdges.size(), {myEdgeSetPart}),
        std::vector<stk::mesh::PartVector>(ownedEdges.size()));
  }

  // Whether the edge lies along one of the original edge_set edges (which are all edges of element
  // 0, known on its owner; other processors gather them).
  bool is_along_original_edge_set_edge(const stk::mesh::Entity edge) const
  {
    const stk::mesh::Entity * nodes = mMesh.begin_nodes(edge);
    return this->is_along_segment(nodes[0], nodes[1], myAllOriginalEdgeSetSegments);
  }

  void check_edge_set_edges_are_along_original_edge_set_edges()
  {
    size_t numEdgeSetEdgesElsewhere = 0, numOtherEdgesAlongEdgeSet = 0;
    for (auto * bucketPtr :
        mMesh.get_buckets(stk::topology::EDGE_RANK, mMesh.mesh_meta_data().locally_owned_part()))
      for (auto edge : *bucketPtr)
      {
        const bool inEdgeSet = bucketPtr->member(*myEdgeSetPart);
        const bool along = is_along_original_edge_set_edge(edge);
        if (inEdgeSet && !along) ++numEdgeSetEdgesElsewhere;
        if (!inEdgeSet && along) ++numOtherEdgesAlongEdgeSet;
      }
    all_reduce_sum(mMesh.parallel(), numEdgeSetEdgesElsewhere, numOtherEdgesAlongEdgeSet);
    EXPECT_EQ(0u, numEdgeSetEdgesElsewhere);
    EXPECT_EQ(0u, numOtherEdgesAlongEdgeSet);
  }

  stk::mesh::Part * myEdgePart{nullptr};
  stk::mesh::Part * myEdgeSetPart{nullptr};
  std::vector<std::array<stk::math::Vector3d, 2>> myAllOriginalEdgeSetSegments;
};

using TetMeshSpecs = ::testing::Types<RightTetSurroundedByFaceTets, RightTet10SurroundedByFaceTets>;
TYPED_TEST_SUITE(RefineEdgesOfTetsFixture, TetMeshSpecs);

TYPED_TEST(RefineEdgesOfTetsFixture,
    everyRefinementOfTheFaceNeighbors_edgesAttachedToAllActiveElements_andUnrefineRestoresEdges)
{
  if (!this->is_valid_proc_size_for_test()) return;

  const size_t numOriginalEdges = this->num_edges_in(*this->myEdgePart);
  const size_t numOriginalEdgeSetEdges = this->num_edges_in(*this->myEdgeSetPart);
  EXPECT_EQ(6u, numOriginalEdgeSetEdges);

  const unsigned numCases = 1 << 4;
  for (unsigned i = 0; i < numCases; ++i)
  {
    std::vector<unsigned> faceElementsToRefine;
    for (unsigned iFace = 0; iFace < 4; ++iFace)
      if (i & (1 << iFace)) faceElementsToRefine.push_back(1 + iFace);

    this->refine_elements_with_given_indices(faceElementsToRefine);
    this->check_edges_of_active_elements();
    this->check_edge_set_edges_are_along_original_edge_set_edges();

    this->unrefine_mesh();
    EXPECT_EQ(numOriginalEdges, this->num_edges_in(*this->myEdgePart)) << "case " << i;
    EXPECT_EQ(numOriginalEdgeSetEdges, this->num_edges_in(*this->myEdgeSetPart)) << "case " << i;
    EXPECT_EQ(0u, this->get_global_num_refined_edge_nodes()) << "case " << i;
  }
}

TYPED_TEST(
    RefineEdgesOfTetsFixture, uniformRefinementTwice_edgesAttachedAndEdgeSetInheritedByChildren)
{
  if (!this->is_valid_proc_size_for_test()) return;

  const size_t numOriginalEdgeSetEdges = this->num_edges_in(*this->myEdgeSetPart);
  this->myRefinement.do_uniform_refinement(1);
  this->check_edges_of_active_elements();
  // Each edge_set edge is split in two.
  EXPECT_EQ(numOriginalEdgeSetEdges + 2 * numOriginalEdgeSetEdges,
      this->num_edges_in(*this->myEdgeSetPart));

  this->myRefinement.do_uniform_refinement(1);
  this->check_edges_of_active_elements();
  EXPECT_EQ(numOriginalEdgeSetEdges * (1 + 2 + 4), this->num_edges_in(*this->myEdgeSetPart));

  this->check_edge_set_edges_are_along_original_edge_set_edges();
}

TYPED_TEST(RefineEdgesOfTetsFixture, refineThenUnrefineWithMarkers_childEdgesRemoved)
{
  if (!this->is_valid_proc_size_for_test()) return;

  const size_t numOriginalEdges = this->num_edges_in(*this->myEdgePart);
  this->refine_elements_with_given_indices({0});
  this->check_edges_of_active_elements();
  EXPECT_LT(numOriginalEdges, this->num_edges_in(*this->myEdgePart));

  this->clear_refinement_marker();
  std::vector<stk::mesh::Entity> leaves;
  for (auto * bucketPtr : this->mMesh.get_buckets(stk::topology::ELEMENT_RANK,
           this->mMesh.mesh_meta_data().locally_owned_part() & this->myRefinement.child_part()))
    for (auto elem : *bucketPtr)
      leaves.push_back(elem);
  this->mark_elements_for_unrefinement(leaves);
  this->do_refinement();
  this->check_edges_of_active_elements();
  EXPECT_EQ(numOriginalEdges, this->num_edges_in(*this->myEdgePart));
  EXPECT_EQ(0u, this->get_global_num_refined_edge_nodes());
}

TYPED_TEST(RefineEdgesOfTetsFixture,
    refineChildOfTransitionElementThenUnrefineWithMarkers_edgesAttachedAndOriginalEdgesRestored)
{
  if (!this->is_valid_proc_size_for_test()) return;

  const size_t numOriginalElements =
      get_global_num_entities(this->mMesh, stk::topology::ELEMENT_RANK);
  const size_t numOriginalEdges = this->num_edges_in(*this->myEdgePart);
  this->refine_elements_with_given_indices({0});
  EXPECT_EQ(29u, get_global_num_entities(this->mMesh, stk::topology::ELEMENT_RANK));

  // Refining a transition element fully refines its parent instead.
  std::vector<stk::mesh::Entity> elementsToRefine;
  const stk::mesh::Entity transitionParent = this->get_assigned_element_for_index(1);
  if (this->mMesh.is_valid(transitionParent) && this->mMesh.bucket(transitionParent).owned())
  {
    const std::vector<stk::mesh::Entity> transitionElements = this->get_children(transitionParent);
    ASSERT_EQ(4u, transitionElements.size());
    elementsToRefine.push_back(transitionElements[0]);
  }
  this->refine_elements(elementsToRefine);
  EXPECT_FALSE(this->myRefinement.have_any_hanging_refined_nodes());
  EXPECT_EQ(33u, get_global_num_entities(this->mMesh, stk::topology::ELEMENT_RANK));
  this->check_edges_of_active_elements();

  this->clear_refinement_marker();
  std::vector<stk::mesh::Entity> leaves;
  for (auto * bucketPtr : this->mMesh.get_buckets(stk::topology::ELEMENT_RANK,
           this->mMesh.mesh_meta_data().locally_owned_part() & this->myRefinement.child_part()))
    for (auto elem : *bucketPtr)
      leaves.push_back(elem);
  this->mark_elements_for_unrefinement(leaves);
  this->do_refinement();
  this->check_edges_of_active_elements();
  EXPECT_EQ(numOriginalElements, get_global_num_entities(this->mMesh, stk::topology::ELEMENT_RANK));
  EXPECT_EQ(numOriginalEdges, this->num_edges_in(*this->myEdgePart));
  EXPECT_EQ(0u, this->get_global_num_refined_edge_nodes());
}

TYPED_TEST(RefineEdgesOfTetsFixture,
    unrefineFullyRefinedElementNextToRefinedNeighbor_becomesTransitionAndEdgesAttached)
{
  if (!this->is_valid_proc_size_for_test()) return;

  this->refine_elements_with_given_indices({0, 1});

  // The refined nodes on the face shared with the refined neighbor are kept, so the central element
  // goes from fully refined to a transition with that face refined.
  this->clear_refinement_marker();
  const stk::mesh::Entity centralElem = this->get_assigned_element_for_index(0);
  const bool ownCentralElem =
      this->mMesh.is_valid(centralElem) && this->mMesh.bucket(centralElem).owned();
  if (ownCentralElem)
  {
    ASSERT_EQ(8u, this->get_num_children(centralElem));
    this->mark_elements_for_unrefinement(this->get_children(centralElem));
  }
  this->do_refinement();
  if (ownCentralElem)
  {
    EXPECT_EQ(4u, this->get_num_children(centralElem));
  }
  this->check_edge_set_edges_are_along_original_edge_set_edges();
}

// The same checks of the edges after every adaptation, through the existing refinement tests of
// other meshes.

template <typename FIXTURE> class WithEdgesCreated : public WithExplicitEdges<FIXTURE>
{
public:
  WithEdgesCreated()
  {
    if (this->is_valid_proc_size_for_test()) this->create_edges();
  }
};

using RightTetSurroundedByEdgeTetsRefinementWithEdges =
    WithEdgesCreated<RightTetSurroundedByEdgeTetsRefinement>;
using RightTet10SurroundedByEdgeTetsRefinementWithEdges =
    WithEdgesCreated<RightTet10SurroundedByEdgeTetsRefinement>;
using UMRRegularTetRefinementWithEdges = WithEdgesCreated<UMRRegularTetRefinement>;
using UMRRegularTet10RefinementWithEdges = WithEdgesCreated<UMRRegularTet10Refinement>;

TEST_F(RightTetSurroundedByEdgeTetsRefinementWithEdges,
    markedAnyTransitionElementForEveryEdgeConfiguration_edgesAttachedToAllActiveElements)
{
  if (!is_valid_proc_size_for_test()) return;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(0);
}

TEST_F(RightTet10SurroundedByEdgeTetsRefinementWithEdges,
    markedAnyTransitionElementForEveryEdgeConfiguration_edgesAttachedToAllActiveElements)
{
  if (!is_valid_proc_size_for_test()) return;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(0);
}

TEST_F(
    UMRRegularTetRefinementWithEdges, refinementThenUnrefinement_edgesAttachedToAllActiveElements)
{
  if (!is_valid_proc_size_for_test()) return;
  test_refinement_of_elements_spanning_x_equal_0_then_unrefinement({72, 512, 2892}, {546, 72, 8});
  EXPECT_EQ(num_original_edges(), num_edges_in(mMesh.mesh_meta_data().universal_part()));
}

TEST_F(
    UMRRegularTet10RefinementWithEdges, refinementThenUnrefinement_edgesAttachedToAllActiveElements)
{
  if (!is_valid_proc_size_for_test()) return;
  test_refinement_of_elements_spanning_x_equal_0_then_unrefinement({72, 512, 2892}, {546, 72, 8});
  EXPECT_EQ(num_original_edges(), num_edges_in(mMesh.mesh_meta_data().universal_part()));
}

TEST_F(UMRRegularTetRefinementWithEdges,
    parentAndChildrenMovedToNewProc_unrefinementReturnsParentToOriginatingProcWithItsEdges)
{
  if (!is_valid_proc_size_for_test() || parallel_size() == 1) return;

  perform_iterations_of_uniform_refinement_with_general_element_marker(1);

  // A parent on processor 0 and its children move to processor 1.
  std::vector<stk::mesh::EntityId> idsOfElemsToMove;
  stk::mesh::EntityId parentId = 0;
  if (0 == mMesh.parallel_rank())
  {
    const stk::mesh::Selector ownedParents =
        mMesh.mesh_meta_data().locally_owned_part() & myRefinement.parent_part();
    const stk::mesh::Entity parent =
        (*mMesh.get_buckets(stk::topology::ELEMENT_RANK, ownedParents)[0])[0];
    parentId = mMesh.identifier(parent);
    idsOfElemsToMove.push_back(parentId);
    for (auto child : get_children(parent))
      idsOfElemsToMove.push_back(mMesh.identifier(child));
  }
  all_reduce_max(mMesh.parallel(), parentId);
  move_owned_elements_with_given_ids_and_owned_attached_entities_to_processor(idsOfElemsToMove, 1);
  check_edges_of_active_elements();

  mark_all_elements_for_unrefinement();
  refine_marked_elements();

  const stk::mesh::Entity parent = mMesh.get_entity(stk::topology::ELEMENT_RANK, parentId);
  if (mMesh.is_valid(parent) && mMesh.bucket(parent).owned())
  {
    EXPECT_EQ(0, mMesh.parallel_rank());
  }
  EXPECT_EQ(8u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  EXPECT_EQ(num_original_edges(), num_edges_in(mMesh.mesh_meta_data().universal_part()));
}

} // namespace krino
