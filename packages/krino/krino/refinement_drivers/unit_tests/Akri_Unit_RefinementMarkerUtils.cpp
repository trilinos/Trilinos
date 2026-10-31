// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MeshSpecs.hpp>
#include <Akri_RefinementManager.hpp>
#include <Akri_RefinementMarkerUtils.hpp>
#include <Akri_RefinementTreeElementParts.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

#include <algorithm>
#include <set>

namespace krino {

TEST(RefinementMarkerUtils, determineRefinementMarker)
{
  constexpr bool selected = true;
  constexpr bool notSelected = false;
  constexpr bool defaultCoarsen = true;
  constexpr bool defaultNothing = false;

  // Selected elements are refined until they reach the level, then left alone.
  EXPECT_EQ(Refinement_Marker::REFINE, determine_refinement_marker(selected, 2, 0, defaultCoarsen));
  EXPECT_EQ(Refinement_Marker::REFINE, determine_refinement_marker(selected, 2, 1, defaultNothing));
  EXPECT_EQ(Refinement_Marker::NOTHING, determine_refinement_marker(selected, 2, 2, defaultCoarsen));
  EXPECT_EQ(Refinement_Marker::NOTHING, determine_refinement_marker(selected, 2, 3, defaultCoarsen));

  // Unselected elements get the default.
  EXPECT_EQ(Refinement_Marker::COARSEN, determine_refinement_marker(notSelected, 2, 1, defaultCoarsen));
  EXPECT_EQ(Refinement_Marker::NOTHING, determine_refinement_marker(notSelected, 2, 1, defaultNothing));
}

TEST(RefinementMarkerUtils, shouldContinueToCoarsen_onlyForTheFirstTargetLevelIterations)
{
  EXPECT_TRUE(should_continue_to_coarsen(0, 2));
  EXPECT_TRUE(should_continue_to_coarsen(1, 2));
  EXPECT_FALSE(should_continue_to_coarsen(2, 2));
  EXPECT_FALSE(should_continue_to_coarsen(0, 0));
}

class FourRightTrisRefinementMarker : public StkMeshTriFixture
{
protected:
  FourRightTrisRefinementMarker()
  {
    build_mesh_from_spec(meshSpec);
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    meta.enable_late_fields();
    myRefinement = &RefinementManager::create(meta);
    meta.disable_late_fields();
  }

  int get_marker(const stk::mesh::Entity elem) const
  {
    return myRefinement->get_marker_field_and_sync_to_host().data<int>().entity_values(elem)();
  }

  std::vector<stk::mesh::Entity> owned_elements_with_indices(const std::vector<unsigned> & indices) const
  {
    std::vector<stk::mesh::Entity> elems;
    for (unsigned index : indices)
      if (is_owned_element(index))
        elems.push_back(get_assigned_element_for_index(index));
    return elems;
  }

  std::vector<stk::mesh::Entity> get_children(const stk::mesh::Entity parent) const
  {
    std::vector<stk::mesh::Entity> children;
    myRefinement->fill_children(parent, children);
    return children;
  }

  FourRightTris meshSpec;
  RefinementManager * myRefinement{nullptr};
};

TEST_F(FourRightTrisRefinementMarker, initializeMarker_fillsWithDefault)
{
  if (!is_valid_proc_size_for_test()) return;

  initialize_marker(mMesh, *myRefinement, true);
  for (auto elem : get_owned_elements())
    EXPECT_EQ(static_cast<int>(Refinement_Marker::COARSEN), get_marker(elem));

  initialize_marker(mMesh, *myRefinement, false);
  for (auto elem : get_owned_elements())
    EXPECT_EQ(static_cast<int>(Refinement_Marker::NOTHING), get_marker(elem));
}

TEST_F(FourRightTrisRefinementMarker, markElementsToReachRefinementLevel_marksOnlyGivenElementsUntilTheyReachTheLevel)
{
  if (!is_valid_proc_size_for_test()) return;

  // Element 0 is marked on its owner; the markers are checked on the owners.
  const stk::mesh::Entity markedElem = get_assigned_element_for_index(0);
  const std::vector<stk::mesh::Entity> elemsToMark = owned_elements_with_indices({0});

  initialize_marker(mMesh, *myRefinement, false);
  mark_elements_to_reach_refinement_level(mMesh, *myRefinement, elemsToMark, 1, false);
  if (is_owned_element(0))
  {
    EXPECT_EQ(static_cast<int>(Refinement_Marker::REFINE), get_marker(markedElem));
  }
  if (is_owned_element(1))
  {
    EXPECT_EQ(static_cast<int>(Refinement_Marker::NOTHING), get_marker(get_assigned_element_for_index(1)));
  }

  myRefinement->do_refinement();

  // The marked element is now refined to level 1: its children are at the level, so are left alone.
  std::vector<stk::mesh::Entity> ownedChildren;
  if (is_owned_element(0))
    for (auto child : get_children(markedElem))
      if (mMesh.bucket(child).owned())
        ownedChildren.push_back(child);
  EXPECT_LT(0u, stk::get_global_sum(mComm, ownedChildren.size()));

  initialize_marker(mMesh, *myRefinement, false);
  mark_elements_to_reach_refinement_level(mMesh, *myRefinement, ownedChildren, 1, false);
  for (auto child : ownedChildren)
  {
    EXPECT_EQ(1, myRefinement->fully_refined_level(child));
    EXPECT_EQ(static_cast<int>(Refinement_Marker::NOTHING), get_marker(child));
  }
}

TEST_F(FourRightTrisRefinementMarker, ownedRefinementParentsAndTheirElementPart_eachParentGetsPartChosenFromItsChildren)
{
  if (!is_valid_proc_size_for_test()) return;

  // Two levels: refine element 0, then one of its children, each marked on its owner.
  const stk::mesh::Entity root = get_assigned_element_for_index(0);
  mark_elements_for_refinement(*myRefinement, owned_elements_with_indices({0}));
  myRefinement->do_refinement();

  std::vector<stk::mesh::Entity> childToRefine;
  if (is_owned_element(0))
  {
    const std::vector<stk::mesh::Entity> children = get_children(root);
    EXPECT_FALSE(children.empty());
    if (!children.empty() && mMesh.bucket(children[0]).owned())
      childToRefine.push_back(children[0]);
  }
  const stk::mesh::EntityId childToRefineId = stk::get_global_max(mComm, childToRefine.empty() ? stk::mesh::EntityId(0) : mMesh.identifier(childToRefine[0]));
  ASSERT_NE(0u, childToRefineId);
  clear_refinement_marker(*myRefinement);
  mark_elements_for_refinement(*myRefinement, childToRefine);
  myRefinement->do_refinement();

  const unsigned blockOrdinal = mMesh.mesh_meta_data().get_part("block_1")->mesh_meta_data_ordinal();
  constexpr unsigned chosenPart = 1000;
  std::vector<std::set<unsigned>> callbackInputs;
  const auto parentsAndPart = get_owned_refinement_parents_and_their_element_part(mMesh, *myRefinement,
      [&](const std::set<unsigned> & childPartOrdinals) { callbackInputs.push_back(childPartOrdinals); return chosenPart; });

  // Every owned parent element appears once, on its owner, with the chosen part.  (Refinement also makes
  // transition parents of the neighbors of refined elements, possibly on other processors.)
  std::vector<stk::mesh::Entity> parents;
  for (auto && [parent, part] : parentsAndPart)
  {
    EXPECT_TRUE(myRefinement->is_parent(parent));
    EXPECT_TRUE(mMesh.bucket(parent).owned());
    EXPECT_EQ(chosenPart, part);
    parents.push_back(parent);
  }
  std::sort(parents.begin(), parents.end());
  EXPECT_EQ(parents.end(), std::adjacent_find(parents.begin(), parents.end()));

  const auto is_in_parents_on_one_proc = [&](const stk::mesh::Entity elem)
  {
    const bool isInParents = mMesh.is_valid(elem) && std::find(parents.begin(), parents.end(), elem) != parents.end();
    return 1u == stk::get_global_sum(mComm, isInParents ? 1u : 0u);
  };
  EXPECT_TRUE(is_in_parents_on_one_proc(root));
  EXPECT_TRUE(is_in_parents_on_one_proc(mMesh.get_entity(stk::topology::ELEMENT_RANK, childToRefineId)));

  size_t numOwnedParents = 0;
  for (auto * bucket : mMesh.get_buckets(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part()))
    for (auto elem : *bucket)
      if (myRefinement->is_parent(elem))
        ++numOwnedParents;
  EXPECT_EQ(numOwnedParents, parents.size());

  // The callback is called once per parent.  Leaf children are all in block_1; a refined child passes on the
  // part chosen for it, so the root's children include the chosen part.
  ASSERT_EQ(parents.size(), callbackInputs.size());
  bool sawChosenPartFromRefinedChild = false;
  for (auto && childPartOrdinals : callbackInputs)
  {
    for (const unsigned ordinal : childPartOrdinals)
      EXPECT_TRUE(ordinal == blockOrdinal || ordinal == chosenPart);
    if (childPartOrdinals.count(chosenPart) > 0)
      sawChosenPartFromRefinedChild = true;
  }
  EXPECT_TRUE(stk::is_true_on_any_proc(mComm, sawChosenPartFromRefinedChild));
}

} // namespace krino
