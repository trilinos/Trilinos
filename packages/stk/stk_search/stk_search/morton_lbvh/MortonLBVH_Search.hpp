// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_SEARCH_MORTON_LBVH_MORTONLBVH_SEARCH_HPP
#define STK_SEARCH_MORTON_LBVH_MORTONLBVH_SEARCH_HPP

#include <stk_search/morton_lbvh/MortonLBVH_Tree.hpp>
#include <stk_search/morton_lbvh/MortonLBVH_TreeManipulationUtils.hpp>
#include <stk_search/morton_lbvh/MortonLBVH_CollisionList.hpp>
#include <stk_search/BoxIdent.hpp>
#include <stk_search/HelperTraits.hpp>
#include <Kokkos_Timer.hpp>
#include <vector>
#include <utility>

#define ALWAYS_SORT_MORTON_TREES 0

namespace stk::search {

constexpr size_t COLLISION_SCALE_FACTOR = 16;

template <typename Tree1Type, typename Tree2Type, typename ExecutionSpace>
inline void determine_mas_calling(const Tree1Type &partialTree1,
                                  const Tree2Type &partialTree2,
                                  bool &sortTree1, bool &sortTree2, bool &flipTrees)
{
  if (ALWAYS_SORT_MORTON_TREES) {
    sortTree1 = sortTree2 = true;
    flipTrees = false;
    return;
  }

  LocalOrdinal tree1NumLeaves = partialTree1.hm_numLeaves();
  LocalOrdinal tree2NumLeaves = partialTree2.hm_numLeaves();

  if ((tree1NumLeaves == 0) || (tree2NumLeaves == 0)) {
    // Search will shortcut out.
    sortTree1 = sortTree2 = flipTrees = false;
    return;
  }

  flipTrees = (tree1NumLeaves < tree2NumLeaves);

  const LocalOrdinal magicThreshold = 64;
  if ((tree1NumLeaves <= magicThreshold) && (tree2NumLeaves <= magicThreshold)) {
    sortTree1 = flipTrees;
    sortTree2 = !flipTrees;
    return;
  }

  const double divergenceFactor = 32.0;
  const double magicRatio = (divergenceFactor - 1) * 0.2;

  if (flipTrees) {
    sortTree1 = true;
    double depthRatio = log(tree2NumLeaves);
    if (tree1NumLeaves > 2) {
      depthRatio = depthRatio / log(tree1NumLeaves);
    }
    sortTree2 = (depthRatio < magicRatio);
  }
  else {
    sortTree2 = true;
    double depthRatio = log(tree1NumLeaves);
    if (tree2NumLeaves > 2) {
      depthRatio = depthRatio / log(tree2NumLeaves);
    }
    sortTree1 = (depthRatio < magicRatio);
  }
}

template <typename TreeType, typename BoxType, typename IdentProcType>
inline void export_from_box_ident_proc_vec_to_morton_tree(
    const std::vector<std::pair<BoxType, IdentProcType>> &boxIdentProcVec,
    TreeType &tree)
{
  int numBoxes = static_cast<int>(boxIdentProcVec.size());
  tree.reset(numBoxes);

  for (int i = 0; i < numBoxes; ++i) {
    const auto & [box, identProc] = boxIdentProcVec[i];
    tree.host_set_box(i, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), box.get_z_min(),
                      box.get_z_max());
  }
}

template <typename ViewType, typename TreeType, typename ExecutionSpace>
inline void export_box_ident_view_to_morton_tree(
  const ViewType& boxIdentProcs,
  TreeType& tree,
  ExecutionSpace execSpace)
{
  check_view_is_usable_from<ExecutionSpace, ViewType>();
  static_assert(is_box_ident_proc_container_v<ViewType> || is_box_ident_container_v<ViewType>,
                "view must be a View<BoxIdent> or View<BoxIdentProc>");
  using BoxType = typename ViewType::value_type::box_type;
  tree.reset(boxIdentProcs.extent(0));

  Kokkos::RangePolicy<ExecutionSpace> policy(execSpace, 0, boxIdentProcs.extent(0));
  auto func = KOKKOS_LAMBDA(int index)
  {
    const BoxType box = boxIdentProcs(index).box;
    tree.device_set_box(index, box.get_x_min(), box.get_y_min(), box.get_z_min(),
                               box.get_x_max(), box.get_y_max(), box.get_z_max());
  };

  Kokkos::parallel_for("export box-ident view to tree", policy, func);
}

template <typename TreeType, typename BoxType, typename IdentType, typename ExecutionSpace>
inline void export_from_box_ident_vector_to_morton_tree(
    const std::vector<std::pair<BoxType, IdentType>> &boxIdentList, TreeType &tree)
{
  static_assert(Kokkos::SpaceAccessibility<typename ExecutionSpace::memory_space, Kokkos::HostSpace>::assignable);
  int numBoxes = static_cast<int>(boxIdentList.size());
  tree.reset(numBoxes);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, numBoxes), KOKKOS_LAMBDA(int index) {
        const auto &box = boxIdentList[index].first;
        tree.device_set_box(index, box.get_x_min(), box.get_y_min(), box.get_z_min(),
                                   box.get_x_max(), box.get_y_max(), box.get_z_max());
      });
}


template <typename TreeType, typename ExecutionSpace, typename BoxType>
inline void export_from_box_vec_to_morton_tree(const std::vector<BoxType> &boxVec,
                                               TreeType &tree)
{
  int numBoxes = static_cast<int>(boxVec.size());
  tree.reset(numBoxes);
  for (int i = 0; i < numBoxes; ++i) {
    const BoxType &box = boxVec[i];
    tree.host_set_box(i, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), box.get_z_min(),
                      box.get_z_max());
  }
}

template <typename TreeType, typename ExecutionSpace, typename BoxView>
inline void export_box_view_to_morton_tree(const BoxView boxes,
                                                TreeType& tree,
                                                ExecutionSpace execSpace)
{
  check_view_is_usable_from<ExecutionSpace, BoxView>();
  using BoxType = typename BoxView::value_type;
  tree.reset(boxes.extent(0));

  Kokkos::RangePolicy<ExecutionSpace> policy(execSpace, 0, boxes.extent(0));
  auto func = KOKKOS_LAMBDA (int index)
  {
    const BoxType box = boxes(index);
    tree.device_set_box(index, box.get_x_min(), box.get_y_min(), box.get_z_min(),
                               box.get_x_max(), box.get_y_max(), box.get_z_max());
  };

  Kokkos::parallel_for(policy, func);
}

template <typename View1Type, typename View2Type, typename ExecutionSpace, typename Callback>
inline void morton_lbvh_search(MortonAabbTree<View1Type, ExecutionSpace> &tree1,
                               MortonAabbTree<View2Type, ExecutionSpace> &tree2,
                               Callback& resultCallback,
                               ExecutionSpace const& execSpace = ExecutionSpace{})
{
  Kokkos::Profiling::pushRegion("Initialize the trees");
  Kokkos::Profiling::pushRegion("Determine need to sort/flip trees");
  bool sortTree1, sortTree2, flipOrder;
  determine_mas_calling<MortonAabbTree<View1Type, ExecutionSpace>,
                        MortonAabbTree<View2Type, ExecutionSpace>,
                        ExecutionSpace>(tree1, tree2, sortTree1, sortTree2, flipOrder);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Get total bounds");
  TotalBoundsFunctor<View1Type, ExecutionSpace>::apply(tree1, execSpace);
  TotalBoundsFunctor<View2Type, ExecutionSpace>::apply(tree2, execSpace);
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Morton encoding of leaves");
  MortonEncoder<View1Type, ExecutionSpace>::apply(tree1, execSpace, sortTree1);
  MortonEncoder<View2Type, ExecutionSpace>::apply(tree2, execSpace, sortTree2);
  Kokkos::Profiling::popRegion();

  execSpace.fence();

  if (sortTree1) {
    Kokkos::Profiling::pushRegion("Sort the trees");
    SortByCode<MortonAabbTree<View1Type,ExecutionSpace>, ExecutionSpace>::apply(tree1, execSpace);
    Kokkos::Profiling::popRegion();
  }

  if (sortTree2) {
    Kokkos::Profiling::pushRegion("Sort the trees");
    SortByCode<MortonAabbTree<View2Type,ExecutionSpace>, ExecutionSpace>::apply(tree2, execSpace);
    Kokkos::Profiling::popRegion();
  }

  execSpace.fence();

  // Build the tree structures, if appropriate, following Karras's algorithm
  Kokkos::Profiling::pushRegion("Build the tree structures");
  bool buildTree1 = (sortTree1 && flipOrder);
  bool buildTree2 = (sortTree2 && !flipOrder);
  using Box1Type = typename View1Type::value_type::box_type;
  using Real1Type = typename Box1Type::value_type;
  using Box2Type = typename View2Type::value_type::box_type;
  using Real2Type = typename Box2Type::value_type;
  if (buildTree1) {
    BuildRadixTree<Real1Type, MortonAabbTree<View1Type, ExecutionSpace>, ExecutionSpace>::apply(tree1, execSpace);
  }
  if (buildTree2) {
    BuildRadixTree<Real2Type, MortonAabbTree<View2Type, ExecutionSpace>, ExecutionSpace>::apply(tree2, execSpace);
  }
  execSpace.fence();
  Kokkos::Profiling::popRegion();

  // Augment the trees to be bounding volume (box) hierarchies
  Kokkos::Profiling::pushRegion("Augment the trees to be bounding volume hierarchies");
  if (buildTree1) {
    UpdateInteriorNodeBVs<View1Type, ExecutionSpace>::apply(tree1, execSpace);
  }
  if (buildTree2) {
    UpdateInteriorNodeBVs<View2Type, ExecutionSpace>::apply(tree2, execSpace);
  }
  execSpace.fence();
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  // Test the boxes from the non-tree against the tree that was built.
  Kokkos::Profiling::pushRegion("Search query");

  if (flipOrder) {
    search_tree(tree2, tree1, resultCallback, execSpace, true);
  } else {
    search_tree(tree1, tree2, resultCallback, execSpace);
  }
  execSpace.fence();
  Kokkos::Profiling::popRegion();
}

template <typename Tree1Type, typename Tree2Type, typename ExecutionSpace>
inline void morton_lbvh_search(Tree1Type &tree1,
                               Tree2Type &tree2,
                               CollisionList<ExecutionSpace> &searchResults,
                               ExecutionSpace const& execSpace = ExecutionSpace{})
{
  if (searchResults.get_capacity() == 0) {
    const int numDomainLeaves = tree1.hm_numLeaves();
    const int numRangeLeaves  = tree2.hm_numLeaves();

    const int collisionEstimate = std::max(numDomainLeaves, numRangeLeaves) * COLLISION_SCALE_FACTOR;
    searchResults.reset(collisionEstimate);
  }

  CollisionListCallback<ExecutionSpace> resultCallback(searchResults);
  morton_lbvh_search(tree1, tree2, resultCallback, execSpace);
  searchResults = resultCallback.get_collision_list();
}

template <typename RealType, class ExecutionSpace, typename BoxType>
inline void morton_lbvh_search(const std::vector<BoxType> &boxA,
                               const std::vector<BoxType> &boxB,
                               CollisionList<ExecutionSpace> &searchResults,
                               ExecutionSpace const& execSpace = ExecutionSpace{})
{
  Kokkos::Profiling::pushRegion("morton_lbvh_search: export boxes to trees");
  using ViewType = Kokkos::View<BoxIdent<BoxType,int>*,ExecutionSpace>;
  using TreeType = MortonAabbTree<ViewType, ExecutionSpace>;
  TreeType mlbvhA("a"), mlbvhB("b");
  export_from_box_vec_to_morton_tree<TreeType,ExecutionSpace,BoxType>(boxA, mlbvhA);
  export_from_box_vec_to_morton_tree<TreeType,ExecutionSpace,BoxType>(boxB, mlbvhB);
  mlbvhA.sync_to_device();
  mlbvhB.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("morton_lbvh_search: execute search");
  morton_lbvh_search<TreeType, TreeType, ExecutionSpace>(mlbvhA, mlbvhB, searchResults, execSpace);
  Kokkos::Profiling::popRegion();
}

}  // namespace stk::search

#endif  // STK_SEARCH_MORTON_LBVH_MORTONLBVH_TREE_HPP
