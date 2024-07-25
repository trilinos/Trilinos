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
#include <Kokkos_Timer.hpp>
#include <vector>
#include <utility>

#define ALWAYS_SORT_MORTON_TREES 0

namespace stk::search {

template <typename RealType, typename ExecutionSpace>
inline void determine_mas_calling(const MortonAabbTree<RealType, ExecutionSpace> &partialTree1,
                                  const MortonAabbTree<RealType, ExecutionSpace> &partialTree2,
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

template <typename RealType, typename ExecutionSpace, typename BoxType, typename IdentProcType>
inline void export_from_box_ident_proc_vec_to_morton_tree(
    const std::vector<std::pair<BoxType, IdentProcType>> &boxIdentProcVec,
    MortonAabbTree<RealType, ExecutionSpace> &tree)
{
  int numBoxes = static_cast<int>(boxIdentProcVec.size());
  tree.reset(numBoxes);

  for (int i = 0; i < numBoxes; ++i) {
    const auto & [box, identProc] = boxIdentProcVec[i];
    tree.host_set_box(i, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), box.get_z_min(),
                      box.get_z_max());
  }
}

template <typename RealType, typename BoxType, typename IdentType, typename ExecutionSpace>
inline void export_from_box_ident_vector_to_morton_tree(
    const std::vector<std::pair<BoxType, IdentType>> &boxIdentList, MortonAabbTree<RealType, ExecutionSpace> &tree)
{
  static_assert(Kokkos::SpaceAccessibility<typename ExecutionSpace::memory_space, Kokkos::HostSpace>::assignable);
  int numBoxes = static_cast<int>(boxIdentList.size());
  tree.reset(numBoxes);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, numBoxes), KOKKOS_LAMBDA(int index) {
        const auto &box = boxIdentList[index].first;
        tree.device_set_box(index, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), box.get_z_min(),
            box.get_z_max());
      });
}

template <typename RealType, typename ExecutionSpace, typename BoxType, typename IdentType>
inline void export_from_box_ident_view_to_morton_tree(
    const Kokkos::View<BoxIdent<BoxType, IdentType>*, ExecutionSpace> & boxIdentList,
    MortonAabbTree<RealType, ExecutionSpace> & tree)
{
  int numBoxes = static_cast<int>(boxIdentList.extent(0));
  tree.reset(numBoxes);

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, numBoxes),
      KOKKOS_LAMBDA(int index) {
        const auto & box = boxIdentList[index].box;
        tree.device_set_box(index, box.get_x_min(), box.get_x_max(),
                            box.get_y_min(), box.get_y_max(),
                            box.get_z_min(), box.get_z_max());
      });
}

template <typename RealType, typename ExecutionSpace, typename BoxType>
inline void export_from_box_vec_to_morton_tree(const std::vector<BoxType> &boxVec,
                                               MortonAabbTree<RealType, ExecutionSpace> &tree)
{
  int numBoxes = static_cast<int>(boxVec.size());
  tree.reset(numBoxes);
  for (int i = 0; i < numBoxes; ++i) {
    const BoxType &box = boxVec[i];
    tree.host_set_box(i, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), box.get_z_min(),
                      box.get_z_max());
  }
}

template <typename RealType, typename ExecutionSpace>
inline void morton_lbvh_search(MortonAabbTree<RealType, ExecutionSpace> &tree1,
                               MortonAabbTree<RealType, ExecutionSpace> &tree2,
                               CollisionList<ExecutionSpace> &searchResults,
                               ExecutionSpace const& execSpace = ExecutionSpace{})
{
  Kokkos::Profiling::pushRegion("Initialization");
  Kokkos::Profiling::pushRegion("Get global bounds");
  // Get total bounds
  TotalBoundsFunctor<RealType, ExecutionSpace>::apply(tree1, execSpace);
  TotalBoundsFunctor<RealType, ExecutionSpace>::apply(tree2, execSpace);
  Kokkos::fence();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Determine need to sort/flip trees");
  bool sortTree1, sortTree2, flipOrder;
  determine_mas_calling(tree1, tree2, sortTree1, sortTree2, flipOrder);
  Kokkos::Profiling::popRegion();

  // Morton encode the centroids of the leaves
  Kokkos::Profiling::pushRegion("Morton encoding of leaves");
  MortonEncoder<RealType, ExecutionSpace>::apply(tree1, execSpace, sortTree1);
  MortonEncoder<RealType, ExecutionSpace>::apply(tree2, execSpace, sortTree2);
  Kokkos::fence();
  Kokkos::Profiling::popRegion();

  // Sort the leaves if appropriate
  Kokkos::Profiling::pushRegion("Sort the trees");
  if (sortTree1) {
    // printf("Sorting tree with %d leaves\n", tree1.hm_numLeaves());
    SortByCode<RealType, ExecutionSpace>::apply(tree1, execSpace);
  }
  if (sortTree2) {
    // printf("Sorting tree with %d leaves\n", tree1.hm_numLeaves());
    SortByCode<RealType, ExecutionSpace>::apply(tree2, execSpace);
  }
  Kokkos::fence();
  Kokkos::Profiling::popRegion();

  // Build the tree structures, if appropriate, following Karras's algorithm
  Kokkos::Profiling::pushRegion("Build the tree structures");
  bool buildTree1 = (sortTree1 && flipOrder);
  bool buildTree2 = (sortTree2 && !flipOrder);
  if (buildTree1) {
    BuildRadixTree<RealType, ExecutionSpace>::apply(tree1, execSpace);
  }
  if (buildTree2) {
    BuildRadixTree<RealType, ExecutionSpace>::apply(tree2, execSpace);
  }
  Kokkos::fence();
  Kokkos::Profiling::popRegion();

  // Augment the trees to be bounding volume (box) hierarchies
  Kokkos::Profiling::pushRegion("Augment the trees to be bounding volume hierarchies");
  if (buildTree1) {
    UpdateInteriorNodeBVs<RealType, ExecutionSpace>::apply(tree1, execSpace);
  }
  if (buildTree2) {
    UpdateInteriorNodeBVs<RealType, ExecutionSpace>::apply(tree2, execSpace);
  }
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();

  // Test the boxes from the non-tree against the tree that was built.
  Kokkos::Profiling::pushRegion("Search query");
  if (flipOrder) {
    Traverse_MASTB_BVH_Functor<RealType, ExecutionSpace>::apply_tree(tree2, tree1, searchResults, execSpace, true);
  }
  else {
    Traverse_MASTB_BVH_Functor<RealType, ExecutionSpace>::apply_tree(tree1, tree2, searchResults, execSpace);
  }
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
}

template <typename RealType, class ExecutionSpace, typename BoxType>
inline void morton_lbvh_search(const std::vector<BoxType> &boxA,
                               const std::vector<BoxType> &boxB,
                               CollisionList<ExecutionSpace> &searchResults,
                               ExecutionSpace const& execSpace = ExecutionSpace{})
{
  Kokkos::Profiling::pushRegion("morton_lbvh_search: export boxes to trees");
  MortonAabbTree<RealType, ExecutionSpace> mlbvhA("a"), mlbvhB("b");
  export_from_box_vec_to_morton_tree(boxA, mlbvhA);
  export_from_box_vec_to_morton_tree(boxB, mlbvhB);
  mlbvhA.sync_to_device();
  mlbvhB.sync_to_device();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("morton_lbvh_search: execute search");
  morton_lbvh_search<RealType, ExecutionSpace>(mlbvhA, mlbvhB, searchResults, execSpace);
  Kokkos::Profiling::popRegion();
}

}  // namespace stk::search

#endif  // STK_SEARCH_MORTON_LBVH_MORTONLBVH_TREE_HPP
