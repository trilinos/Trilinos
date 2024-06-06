#ifndef MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP
#define MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP

#include "stk_util/parallel/Parallel.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Tree.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Search.hpp"
#include "stk_search/Box.hpp"
#include "stk_search/BoundingBox.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include <vector>
#include <utility>

namespace stk::search {

template <typename DomainBoxType, typename DomainIdentProcType>
std::vector<Box<typename DomainBoxType::value_type>>
gather_all_processor_superset_domain_boxes(const std::vector<std::pair<DomainBoxType, DomainIdentProcType>> & localDomain,
                                           MPI_Comm & comm)
{
  std::vector<Box<typename DomainBoxType::value_type>> globalSupersetBoxes;

  Box<typename DomainBoxType::value_type> localSupersetBox;
  for (const auto & [box, ident] : localDomain) {
    stk::search::add_to_box(localSupersetBox, box);
  }

  stk::search::all_gather_helper(localSupersetBox, globalSupersetBoxes, comm);

  return globalSupersetBoxes;
}


template<typename ExecutionSpace, typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
std::pair<std::vector<RangeBoxType>, std::vector<RangeIdentProcType>>
morton_extend_local_range_with_remote_boxes_that_might_intersect(
    const std::vector<std::pair<DomainBoxType, DomainIdentProcType>> & localDomain,
    const std::vector<std::pair<RangeBoxType, RangeIdentProcType>> & localRange,
    MPI_Comm comm,
    ExecutionSpace const& execSpace)
{
  using DomainValueType = typename DomainBoxType::value_type;

  const int numProcs = stk::parallel_machine_size(comm);
  const int procId = stk::parallel_machine_rank(comm);

  const auto globalSupersetBoxes = gather_all_processor_superset_domain_boxes(localDomain, comm);

  stk::search::MortonAabbTree<DomainValueType, ExecutionSpace> domainTree("Proc Domain Tree",
                                                                                        localRange.size());
  stk::search::MortonAabbTree<DomainValueType, ExecutionSpace> rangeTree("Proc Range Tree",
                                                                                       globalSupersetBoxes.size());

  export_from_box_ident_proc_vec_to_morton_tree(localRange, domainTree);
  export_from_box_vec_to_morton_tree(globalSupersetBoxes, rangeTree);
  domainTree.sync_to_device();
  rangeTree.sync_to_device();

  stk::search::CollisionList<ExecutionSpace> collisionList("Proc Collision List");
  stk::search::morton_lbvh_search<DomainValueType, ExecutionSpace>(domainTree, rangeTree, collisionList, execSpace);
  collisionList.sync_from_device();

  using GlobalIdType = typename RangeIdentProcType::ident_type;
  using BoxIdPair = std::pair<RangeBoxType, GlobalIdType>;
  std::vector<std::vector<BoxIdPair>> sendList(numProcs);
  std::vector<std::vector<BoxIdPair>> recvList(numProcs);

  const unsigned numCollisions = collisionList.hm_idx();

  for (unsigned i = 0; i < numCollisions; ++i) {
    const int entityIndex = collisionList.hm_data(i, 0);
    const int remoteProcId = collisionList.hm_data(i, 1);
    const auto & [localRangeBox, localRangeIdentProc] = localRange[entityIndex];
    if (remoteProcId != procId) {
      sendList[remoteProcId].emplace_back(localRangeBox, localRangeIdentProc.id());
    }
  }

  stk::parallel_data_exchange_t(sendList, recvList, comm);

  std::pair<std::vector<RangeBoxType>, std::vector<RangeIdentProcType>> result;
  auto & [extendedRangeBoxes, remoteRangeIdentProcs] = result;

  extendedRangeBoxes.reserve(localRange.size());
  for (const auto & [box, identProc] : localRange) {
    extendedRangeBoxes.push_back(box);
  }

  for (size_t proc = 0; proc < recvList.size(); ++proc) {
    for (const auto & [box, id] : recvList[proc]) {
      extendedRangeBoxes.push_back(box);
      remoteRangeIdentProcs.emplace_back(id, proc);
    }
  }

  return result;
}

}

#endif // MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP
