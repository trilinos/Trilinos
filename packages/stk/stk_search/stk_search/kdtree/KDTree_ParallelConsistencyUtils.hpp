#ifndef KDTREE_PARALLELCONSISTENCYUTILS_HPP
#define KDTREE_PARALLELCONSISTENCYUTILS_HPP

#include "stk_util/parallel/Parallel.hpp"
#include "stk_search/kdtree/KDTree_BoundingBox.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include <vector>
#include <utility>

namespace stk::search {

template <typename ObjType, typename IdentifierType, typename BaseBoxType>
void ParallelComputeProcObjectBoxes(const std::vector<std::pair<ObjType, IdentifierType>>& local_objsWithIdents,
                                    std::vector<stk::search::ObjectBoundingBox_T<BaseBoxType>> &objBB_proc_box_array,
                                    MPI_Comm &comm)
{
  const unsigned numBoxDomain = local_objsWithIdents.size();
  stk::search::ObjectBoundingBox_T<BaseBoxType> objBB_proc;

#ifdef _OPENMP
  std::vector<stk::search::ObjectBoundingBox_T<BaseBoxType> > threadBoxes( omp_get_max_threads() );
#endif

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP
    stk::search::ObjectBoundingBox_T<BaseBoxType>& curBox = threadBoxes[omp_get_thread_num()];
#else
    stk::search::ObjectBoundingBox_T<BaseBoxType>& curBox = objBB_proc;
#endif
#ifdef _OPENMP
#pragma omp for
#endif
    for(unsigned ibox = 0; ibox < numBoxDomain; ++ibox) {
      stk::search::add_to_box(curBox.GetBox(), local_objsWithIdents[ibox].first);
    }
  }

#ifdef _OPENMP
  for(unsigned i=0; i<threadBoxes.size(); ++i) {
    stk::search::add_to_box(objBB_proc.GetBox(), threadBoxes[i].GetBox());
  }
#endif
  //
  //  Do a global communication to communicate all processor boxA bounding boxes
  //  to all processors in the group
  //
  stk::search::all_gather_helper(objBB_proc, objBB_proc_box_array, comm);

  int numProcs;
  MPI_Comm_size(comm, &numProcs);
  for (int iproc = 0; iproc < numProcs; ++iproc) {
    objBB_proc_box_array[iproc].set_object_number(iproc);
  }
}


// Ghost the range boxes needed for each processor to search against its local domain boxes
// in distributed AABB overlap search ("coarse search")..
template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void
ComputeRangeWithGhostsForCoarseSearch(const std::vector<std::pair<DomainObjType, DomainIdentifier>>& local_domain,
                                      const std::vector<std::pair<RangeBoxType, RangeIdentifier>>& local_range,
                                      int num_procs,
                                      std::vector<RangeBoxType>& rangeBoxes,
                                      std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm)
{
  const unsigned numBoxRange  = local_range.size();

  using domainValueType = typename DomainObjType::value_type;
  using DomainBox = stk::search::Box<domainValueType>;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < numBoxRange; ++i) {
    rangeBoxes[i] = local_range[i].first;
  }

  //
  //  Determine the total number of processors involved in the communication and the current processor number
  //
  int current_proc(0);
  MPI_Comm_rank(comm, &current_proc);
  if(num_procs == 0) {
    return;
  }

  //
  //  Compute the processor local bounding boxes for the box sets
  //  Store the boxes in unique entries in a global processor bounding box array.
  //
  std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array;
  ParallelComputeProcObjectBoxes(local_domain, boxA_proc_box_array, comm);

  //
  //  Create a hierarchy of boxA processor bounding boxes.
  //  This hierarchy will be used to search for overlaps between processors and
  //  objects.
  //
  stk::search::ProximitySearchTree_T<DomainBox> boxA_box_hierarchy(boxA_proc_box_array);

  //
  //  Determine what to ghost.  If a boxB box from this processor overlaps another processor's
  //  processor all-boxA box, then we need to ghost the boxB's data to that other processor
  //  to include in the range boxes (local + global) to search against its local domain boxes.
  //
  typedef typename RangeIdentifier::ident_type GlobalIdType;
  typedef std::pair<RangeBoxType, GlobalIdType> BoxIdPair;
  std::vector<std::vector<BoxIdPair>> send_list(num_procs);
  std::vector<std::vector<BoxIdPair>> recv_list(num_procs);

  std::vector<int> proc_list(num_procs);
  for (unsigned int iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    boxA_box_hierarchy.SearchForOverlap(local_range[iboxB].first, proc_list);
    for (auto&& overlapping_proc : proc_list) {
      if (overlapping_proc == current_proc) continue;
      GlobalIdType id = local_range[iboxB].second.id();
      send_list[overlapping_proc].push_back(BoxIdPair(rangeBoxes[iboxB], id));
    }
  }

  stk::parallel_data_exchange_t(send_list, recv_list, comm);
  rangeGhostIdentifiers.clear();
  for (size_t i = 0; i < recv_list.size(); ++i) {
    for (size_t j = 0; j < recv_list[i].size(); ++j) {
      const BoxIdPair& recvd_boxIdPair = recv_list[i][j];
      rangeBoxes.push_back(recvd_boxIdPair.first);
      rangeGhostIdentifiers.push_back(RangeIdentifier(recvd_boxIdPair.second, i));
    }
  }
}

}

#endif // KDTREE_PARALLELCONSISTENCYUTILS_HPP
