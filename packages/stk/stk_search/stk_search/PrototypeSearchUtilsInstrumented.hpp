#ifndef PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_
#define PROTOTYPE_SEARCH_UTILS_INSTURMENTED_H_

namespace stk {
namespace search {
namespace experimental {


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
class GhostingSearcher {

  typedef typename DomainObjType::value_type domainValueType;
  typedef stk::search::Box<domainValueType>  DomainBox;
  typedef typename RangeIdentifier::ident_type GlobalIdType;
  typedef std::pair<RangeBoxType, GlobalIdType> BoxIdPair;

 public:
  GhostingSearcher(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm,
      instrumented::GhostingSearchTimeBreakdown &timeBreakdown);

  void searchFromScratch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      instrumented::GhostingSearchTimeBreakdown &timeBreakdown);

  void updateAndSearch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      instrumented::GhostingSearchTimeBreakdown &timeBreakdown);

 private:
  MPI_Comm mpiComm;
  int      pSize;
  int      pRank;

  void findNeighborsFromScratch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
      instrumented::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void intersectLocalRangeBoxesWithNeighborBVs(
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array,
      const unsigned numBoxRange,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<std::vector<BoxIdPair> >& send_list,
      std::vector<RangeBoxType>& rangeBoxes,
      instrumented::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void commIntersectingRangeBoxesToNeighbors(
      std::vector<std::vector<BoxIdPair> > send_list,
      std::vector<std::vector<BoxIdPair> >& recv_list,
      instrumented::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void writeRangeBoxesAndGhostIdentifiers(
      const std::vector<std::vector<BoxIdPair> >& recv_list,
      SplitTimer stopwatch, std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      std::vector<RangeBoxType>& rangeBoxes,
      instrumented::GhostingSearchTimeBreakdown& timeBreakdown);
};


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::GhostingSearcher(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm,
      instrumented::GhostingSearchTimeBreakdown &timeBreakdown)
      : mpiComm(comm), pSize(-1), pRank(-1)
{
  MPI_Comm_rank(mpiComm, &pRank);
  MPI_Comm_size(mpiComm, &pSize);
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
    instrumented::GhostingSearchTimeBreakdown &timeBreakdown,
    SplitTimer &stopwatch)
{
  //
  //  Compute the processor local bounding boxes for the box sets
  //  Store the boxes in unique entries in a global processor bounding box array.
  //

  stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc;
#ifdef _OPENMP
  std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > threadBoxes( omp_get_max_threads() );
#endif

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP
    stk::search::ObjectBoundingBox_T<DomainBox>& curBox = threadBoxes[omp_get_thread_num()];
#else
    stk::search::ObjectBoundingBox_T<DomainBox>& curBox = boxA_proc;
#endif
#ifdef _OPENMP
#pragma omp for
#endif
    const unsigned numBoxDomain = local_domain.size();
    for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
      stk::search::add_to_box(curBox.GetBox(), local_domain[iboxA].first);
    }
  }

#ifdef _OPENMP
  for(unsigned i=0; i<threadBoxes.size(); ++i) {
    stk::search::add_to_box(boxA_proc.GetBox(), threadBoxes[i].GetBox());
  }
#endif

  boxA_proc_box_array[pRank] = boxA_proc;

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  //
  //  Do a global communication to communicate all processor boxA bounding boxes
  //  to all processors in the group
  //
  instrumented::GlobalBoxCombine(boxA_proc_box_array, mpiComm);
  for(int iproc = 0; iproc < pSize; ++iproc) {
    boxA_proc_box_array[iproc].set_object_number(iproc);
  }
  timeBreakdown.findNeighbors += stopwatch.split();

}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
  RangeBoxType>::intersectLocalRangeBoxesWithNeighborBVs(
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array,
    const unsigned numBoxRange,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<std::vector<BoxIdPair> >& send_list,
    std::vector<RangeBoxType>& rangeBoxes,
    instrumented::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch)
{
  //  Create a hierarchy of boxA processor bounding boxes.
  //  This hierarchy will be used to search for overlaps between processors and
  //  objects.
  stk::search::ProximitySearchTree_T<DomainBox> boxA_box_hierarchy(
      boxA_proc_box_array);

  //  Determine what to ghost.  If a boxB box from this processor overlaps another processor's
  //  processor all-boxA box, then we need to ghost the boxB's data to that other processor.
  //  (The stricter criteria used by the full BoxA_BoxB_Ghost function would make sure that
  //  the individual boxB box overlaps some individual boxA box from the other processor.)
  std::vector<int> proc_list(pSize);
  for (unsigned int iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    boxA_box_hierarchy.SearchForOverlap(local_range[iboxB].first, proc_list);
    for (unsigned i = 0; i < proc_list.size(); ++i) {
      int overlapping_proc = proc_list[i];
      if (overlapping_proc == pRank)
        continue;

      GlobalIdType id = local_range[iboxB].second.id();
      send_list[overlapping_proc].push_back(BoxIdPair(rangeBoxes[iboxB], id));
    }
  }
  timeBreakdown.intersectLocalRangeBoxesWithNeighborBVs += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier,
    typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
    RangeBoxType>::commIntersectingRangeBoxesToNeighbors(
    std::vector<std::vector<BoxIdPair> > send_list,
    std::vector<std::vector<BoxIdPair> >& recv_list,
    instrumented::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch) {
  stk::parallel_data_exchange_t(send_list, recv_list, mpiComm);
  timeBreakdown.communicateRangeBoxesToNeighbors += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier,
    typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
    RangeBoxType>::writeRangeBoxesAndGhostIdentifiers(
    const std::vector<std::vector<BoxIdPair> >& recv_list, SplitTimer stopwatch,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    std::vector<RangeBoxType>& rangeBoxes,
    instrumented::GhostingSearchTimeBreakdown& timeBreakdown) {
  rangeGhostIdentifiers.clear();
  for (size_t i = 0; i < recv_list.size(); i++) {
    for (size_t j = 0; j < recv_list[i].size(); j++) {
      const BoxIdPair& recvd_boxIdPair = recv_list[i][j];
      rangeBoxes.push_back(recvd_boxIdPair.first);
      rangeGhostIdentifiers.push_back(
          RangeIdentifier(recvd_boxIdPair.second, i));
    }
  }
  timeBreakdown.writeRangeBoxesAndGhostIdentifiers += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::searchFromScratch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<RangeBoxType >& rangeBoxes,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    instrumented::GhostingSearchTimeBreakdown &timeBreakdown)
{
  // Pretend-o-type for TDD.

  const unsigned numBoxRange  = local_range.size();

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < numBoxRange; i++) {
    rangeBoxes[i] = local_range[i].first;
  }
  if(pSize == 0) {
    return;
  }

  SplitTimer stopwatch;

  std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array(pSize);
  findNeighborsFromScratch(local_domain, boxA_proc_box_array, timeBreakdown, stopwatch);

  std::vector<std::vector<BoxIdPair> > send_list(pSize);
  intersectLocalRangeBoxesWithNeighborBVs(boxA_proc_box_array, numBoxRange,
                                          local_range, send_list, rangeBoxes,
                                          timeBreakdown, stopwatch);

  std::vector<std::vector<BoxIdPair> > recv_list(pSize);
  commIntersectingRangeBoxesToNeighbors(send_list, recv_list, timeBreakdown,
                                        stopwatch);

  writeRangeBoxesAndGhostIdentifiers(recv_list, stopwatch,
                                     rangeGhostIdentifiers, rangeBoxes,
                                     timeBreakdown);
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::updateAndSearch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<RangeBoxType >& rangeBoxes,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    instrumented::GhostingSearchTimeBreakdown &timeBreakdown)
{
  // Pretend-o-type!
  searchFromScratch(local_domain, local_range, rangeBoxes, rangeGhostIdentifiers, timeBreakdown);
}


}}}

#endif // PROTOTYPE_SEARCH_UTILS_INSTURMENTED_H_
