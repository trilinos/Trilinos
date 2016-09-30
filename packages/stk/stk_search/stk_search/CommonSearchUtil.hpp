#ifndef COMMON_SEARCH_UTIL_H_
#define COMMON_SEARCH_UTIL_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stk_search/OctTreeOps.hpp>
#include "stk_util/environment/WallTime.hpp"
#include "stk_search/KDTree_BoundingBox.hpp"
#include "stk_search/KDTree.hpp"

namespace stk {
  namespace search {



    template <typename DomainBox>
    inline void GlobalBoxCombine(DomainBox &box_array, MPI_Comm &communicator)
    {
      typedef typename DomainBox::value_type::coordinate_t coordinate_t;

      int num_boxes = static_cast<int>(box_array.size());
      //
      //  Allocate a common set of arrays to perform the reductions on
      //
      int array_length = num_boxes * 3;
      std::vector<coordinate_t> all_box_min_local (array_length);
      std::vector<coordinate_t> all_box_max_local (array_length);
      std::vector<coordinate_t> all_box_min_global(array_length);
      std::vector<coordinate_t> all_box_max_global(array_length);
      //
      //  Fill the local arrays
      //
      for(int ibox = 0; ibox < num_boxes; ++ibox) {
        all_box_min_local[ibox * 3 + 0] = box_array[ibox].GetBox().get_x_min();
        all_box_min_local[ibox * 3 + 1] = box_array[ibox].GetBox().get_y_min();
        all_box_min_local[ibox * 3 + 2] = box_array[ibox].GetBox().get_z_min();
        all_box_max_local[ibox * 3 + 0] = box_array[ibox].GetBox().get_x_max();
        all_box_max_local[ibox * 3 + 1] = box_array[ibox].GetBox().get_y_max();
        all_box_max_local[ibox * 3 + 2] = box_array[ibox].GetBox().get_z_max();
      }
      //
      //  Perform the global MPI reductions
      //
      MPI_Datatype floatType;
      if(sizeof(coordinate_t) == sizeof(float)) {
        floatType = MPI_FLOAT;
      } else if (sizeof(coordinate_t) == sizeof(double)) {
        floatType = MPI_DOUBLE;
      } else {
        floatType = MPI_DOUBLE;
      }
      MPI_Allreduce( all_box_min_local.data(), all_box_min_global.data(), array_length, floatType, MPI_MIN, communicator );
      MPI_Allreduce( all_box_max_local.data(), all_box_max_global.data(), array_length, floatType, MPI_MAX, communicator );
      //
      //  Scatter the local arrays back to the boxes
      //
      for(int ibox = 0; ibox < num_boxes; ++ibox) {
        box_array[ibox].GetBox().set_box(all_box_min_global[ibox * 3 + 0],
                                         all_box_min_global[ibox * 3 + 1],
                                         all_box_min_global[ibox * 3 + 2],
                                         all_box_max_global[ibox * 3 + 0],
                                         all_box_max_global[ibox * 3 + 1],
                                         all_box_max_global[ibox * 3 + 2]);
      }
    }


    template<typename DomainIdentifier, typename RangeIdentifier, typename AxisAlignedBBType>
      void
      ComputeRangeWithGhostsForCoarseSearch(
                                            const std::vector<std::pair<AxisAlignedBBType, DomainIdentifier> >& local_domain,
                                            const std::vector<std::pair<AxisAlignedBBType, RangeIdentifier> >& local_range,
                                            int num_procs,
                                            std::vector<AxisAlignedBBType>& rangeBoxes, std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm)
    {

      const unsigned numBoxRange  = local_range.size();
      const unsigned numBoxDomain = local_domain.size();

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
      for (size_t i = 0; i < numBoxRange; i++) {
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

      stk::search::ObjectBoundingBox_T<AxisAlignedBBType> boxA_proc;
#ifdef _OPENMP
      std::vector<stk::search::ObjectBoundingBox_T<AxisAlignedBBType> > threadBoxes( omp_get_max_threads() );
#endif

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef _OPENMP
        stk::search::ObjectBoundingBox_T<AxisAlignedBBType>& curBox = threadBoxes[omp_get_thread_num()];
#else
        stk::search::ObjectBoundingBox_T<AxisAlignedBBType>& curBox = boxA_proc;
#endif
#ifdef _OPENMP
#pragma omp for
#endif
        for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
          stk::search::add_to_box(curBox.GetBox(), local_domain[iboxA].first);
        }

      }

#ifdef _OPENMP
      for(unsigned i=0; i<threadBoxes.size(); ++i) {
        stk::search::add_to_box(boxA_proc.GetBox(), threadBoxes[i].GetBox());
      }
#endif

      std::vector<stk::search::ObjectBoundingBox_T<AxisAlignedBBType>> boxA_proc_box_array(num_procs);
      boxA_proc_box_array[current_proc] = boxA_proc;

      //
      //  Do a global communication to communicate all processor boxA bounding boxes
      //  to all processors in the group
      //
      GlobalBoxCombine(boxA_proc_box_array, comm);
      for(int iproc = 0; iproc < num_procs; ++iproc) {
        boxA_proc_box_array[iproc].set_object_number(iproc);
      }

      //
      //  Create a hierarchy of boxA processor bounding boxes.
      //  This hierarchy will be used to search for overlaps between processors and
      //  objects.
      //
      stk::search::ProximitySearchTree_T<AxisAlignedBBType> boxA_box_hierarchy(boxA_proc_box_array);

      //
      //  Determine what to ghost.  If a boxB box from this processor overlaps another processor's
      //  processor all-boxA box, then we need to ghost the boxB's data to that other processor.
      //  (The stricter criteria used by the full BoxA_BoxB_Ghost function would make sure that
      //  the individual boxB box overlaps some individual boxA box from the other processor.)
      //

      typedef typename RangeIdentifier::ident_type GlobalIdType;
      typedef std::pair<AxisAlignedBBType, GlobalIdType> BoxIdPair;
      std::vector<std::vector<BoxIdPair> > send_list(num_procs);
      std::vector<std::vector<BoxIdPair> > recv_list(num_procs);


      std::vector<int> proc_list(num_procs);



      for(unsigned int iboxB = 0; iboxB < numBoxRange; ++iboxB) {
        boxA_box_hierarchy.SearchForOverlap(local_range[iboxB].first, proc_list);
        for(unsigned i = 0; i < proc_list.size(); ++i) {
          int overlapping_proc = proc_list[i];
          if(overlapping_proc == current_proc) continue;
          GlobalIdType id = local_range[iboxB].second.id();
          send_list[overlapping_proc].push_back(BoxIdPair(rangeBoxes[iboxB], id));
        }
      }

      stk::parallel_data_exchange_t(send_list, recv_list, comm);
      rangeGhostIdentifiers.clear();
      for (size_t i = 0; i < recv_list.size(); i++) {
        for (size_t j = 0; j < recv_list[i].size(); j++) {
          const BoxIdPair& recvd_boxIdPair = recv_list[i][j];
          rangeBoxes.push_back(recvd_boxIdPair.first);
          rangeGhostIdentifiers.push_back(RangeIdentifier(recvd_boxIdPair.second, i));
        }
      }

    }

    template <typename DataType>
      inline void ConcatenateThreadLists(const std::vector<std::vector<DataType> > &vectorIn,
                                         std::vector<DataType> &vectorOut)
    {
      const unsigned numThreadLists = vectorIn.size();
      std::vector<unsigned> offsets(numThreadLists);
      unsigned totSize = 0;
      for (unsigned i = 0; i < numThreadLists; ++i) {
        offsets[i] = totSize;
        totSize += vectorIn[i].size();
      }

      vectorOut.resize(totSize);

#ifdef _OPENMP
#pragma omp parallel default(shared)
      {
        const unsigned ithread = omp_get_thread_num();
        const std::vector<DataType> &data = vectorIn[ithread];
        std::copy(data.begin(), data.end(), &vectorOut[offsets[ithread]]);
      }
#else
      for (unsigned ithread = 0; ithread < numThreadLists; ++ithread) {
        const std::vector<DataType> &data = vectorIn[ithread];
        std::copy(data.begin(), data.end(), &vectorOut[offsets[ithread]]);
      }
#endif

    }
  }


}

#endif /* GEOMETRY_TOOLKIT_INCLUDE_GEOM_COARSESEARCHMULTIALG_IMPL_H_ */
