// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef COMMON_SEARCH_UTIL_H_
#define COMMON_SEARCH_UTIL_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stk_util/environment/WallTime.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_search/KDTree_BoundingBox.hpp"
#include "stk_search/KDTree.hpp"
#include "mpi.h"

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


    //
    //  Exchange boxes so that current proc local box is sent to the global box on all processors.
    //
    template <typename DomainBox>
    inline void AllGatherHelper(const DomainBox& localBox, std::vector<DomainBox> &global_box_array, MPI_Comm &comm)
    {
      int numProc;
      MPI_Comm_size(comm, &numProc); 
      global_box_array.clear();
      global_box_array.resize(numProc);
      const char* localDataConst = reinterpret_cast<const char*>(&localBox);
      //  NKC, hack to support old MPI version used by Goodyear
      char* localData = const_cast<char*>(localDataConst);

      MPI_Allgather(localData, sizeof(DomainBox), MPI_CHAR, global_box_array.data(), sizeof(DomainBox), MPI_CHAR, comm);
    }



    template <typename ObjType, typename IdentifierType, typename BaseBoxType>
    void ParallelComputeProcObjectBoxes(const std::vector<std::pair<ObjType, IdentifierType> >& local_objsWithIdents,
                                        std::vector<stk::search::ObjectBoundingBox_T<BaseBoxType> > &objBB_proc_box_array,
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
      stk::search::AllGatherHelper(objBB_proc, objBB_proc_box_array, comm);

      int numProcs;
      MPI_Comm_size(comm, &numProcs);
      for(int iproc = 0; iproc < numProcs; ++iproc) {
        objBB_proc_box_array[iproc].set_object_number(iproc);
      }
    }


    // Ghost the range boxes needed for each processor to search against its local domain boxes
    // in distributed AABB overlap search ("coarse search")..
    template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
      void
      ComputeRangeWithGhostsForCoarseSearch(
                                            const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
                                            const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
                                            int num_procs,
                                            std::vector<RangeBoxType >& rangeBoxes, 
                                            std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm)
    {

      const unsigned numBoxRange  = local_range.size();

      using domainValueType = typename DomainObjType::value_type;
      using DomainBox       = stk::search::Box<domainValueType>;


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

template <typename DomainKey, typename RangeKey>
void communicateVector(
  stk::ParallelMachine arg_comm ,
  const std::vector< std::pair< DomainKey, RangeKey> > & send_relation ,
        std::vector< std::pair< DomainKey, RangeKey> > & recv_relation ,
        bool communicateRangeBoxInfo = true )
{
  typedef std::pair<DomainKey, RangeKey> ValueType ;

  CommSparse commSparse( arg_comm );

  const int p_rank = commSparse.parallel_rank();
  const int p_size = commSparse.parallel_size();

  typename std::vector< ValueType >::const_iterator i ; 

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) { 
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) == p_rank || ( communicateRangeBoxInfo && static_cast<int>(val.second.proc()) == p_rank) )
    {   
      recv_relation.push_back( val );
    }   
    if ( static_cast<int>(val.first.proc()) != p_rank ) { 
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.skip<ValueType>( 1 );
    }   
    if ( communicateRangeBoxInfo )
    {   
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) { 
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.skip<ValueType>( 1 );
        }
    }   
  }

  commSparse.allocate_buffers();

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) { 
    const ValueType & val = *i ;
    if ( static_cast<int>(val.first.proc()) != p_rank ) { 
      CommBuffer & buf = commSparse.send_buffer( val.first.proc() );
      buf.pack<ValueType>( val );
    }   
    if ( communicateRangeBoxInfo )
    {   
        if ( static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc() ) { 
          CommBuffer & buf = commSparse.send_buffer( val.second.proc() );
          buf.pack<ValueType>( val );
        }
    }   
  }

  commSparse.communicate();

  for ( int p = 0 ; p < p_size ; ++p ) { 
    CommBuffer & buf = commSparse.recv_buffer( p );
    while ( buf.remaining() ) { 
      ValueType val ;
      buf.unpack<ValueType>( val );
      recv_relation.push_back( val );
    }   
  }
}

  } // end namespace search
} // end namespace stk

#endif
