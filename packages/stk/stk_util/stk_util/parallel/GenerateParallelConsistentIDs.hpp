
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
// 

#ifndef stk_util_parallel_GenerateParallelConsistentIDs_hpp
#define stk_util_parallel_GenerateParallelConsistentIDs_hpp

#include "stk_util/parallel/ParallelVectorConcat.hpp"
#include "stk_util/parallel/MPI.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelIndexGapFinder.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>                       // for string
#include <sstream>
#include "mpi.h"
#include <assert.h>

namespace stk {

  //--------------------------------------------------------------------------------------------
  //
  //  Given an input set of existing ids, generate the requested number of new ids.  Subject to the following
  //  requirements:
  //
  //  1) The new ids will not overlap any existing ids, this will hold globally in parallel
  //  2) All new ids will be less than or equal to the provided maxAllowedId
  //  3) Provided the 'orderArray' given is always the same irrespective of parallel decomposition
  //     the numerical values of the returned as ids should always be the same and the same id 
  //     should be returned for the same orderArray value irrespective of parallel decomposition.
  //
  //  Arguments:
  //  existingIds  : input  : List of ids already in use.  Note, the id list may be unique to a processor 
  //                          or different processor lists may overlap
  //  orderArray   : input  : Defines a unique ordering number of objects.  Node OrderType objects must be 
  //                          sortable and have a '<' operator defined.  Additionally orderArray objects
  //                          should be 'plain old data' without any pointers or allocated memory.  Finally 
  //                          the orderArray should contain unique keys only.
  //
  //                          The length of the order array also defines the number of ids needed
  //                          on the current processor.
  //  maxAllowedId : input  : Max id value allowed to be returned.  If a valid set of ids cannot be found
  //                          below maxAllowedId routine will not return any ids and throw a failure.
  //
  //  Returns:  Array of new ids created.  Returned array will have the same size as the input order array.  
  //            Will generally throw if encountering any internal errors.  Note, the error
  //            throws are NOT guaranteed to be parallel consistent and usually are not recoverable.
  //
  //  Assumptions and Usage Guidelines
  //    This algorithm assumes current existingIds are relatively dense packed.  The specific requirement for
  //    success is 'maxAllowedId - global_max(existingIds) > global_sum(localOrderArray.size())'
  //


  template <typename OrderType> 
  std::vector<uint64_t> generate_parallel_consistent_ids(const uint64_t maxAllowedId,
                                                         const std::vector<uint64_t>& existingIds,
                                                         const std::vector<OrderType>& localOrderArray,
                                                         const ParallelMachine comm) {

    std::vector<uint64_t> newIds;
    //
    //  Extract global max existing id.  For basic use case just start generating ids starting
    //  at the previous max_id + 1.  

    uint64_t localMaxId = 0;
    uint64_t globalMaxId;    
    for(uint64_t i=0; i<existingIds.size(); ++i) {
      if(existingIds[i] > localMaxId) {
        localMaxId = existingIds[i];
      }
    }

    int mpiResult = MPI_Allreduce(&localMaxId, &globalMaxId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_MAX, comm);
    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
    }
    //
    //  Communicate the entire ordering function to every single other processor.
    //    Expensive - YES
    //    Nessecary - YES
    //    
    std::vector<OrderType> globalOrderArray;
    stk::parallel_vector_concat(comm, localOrderArray, globalOrderArray);

    uint64_t globalNumIdsRequested = globalOrderArray.size();
    if(globalNumIdsRequested == 0) return newIds;

    //
    //  Sort the ordering array and confirm it has no duplicated keys
    //
    std::sort(globalOrderArray.begin(), globalOrderArray.end());
    for(uint64_t i=0; i<globalNumIdsRequested-1; ++i) {
      if(globalOrderArray[i] == globalOrderArray[i+1]) {
        std::ostringstream msg;
        msg << "In generate_parallel_consistent_ids, an ordering array with non unique keys provided.\nThis is not allowed.\n";
        throw std::runtime_error(msg.str());
      }
    }


    //
    //  Check if sufficent extra ids exist to run this algorithm
    //    NKC, for now fail if not, later maybe switch to a more robust but more
    //    expensive algorithm.
    //  Note, below computations organized to avoid potential overflow
    //
    newIds.reserve(localOrderArray.size());
    uint64_t availableIds = maxAllowedId - globalMaxId;
    if(availableIds < globalNumIdsRequested) {
      //
      // Run fallback fill algorithm
      //
      uint64_t numNewIdsLocal = localOrderArray.size();
      uint64_t myFirstNewId;
      mpiResult = MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1,  sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, comm);
      myFirstNewId -= numNewIdsLocal;
      std::vector<uint64_t> allowedIdsLocal;
      std::vector<uint64_t> allowedIdsGlobal;
      int returnCode = parallel_index_gap_finder_global(comm, 0, maxAllowedId, 
                                       existingIds, numNewIdsLocal, globalNumIdsRequested, 
                                       myFirstNewId, allowedIdsLocal);
      if(returnCode != 0) {
        std::ostringstream msg;
        msg << "In generate_parallel_consistent_ids, failure in id allocation \n";
        throw std::runtime_error(msg.str());
      }

      parallel_vector_concat(comm, allowedIdsLocal, allowedIdsGlobal);
     
      for(uint64_t i=0, n=localOrderArray.size(); i<n; ++i) {
        typename std::vector<OrderType>::iterator index = std::lower_bound(globalOrderArray.begin(), globalOrderArray.end(), localOrderArray[i]);
        assert(index != globalOrderArray.end());
        int offset = std::distance(globalOrderArray.begin(), index);
        newIds.push_back(allowedIdsGlobal[offset]);
      }

    } else {

      for(uint64_t i=0, n=localOrderArray.size(); i<n; ++i) {
        typename std::vector<OrderType>::iterator index = std::lower_bound(globalOrderArray.begin(), globalOrderArray.end(), localOrderArray[i]);
        assert(index != globalOrderArray.end());
        int offset = std::distance(globalOrderArray.begin(), index);
        newIds.push_back(globalMaxId+1+offset);
      }
    }

    return newIds;
  }

}

#endif

