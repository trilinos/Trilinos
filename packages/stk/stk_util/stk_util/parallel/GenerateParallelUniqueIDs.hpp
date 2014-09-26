
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
// 

#ifndef stk_util_parallel_GenerateParallelUniqueIDs_hpp
#define stk_util_parallel_GenerateParallelUniqueIDs_hpp

#include "stk_util/parallel/ParallelVectorConcat.hpp" 
#include "stk_util/parallel/Parallel.hpp" 
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
  //
  //  Note this if this routine is called on different parallel decompositions it may return different id sets both in 
  //  order and in content.
  //
  //  Arguments:
  //  existingIds    : input  : List of ids already in use.  Note, the id list may be unique to a processor 
  //                            or different processor lists may overlap
  //  numNewIdsLocal : input  : Number of new ids requested on this processor
  //  maxAllowedId   : input  : Max id value allowed to be returned.  If a valid set of ids cannot be found
  //                            below maxAllowedId routine will return a failure code.
  //
  //  Returns:  Vector of ids created for the local processor.  Will generally throw if encountering any internal 
  //            errors.  Note, the error throws are NOT guaranteed to be parallel consistent and are usually unrecoverable.
  //
  //  Assumptions and Usage Guidelines
  //    This algorithm assumes current existingIds are relatively dense packed.  The specific requirement for
  //    success is 'maxAllowedId - global_max(existingIds) > orderArray.size()'
  //

  std::vector<unsigned> generate_parallel_unique_ids(const unsigned maxAllowedId,
                                                     const std::vector<unsigned>& existingIds,
                                                     unsigned numNewIdsLocal,
                                                     ParallelMachine comm) {
    std::vector<unsigned> newIds;
    newIds.reserve(numNewIdsLocal);
    //
    //  Extract global max existing id.  For basic use case just start generating ids starting
    //  at the previous max_id + 1.  
    //
    unsigned localMaxId = 0;
    unsigned globalMaxId;    
    for(unsigned i=0; i<existingIds.size(); ++i) {
      if(existingIds[i] > localMaxId) {
        localMaxId = existingIds[i];
      }
    }

    int mpiResult = MPI_SUCCESS ;  
    mpiResult = MPI_Allreduce(&localMaxId, &globalMaxId, 1, MPI_UNSIGNED, MPI_MAX, comm);
    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
      return newIds;
    }
    //
    //  Compute the total number of ids requested
    //
    unsigned globalNumIdsRequested;
    mpiResult = MPI_Allreduce(&numNewIdsLocal, &globalNumIdsRequested, 1, MPI_UNSIGNED, MPI_SUM, comm);
    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
      return newIds;
    }
    if(globalNumIdsRequested == 0) {
      return newIds;
    }


    unsigned myFirstNewId;
    mpiResult = MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1, MPI_UNSIGNED, MPI_SUM, comm);


    //
    //  Check if sufficent extra ids exist to run this algorithm
    //    NKC, for now fail if not, later maybe switch to a more robust but more
    //    expensive algorithm.
    //  Note, below computations organized to avoid potential overflow
    //
    unsigned availableIds = maxAllowedId - globalMaxId;
    if(availableIds < globalNumIdsRequested) {
      std::ostringstream msg;
      msg << "In generate_parallel_unique_ids, insufficent ids available \n";
      msg << "  Max current id:       "<<globalMaxId<<"\n";
      msg << "  Max allowed id:       "<<maxAllowedId<<"\n";
      msg << "  Number ids requested: "<<globalNumIdsRequested<<"\n";
      throw std::runtime_error(msg.str());
      return newIds;
    }
    //
    //  Create new ids starting at 'globalMaxId+1' and ending at 'globalMaxId+globalNumIdsRequested'
    //  Compute the starting offset for ids in each processor.
    //
    myFirstNewId = (myFirstNewId - numNewIdsLocal) + globalMaxId + 1;
    for(unsigned i=0; i<numNewIdsLocal; ++i) {
      newIds.push_back(myFirstNewId+i);
    }
    return newIds;
  }

}

#endif

