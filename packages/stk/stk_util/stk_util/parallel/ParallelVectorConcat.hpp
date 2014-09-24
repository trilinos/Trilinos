
// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#ifndef stk_util_parallel_ParallelVectorConcat_hpp
#define stk_util_parallel_ParallelVectorConcat_hpp

#include "stk_util/parallel/Parallel.hpp" 
#include <vector> 
#include "mpi.h"  

namespace stk {

  //------------------------------------------------------------------------
  //
  //  Take a list vector of T's on each processor.  Sum it into a single list that will be placed on all processors.
  //  The list contents and order will be guaranteed identical on every processor and formed by concatenation of the list
  //  fragments in processor order.
  //
  //  Return Code:  An MPI error code, MPI_SUCESS if correct
  //
  //  Example:
  //    Processor 1: localVec = {1, 2}
  //    Processor 2: localVec = {30, 40}
  //    Processor 3: localVec = {500, 600}
  //    Result on all processors: globalVec = {1, 2, 30, 40, 500, 600}
  // 
  template <typename T> int parallel_vector_concat(ParallelMachine comm, std::vector<T>& localVec, std::vector<T>& globalVec ) {
    globalVec.clear();
  
    const unsigned p_size = parallel_machine_size( comm );

    unsigned int sizeT     = sizeof(T);
    unsigned int localSize = sizeT * localVec.size();

    //
    //  Determine the total number of bytes being sent by each other processor.
    //
    std::vector<int> messageSizes(p_size);
    int mpiResult = MPI_SUCCESS ;  
    mpiResult = MPI_Allgather(&localSize, 1, MPI_INT, &messageSizes[0], 1, MPI_INT, comm);
    if(mpiResult != MPI_SUCCESS) {
      // Unknown failure, pass error code up the chain
      return mpiResult;
    }
    //
    //  Compute the offsets into the resultant array
    //
    std::vector<int> offsets(p_size+1);
    offsets[0] = 0;
    for(unsigned iproc=1; iproc<p_size+1; ++iproc) {
      offsets[iproc] = offsets[iproc-1] + messageSizes[iproc-1];
    }
    unsigned int totalSize = (offsets[p_size])/sizeT;
    globalVec.resize(totalSize);

    //
    //  Do the all gather to copy the actual array data and propogate to all processors
    //
    mpiResult = MPI_Allgatherv(&localVec[0], localSize, MPI_CHAR, &globalVec[0], &messageSizes[0], &offsets[0], MPI_CHAR, comm);
    return mpiResult;
  }

}

#endif

