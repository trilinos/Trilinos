
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

#ifndef stk_util_parallel_ParallelIndexGapFinder_hpp
#define stk_util_parallel_ParallelIndexGapFinder_hpp

#include "stk_util/parallel/Parallel.hpp" 
#include <vector> 
#include <algorithm> 
#include "mpi.h"  
#include <assert.h>
#include "stk_util/diag/String.hpp"
#include <stdexcept>            
#include <math.h>
#include <stdint.h>


namespace stk {


  
  //------------------------------------------------------------------------
  //
  //  Take a list vector of 'in use' indexes on each processor.  Purpose of this routine is to find a 
  //  set of indexes of the given size that are not used by any current process
  //
  //  Arguments:
  //    minId              : input  : Minimum available id, all returnIds will be greater than or equal to this
  //    maxId              : input  : Maximum number of available id, all returnIds will be less than this number
  //    localIdsInUse      : input  : List of ids currently in use by this processor.  Note, overlap 
  //                                  of this array between processors is allowed, but optimal performance
  //                                  will occur with no overlap
  //    globalNumIdsNeeded : input  : Num of ids required ids on all processors
  //    returnIds          : output : Global consistent list of available ids generated
  //
  //  Returns:  An MPI error code, 0 if correct
  //
  int parallel_index_gap_finder_global(ParallelMachine comm, const uint64_t minId, 
                                       const uint64_t maxId, 
                                       const std::vector<uint64_t>& localIdsInUse, 
                                       const uint64_t localNumIdsNeeded,            
                                       const uint64_t globalNumIdsNeeded,
                                       const uint64_t prefixSumOfLocalIdsNeeded,
                                       std::vector<uint64_t>& returnIds );
}

#endif

