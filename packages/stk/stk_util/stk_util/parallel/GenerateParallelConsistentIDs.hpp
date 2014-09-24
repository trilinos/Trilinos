
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

#ifndef stk_util_parallel_GenerateParallelConsistentIDs_hpp
#define stk_util_parallel_GenerateParallelConsistentIDs_hpp

#include "stk_util/parallel/ParallelVectorConcat.hpp" 
#include "stk_util/parallel/Parallel.hpp" 
#include <vector> 
#include "mpi.h"  

namespace stk {

  //--------------------------------------------------------------------------------------------
  //
  //  Given an input set of existing ids, generate the requested number of new ids.  Subject to the following
  //  requirements:
  //
  //  1) The new ids will not overlap any existing ids, this will hold globally in parallel
  //  2) All new ids will be less than or equal to the provided maxId
  //  3) Provided the 'orderArray' given is always the same irrespective of parallel decomposition
  //     the numerical values of the returned as ids should always be the same and the same id 
  //     should be returned for the same orderArray value.
  //
  //  Arguments:
  //  existingIds : input  : List of ids already in use.  Note, the id list may be unique to a processor 
  //                         or different processor lists may overlap
  //  orderArray  : input  : Defines a unique ordering number of objects.  Node OrderType objects must be 
  //                         sortable and have a '<' operator defined.  The length of the order array also
  //                         defines the number of new objects needed on the current processor.
  //  newIds      : output : Returns orderArray.size() new ids that match one to one with the passed order
  //                         keys. 
  //  maxId       : input  : Max id value allowed to be returned.  If a valid set of ids cannot be found
  //                         below maxId routine will return a failure code.
  //
  //  Returns:  Status enum indicating success or failure, note, some sub lower level functions could 
  //            potentially also throw on other fatal errors.
  //     
  //
  //

  void generate_parallel_consistent_ids(ParallelMachine comm) {


  }

}

#endif

