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

#include "stk_util/parallel/GenerateParallelUniqueIDs.hpp"
#include "stk_util/parallel/MPI.hpp"                     // for Datatype
#include "stk_util/parallel/ParallelIndexGapFinder.hpp"  // for parallel_index_gap_finder_global
#include "stk_util/parallel/ParallelReduce.hpp"          // for all_reduce_max
#include "stk_util/util/ReportHandler.hpp"               // for ThrowRequireMsg
#include <algorithm>                                     // for max, max_element
#include <cstdint>                                       // for uint64_t
#include <sstream>                                       // for operator<<, ostringstream, basic...
#include <stdexcept>                                     // for runtime_error

namespace stk {

  //
  //  MPI operation to compute the sum and max of a value simultaneously
  //
  void MpiSumMax(uint64_t* in, uint64_t* inout, int* len, MPI_Datatype * /*dptr*/) {
    if(*len != 2) {
      std::ostringstream msg;
      msg << "MpiSumMax MPI operation can only be given two arguements \n";
      throw std::runtime_error(msg.str());
    }
    inout[0] += in[0];
    inout[1] = std::max(inout[1], in[1]);
  }

void compute_global_sum_and_max(ParallelMachine comm,
                                uint64_t numLocalIds,
                                uint64_t& globalSumNumIds,
                                uint64_t& globalMaxNumIds)
{
  uint64_t numLocalReduce[2] = {numLocalIds, numLocalIds};

  MPI_Op myOp;
  MPI_Op_create((MPI_User_function *)MpiSumMax, true, &myOp);

  uint64_t numGlobalReduce[2] = {0u, 0u};
  STK_ThrowRequireMsg(MPI_SUCCESS == MPI_Allreduce(numLocalReduce, numGlobalReduce, 2, sierra::MPI::Datatype<uint64_t>::type(), myOp, comm), "MPI_Allreduce failed");
  MPI_Op_free(&myOp);

  globalSumNumIds = numGlobalReduce[0];
  globalMaxNumIds = numGlobalReduce[1];
}

void generate_parallel_ids_in_gap(ParallelMachine comm,
                                  const std::vector<uint64_t>& existingIds,
                                  uint64_t maxAllowedId,
                                  uint64_t numNewIdsLocal,
                                  uint64_t globalNumIdsRequested,
                                  std::vector<uint64_t>& newIds)
{
  //
  // Run fallback sparse fill algorithm
  //

  newIds.clear();
  newIds.reserve(numNewIdsLocal);
  uint64_t myFirstNewId;
  STK_ThrowRequireMsg(MPI_SUCCESS == MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, comm), "MPI_Scan failed");
  myFirstNewId -= numNewIdsLocal;

  STK_ThrowRequireMsg(0 == parallel_index_gap_finder_global(comm, 0, maxAllowedId,
                                                    existingIds, numNewIdsLocal,
                                                    globalNumIdsRequested,
                                                    myFirstNewId, newIds),
                  "Id allocation failure beneath generate_parallel_ids_in_gap.");
}

void generate_parallel_ids_above_existing_max(ParallelMachine comm,
                                              uint64_t numNewIdsLocal,
                                              uint64_t globalNumIdsRequested,
                                              uint64_t maxIdsRequested,
                                              uint64_t availableIds,
                                              uint64_t globalMaxId,
                                              std::vector<uint64_t>& newIds)
{
  newIds.clear();
  newIds.reserve(numNewIdsLocal);
  int mpi_size = stk::parallel_machine_size(comm);
  uint64_t wastfullFillSize = maxIdsRequested*mpi_size;
  if(availableIds >= wastfullFillSize && wastfullFillSize/4 < globalNumIdsRequested) {
    //
    //  Check if can use the super-cheap algorithm, if space wastage less than 400%, define an equal number of ids per processor
    //  and avoid the scan
    //
    int mpi_rank = stk::parallel_machine_rank(comm);

    uint64_t myFirstNewId = globalMaxId + mpi_rank*maxIdsRequested + 1;
    for(uint64_t i=0; i<numNewIdsLocal; ++i) {
      newIds.push_back(myFirstNewId+i);
    }
  } else {
    //
    //  Otherwise still use a very cheap algorithm, densely pack the resulting ids
    //
    uint64_t myFirstNewId;
    STK_ThrowRequireMsg(MPI_SUCCESS == MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, comm), "MPI_Scan failed");
    myFirstNewId -= numNewIdsLocal;
    //
    //  Run basic cheap algorithm.  Start counting new ids at end.
    //  Create new ids starting at 'globalMaxId+1' and ending at 'globalMaxId+globalNumIdsRequested'
    //  Compute the starting offset for ids in each processor.
    //
    myFirstNewId = (myFirstNewId) + globalMaxId + 1;
    for(uint64_t i=0; i<numNewIdsLocal; ++i) {
      newIds.push_back(myFirstNewId+i);
    }
  }
}

std::vector<uint64_t> generate_parallel_unique_ids(const uint64_t maxAllowedId,
                                                   const std::vector<uint64_t>& existingIds,
                                                   uint64_t numNewIdsLocal,
                                                   ParallelMachine comm,
                                                   uint64_t localMaxExistingId)
{
  std::vector<uint64_t> newIds;
  //
  //  Extract global max existing id.  For basic use case just start generating ids starting
  //  at the previous max_id + 1.
  //
  uint64_t globalMaxId     = 0;
  uint64_t maxIdsRequested = 0;
  uint64_t localMaxId = localMaxExistingId;
  if (localMaxId == 0) {
    localMaxId = existingIds.empty() ? 0 : *std::max_element(existingIds.begin(), existingIds.end());
  }

  const unsigned numMax = 2;
  uint64_t locals[numMax] = {localMaxId, numNewIdsLocal};
  uint64_t globals[numMax] = {0, 0};
  stk::all_reduce_max(comm, locals, globals, numMax);
  globalMaxId = globals[0];
  maxIdsRequested = globals[1];

  uint64_t globalNumIdsRequested = stk::get_global_sum(comm, numNewIdsLocal);

  if(globalNumIdsRequested == 0) {
    return newIds;
  }

  //
  //  Check if sufficent extra ids exist to run this algorithm
  //
  //  Note, below computations organized to avoid potential overflow
  //
  uint64_t availableIds = maxAllowedId - globalMaxId;
  if(availableIds < globalNumIdsRequested) {
    generate_parallel_ids_in_gap(comm, existingIds, maxAllowedId, numNewIdsLocal,
                                 globalNumIdsRequested, newIds);
  }
  else {
    generate_parallel_ids_above_existing_max(comm, numNewIdsLocal,
                                      globalNumIdsRequested, maxIdsRequested,
                                      availableIds, globalMaxId, newIds);
  }
  return newIds;
}

}

