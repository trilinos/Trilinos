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

#include <stk_util/parallel/GenerateParallelUniqueIDs.hpp>
#include <stk_util/parallel/MPI.hpp>

namespace stk {

  //
  //  MPI operation to compute the sum and max of a value simultaneously
  //
  void MpiSumMax(uint64_t* in, uint64_t* inout, int* len, MPI_Datatype * dptr) {
    if(*len != 2) {
      std::ostringstream msg;
      msg << "MpiSumMax MPI operation can only be given two arguements \n";
      throw std::runtime_error(msg.str());
    }
    inout[0] += in[0];
    inout[1] = std::max(inout[1], in[1]);
  }


  std::vector<uint64_t> generate_parallel_unique_ids(const uint64_t maxAllowedId,
                                                     const std::vector<uint64_t>& existingIds,
                                                     uint64_t numNewIdsLocal,
                                                     ParallelMachine comm) {

    int mpi_rank = stk::parallel_machine_rank(comm);

    std::vector<uint64_t> newIds;
    newIds.reserve(numNewIdsLocal);
    //
    //  Extract global max existing id.  For basic use case just start generating ids starting
    //  at the previous max_id + 1.
    //
    uint64_t localMaxId = 0;
    uint64_t globalMaxId;
    for(uint64_t i=0; i<existingIds.size(); ++i) {
      if(existingIds[i] > localMaxId) {
        localMaxId = existingIds[i];
      }
    }

    int mpiResult = MPI_SUCCESS ;
    mpiResult = MPI_Allreduce(&localMaxId, &globalMaxId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_MAX, comm);
    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
    }
    //
    //  Compute the total number of ids requested
    //
    uint64_t numLocalReduce[2];
    numLocalReduce[0] = numNewIdsLocal;
    numLocalReduce[1] = numNewIdsLocal;

    MPI_Op myOp;
    MPI_Op_create((MPI_User_function *)MpiSumMax, true, &myOp);

    uint64_t numGlobalReduce[2] = {0u, 0u};
    mpiResult = MPI_Allreduce(numLocalReduce, numGlobalReduce, 2, sierra::MPI::Datatype<uint64_t>::type(), myOp, comm);
    MPI_Op_free(&myOp);

    uint64_t globalNumIdsRequested = numGlobalReduce[0];
    uint64_t maxIdsRequested       = numGlobalReduce[1];

    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
    }
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
      //
      // Run fallback sparse fill algorithm
      //

      uint64_t myFirstNewId;
      mpiResult = MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, comm);
      myFirstNewId -= numNewIdsLocal;

      int returnCode = parallel_index_gap_finder_global(comm, 0, maxAllowedId,
                                                        existingIds, numNewIdsLocal, globalNumIdsRequested,
                                                        myFirstNewId, newIds);
      if(returnCode != 0) {
        std::ostringstream msg;
        msg << "In generate_parallel_unique_ids, failure in id allocation \n";
        throw std::runtime_error(msg.str());
      }
      return newIds;
    } else {
      int mpi_size = stk::parallel_machine_size(comm);
      uint64_t wastfullFillSize = maxIdsRequested*mpi_size;
      if(availableIds >= wastfullFillSize && wastfullFillSize/4 < globalNumIdsRequested) {
        //
        //  Check if can use the super-cheap algorithm, if space wastage less than 400%, define an equal number of ids per processor
        //  and avoid the scan
        //
        uint64_t myFirstNewId = globalMaxId + mpi_rank*maxIdsRequested + 1;
        for(uint64_t i=0; i<numNewIdsLocal; ++i) {
          newIds.push_back(myFirstNewId+i);
        }
      } else {
        //
        //  Otherwise still use a very cheap algorithm, densely pack the reultant ids
        //
        uint64_t myFirstNewId;
        mpiResult = MPI_Scan(&numNewIdsLocal, &myFirstNewId, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, comm);
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
    return newIds;
  }

}

