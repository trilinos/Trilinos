
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

#ifndef stk_util_parallel_ParallelVectorConcat_hpp
#define stk_util_parallel_ParallelVectorConcat_hpp

#include "stk_util/parallel/Parallel.hpp" 
#include "stk_util/parallel/MPIDatatypeGenerator.hpp"
#include <limits>
#include <vector>
#include <type_traits>
#include "stk_util/util/ReportHandler.hpp"

namespace stk {

#if defined( STK_HAS_MPI )
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
  //  Usage Guidelines:
  //    Generally type T must be a plain data type with no pointers or allocated memory.  
  //    Thus T could be standard types such as int or double or structs or classes that contain only ints
  //    and double (such as mtk::Vec3<double>).
  //    T should NOT be a general structure that contains pointers, strings, or vectors as these structures cannot be
  //    properly transfered between processors.  
  //    A handful of specializations are available to handle certain more complex types
  //
  //  Specializations for non-PODs.
  //    std::string
  //
  template <typename T> inline int parallel_vector_concat(ParallelMachine comm, const std::vector<T>& localVec, std::vector<T>& globalVec )
  {
    static_assert(std::is_trivially_destructible_v<T>,
      "parallel_vector_concat's generic template byte-copies T across MPI; a type that "
      "owns heap memory (e.g. std::string, sierra::String, std::vector) is not trivially "
      "destructible and needs an explicit specialization whose declaring header is included "
      "at the call site.");
    const unsigned p_size = parallel_machine_size( comm );

    //  Check for serial simplified early out condition
    if(p_size == 1) {
      globalVec = localVec;
      return MPI_SUCCESS;
    }

    globalVec.clear();

    STK_ThrowRequireMsg(localVec.size() <= std::numeric_limits<int>::max(), "input vector length must fit in an int");
    int localSize = localVec.size();
    const int my_rank = parallel_machine_rank( comm );

    //
    //  Determine the total number of bytes being sent by each other processor.
    //
    std::vector<int> messageSizes(p_size);
    int mpiResult = MPI_SUCCESS ;
    mpiResult = MPI_Allgather(&localSize, 1, MPI_INT, messageSizes.data(), 1, MPI_INT, comm);
    if(mpiResult != MPI_SUCCESS) {
      // Unknown failure, pass error code up the chain
      return mpiResult;
    }

    size_t totalSize = 0;
    for (auto& size : messageSizes)
    {
      totalSize += size;
    }

    globalVec.resize(totalSize);

    //
    //  Compute the offsets into the resultant array.  Use size_t here since the total
    //  concatenated size may exceed the range of an int (the per-processor size is
    //  guaranteed to fit in an int by the check above).
    //
    std::vector<size_t> globalOffsets(p_size);
    globalOffsets[0] = 0;
    for(unsigned iproc=1; iproc<p_size; ++iproc) {
      globalOffsets[iproc] = globalOffsets[iproc-1] + messageSizes[iproc-1];
    }

    //
    //  Do the all gather to copy the actual array data and propogate to all processors
    //  Note, localVec should not be modified by the MPI call, but MPI does not guarntee the const in the
    //  interface argument.
    //
    T* ptrNonConst = const_cast<T*>(localVec.data());
    MPI_Datatype datatype = stk::generate_mpi_datatype<T>();

    //
    //  MPI_Allgatherv expresses its receive counts and displacements as ints, so the
    //  cumulative displacement into the receive buffer must fit in an int.  When the total
    //  concatenated size exceeds the 32-bit int limit, split the transfer into multiple
    //  MPI_Allgatherv calls over contiguous groups of processors, where each group's
    //  combined size fits in an int.  For a given call, processors outside the current
    //  group send and receive nothing.  Since messageSizes is identical on every processor,
    //  the group decomposition (and thus the sequence of collective calls) is identical on
    //  all processors.
    //
    const size_t maxCountPerCall = size_t(std::numeric_limits<int>::max());
    std::vector<int> groupCounts(p_size);
    std::vector<int> groupDispls(p_size);

    unsigned groupStart = 0;
    while(groupStart < p_size) {
      //  Grow the group while its combined size stays within the int limit.  Each single
      //  processor's size is <= int max (checked above), so every group holds at least one
      //  processor and the loop always makes progress.
      unsigned groupEnd = groupStart;
      size_t groupSize = 0;
      while(groupEnd < p_size && groupSize + size_t(messageSizes[groupEnd]) <= maxCountPerCall) {
        groupSize += size_t(messageSizes[groupEnd]);
        ++groupEnd;
      }

      const size_t groupBaseOffset = globalOffsets[groupStart];
      for(unsigned iproc=0; iproc<p_size; ++iproc) {
        const bool inGroup = (iproc >= groupStart && iproc < groupEnd);
        groupCounts[iproc] = inGroup ? messageSizes[iproc] : 0;
        groupDispls[iproc] = inGroup ? int(globalOffsets[iproc] - groupBaseOffset) : 0;
      }

      const bool selfInGroup = (my_rank >= int(groupStart) && my_rank < int(groupEnd));
      const int sendCount = selfInGroup ? localSize : 0;

      mpiResult = MPI_Allgatherv(ptrNonConst, sendCount, datatype,
                                 globalVec.data() + groupBaseOffset,
                                 groupCounts.data(), groupDispls.data(), datatype, comm);
      if(mpiResult != MPI_SUCCESS) {
        // Unknown failure, pass error code up the chain
        return mpiResult;
      }

      groupStart = groupEnd;
    }

    return MPI_SUCCESS;
  }


  //------------------------------------------------------------------------
  //
  //  std::string specializations for parallel_vector_concat.  As strings are not PODs they need special handling
  //  to concat correctly.  String concatentaion is a common use case particularly for generating
  //  parallel consistent error messages.
  //
  template<>
  inline int parallel_vector_concat(ParallelMachine comm, const std::vector<std::string>& localList, std::vector<std::string>& globalList ) {
    //
    //  Convert the local vector of strings into a single vector of null seperated char bits 
    //  so that standardized list concact can be used.  
    //
    std::vector<char> charLocalList;
    for(unsigned istring=0; istring<localList.size(); ++istring) {
      unsigned numChar = localList[istring].size();
      const char* str = localList[istring].c_str();
      for(unsigned ichar=0; ichar<numChar; ++ichar) {
        charLocalList.push_back(str[ichar]);
      }
      charLocalList.push_back(0);
    }
    //
    //  Parallel concat the character lists
    //
    std::vector<char> charGlobalList;
    int mpiResult = stk::parallel_vector_concat<char>(comm, charLocalList, charGlobalList);
    if(mpiResult != MPI_SUCCESS) {
      // Unknown failure, pass error code up the chain
      return mpiResult;
    }
    //
    //  Turn the character arrays back into strings for output
    //
    unsigned curCharIndex = 0;
    unsigned charGlobalListLen = charGlobalList.size();
    std::vector<char> nextString;
    while(curCharIndex < charGlobalListLen) {
      char curChar = charGlobalList[curCharIndex];
      nextString.push_back(curChar);
      if(curChar == 0) {
        globalList.emplace_back(nextString.data());
        nextString.clear();
      }
      curCharIndex++;
    }
    return MPI_SUCCESS;
  }

  //------------------------------------------------------------------------
  //
  //  bool specialization for parallel_vector_concat.
  //
  template<>
  inline int parallel_vector_concat(ParallelMachine comm, const std::vector<bool>& localVec, std::vector<bool>& globalVec ) {
    //it turns out that std::vector<bool> is a weird beast, it doesn't have a .data() method.
    //In general, its contents can't be treated like a 'bool*'.
    //Thus the best approach here is to copy to a vector of chars and
    //call the general implementation of parallel_vector_concat.

    std::vector<unsigned char> localChars(localVec.size());
    for(unsigned i=0; i<localVec.size(); ++i) {
      localChars[i] = localVec[i] ? 1 : 0;
    }

    std::vector<unsigned char> globalChars;

    int returnValue = parallel_vector_concat(comm, localChars, globalChars);

    globalVec.resize(globalChars.size());
    for(unsigned i=0; i<globalVec.size(); ++i) {
      globalVec[i] = globalChars[i] == 1 ? true : false;
    }

    return returnValue;
  }

#else
  template <typename T> inline int parallel_vector_concat(ParallelMachine comm, const std::vector<T>& localVec, std::vector<T>& globalVec ) {
    static_assert(std::is_trivially_destructible_v<T>,
      "parallel_vector_concat's generic template byte-copies T across MPI; a type that "
      "owns heap memory (e.g. std::string, sierra::String, std::vector) is not trivially "
      "destructible and needs an explicit specialization whose declaring header is included "
      "at the call site.");
    globalVec = localVec;
    return 0;
}
#endif

}

#endif

