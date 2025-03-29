
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
    const unsigned p_size = parallel_machine_size( comm );

    //  Check for serial simplified early out condition
    if(p_size == 1) {
      globalVec = localVec;
      return MPI_SUCCESS;
    }

    globalVec.clear();

    STK_ThrowRequireMsg(localVec.size() <= std::numeric_limits<int>::max(), "input vector length must fit in an int");
    int localSize = localVec.size();

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

    STK_ThrowRequireMsg(totalSize <= size_t(std::numeric_limits<int>::max()), "input vector length must fit in an int");
    globalVec.resize(totalSize);

    //
    //  Compute the offsets into the resultant array
    //
    std::vector<int> offsets(p_size);
    offsets[0] = 0;
    for(unsigned iproc=1; iproc<p_size; ++iproc) {
      offsets[iproc] = offsets[iproc-1] + messageSizes[iproc-1];
    }

    //
    //  Do the all gather to copy the actual array data and propogate to all processors
    //  Note, localVec should not be modified by the MPI call, but MPI does not guarntee the const in the 
    //  interface argument.
    //
    T* ptrNonConst = const_cast<T*>(localVec.data());
    MPI_Datatype datatype = stk::generate_mpi_datatype<T>();

    mpiResult = MPI_Allgatherv(ptrNonConst, localSize, datatype, globalVec.data(), messageSizes.data(), offsets.data(), datatype, comm);
    return mpiResult;
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
    globalVec = localVec;
    return 0;
}
#endif

}

#endif

