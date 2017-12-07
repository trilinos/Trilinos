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

#ifndef THREAD_SORT_H
#define THREAD_SORT_H

#include <vector>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
 #include <omp.h>
#if !defined(__INTEL_COMPILER) && !defined(SALINAS_LINUX)
#if defined(GCC_VERSION) && GCC_VERSION >= 500
  #define USE_NATIVE_THREADED_SORT 1
  #include <parallel/algorithm>
#endif
#endif
#endif


namespace stk {
namespace search {

#ifdef _OPENMP
  template<typename T>
  void threadedMergeSort(T* inputVec, T* scratchVec, unsigned vecSize, const unsigned numThreadsToUse) {
    //
    //  Reached the serial limit, just sort it
    //
    if(numThreadsToUse == 1) {
      std::sort(inputVec, inputVec+vecSize);
      return;
    }

    unsigned rightNumThreads = numThreadsToUse/2;
    unsigned leftNumThreads  = numThreadsToUse-rightNumThreads;
    unsigned leftSize        = vecSize/2;
    unsigned rightSize       = vecSize - leftSize;
    T*  leftList         = inputVec;
    T*  rightList        = inputVec+leftSize;
    T*  leftScratch      = scratchVec;
    T*  rightScratch     = scratchVec+leftSize;

#pragma omp parallel sections num_threads(2) default(none), shared(leftNumThreads, rightNumThreads, leftSize, rightSize, leftScratch, rightScratch, leftList, rightList)
    {
#pragma omp section
      {
        threadedMergeSort(leftList,  leftScratch,  leftSize,  leftNumThreads );
      }
#pragma omp section
      {
        threadedMergeSort(rightList, rightScratch, rightSize, rightNumThreads);
      }
    }

    //
    //  Merge operation, not uses just two threads irrespective of the number of available threads, could prove to be some bottle neck
    //  when a large number of threads are available this is a bidirection meet at the middle merge with no scan overhead.  
    //


#pragma omp parallel sections num_threads(2) default(none), shared(leftSize, rightSize, scratchVec, rightScratch, leftList, rightList, vecSize)
    {
#pragma omp section
      {
        unsigned leftSubListPos  = 0;
        unsigned rightSubListPos = 0;
        unsigned mergedListPos = 0;       
        for(unsigned ilist=0; ilist<leftSize; ++ilist) {
          if(leftList[leftSubListPos] < rightList[rightSubListPos]) {
            scratchVec[mergedListPos++] = leftList[leftSubListPos++];
          } else {
            scratchVec[mergedListPos++] = rightList[rightSubListPos++];
          }
        }
      }
#pragma omp section
      {
        unsigned leftSubListPos  = leftSize-1;
        unsigned rightSubListPos = rightSize-1;
        unsigned mergedListPos = vecSize-1;       
        for(unsigned ilist=0; ilist<rightSize-1; ++ilist) {
          if(leftList[leftSubListPos] > rightList[rightSubListPos]) {
            scratchVec[mergedListPos--] = leftList[leftSubListPos--];
          } else {
            scratchVec[mergedListPos--] = rightList[rightSubListPos--];
          }
        }

        if(leftSubListPos<leftSize && leftList[leftSubListPos] > rightList[rightSubListPos]) {
          scratchVec[mergedListPos--] = leftList[leftSubListPos--];
        } else {
          scratchVec[mergedListPos--] = rightList[rightSubListPos--];
        }
      }
    }
    //
    //  Copy the scratch data back to the main data array
    //    NKC, could likley improve this somewhat with a scratch state swap, though this may require isolating to 
    //    a power^2 number of threads
    //   ?????? !!!!! do swaps!!!!
    //
#pragma omp parallel for num_threads(numThreadsToUse) default(none) shared(inputVec, scratchVec, vecSize)
    for(unsigned i=0; i<vecSize; ++i) {
      inputVec[i] = scratchVec[i];
    }
  }
#endif


  template<typename T>
  void MyThreadedSort(std::vector<T>& inputVec) {
#ifdef _OPENMP
    unsigned maxNumThread = omp_get_max_threads();
    if(maxNumThread == 1) {
      std::sort(inputVec.begin(), inputVec.end());
      return;
    }
    unsigned inputVecSize = inputVec.size();

    //
    //  Need to tune this, though at some point using too many threads as compared to list size will definitely 
    //  cause more harm than good.  "5" below picked somewhat arbitrarly.
    //
    if(inputVecSize < 5) {
      std::sort(inputVec.begin(), inputVec.end());
      return;
    }
    unsigned numThreadToUse = std::min(maxNumThread, inputVecSize/5);

    int origNesting = omp_get_nested();
    omp_set_nested(1);

    std::vector<T> scratchSpace(inputVecSize);

    threadedMergeSort(inputVec.data(), scratchSpace.data(), inputVecSize, numThreadToUse);
  
    omp_set_nested(origNesting);
#else
    std::sort(inputVec.begin(), inputVec.end());
#endif
  }

  //
  //  General threaded sort capability
  //  Just wrapping gnu parallel sort with a few heuristics for now.  May need a fallback implementation
  //  if hit any platform specific support issues with the gnu code,  will see.
  //
  template<typename T>
  void threadedSort(std::vector<T>& inputVec) {
#ifdef _OPENMP
    typename std::vector<T>::iterator beginIter = inputVec.begin();
    typename std::vector<T>::iterator endIter   = inputVec.end();
    unsigned maxNumThread = omp_get_max_threads();
    if(maxNumThread == 1 || inputVec.size()<10) {
      std::sort(inputVec.begin(), inputVec.end());
    } else {
#if defined(USE_NATIVE_THREADED_SORT)
      //  Use a gnu implementation if at all possible
      unsigned targetNumThread = ceil(double(inputVec.size())/10.0);
      if(targetNumThread < maxNumThread) {
        omp_set_num_threads(targetNumThread);
        __gnu_parallel::sort(inputVec.begin(), inputVec.end());
        omp_set_num_threads(maxNumThread);
      } else {
        __gnu_parallel::sort(inputVec.begin(), inputVec.end());
      }
#else
      //  Home grown fallback version, needed at least on morgan
      MyThreadedSort(inputVec);
#endif
    }
#else
    std::sort(inputVec.begin(), inputVec.end());
#endif
  }

}
}

#endif // THREADED_SORT_H
