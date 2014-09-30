
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



namespace stk {


  

  //
  //  Recursive driver helper routine, Try to find at least numIdsRequested free ids in the range rangeStart to (rangeEnd-1). 
  //  Return success(0) or failure code if the ids could not be found.  If sucessful append the found ids to the returnIds vector.
  //  Note, existing ids must be sorted.
  //
  inline int ParallelFindFreeIdsInRangeGlobal(const ParallelMachine        comm, 
                                       const unsigned               rangeStart, 
                                       const unsigned               rangeEnd, 
                                       const unsigned               numIdsRequested,
                                       const std::vector<unsigned>::iterator existingIdStart,
                                       const std::vector<unsigned>::iterator existingIdEnd,
                                       const unsigned               myProcFirstReturnIndex,
                                       const unsigned               myProcLastReturnIndex,
                                       unsigned&                    curReturnIndex,
                                       std::vector<unsigned>&       returnIds) {
    assert(rangeEnd >= rangeStart);

    if(numIdsRequested == 0) {
      return 0;
    }

    //  Defines length of arrays used in MPI collectives.  Want to strike a balance where longer message takes somewhat longer to
    //  communicate but more finely subdivides the domain.  Ideally pick this number to be roughly where the mpi collective starts
    //  hitting linear scaling in N.  Try 64 for now.
    const unsigned MAX_SUBDIV = 64;
    //
    //  Check for 'leaf' case where we are effectively checking individual id blocks rather than ranges
    //
    unsigned numSubdivision;  
    unsigned rangeSize = (rangeEnd-rangeStart)+1;
    if(rangeSize < MAX_SUBDIV) {
      numSubdivision = rangeSize;
    } else if (rangeSize <= 0) {
      return 1;
    } else {
      numSubdivision = MAX_SUBDIV;
    }

    //
    //  Divide the range into 'numSubdivision' equal segments.  Determine locally the segment each existing id falls in.
    //  Note, last subdivision could be somewhat large due to rounding effects
    //
    unsigned subdivisionSize = rangeSize/numSubdivision;

    unsigned localSubDivCount[MAX_SUBDIV];
    unsigned subDivStart     [MAX_SUBDIV];
    unsigned subDivEnd       [MAX_SUBDIV];
    for(unsigned i=0; i < numSubdivision-1; ++i) {
      localSubDivCount[i] = 0;
      subDivStart[i] = subdivisionSize*i        + rangeStart;
      subDivEnd[i]   = (subdivisionSize*(i+1)-1 + rangeStart);
    }

    localSubDivCount[numSubdivision-1] = 0;
    subDivStart     [numSubdivision-1] = subdivisionSize*(numSubdivision-1) + rangeStart;
    subDivEnd       [numSubdivision-1] = rangeEnd;

    for(std::vector<unsigned>::iterator i=existingIdStart; i != existingIdEnd; ++i) {
      unsigned curId = (*i);
      if(curId < rangeStart || curId > rangeEnd) {
        continue;
      }      
      unsigned curSubdivision = (curId-rangeStart)/subdivisionSize;
      if(curSubdivision >= numSubdivision) {
        curSubdivision = numSubdivision-1;
      }
      assert(curSubdivision < numSubdivision);
      localSubDivCount[curSubdivision]++;
    }
    //
    //  Get a global consistent sum of the subdomain usage
    //
    unsigned globalSubDivCount[MAX_SUBDIV];
    int mpiResult = MPI_Allreduce(localSubDivCount, globalSubDivCount, numSubdivision, MPI_UNSIGNED, MPI_SUM, comm);
    if(mpiResult != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed");
      return 2;
    }
    //
    //  Extract Ids, if one or more divisions are empty grab up all ids for that division until
    //  either run out of empty divisions or all ids are alloacated.  Also extract the array
    //  of subdivision fill ammount for potential later use.
    //

    unsigned subFillLength = 0;
    std::pair<unsigned, unsigned> subFill[MAX_SUBDIV];

    unsigned numNeededIds = numIdsRequested;    
    unsigned totalIdsAvailable = 0;


    for(unsigned idiv = 0; idiv<numSubdivision; ++idiv) {

      if(globalSubDivCount[idiv] == 0 && numNeededIds > 0) {
        /*
        //
        //  Empty subdivision, Allocate its indexes to needed processors
        //
        unsigned rangeSize = (subDivEnd[idiv]-subDivStart[idiv])+1;
        unsigned numIdsToUse = min(numNeededIDs, rangeSize);
        //
        //  Case 1, return range above proc range
        //
        if(curReturnIndex >= myProcLastReturnIndex) {
          numNeededIds -= numIdsToUse;
          curReturnIndex += numNeededIds;
        }
        //
        //  Case 2, return range below proc range
        //
        */

        for(unsigned index=subDivStart[idiv]; index<=subDivEnd[idiv]; ++index) {
          if(numNeededIds == 0) break;

          

          if(curReturnIndex >= myProcFirstReturnIndex && curReturnIndex < myProcLastReturnIndex) {
            returnIds.push_back(index);
          }
          curReturnIndex++;
          --numNeededIds;
        }
      } else {
        unsigned curRangeSize = (subDivEnd[idiv]-subDivStart[idiv])+1;
        if(curRangeSize > globalSubDivCount[idiv]) {  
          unsigned curNumAvail = curRangeSize-globalSubDivCount[idiv];
          subFill[subFillLength] = std::pair<unsigned, unsigned>(curNumAvail, idiv);
          subFillLength++;
          totalIdsAvailable += curNumAvail;
        }
      }
    }
    if(numNeededIds == 0) {
      return 0; //SUCCESS!
    }
    if(numNeededIds > totalIdsAvailable) {
      throw std::runtime_error("In ParallelFindFreeIdsInRange, cannot allocate sufficent ids, aborting");
      return 6;
    }


    //
    //  Sufficent ids not found.  Sort subdivisions to find the most promising chunks and then fill those
    //
    std::sort(subFill, subFill+subFillLength);

    unsigned densityFactor = totalIdsAvailable/numNeededIds;


    

    for(unsigned  idiv = subFillLength-1; idiv>=0; --idiv) {

      unsigned curDiv        = subFill[idiv].second;

      unsigned curRangeStart = subDivStart[curDiv];
      unsigned curRangeEnd   = subDivEnd[curDiv];

      unsigned curRangeSize  = (curRangeEnd-curRangeStart)+1;
      unsigned numFree       = curRangeSize-globalSubDivCount[curDiv];

      if(numFree > 0) {
        unsigned maxNumToFill;
        if(densityFactor > 5) {  
          //  Sparse fill.  Fill no more than half of each range
          maxNumToFill = ceil(double(numFree)/2.0);
        } else {
          // Dense fill.  Fill entire range
          maxNumToFill = numFree;
        }

        unsigned numToFill = std::min(maxNumToFill, numNeededIds);
        //
        //  Extract the existing id range that overlaps this subdivision
        //
        const std::vector<unsigned>::iterator newExistingStart = std::lower_bound(existingIdStart, existingIdEnd, curRangeStart);
        const std::vector<unsigned>::iterator newExistingEnd   = std::upper_bound(existingIdStart, existingIdEnd, curRangeEnd-1);

        int returnCode = ParallelFindFreeIdsInRangeGlobal(comm, curRangeStart, curRangeEnd, numToFill, newExistingStart, newExistingEnd, 
                                                          myProcFirstReturnIndex, myProcLastReturnIndex, curReturnIndex, returnIds);

        if(returnCode != 0) {
          return returnCode;
        }
        numNeededIds -= numToFill;

        if(numNeededIds == 0) {
          return 0; //SUCCESS
        }
      }
    }
    //  Something bad has happened, if ids where successfully filled should not get here
    return 5;
  }




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
  inline int parallel_index_gap_finder_global(ParallelMachine comm, const unsigned minId, 
                                       const unsigned maxId, 
                                       const std::vector<unsigned>& localIdsInUse, 
                                       const unsigned localNumIdsNeeded,            
                                       const unsigned globalNumIdsNeeded,
                                       const unsigned prefixSumOfLocalIdsNeeded,
                                       std::vector<unsigned>& returnIds ) {

    if(minId > maxId) {
      throw std::runtime_error("In parallel_index_gap_finder_global, invalid max and min id range");
    }

    //const unsigned p_size = parallel_machine_size( comm );
    returnIds.clear();
    returnIds.reserve(globalNumIdsNeeded);
    //
    //  Create a sorted id vector needed by the id finding algorithm
    //
    std::vector<unsigned> sortedIdsInUse = localIdsInUse;
    std::sort(sortedIdsInUse.begin(), sortedIdsInUse.end());
    unsigned curIndex = 0;
    return ParallelFindFreeIdsInRangeGlobal(comm, minId, maxId-1, globalNumIdsNeeded, sortedIdsInUse.begin(), sortedIdsInUse.end(), 
                                            prefixSumOfLocalIdsNeeded, prefixSumOfLocalIdsNeeded+localNumIdsNeeded, curIndex, returnIds);
  }



}

#endif

