/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _testData_h_
#define _testData_h_

#include <fei_macros.hpp>

#include <vector>

/** Simple container of arbitrarily chosen test data.
 */
class testData {
public:
  testData(int localProc, int numProcs)
    : fieldIDs(2),
    fieldSizes(2),
    idTypes(2),
    ids(4),
    sharedIDs(0),
    numSharingProcsPerID(0),
    sharingProcs(0)
    {
      //this testData object contains the following:
      //
      //fieldIDs         3   9
      //fieldSizes       1   3
      //idTypes          0   5
      //ids  length 4, first 2 ids shared with localProc-1,
      //                last 2 ids shared with localProc+1
      //ids[i] = localProc*2 + i
      //sharedIDs, numSharingProcsPerID, sharingProcs
      //
      fieldIDs[0] = 3; fieldIDs[1] = 9;
      fieldSizes[0] = 1; fieldSizes[1] = 3;
      idTypes[0] = 0; idTypes[1] = 5;
      for(int i=0; i<4; ++i) {
	ids[i] = localProc*2 + i;
      }

      if (localProc > 0) {
	sharedIDs.push_back(ids[0]);
	sharedIDs.push_back(ids[1]);
	numSharingProcsPerID.push_back(1);
	sharingProcs.push_back(localProc-1);
	numSharingProcsPerID.push_back(1);
	sharingProcs.push_back(localProc-1);
      }

      if (localProc < numProcs-1) {
	sharedIDs.push_back(ids[2]);
	sharedIDs.push_back(ids[3]);
	numSharingProcsPerID.push_back(1);
	sharingProcs.push_back(localProc+1);
	numSharingProcsPerID.push_back(1);
	sharingProcs.push_back(localProc+1);
      }
    }

  virtual ~testData()
    {
    }

  std::vector<int> fieldIDs;
  std::vector<int> fieldSizes;
  std::vector<int> idTypes;
  std::vector<int> ids;
  std::vector<int> sharedIDs;
  std::vector<int> numSharingProcsPerID;
  std::vector<int> sharingProcs;
};

#endif // _testData_h_

