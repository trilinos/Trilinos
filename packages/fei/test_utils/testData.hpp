/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

