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

#ifndef _InputData_h_
#define _InputData_h_


#include <fei_macros.hpp>

#include <vector>

class ElemContribution {
 public:
  ElemContribution(){}
  ElemContribution(const ElemContribution& src)
    {
      rhsContributions = src.rhsContributions;
    }

  ~ElemContribution(){}

  bool operator==(const ElemContribution& rhs)
    {
      if (matrixContributions != rhs.matrixContributions) {
	cout << "matrixContributions don't match." << endl;
	return(false);
      }

      if ( rhsContributions != rhs.rhsContributions ) {
	cout << "rhsContributions don't match." << endl;
	return(false);
      }

      return(true);
    }

  bool operator!=(const ElemContribution& rhs)
    {
      return( !( *this == rhs) );
    }

  std::vector<double> matrixContributions;
  std::vector<double> rhsContributions;
};

class InputData {
 public:
  InputData(){}
  ~InputData()
    {
      for(int i=0; i<elemIDs.length(); i++) delete elemIDs[i];
    }

  std::vector<int> elemBlockIDs;
  std::vector<std::vector<int>*> elemIDs;
  std::vector<std::vector<ElemContribution>*> elemContributions;

  bool operator==(const InputData& rhs)
    {
      if (elemBlockIDs != rhs.elemBlockIDs) {
	cout << "elemBlockIDs don't match." << endl;
	return(false);
      }

      for(int i=0; i<elemIDs.length(); i++) {
	std::vector<ElemContribution>& elems = *(elemContributions[i]);
	std::vector<ElemContribution>& rhsElems = *(rhs.elemContributions[i]);

	for(int j=0; j<elemIDs[i]->length(); j++) {
	  int id1 = (*(elemIDs[i]))[j];
	  int id2 = (*(rhs.elemIDs[i]))[j];

	  if ( id1 != id2 ) {
	    cout << "elemIDs don't match. element-block " << elemBlockIDs[i]
		 << ", elemID in position " << j << " is " << id1
		 << ", doesn't match " << id2 << "." << endl;
	    return(false);
	  }

	  if (elems[j] != rhsElems[j]) {
	    cout << "element-block " << elemBlockIDs[i] << ", elemID " << id1
	      << "'s element-contributions don't match." << endl;
	    return(false);
	  }
	}
      }

      return(true);
    }

  bool operator!=(const InputData& rhs)
    {
      return( !( (*this) == rhs) );
    }

  int addElemID(int elemBlockID, int elemID)
    {
      //add and elemBlockID/elemID pair to the internal arrays if not already
      //present.

      int err, insertPoint = -1;
      int blkInd = elemBlockIDs.binarySearch(elemBlockID, insertPoint);
      if (blkInd < 0) {
	err = elemBlockIDs.insert(elemBlockIDs.begin()+insertPoint, elemBlockID);
	err += elemIDs.insert(elemIDs.begin()+insertPoint, new std::vector<int>);
	err += elemContributions.insert(elemContributions.begin()+insertPoint, new std::vector<ElemContribution>);
	if (err != 0) return(err);
	blkInd = insertPoint;
      }

      std::vector<int>& IDs = *(elemIDs[blkInd]);
      std::vector<ElemContribution>& ec = *(elemContributions[blkInd]);

      err = IDs.insertSorted(elemID);      
      if (err == -2) return(err);

      ElemContribution dummy;
      if (err >= 0) err = ec.insert(ec.begin()+err, dummy);
      if (err == -2) return(err);

      return(0);
    }

  int addElemMatrix(int elemBlockID, int elemID, std::vector<double>& matrixData)
    {
      int insertPoint = -1;
      int blkInd = elemBlockIDs.binarySearch(elemBlockID, insertPoint);
      if (blkInd < 0) {
	cerr << " addElemMatrix ERROR, elemBlockID " << (int)elemBlockID
	     << " not found" << endl;
	return(-1);
      }

      int elemIdx = elemIDs[blkInd]->binarySearch(elemID);
      if (elemIdx < 0) {
	cerr << "addElemMatrix ERROR, elemID " << (int)elemID << " not found."
	     <<endl;
	return(-1);
      }

      ElemContribution& elemContr = (*(elemContributions[blkInd]))[elemIdx];

      std::vector<double>& elemContrMatrix = elemContr.matrixContributions;
      int len = matrixData.length();
      int oldLen = elemContrMatrix.length();
      if (oldLen < len) {
	elemContrMatrix.resize(len);
	for(int i=oldLen; i<len; i++) elemContrMatrix[i] = 0.0;
      }

      for(int i=0; i<matrixData.length(); i++) {
	elemContrMatrix[i] += matrixData[i];
      }

      return(0);
    }

  int addElemRHS(int elemBlockID, int elemID, std::vector<double>& rhsData)
    {
      int insertPoint = -1;
      int blkInd = elemBlockIDs.binarySearch(elemBlockID, insertPoint);
      if (blkInd < 0) {
	cerr << " addElemRHS ERROR, elemBlockID " << (int)elemBlockID
	     << " not found" << endl;
	return(-1);
      }

      int elemIdx = elemIDs[blkInd]->binarySearch(elemID);
      if (elemIdx < 0) {
	cerr << "addElemRHS ERROR, elemID " << (int)elemID << " not found."<<endl;
	return(-1);
      }

      ElemContribution& elemContr = (*(elemContributions[blkInd]))[elemIdx];

      std::vector<double>& elemContrRHS = elemContr.rhsContributions;
      int len = rhsData.length();
      int oldLen = elemContrRHS.length();
      if (oldLen < len) {
	elemContrRHS.resize(len);
	for(int i=oldLen; i<len; i++) elemContrRHS[i] = 0.0;
      }

      for(int i=0; i<rhsData.length(); i++) {
	elemContrRHS[i] += rhsData[i];
      }

      return(0);
    }
};

#endif // _InputData_h_
