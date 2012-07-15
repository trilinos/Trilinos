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


#include <fei_macros.hpp>

#include <fei_defs.h>

#include <fei_Pattern.hpp>

#include <fei_ConnectivityBlock.hpp>

#undef fei_file
#define fei_file "fei_ConnectivityBlock.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
fei::ConnectivityBlock::ConnectivityBlock(int blockID,
				      fei::Pattern* pattern,
				      int numConnectivities)
  : blockID_(blockID),
    pattern_(pattern),
    colPattern_(NULL),
    isSymmetric_(true),
    isDiagonal_(false),
    connIDsOffsetMap_(),
    connectivityOffsets_(),
    numRecordsPerConnectivity_(pattern->getNumIDs()),
    connectivities_(pattern->getNumIDs()*numConnectivities),
    numRecordsPerColConnectivity_(0),
    colConnectivities_(),
    fieldID_(-99),
    haveFieldID_(false)
{
}

//----------------------------------------------------------------------------
fei::ConnectivityBlock::ConnectivityBlock(int blockID,
				      fei::Pattern* rowpattern,
				      fei::Pattern* colpattern,
				      int numConnectivities)
  : blockID_(blockID),
    pattern_(rowpattern),
    colPattern_(colpattern),
    isSymmetric_(false),
    isDiagonal_(false),
    connIDsOffsetMap_(),
    connectivityOffsets_(),
    numRecordsPerConnectivity_(rowpattern->getNumIDs()),
    connectivities_(rowpattern->getNumIDs()*numConnectivities),
    numRecordsPerColConnectivity_(colpattern->getNumIDs()),
    colConnectivities_(colpattern->getNumIDs()*numConnectivities),
    fieldID_(-99),
    haveFieldID_(false)
{
}

//----------------------------------------------------------------------------
fei::ConnectivityBlock::ConnectivityBlock(int numRowIDs,
				      const int* rowIDs,
				      const int* rowOffsets,
				      bool offsets_are_lengths)
  : blockID_(-1),
    pattern_(NULL),
    colPattern_(NULL),
    isSymmetric_(false),
    isDiagonal_(false),
    connIDsOffsetMap_(),
    connectivityOffsets_(),
    numRecordsPerConnectivity_(0),
    connectivities_(),
    numRecordsPerColConnectivity_(0),
    colConnectivities_(),
    fieldID_(-99),
    haveFieldID_(false)
{
  connectivities_.resize(numRowIDs);
  connectivityOffsets_.resize(numRowIDs+1);

  int clen = 0;
  if (offsets_are_lengths) {
    int sum = 0;
    for(int ii=0; ii<numRowIDs; ++ii) {
      sum += rowOffsets[ii];
    }
    clen = sum;
  }
  else clen = rowOffsets[numRowIDs];

  colConnectivities_.resize(clen);

  int i;
  if (offsets_are_lengths) {
    int offset = 0;
    for(i=0; i<numRowIDs; ++i) {
      connIDsOffsetMap_[rowIDs[i]] = i;
      connectivityOffsets_[i] = offset;
      offset += rowOffsets[i];
    }
    connectivityOffsets_[numRowIDs] = offset;
  }
  else {
    for(i=0; i<numRowIDs; ++i) {
      connIDsOffsetMap_[rowIDs[i]] = i;
      connectivityOffsets_[i] = rowOffsets[i];
    }
    connectivityOffsets_[numRowIDs] = rowOffsets[numRowIDs];
  }
}

//----------------------------------------------------------------------------
fei::ConnectivityBlock::ConnectivityBlock(int fldID,
				      int numRowIDs,
				      const int* rowIDs,
				      const int* rowOffsets,
				      bool offsets_are_lengths)
  : blockID_(-1),
    pattern_(NULL),
    colPattern_(NULL),
    isSymmetric_(false),
    isDiagonal_(false),
    connIDsOffsetMap_(),
    connectivityOffsets_(),
    numRecordsPerConnectivity_(0),
    connectivities_(),
    numRecordsPerColConnectivity_(0),
    colConnectivities_(),
    fieldID_(fldID),
    haveFieldID_(true)
{
  connectivities_.resize(numRowIDs);
  connectivityOffsets_.resize(numRowIDs+1);

  int clen = 0;
  if (offsets_are_lengths) {
    int sum = 0;
    for(int ii=0; ii<numRowIDs; ++ii) {
      sum += rowOffsets[ii];
    }
    clen = sum;
  }
  else clen = rowOffsets[numRowIDs];

  colConnectivities_.resize(clen);

  int i;
  if (offsets_are_lengths) {
    int offset = 0;
    for(i=0; i<numRowIDs; ++i) {
      connIDsOffsetMap_[rowIDs[i]] = i;
      connectivityOffsets_[i] = offset;
      offset += rowOffsets[i];
    }
    connectivityOffsets_[numRowIDs] = offset;
  }
  else {
    for(i=0; i<numRowIDs+1; ++i) {
      connIDsOffsetMap_[rowIDs[i]] = i;
      connectivityOffsets_[i] = rowOffsets[i];
    }
    connectivityOffsets_[numRowIDs] = rowOffsets[numRowIDs];
  }
}

//----------------------------------------------------------------------------
fei::ConnectivityBlock::~ConnectivityBlock()
{
}

//----------------------------------------------------------------------------
const int* fei::ConnectivityBlock::getRowConnectivity(int ID) const
{

  MapIntInt::const_iterator
    iter = connIDsOffsetMap_.find(ID);
  if (iter == connIDsOffsetMap_.end()) {
    return(NULL);
  }

  int ind = iter->second;
  const int* ptr = &connectivities_[0];
  return( ptr + ind*numRecordsPerConnectivity_);
}

//----------------------------------------------------------------------------
int* fei::ConnectivityBlock::getRowConnectivity(int ID)
{
  MapIntInt::const_iterator
    iter = connIDsOffsetMap_.find(ID);
  if (iter == connIDsOffsetMap_.end()) {
    return(NULL);
  }

  int ind = iter->second;
  int* ptr = &connectivities_[0];
  return( ptr + ind*numRecordsPerConnectivity_);
}

//----------------------------------------------------------------------------
const int* fei::ConnectivityBlock::getColConnectivity(int ID) const
{
  MapIntInt::const_iterator
    iter = connIDsOffsetMap_.find(ID);
  if (iter == connIDsOffsetMap_.end()) {
    return(NULL);
  }

  int ind = iter->second;
  const int* ptr = &colConnectivities_[0];
  return(ptr+ind*numRecordsPerColConnectivity_);
}

//----------------------------------------------------------------------------
int* fei::ConnectivityBlock::getColConnectivity(int ID)
{
  MapIntInt::const_iterator
    iter = connIDsOffsetMap_.find(ID);
  if (iter == connIDsOffsetMap_.end()) {
    return(NULL);
  }

  int ind = iter->second;
  int* ptr = &colConnectivities_[0];
  return(ptr+ind*numRecordsPerColConnectivity_);
}

