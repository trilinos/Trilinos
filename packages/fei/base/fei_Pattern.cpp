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


#include "fei_Pattern.hpp"

//-----------------------------------------------------------------------------
fei::Pattern::Pattern(int numIDs, int idType, snl_fei::RecordCollection* recordCollection)
  : type_(Pattern::NO_FIELD),
    numIDs_(numIDs),
    totalNumFields_(0),
    numIndices_(numIDs),
    data_(),
    recordCollections_(numIDs, recordCollection)
{
  int i, len = numIDs_*4;
  data_.resize(len);
  int offset = 0;

  //set idTypes
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = idType;
  }

  //set numFieldsPerID
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = 0;
  }

  //set numIndicesPerID
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = 1;
  }

  //set fieldIDs
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = -1;
  }

  idTypes_        = &(data_[0]);
  numFieldsPerID_ = idTypes_+numIDs_;
  numIndicesPerID_= idTypes_+2*numIDs_;
  fieldIDs_       = idTypes_+3*numIDs_;
}

//-----------------------------------------------------------------------------
fei::Pattern::Pattern(int numIDs, int idType, snl_fei::RecordCollection* recordCollection,
			  int fieldID, int fieldSize)
  : type_(Pattern::SIMPLE),
    numIDs_(numIDs),
    totalNumFields_(numIDs),
    numIndices_(0), //numIndices_ to be calculated
    data_(),
    recordCollections_(numIDs, recordCollection)
{
  int i, len = numIDs_*4;
  data_.resize(len);
  int offset = 0;

  //set idTypes
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = idType;
  }

  //set numFieldsPerID
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = 1;
  }

  //set numIndicesPerID
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = fieldSize;
  }

  //set fieldIDs
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = fieldID;
  }

  numIndices_ = numIDs_*fieldSize;

  idTypes_        = &(data_[0]);
  numFieldsPerID_ = idTypes_+numIDs_;
  numIndicesPerID_= idTypes_+2*numIDs_;
  fieldIDs_       = idTypes_+3*numIDs_;
}

//-----------------------------------------------------------------------------
fei::Pattern::Pattern(int numIDs, int idType, snl_fei::RecordCollection* recordCollection,
			  const int* numFieldsPerID,
			  const int* fieldIDs,
			  const int* fieldSizes)
  : type_(Pattern::SINGLE_IDTYPE),
    numIDs_(numIDs),
    totalNumFields_(0), //totalNumFields_ to be calculated
    numIndices_(0), //numIndices_ to be calculated
    data_(),
    recordCollections_(numIDs, recordCollection)
{
  int i, len = numIDs_*3;
  int maxNumFieldsPerID = 0;
  bool oneDistinctFieldID = true;

  for(i=0; i<numIDs; ++i) len += numFieldsPerID[i];
  data_.resize(len);

  int offset = 0;
  //set idTypes
  for(i=0; i<numIDs_; ++i) {
    data_[offset++] = idType;
  }

  //set numFieldsPerID
  for(i=0; i<numIDs; ++i) {
    data_[offset++] = numFieldsPerID[i];
    if (numFieldsPerID[i] > maxNumFieldsPerID) {
      maxNumFieldsPerID = numFieldsPerID[i];
    }
  }

  //next set numIndicesPerID and fieldIDs
  int firstFieldID = 0;
  if (numIDs > 0) firstFieldID = fieldIDs[0];

  int fieldIDOffset = offset + numIDs;
  for(i=0; i<numIDs; ++i) {
    int thisNumIndices = 0;
    for(int j=0; j<numFieldsPerID[i]; ++j) {
      int fieldSize = fieldSizes[totalNumFields_];
      int fieldID = fieldIDs[totalNumFields_++];
      if (fieldID != firstFieldID) oneDistinctFieldID = false;
      data_[fieldIDOffset++] = fieldID;
      numIndices_ += fieldSize;
      thisNumIndices += fieldSize;
    }
    data_[offset+i] = thisNumIndices;
  }

  if (oneDistinctFieldID == true && maxNumFieldsPerID < 2) {
    type_ = Pattern::SIMPLE;
  }

  idTypes_        = &(data_[0]);
  numFieldsPerID_ = idTypes_+numIDs_;
  numIndicesPerID_= idTypes_+2*numIDs_;
  fieldIDs_       = idTypes_+3*numIDs_;
}

//-----------------------------------------------------------------------------
fei::Pattern::Pattern(int numIDs, const int* idTypes, snl_fei::RecordCollection*const* recordCollections,
			  const int* numFieldsPerID,
			  const int* fieldIDs,
			  const int* fieldSizes)
  : type_(Pattern::GENERAL),
    numIDs_(numIDs),
    totalNumFields_(0), //totalNumFields_ to be calculated
    numIndices_(0), //numIndices_ to be calculated
    data_(),
    recordCollections_(recordCollections, recordCollections+numIDs)
{
  int i, len = numIDs*3;
  int maxNumFieldsPerID = 0;
  bool oneDistinctFieldID = true;
  bool oneDistinctIDType = true;

  for(i=0; i<numIDs; ++i) len += numFieldsPerID[i];
  data_.resize(len);

  int firstIDType = 0;
  if (numIDs > 0) firstIDType = idTypes[0];

  int offset = 0;
  //set idTypes
  for(i=0; i<numIDs; ++i) {
    data_[offset++] = idTypes[i];
    if (idTypes[i] != firstIDType) oneDistinctIDType = false;
  }

  //set numFieldsPerID
  for(i=0; i<numIDs; ++i) {
    data_[offset++] = numFieldsPerID[i];
    if (numFieldsPerID[i] > maxNumFieldsPerID) {
      maxNumFieldsPerID = numFieldsPerID[i];
    }
  }

  //next set numIndicesPerID and fieldIDs
  int firstFieldID = 0;
  if (numIDs > 0) firstFieldID = fieldIDs[0];

  int fieldIDOffset = offset + numIDs;
  for(i=0; i<numIDs; ++i) {
    int thisNumIndices = 0;
    for(int j=0; j<numFieldsPerID[i]; ++j) {
      int fieldSize = fieldSizes[totalNumFields_];
      int fieldID = fieldIDs[totalNumFields_++];
      if (fieldID != firstFieldID) oneDistinctFieldID = false;
      data_[fieldIDOffset++] = fieldID;
      numIndices_ += fieldSize;
      thisNumIndices += fieldSize;
    }
    data_[offset+i] = thisNumIndices;
  }

  if (oneDistinctFieldID == true && maxNumFieldsPerID < 2 &&
      oneDistinctIDType == true) {
    type_ = Pattern::SIMPLE;
  }
  else if (oneDistinctIDType == true) {
    type_ = Pattern::SINGLE_IDTYPE;
  }

  idTypes_        = &(data_[0]);
  numFieldsPerID_ = idTypes_+numIDs_;
  numIndicesPerID_= idTypes_+2*numIDs_;
  fieldIDs_       = idTypes_+3*numIDs_;
}

//-----------------------------------------------------------------------------
fei::Pattern::~Pattern()
{
}

//-----------------------------------------------------------------------------
bool fei::Pattern::operator==(const fei::Pattern& rhs) const
{
  return type_ == rhs.type_ &&
         numIDs_ == rhs.numIDs_ &&
         totalNumFields_ == rhs.totalNumFields_ &&
         numIndices_ == rhs.numIndices_ &&
         data_ == rhs.data_;
}

//-----------------------------------------------------------------------------
bool fei::Pattern::operator!=(const fei::Pattern& rhs) const
{
  return !(*this == rhs);
}

