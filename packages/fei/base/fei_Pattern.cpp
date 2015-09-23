/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

