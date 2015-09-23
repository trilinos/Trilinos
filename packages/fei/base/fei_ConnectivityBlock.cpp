/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
  std::map<int,int>::const_iterator
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
  std::map<int,int>::const_iterator
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
  std::map<int,int>::const_iterator
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
  std::map<int,int>::const_iterator
    iter = connIDsOffsetMap_.find(ID);
  if (iter == connIDsOffsetMap_.end()) {
    return(NULL);
  }

  int ind = iter->second;
  int* ptr = &colConnectivities_[0];
  return(ptr+ind*numRecordsPerColConnectivity_);
}

