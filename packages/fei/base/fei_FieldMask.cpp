/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"

#include <string>
#include <exception>
#include <stdexcept>

#include "fei_ArrayUtils.hpp"
#include "fei_FieldMask.hpp"


fei::FieldMask::FieldMask()
  : maskID_(0),
    fieldIDs_(),
    fieldSizes_(),
    fieldInstances_(),
    fieldEqnOffsets_(),
    numFields_(0),
    numIndices_(1)
{
}

fei::FieldMask::FieldMask(const FieldMask& fm)
  : maskID_(fm.maskID_),
    fieldIDs_(fm.fieldIDs_),
    fieldSizes_(fm.fieldSizes_),
    fieldInstances_(fm.fieldInstances_),
    fieldEqnOffsets_(fm.fieldEqnOffsets_)
{
  numFields_ = fieldIDs_.size();
  numIndices_ = fm.numIndices_;
}

fei::FieldMask::FieldMask(int numFields,
			      const int* fieldIDs,
			      const int* fieldSizes,
			      const int* numInstancesOfThisFieldPerID)
  : maskID_(0),
    fieldIDs_(0, 4),
    fieldSizes_(0, 4),
    fieldInstances_(0, 4),
    fieldEqnOffsets_(0, 4)
{
  for(int i=0; i<numFields; ++i) {
    addField(fieldIDs[i], fieldSizes[i],
	     numInstancesOfThisFieldPerID[i]);
  }
}

fei::FieldMask::~FieldMask()
{
}

bool fei::FieldMask::hasFieldID(int fieldID) const
{
  return(
    fei::binarySearch(fieldID,
      fieldIDs_.size() ? &fieldIDs_[0] : 0,
      fieldIDs_.size())
    >= 0
    );
}

void fei::FieldMask::getFieldEqnOffset(int fieldID,
					   int& offset,
					   int& numInstancesPerID) const
{
  int idindex = 0;
  if (numFields_ < 2) {
    if (numFields_ < 1) {
      offset = 0;
      numInstancesPerID = 1;
      return;
    }

    if (fieldIDs_[0] != fieldID) {
      throw std::runtime_error("fei::FieldMask::getFieldEqnOffset: fieldID not found");
    }
  }
  else idindex = fei::binarySearch(fieldID, &fieldIDs_[0], fieldIDs_.size());

  if (idindex < 0) {
     throw std::runtime_error("fei::FieldMask::getFieldEqnOffset: fieldID not found");
  }

  offset = fieldEqnOffsets_[idindex];
  numInstancesPerID = fieldInstances_[idindex];
}

void fei::FieldMask::addField(int fieldID, int fieldSize, int numInstances)
{
  if (fieldID < 0) {
    throw std::runtime_error("fei::FieldMask ERROR, fieldID should be >= 0.");
  }

  int insertPoint = -1;
  int idindex = fei::binarySearch(fieldID, fieldIDs_, insertPoint);
  if (idindex >= 0) {
    fieldInstances_[idindex] += numInstances;
    for(unsigned i=idindex+1; i<fieldEqnOffsets_.size(); ++i) {
      fieldEqnOffsets_[i] += numInstances*fieldSize;
    }
  }
  else {
    fieldIDs_.insert(fieldIDs_.begin()+insertPoint, fieldID);

    fieldSizes_.insert(fieldSizes_.begin()+insertPoint, fieldSize);

    fieldInstances_.insert(fieldInstances_.begin()+insertPoint, numInstances);

    fieldEqnOffsets_.push_back(numInstances);

    int eqnOffset = 0;
    numIndices_ = 0;
    for(unsigned i=0; i<fieldIDs_.size(); ++i) {
      fieldEqnOffsets_[i] = eqnOffset;
      eqnOffset += fieldInstances_[i]*fieldSizes_[i];
      numIndices_ += fieldInstances_[i]*fieldSizes_[i];
    }

    numFields_ = fieldIDs_.size();
  }

  maskID_ = calculateMaskID();
}

int fei::FieldMask::calculateMaskID()
{
  return( calculateMaskID(fieldIDs_.size(), &fieldIDs_[0],
			  &fieldInstances_[0]) );
}

int fei::FieldMask::calculateMaskID(int numFields, const int* fieldIDs,
					const int* numInstancesPerID)
{
  int maskID = 0;
  for(int i=0; i<numFields; ++i) {
    maskID += (fieldIDs[i]+1)*numInstancesPerID[i] +(i+1)*1000;
  }

  return(maskID);
}

int fei::FieldMask::calculateMaskID(const FieldMask& fm,
					int fieldID, int numInstances)
{
  return( fm.maskID_ + (fieldID+1)*numInstances + (fm.numFields_+1)*1000 );
}

