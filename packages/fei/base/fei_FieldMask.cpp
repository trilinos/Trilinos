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
    fieldEqnOffsets_(),
    numFields_(0),
    numIndices_(1)
{
}

fei::FieldMask::FieldMask(const FieldMask& fm)
  : maskID_(fm.maskID_),
    fieldIDs_(fm.fieldIDs_),
    fieldSizes_(fm.fieldSizes_),
    fieldEqnOffsets_(fm.fieldEqnOffsets_)
{
  numFields_ = fieldIDs_.size();
  numIndices_ = fm.numIndices_;
}

fei::FieldMask::FieldMask(int numFields,
			      const int* fieldIDs,
			      const int* fieldSizes)
  : maskID_(0),
    fieldIDs_(0, 4),
    fieldSizes_(0, 4),
    fieldEqnOffsets_(0, 4)
{
  for(int i=0; i<numFields; ++i) {
    addField(fieldIDs[i], fieldSizes[i]);
  }
}

fei::FieldMask::~FieldMask()
{
}

void fei::FieldMask::getFieldEqnOffset(int fieldID,
					   int& offset) const
{
  int idindex = 0;
  if (numFields_ < 2) {
    if (numFields_ < 1) {
      offset = 0;
      return;
    }

    if (fieldIDs_[0] != fieldID) {
      throw std::runtime_error("fei::FieldMask::getFieldEqnOffset: fieldID not found");
    }
  }
  else {
    idindex = -1;
    for(size_t i=0; i<fieldIDs_.size(); ++i) {
      if (fieldIDs_[i] == fieldID) {
        idindex = i; break;
      }
    }
  }

  if (idindex < 0) {
    throw std::runtime_error("fei::FieldMask::getFieldEqnOffset: fieldID not found");
  }

  offset = fieldEqnOffsets_[idindex];
}

void fei::FieldMask::addField(int fieldID, int fieldSize)
{
  if (fieldID < 0) {
    throw std::runtime_error("fei::FieldMask ERROR, fieldID should be >= 0.");
  }

  int insertPoint = -1;
  int idindex = fei::binarySearch(fieldID, fieldIDs_, insertPoint);
  if (idindex >= 0) {
    for(unsigned i=idindex+1; i<fieldEqnOffsets_.size(); ++i) {
      fieldEqnOffsets_[i] += fieldSize;
    }
  }
  else {
    fieldIDs_.insert(fieldIDs_.begin()+insertPoint, fieldID);

    fieldSizes_.insert(fieldSizes_.begin()+insertPoint, fieldSize);

    fieldEqnOffsets_.push_back(1);

    int eqnOffset = 0;
    numIndices_ = 0;
    for(unsigned i=0; i<fieldIDs_.size(); ++i) {
      fieldEqnOffsets_[i] = eqnOffset;
      eqnOffset += fieldSizes_[i];
      numIndices_ += fieldSizes_[i];
    }

    numFields_ = fieldIDs_.size();
  }

  maskID_ = calculateMaskID();
}

int fei::FieldMask::calculateMaskID()
{
  return( calculateMaskID(fieldIDs_.size(), &fieldIDs_[0]));
}

int fei::FieldMask::calculateMaskID(int numFields, const int* fieldIDs)
{
  int maskID = 0;
  for(int i=0; i<numFields; ++i) {
    maskID += (fieldIDs[i]+1) +(i+1)*1000;
  }

  return(maskID);
}

int fei::FieldMask::calculateMaskID(const FieldMask& fm, int fieldID)
{
//  if (fm.hasFieldID(fieldID)) return fm.maskID_;
  return( fm.maskID_ + (fieldID+1) + (fm.numFields_+1)*1000 );
}

