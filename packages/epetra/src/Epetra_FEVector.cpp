
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Epetra_FEVector.h>

#include <Epetra_LocalMap.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Util.h>
#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_SerialDenseVector.h>

#include <algorithm>

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(const Epetra_BlockMap& map,
                                 int numVectors,
         bool ignoreNonLocalEntries)
  : Epetra_MultiVector(map, numVectors),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_int_(),
    nonlocalIDs_LL_(),
    nonlocalElementSize_(),
    nonlocalCoefs_(),
    nonlocalMap_(0),
    exporter_(0),
    nonlocalVector_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstID_ = map.MinMyGID64();
  myNumIDs_ = map.NumMyElements();
  nonlocalCoefs_.resize(numVectors);
}

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
                                 double *A, int MyLDA, int NumVectors,
                                 bool ignoreNonLocalEntries)
 : Epetra_MultiVector(CV, Map, A, MyLDA, NumVectors),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_int_(),
    nonlocalIDs_LL_(),
    nonlocalElementSize_(),
    nonlocalCoefs_(),
    nonlocalMap_(0),
    exporter_(0),
    nonlocalVector_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstID_ = Map.MinMyGID64();
  myNumIDs_ = Map.NumMyElements();
  nonlocalCoefs_.resize(NumVectors);
}

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
                                 double **ArrayOfPointers, int NumVectors,
                                 bool ignoreNonLocalEntries)
 : Epetra_MultiVector(CV, Map, ArrayOfPointers, NumVectors),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_int_(),
    nonlocalIDs_LL_(),
    nonlocalElementSize_(),
    nonlocalCoefs_(),
    nonlocalMap_(0),
    exporter_(0),
    nonlocalVector_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstID_ = Map.MinMyGID64();
  myNumIDs_ = Map.NumMyElements();
  nonlocalCoefs_.resize(NumVectors);
}

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(const Epetra_FEVector& source)
  : Epetra_MultiVector(source),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_int_(),
    nonlocalIDs_LL_(),
    nonlocalElementSize_(),
    nonlocalCoefs_(),
    nonlocalMap_(0),
    exporter_(0),
    nonlocalVector_(0),
    ignoreNonLocalEntries_(source.ignoreNonLocalEntries_)
{
  *this = source;
}

//----------------------------------------------------------------------------
Epetra_FEVector::~Epetra_FEVector()
{
  destroyNonlocalData();
  destroyNonlocalMapAndExporter();

  nonlocalCoefs_.clear();
}

//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, true, vectorIndex) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const long long* GIDs,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, true, vectorIndex) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
                       const Epetra_SerialDenseVector& values,
                                         int vectorIndex)
{
  if (GIDs.Length() != values.Length()) {
    return(-1);
  }

  return( inputValues( GIDs.Length(), GIDs.Values(), values.Values(), true,
                       vectorIndex ) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& GIDs,
                       const Epetra_SerialDenseVector& values,
                                         int vectorIndex)
{
  if (GIDs.Length() != values.Length()) {
    return(-1);
  }

  return( inputValues( GIDs.Length(), GIDs.Values(), values.Values(), true,
                       vectorIndex ) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
           const int* numValuesPerID,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, true,
                       vectorIndex) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const long long* GIDs,
           const int* numValuesPerID,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, true,
                       vectorIndex) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, false,
                       vectorIndex) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const long long* GIDs,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, false,
                       vectorIndex) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
                       const Epetra_SerialDenseVector& values,
                                         int vectorIndex)
{
  if (GIDs.Length() != values.Length()) {
    return(-1);
  }

  return( inputValues( GIDs.Length(), GIDs.Values(), values.Values(), false,
                       vectorIndex) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& GIDs,
                       const Epetra_SerialDenseVector& values,
                                         int vectorIndex)
{
  if (GIDs.Length() != values.Length()) {
    return(-1);
  }

  return( inputValues( GIDs.Length(), GIDs.Values(), values.Values(), false,
                       vectorIndex) );
}
#endif
//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FEVector::inputValues(int numIDs,
                                 const int_type* GIDs,
                                 const double* values,
                                 bool suminto,
                                 int vectorIndex)
{
 //Important note!! This method assumes that there is only 1 point
 //associated with each element (GID), and writes to offset 0 in that
 //GID's block.

  for(int i=0; i<numIDs; ++i) {
    if (Map().MyGID(GIDs[i])) {
      if (suminto) {
        SumIntoGlobalValue(GIDs[i], 0, vectorIndex, values[i]);
      }
      else {
        ReplaceGlobalValue(GIDs[i], 0, vectorIndex, values[i]);
      }
    }
    else {
      if (!ignoreNonLocalEntries_) {
        EPETRA_CHK_ERR( inputNonlocalValue(GIDs[i], values[i], suminto,
              vectorIndex) );
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
           const int* numValuesPerID,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, false,
                       vectorIndex) );
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const long long* GIDs,
           const int* numValuesPerID,
                       const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, false,
                       vectorIndex) );
}
#endif
//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FEVector::inputValues(int numIDs,
                                 const int_type* GIDs,
                                 const int* numValuesPerID,
                                 const double* values,
                                 bool suminto,
                                 int vectorIndex)
{
  if(!Map().GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FEVector::inputValues mismatch between argument types (int/long long) and map type.", -1);

  int offset=0;
  for(int i=0; i<numIDs; ++i) {
    int numValues = numValuesPerID[i];
    if (Map().MyGID(GIDs[i])) {
      if (suminto) {
        for(int j=0; j<numValues; ++j) {
          SumIntoGlobalValue(GIDs[i], j, vectorIndex, values[offset+j]);
        }
      }
      else {
        for(int j=0; j<numValues; ++j) {
          ReplaceGlobalValue(GIDs[i], j, vectorIndex, values[offset+j]);
        }
      }
    }
    else {
      if (!ignoreNonLocalEntries_) {
        EPETRA_CHK_ERR( inputNonlocalValues(GIDs[i], numValues,
              &(values[offset]), suminto,
              vectorIndex) );
      }
    }
    offset += numValues;
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FEVector::inputNonlocalValue(int_type GID, double value, bool suminto,
                                        int vectorIndex)
{
  return inputNonlocalValues(GID, 1, &value, suminto, vectorIndex);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FEVector::inputNonlocalValues(int_type GID, int numValues,
                                         const double* values, bool suminto,
                                         int vectorIndex)
{
  if(!Map().GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FEVector::inputValues mismatch between argument types (int/long long) and map type.", -1);

  
  //find offset of GID in nonlocalIDs_var

  std::vector<int_type>& nonlocalIDs_var = nonlocalIDs<int_type>();

  typename std::vector<int_type>::iterator it = std::lower_bound(nonlocalIDs_var.begin(), nonlocalIDs_var.end(), GID);
  int offset = (int) (it - nonlocalIDs_var.begin());
  int insertPoint = offset;
  if (it == nonlocalIDs_var.end() || *it != GID) {
    offset = -1;
  }

  int elemSize = Map().MaxElementSize();
  if (offset >= 0) {
    //if offset >= 0 (meaning GID was found)
    //  put value in nonlocalCoefs_[vectorIndex][offset*elemSize]

    if (numValues != nonlocalElementSize_[offset]) {
      cerr << "Epetra_FEVector ERROR: block-size for GID " << GID << " is "
     << numValues<<" which doesn't match previously set block-size of "
     << nonlocalElementSize_[offset] << endl;
      return(-1);
    }

    offset = offset*elemSize;

    if (suminto) {
      for(int j=0; j<numValues; ++j) {
        nonlocalCoefs_[vectorIndex][offset+j] += values[j];
      }
    }
    else {
      for(int j=0; j<numValues; ++j) {
        nonlocalCoefs_[vectorIndex][offset+j] = values[j];
      }
    }
  }
  else {
    //else
    //  insert GID in nonlocalIDs_
    //  insert numValues   in nonlocalElementSize_
    //  insert values in nonlocalCoefs_

    nonlocalIDs_var.insert(it, GID);
    nonlocalElementSize_.insert(nonlocalElementSize_.begin()+insertPoint, numValues);

    //to keep nonlocalCoefs_[i] the same length for each vector in the multi-
    //vector, we'll insert positions for each vector even though values are
    //only being set for one of them...
    for(int i=0; i<NumVectors(); ++i) {
      for(int ii=0; ii<elemSize; ++ii) {
        nonlocalCoefs_[i].insert(nonlocalCoefs_[i].begin()+insertPoint*elemSize+ii, 0.0);
      }
    }

    for(int j=0; j<numValues; ++j) {
      nonlocalCoefs_[vectorIndex][insertPoint*elemSize+j] = values[j];
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FEVector::GlobalAssemble(Epetra_CombineMode mode,
                                    bool reuse_map_and_exporter)
{
  //In this method we need to gather all the non-local (overlapping) data
  //that's been input on each processor, into the (probably) non-overlapping
  //distribution defined by the map that 'this' vector was constructed with.

  //We don't need to do anything if there's only one processor or if
  //ignoreNonLocalEntries_ is true.
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    return(0);
  }

  if (nonlocalMap_ == 0 || !reuse_map_and_exporter) {
    createNonlocalMapAndExporter<int_type>();
  }

  Epetra_MultiVector& nonlocalVector = *nonlocalVector_;
  nonlocalVector.PutScalar(0.0);

  int elemSize = Map().MaxElementSize();
  for(int vi=0; vi<NumVectors(); ++vi) {
    for(size_t i=0; i<nonlocalIDs<int_type>().size(); ++i) {
      for(int j=0; j<nonlocalElementSize_[i]; ++j) {
        nonlocalVector.ReplaceGlobalValue(nonlocalIDs<int_type>()[i], j, vi,
                                          nonlocalCoefs_[vi][i*elemSize+j]);
      }
    }
  }

  EPETRA_CHK_ERR( Export(nonlocalVector, *exporter_, mode) );

  if (reuse_map_and_exporter) {
    zeroNonlocalData<int_type>();
  }
  else {
    destroyNonlocalData();
  }

  return(0);
}

int Epetra_FEVector::GlobalAssemble(Epetra_CombineMode mode,
                                    bool reuse_map_and_exporter)
{
  if(Map().GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    return GlobalAssemble<int>(mode, reuse_map_and_exporter);
#else
    throw ReportError("Epetra_FEVector::GlobalAssemble: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif

  if(Map().GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    return GlobalAssemble<long long>(mode, reuse_map_and_exporter);
#else
    throw ReportError("Epetra_FEVector::GlobalAssemble: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif

  throw ReportError("Epetra_FEVector::GlobalAssemble: Internal error, unable to determine global index type of maps", -1);
}
//----------------------------------------------------------------------------
template<typename int_type>
void Epetra_FEVector::createNonlocalMapAndExporter()
{
  std::vector<int_type>& nonlocalIDs_var = nonlocalIDs<int_type>();
  delete nonlocalMap_;
  int_type* nlIDptr = nonlocalIDs_var.size()>0 ? &nonlocalIDs_var[0] : NULL;
  int* nlElSzptr = nonlocalElementSize_.size()>0 ? &nonlocalElementSize_[0] : NULL;
  nonlocalMap_ = new Epetra_BlockMap ((int_type) -1, (int) nonlocalIDs_var.size(), nlIDptr,
                                      nlElSzptr, Map().IndexBase(), Map().Comm());
  delete exporter_;
  exporter_ = new Epetra_Export (*nonlocalMap_, Map());

  delete nonlocalVector_;
  nonlocalVector_ = new Epetra_MultiVector (*nonlocalMap_, NumVectors());
}

//----------------------------------------------------------------------------
void Epetra_FEVector::destroyNonlocalMapAndExporter()
{
  delete nonlocalMap_;    nonlocalMap_ = 0;
  delete exporter_;       exporter_ = 0;
  delete nonlocalVector_; nonlocalVector_ = 0;
}

//----------------------------------------------------------------------------
Epetra_FEVector& Epetra_FEVector::operator=(const Epetra_FEVector& source)
{
  if (this == &source) {
    // Don't allow self-assignment, since the allocations and
    // deallocations in the code below assume that source is a
    // different object than *this.
    return *this; 
  }
  // This redundantly checks for self-assignment, but the check is
  // inexpensive (just a pointer comparison).
  Epetra_MultiVector::Assign(source);

  nonlocalIDs_int_ = source.nonlocalIDs_int_;
  nonlocalIDs_LL_ = source.nonlocalIDs_LL_;
  nonlocalElementSize_ = source.nonlocalElementSize_;
  nonlocalCoefs_ = source.nonlocalCoefs_;

  return(*this);
}

//----------------------------------------------------------------------------
template<typename int_type>
void Epetra_FEVector::zeroNonlocalData()
{
  if (nonlocalIDs<int_type>().size() > 0) {
    int maxelemSize = Map().MaxElementSize();
    for(int vi=0; vi<NumVectors(); ++vi) {
      for(size_t i=0; i<nonlocalIDs<int_type>().size(); ++i) {
        int elemSize = nonlocalElementSize_[i];
        for(int j=0; j<elemSize; ++j) {
          nonlocalCoefs_[vi][i*maxelemSize+j] = 0.0;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void Epetra_FEVector::destroyNonlocalData()
{
  nonlocalIDs_int_.clear();
  nonlocalIDs_LL_.clear();
  nonlocalElementSize_.clear();

  if (nonlocalCoefs_.size() > 0) {
    for(int i=0; i<NumVectors(); ++i) {
      nonlocalCoefs_[i].clear();
    }
  }
}

