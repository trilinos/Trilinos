
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(const Epetra_BlockMap& Map,
                                 int numVectors,
				 bool ignoreNonLocalEntries)
  : Epetra_MultiVector(Map, numVectors),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_(NULL),
    nonlocalElementSize_(NULL),
    numNonlocalIDs_(0),
    numNonlocalIDsAlloc_(0),
    nonlocalCoefs_(NULL),
    numNonlocalCoefs_(0),
    numNonlocalCoefsAlloc_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstID_ = Map.MinMyGID();
  myNumIDs_ = Map.NumMyElements();
  nonlocalCoefs_ = new double*[numVectors];
  for(int i=0; i<numVectors; ++i) nonlocalCoefs_[i] = NULL;
}

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(const Epetra_FEVector& source)
  : Epetra_MultiVector(source),
    myFirstID_(0),
    myNumIDs_(0),
    nonlocalIDs_(NULL),
    nonlocalElementSize_(NULL),
    numNonlocalIDs_(0),
    numNonlocalIDsAlloc_(0),
    nonlocalCoefs_(NULL),
    numNonlocalCoefsAlloc_(0),
    ignoreNonLocalEntries_(source.ignoreNonLocalEntries_)
{
  *this = source;
}

//----------------------------------------------------------------------------
Epetra_FEVector::~Epetra_FEVector()
{
  destroyNonlocalData();

  delete [] nonlocalCoefs_;
  nonlocalCoefs_ = NULL;
}

//----------------------------------------------------------------------------
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
			                 const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, true, vectorIndex) );
}

//----------------------------------------------------------------------------
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

//----------------------------------------------------------------------------
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
					 const int* numValuesPerID,
			                 const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, true,
                       vectorIndex) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
			                 const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, values, false,
                       vectorIndex) );
}

//----------------------------------------------------------------------------
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

//----------------------------------------------------------------------------
int Epetra_FEVector::inputValues(int numIDs,
                                 const int* GIDs,
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
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
					 const int* numValuesPerID,
			                 const double* values,
                                         int vectorIndex)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, false,
                       vectorIndex) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputValues(int numIDs,
                                 const int* GIDs,
				 const int* numValuesPerID,
                                 const double* values,
                                 bool suminto,
                                 int vectorIndex)
{
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
int Epetra_FEVector::inputNonlocalValue(int GID, double value, bool suminto,
                                        int vectorIndex)
{
  int insertPoint = -1;

  //find offset of GID in nonlocalIDs_
  int offset = Epetra_Util_binary_search(GID, nonlocalIDs_, numNonlocalIDs_,
					 insertPoint);
  int elemSize = Map().MaxElementSize();
  if (offset >= 0) {
    //if offset >= 0
    //  put value in nonlocalCoefs_[vectorIndex][offset*elemSize]

    offset = offset*elemSize;
    if (suminto) {
      nonlocalCoefs_[vectorIndex][offset] += value;
    }
    else {
      nonlocalCoefs_[vectorIndex][offset] = value;
    }
  }
  else {
    //else
    //  insert GID in nonlocalIDs_
    //  insert 1   in nonlocalElementSize_
    //  insert value in nonlocalCoefs_[vectorIndex]

    int tmp1 = numNonlocalIDs_;
    int tmp2 = numNonlocalIDsAlloc_;
    EPETRA_CHK_ERR( Epetra_Util_insert(GID, insertPoint, nonlocalIDs_,
				       tmp1, tmp2) );
    --tmp1;
    EPETRA_CHK_ERR( Epetra_Util_insert(1, insertPoint, nonlocalElementSize_,
				       tmp1, numNonlocalIDsAlloc_) );

    numNonlocalIDs_ = tmp1;

    //to keep nonlocalCoefs_[i] the same length for each vector in the multi-
    //vector, we'll insert positions for each vector even though values are
    //only being set for one of them...
    for(int i=0; i<NumVectors(); ++i) {
      tmp1 = numNonlocalCoefs_;
      tmp2 = numNonlocalCoefsAlloc_;
      EPETRA_CHK_ERR( Epetra_Util_insert_empty_positions(nonlocalCoefs_[i],
                                       tmp1, tmp2,
                                       insertPoint*elemSize, elemSize));
      for(int ii=0; ii<elemSize; ++ii) {
        nonlocalCoefs_[i][insertPoint*elemSize+ii] = 0.0;
      }
    }
    numNonlocalCoefs_ = tmp1;
    numNonlocalCoefsAlloc_ = tmp2;

    nonlocalCoefs_[vectorIndex][insertPoint*elemSize] = value;
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputNonlocalValues(int GID, int numValues,
					 const double* values, bool suminto,
                                         int vectorIndex)
{
  int insertPoint = -1;

  //find offset of GID in nonlocalIDs_
  int offset = Epetra_Util_binary_search(GID, nonlocalIDs_, numNonlocalIDs_,
					 insertPoint);
  int elemSize = Map().MaxElementSize();
  if (offset >= 0) {
    //if offset >= 0
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

    int tmp1 = numNonlocalIDs_;
    int tmp2 = numNonlocalIDsAlloc_;
    EPETRA_CHK_ERR( Epetra_Util_insert(GID, insertPoint, nonlocalIDs_,
				       tmp1, tmp2) );
    --tmp1;
    EPETRA_CHK_ERR( Epetra_Util_insert(numValues, insertPoint, nonlocalElementSize_,
				       tmp1, numNonlocalIDsAlloc_) );

    numNonlocalIDs_ = tmp1;

    //to keep nonlocalCoefs_[i] the same length for each vector in the multi-
    //vector, we'll insert positions for each vector even though values are
    //only being set for one of them...
    for(int i=0; i<NumVectors(); ++i) {
      tmp1 = numNonlocalCoefs_;
      tmp2 = numNonlocalCoefsAlloc_;
      EPETRA_CHK_ERR( Epetra_Util_insert_empty_positions(nonlocalCoefs_[i],
                                       tmp1, tmp2,
				       insertPoint*elemSize, elemSize));
      for(int ii=0; ii<elemSize; ++ii) {
        nonlocalCoefs_[i][insertPoint*elemSize+ii] = 0.0;
      }
    }
    numNonlocalCoefs_ = tmp1;
    numNonlocalCoefsAlloc_ = tmp2;

    for(int j=0; j<numValues; ++j) {
      nonlocalCoefs_[vectorIndex][insertPoint*elemSize+j] = values[j];
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::GlobalAssemble(Epetra_CombineMode mode)
{
  //In this method we need to gather all the non-local (overlapping) data
  //that's been input on each processor, into the (probably) non-overlapping
  //distribution defined by the map that 'this' vector was constructed with.

  //We don't need to do anything if there's only one processor or if
  //ignoreNonLocalEntries_ is true.
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    return(0);
  }

  //First build a map that describes the data in nonlocalIDs_/nonlocalCoefs_.
  //We'll use the arbitrary distribution constructor of Map.

  Epetra_BlockMap sourceMap(-1, numNonlocalIDs_,
                            nonlocalIDs_, nonlocalElementSize_,
			    Map().IndexBase(), Map().Comm());

  //Now build a vector to hold our nonlocalCoefs_, and to act as the source-
  //vector for our import operation.
  Epetra_MultiVector nonlocalVector(sourceMap, NumVectors());

  int i,j;
  int elemSize = Map().MaxElementSize();
  for(int vi=0; vi<NumVectors(); ++vi) {
    for(i=0; i<numNonlocalIDs_; ++i) {
      for(j=0; j<nonlocalElementSize_[i]; ++j) {
        nonlocalVector.ReplaceGlobalValue(nonlocalIDs_[i], j, vi,
                                          nonlocalCoefs_[vi][i*elemSize+j]);
      }
    }
  }

  Epetra_Export exporter(sourceMap, Map());

  EPETRA_CHK_ERR( Export(nonlocalVector, exporter, mode) );

  destroyNonlocalData();

  return(0);
}

//----------------------------------------------------------------------------
Epetra_FEVector& Epetra_FEVector::operator=(const Epetra_FEVector& source)
{
  Epetra_MultiVector::Assign(source);

  destroyNonlocalData();

  delete [] nonlocalCoefs_;

  if (source.numNonlocalIDsAlloc_ > 0) {
    numNonlocalIDsAlloc_ = source.numNonlocalIDsAlloc_;
    numNonlocalIDs_ = source.numNonlocalIDs_;
    nonlocalIDs_ = new int[numNonlocalIDsAlloc_];
    nonlocalElementSize_ = new int[numNonlocalIDsAlloc_];
    for(int i=0; i<numNonlocalIDs_; ++i) {
      int elemSize = source.nonlocalElementSize_[i];
      nonlocalIDs_[i] = source.nonlocalIDs_[i];
      nonlocalElementSize_[i] = elemSize;
    }
  }

  nonlocalCoefs_ = new double*[NumVectors()];
  for(int i=0; i<NumVectors(); ++i) nonlocalCoefs_[i] = NULL;

  numNonlocalCoefs_ = source.numNonlocalCoefs_;
  numNonlocalCoefsAlloc_ = source.numNonlocalCoefsAlloc_;

  if (numNonlocalCoefsAlloc_ > 0) {
    for(int vi=0; vi<NumVectors(); ++vi) {
      nonlocalCoefs_[vi] = new double[numNonlocalCoefsAlloc_];
      int maxelemSize = Map().MaxElementSize();
      for(int i=0; i<numNonlocalIDs_; ++i) {
        int elemSize = source.nonlocalElementSize_[i];
        for(int j=0; j<elemSize; ++j) {
          nonlocalCoefs_[vi][i*maxelemSize+j] = source.nonlocalCoefs_[vi][i*maxelemSize+j];
        }
      }
    }
  }

  return(*this);
}

//----------------------------------------------------------------------------
void Epetra_FEVector::destroyNonlocalData()
{
  if (numNonlocalIDsAlloc_ > 0) {
    delete [] nonlocalIDs_;
    delete [] nonlocalElementSize_;
    nonlocalIDs_ = NULL;
    nonlocalElementSize_ = NULL;
    numNonlocalIDs_ = 0;
    numNonlocalIDsAlloc_ = 0;
  }

  if (nonlocalCoefs_ != NULL && numNonlocalCoefsAlloc_ > 0) {
    for(int i=0; i<NumVectors(); ++i) {
      delete [] nonlocalCoefs_[i];
      nonlocalCoefs_[i] = NULL;
    }

    numNonlocalCoefs_ = 0;
    numNonlocalCoefsAlloc_ = 0;
  }

  return;
}

