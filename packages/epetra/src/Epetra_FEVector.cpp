
#include <Epetra_FEVector.h>

#include <Epetra_LocalMap.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Util.h>

//
//Initial implementation of Epetra_FEVector. At this point it is
//capable of accepting overlapping data and then gathering that data onto
//the "owning" processors when GlobalAssemble() is called.
//
//Issues to be resolved:
//1. Assumptions are currently made implicitly, that ElementSize==1 for the
//Epetra_BlockMap that's provided at construction. Code needs to be
//generalized for varying element-sizes.
//
//2. Implementation of SumIntoGlobalValues() and ReplaceGlobalValues() is
//currently inefficient (probably). Code needs to be optimized.
//

//----------------------------------------------------------------------------
Epetra_FEVector::Epetra_FEVector(const Epetra_BlockMap& Map,
				 bool ignoreNonLocalEntries)
  : Epetra_MultiVector(Map, 1),
    myFirstID_(0),
    myNumIDs_(0),
    myCoefs_(NULL),
    nonlocalIDs_(NULL),
    nonlocalElementSize_(NULL),
    numNonlocalIDs_(0),
    allocatedNonlocalLength_(0),
    nonlocalCoefs_(NULL),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstID_ = Map.MinMyGID();
  myNumIDs_ = Map.NumMyElements();

  //Currently we impose the restriction that NumVectors==1, so we won't
  //need the LDA argument when calling ExtractView. Hence the "dummy" arg.
  int dummy;
  ExtractView(&myCoefs_, &dummy);
}

//----------------------------------------------------------------------------
Epetra_FEVector::~Epetra_FEVector()
{
  destroyNonlocalData();
}

//----------------------------------------------------------------------------
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, values, true) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
					 const int* numValuesPerID,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, true) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, values, false) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
					 const int* numValuesPerID,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, numValuesPerID, values, false) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputValues(int numIDs,
                                 const int* GIDs,
                                 const double* values,
                                 bool accumulate)
{
 //Important note!! This method assumes that there is only 1 point
 //associated with each element.

  for(int i=0; i<numIDs; ++i) {
    if (Map().MyGID(GIDs[i])) {
      if (accumulate) {
        SumIntoGlobalValue(GIDs[i], 0, 0, values[i]);
      }
      else {
        ReplaceGlobalValue(GIDs[i], 0, 0, values[i]);
      }
    }
    else {
      if (!ignoreNonLocalEntries_) {
	EPETRA_CHK_ERR( inputNonlocalValue(GIDs[i], values[i], accumulate) );
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputValues(int numIDs,
                                 const int* GIDs,
				 const int* numValuesPerID,
                                 const double* values,
                                 bool accumulate)
{
  int offset=0;
  for(int i=0; i<numIDs; ++i) {
    int numValues = numValuesPerID[i];
    if (Map().MyGID(GIDs[i])) {
      if (accumulate) {
	for(int j=0; j<numValues; ++j) {
	  SumIntoGlobalValue(GIDs[i], j, 0, values[offset+j]);
	}
      }
      else {
	for(int j=0; j<numValues; ++j) {
	  ReplaceGlobalValue(GIDs[i], j, 0, values[offset+j]);
	}
      }
    }
    else {
      if (!ignoreNonLocalEntries_) {
	EPETRA_CHK_ERR( inputNonlocalValues(GIDs[i], numValues,
					    &(values[offset]), accumulate) );
      }
    }
    offset += numValues;
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputNonlocalValue(int GID, double value, bool accumulate)
{
  int insertPoint = -1;

  //find offset of GID in nonlocalIDs_
  int offset = Epetra_Util_binary_search(GID, nonlocalIDs_, numNonlocalIDs_,
					 insertPoint);
  if (offset >= 0) {
    //if offset >= 0
    //  put value in nonlocalCoefs_[offset][0]

    if (accumulate) {
      nonlocalCoefs_[offset][0] += value;
    }
    else {
      nonlocalCoefs_[offset][0] = value;
    }
  }
  else {
    //else
    //  insert GID in nonlocalIDs_
    //  insert 1   in nonlocalElementSize_
    //  insert value in nonlocalCoefs_

    int tmp1 = numNonlocalIDs_;
    int tmp2 = allocatedNonlocalLength_;
    int tmp3 = allocatedNonlocalLength_;
    EPETRA_CHK_ERR( Epetra_Util_insert(GID, insertPoint, nonlocalIDs_,
				       tmp1, tmp2) );
    --tmp1;
    EPETRA_CHK_ERR( Epetra_Util_insert(1, insertPoint, nonlocalElementSize_,
				       tmp1, tmp3) );
    double* values = new double[1];
    values[0] = value;
    EPETRA_CHK_ERR( Epetra_Util_insert(values, insertPoint, nonlocalCoefs_,
				       numNonlocalIDs_, allocatedNonlocalLength_) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputNonlocalValues(int GID, int numValues,
					 const double* values, bool accumulate)
{
  int insertPoint = -1;

  //find offset of GID in nonlocalIDs_
  int offset = Epetra_Util_binary_search(GID, nonlocalIDs_, numNonlocalIDs_,
					 insertPoint);
  if (offset >= 0) {
    //if offset >= 0
    //  put value in nonlocalCoefs_[offset][0]

    if (numValues != nonlocalElementSize_[offset]) {
      cerr << "Epetra_FEVector ERROR: block-size for GID " << GID << " is "
	   << numValues<<" which doesn't match previously set block-size of "
	   << nonlocalElementSize_[offset] << endl;
      return(-1);
    }

    if (accumulate) {
      for(int j=0; j<numValues; ++j) {
	nonlocalCoefs_[offset][j] += values[j];
      }
    }
    else {
      for(int j=0; j<numValues; ++j) {
	nonlocalCoefs_[offset][j] = values[j];
      }
    }
  }
  else {
    //else
    //  insert GID in nonlocalIDs_
    //  insert numValues   in nonlocalElementSize_
    //  insert values in nonlocalCoefs_

    int tmp1 = numNonlocalIDs_;
    int tmp2 = allocatedNonlocalLength_;
    int tmp3 = allocatedNonlocalLength_;
    EPETRA_CHK_ERR( Epetra_Util_insert(GID, insertPoint, nonlocalIDs_,
				       tmp1, tmp2) );
    --tmp1;
    EPETRA_CHK_ERR( Epetra_Util_insert(numValues, insertPoint, nonlocalElementSize_,
				       tmp1, tmp3) );
    double* newvalues = new double[numValues];
    for(int j=0; j<numValues; ++j) {
      newvalues[j] = values[j];
    }
    EPETRA_CHK_ERR( Epetra_Util_insert(newvalues, insertPoint, nonlocalCoefs_,
				       numNonlocalIDs_, allocatedNonlocalLength_) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVector::GlobalAssemble()
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
  Epetra_MultiVector nonlocalVector(sourceMap, 1);

  int i,j;
  for(i=0; i<numNonlocalIDs_; ++i) {
    for(j=0; j<nonlocalElementSize_[i]; ++j) {
      nonlocalVector.ReplaceGlobalValue(nonlocalIDs_[i], j, 0,
					nonlocalCoefs_[i][j]);
    }
  }

  Epetra_Export exporter(sourceMap, Map());

  EPETRA_CHK_ERR( Export(nonlocalVector, exporter, Add) );

  //and finally, reset the nonlocalCoefs_ values to 0.0, in case the user
  //does some more data input followed by another globalAssemble.
  for(i=0; i<numNonlocalIDs_; ++i) {
    for(j=0; j<nonlocalElementSize_[i]; ++j) {
      nonlocalCoefs_[i][j] = 0.0;
    }
  }

  destroyNonlocalData();

  return(0);
}

//----------------------------------------------------------------------------
Epetra_FEVector& Epetra_FEVector::operator=(const Epetra_FEVector& source)
{
  Epetra_MultiVector::Assign(source);

  destroyNonlocalData();

  if (source.allocatedNonlocalLength_ > 0) {
    allocatedNonlocalLength_ = source.allocatedNonlocalLength_;
    numNonlocalIDs_ = source.numNonlocalIDs_;
    nonlocalIDs_ = new int[allocatedNonlocalLength_];
    nonlocalElementSize_ = new int[allocatedNonlocalLength_];
    nonlocalCoefs_ = new double*[allocatedNonlocalLength_];
    for(int i=0; i<numNonlocalIDs_; ++i) {
      int elemSize = source.nonlocalElementSize_[i];
      nonlocalCoefs_[i] = new double[elemSize];
      nonlocalIDs_[i] = source.nonlocalIDs_[i];
      nonlocalElementSize_[i] = elemSize;
      for(int j=0; j<elemSize; ++j) {
	nonlocalCoefs_[i][j] = source.nonlocalCoefs_[i][j];
      }
    }
  }

  return(*this);
}

//----------------------------------------------------------------------------
void Epetra_FEVector::destroyNonlocalData()
{
  if (allocatedNonlocalLength_ > 0) {
    delete [] nonlocalIDs_;
    delete [] nonlocalElementSize_;
    nonlocalIDs_ = NULL;
    nonlocalElementSize_ = NULL;
    for(int i=0; i<numNonlocalIDs_; ++i) {
      delete [] nonlocalCoefs_[i];
    }
    delete [] nonlocalCoefs_;
    nonlocalCoefs_ = NULL;
    numNonlocalIDs_ = 0;
    allocatedNonlocalLength_ = 0;
  }
  return;
}
