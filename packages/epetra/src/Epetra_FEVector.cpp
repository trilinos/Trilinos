
#include "Epetra_FEVector.h"

#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Util.h"

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
  if (allocatedNonlocalLength_ > 0) {
    delete [] nonlocalIDs_;
    delete [] nonlocalCoefs_;
  }
}

//----------------------------------------------------------------------------
int Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int* GIDs,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, values, true) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int* GIDs,
			                 const double* values)
{
  return( inputValues( numIDs, GIDs, values, false) );
}

//----------------------------------------------------------------------------
int Epetra_FEVector::inputValues(int numIDs,
                                 const int* GIDs,
                                 const double* values,
                                 bool accumulate)
{
 //Important note!! This method currently assumes that there is only 1 point
 //associated with each element. FIX THIS

  for(int i=0; i<numIDs; ++i) {
    if (Map().MyGID(GIDs[i])) {
      if (accumulate) {
        SumIntoGlobalValue(GIDs[i], 0, values[i]);
      }
      else {
        ReplaceGlobalValue(GIDs[i], 0, values[i]);
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
int Epetra_FEVector::inputNonlocalValue(int GID, double value, bool accumulate)
{
  int insertPoint = -1;

  //find offset of GID in nonlocalIDs_
  int offset = Epetra_Util_binary_search(GID, nonlocalIDs_, numNonlocalIDs_,
					 insertPoint);
  if (offset >= 0) {
    //if offset >= 0
    //  put value in nonlocalCoefs_[offset]

    if (accumulate) {
      nonlocalCoefs_[offset] += value;
    }
    else {
      nonlocalCoefs_[offset] = value;
    }
  }
  else {
    //else
    //  insert GID in nonlocalIDs_
    //  insert value in nonlocalCoefs_

    int tmp1 = numNonlocalIDs_;
    int tmp2 = allocatedNonlocalLength_;
    EPETRA_CHK_ERR( Epetra_Util_insert(GID, insertPoint, nonlocalIDs_,
				       tmp1, tmp2) );
    EPETRA_CHK_ERR( Epetra_Util_insert(value, insertPoint, nonlocalCoefs_,
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

  //First build a map that describes the data in nonlocalIDs_/nonlocalCoefs_.
  //We'll use the arbitrary distribution constructor of Map.

  //(We don't need to do anything if there's only one processor.)
  if (Map().Comm().NumProc() < 2) {
    return(0);
  }

  if (ignoreNonLocalEntries_) {
    return(0);
  }

  Epetra_BlockMap sourceMap(-1, numNonlocalIDs_,
                            nonlocalIDs_, 1, Map().IndexBase(), Map().Comm());

  //Now build a vector to hold our nonlocalCoefs_, and to act as the source-
  //vector for our import operation.
  Epetra_MultiVector nonlocalVector(View, sourceMap, nonlocalCoefs_,
				    numNonlocalIDs_, 1);

  Epetra_Export exporter(sourceMap, Map());

  EPETRA_CHK_ERR( Export(nonlocalVector, exporter, Add) );

  //and finally, reset the nonlocalCoefs_ values to 0.0, in case the user
  //does some more data input followed by another globalAssemble.
  for(int i=0; i<numNonlocalIDs_; ++i) {
    nonlocalCoefs_[i] = 0.0;
  }

  return(0);
}

