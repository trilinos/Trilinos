/*
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
*/

#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"
#include "Epetra_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_HashTable.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LongLongSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_OffsetIndex.h"

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
        const Epetra_BlockMap& rowMap, 
        const int* numIndicesPerRow, bool staticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, rowMap, staticProfile))
{
  Allocate(numIndicesPerRow, 1, staticProfile);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
        const Epetra_BlockMap& rowMap, 
        int numIndicesPerRow, bool staticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, rowMap, staticProfile))
{
  Allocate(&numIndicesPerRow, 0, staticProfile);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
				 const Epetra_BlockMap& rowMap, 
				 const Epetra_BlockMap& colMap, 
				 const int* numIndicesPerRow, bool staticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, rowMap, colMap, staticProfile))
{
	if(!rowMap.GlobalIndicesMatch(colMap))
	   throw ReportError("Epetra_CrsGraph::Epetra_CrsGraph: cannot be called with different indices types for rowMap and colMap", -1);

  Allocate(numIndicesPerRow, 1, staticProfile);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
				 const Epetra_BlockMap& rowMap, 
				 const Epetra_BlockMap& colMap, 
				 int numIndicesPerRow, bool staticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, rowMap, colMap, staticProfile))
{
	if(!rowMap.GlobalIndicesMatch(colMap))
	   throw ReportError("Epetra_CrsGraph::Epetra_CrsGraph: cannot be called with different indices types for rowMap and colMap", -1);

  Allocate(&numIndicesPerRow, 0, staticProfile);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph& Graph)
  : Epetra_DistObject(Graph),
    CrsGraphData_(Graph.CrsGraphData_)
{
  CrsGraphData_->IncrementReferenceCount();
}

// private =====================================================================
template<typename int_type>
int Epetra_CrsGraph::TAllocate(const int* numIndicesPerRow, int Inc, bool staticProfile) {
  int i;
  const int numMyBlockRows = CrsGraphData_->NumMyBlockRows_;
  
  // Portions specific to ISDVs
  // Note that NumIndicesPerRow_ will be 1 element longer than needed.
  // This is because, if OptimizeStorage() is called, the storage associated with
  // NumIndicesPerRow_ will be reused implicitly for the IndexOffset_ array.
  // Although a bit fragile, for users who care about efficient memory allocation,
  // the time and memory fragmentation are important: Mike Heroux Feb 2005.
  int errorcode = CrsGraphData_->NumIndicesPerRow_.Size(numMyBlockRows+1); 
  if(errorcode != 0) 
    throw ReportError("Error with NumIndicesPerRow_ allocation.", -99);

  errorcode = CrsGraphData_->NumAllocatedIndicesPerRow_.Size(numMyBlockRows);
  if(errorcode != 0) 
    throw ReportError("Error with NumAllocatedIndicesPerRow_ allocation.", -99);

  int nnz = 0;
  if(CrsGraphData_->CV_ == Copy) {
    if (numIndicesPerRow != 0) {
      for(i = 0; i < numMyBlockRows; i++) {
	int nnzr = numIndicesPerRow[i*Inc];
	nnz += nnzr;
      }
    }
  }
  CrsGraphData_->NumMyEntries_ = nnz; // Define this now since we have it and might want to use it
  //***

  CrsGraphData_->MaxNumIndices_ = 0;

  Epetra_CrsGraphData::IndexData<int_type>& Data = CrsGraphData_->Data<int_type>();
  
  // Allocate and initialize entries if we are copying data
  if(CrsGraphData_->CV_ == Copy) {
    if (staticProfile) Data.All_Indices_.Size(nnz);
    int_type * all_indices = Data.All_Indices_.Values(); // First address of contiguous buffer
    for(i = 0; i < numMyBlockRows; i++) {
      const int NumIndices = numIndicesPerRow==0 ? 0 :numIndicesPerRow[i*Inc];
      const int indexBaseMinusOne = IndexBase() - 1;

      if(NumIndices > 0) {
	if (staticProfile) {
	  Data.Indices_[i] = all_indices;
	  all_indices += NumIndices;
    int_type* ColIndices = Data.Indices_[i];
    for(int j = 0; j < NumIndices; j++) 
      ColIndices[j] = indexBaseMinusOne; // Fill column indices with out-of-range values
	}
	else {
    // reserve memory in the STL vector, and then resize it to zero
    // again in order to signal the program that no data is in there
    // yet.
    Data.SortedEntries_[i].entries_.resize(NumIndices,
                 indexBaseMinusOne);
    Data.Indices_[i] = NumIndices > 0 ? &Data.SortedEntries_[i].entries_[0]: NULL;
    Data.SortedEntries_[i].entries_.resize(0);
	}
      }
      else {
	Data.Indices_[i] = 0;
      }

      CrsGraphData_->NumAllocatedIndicesPerRow_[i] = NumIndices;
    }
    if (staticProfile) assert(Data.All_Indices_.Values()+nnz==all_indices); // Sanity check
  }	 
  else { // CV_ == View
    for(i = 0; i < numMyBlockRows; i++) {
      Data.Indices_[i] = 0;
    }
  }

  SetAllocated(true);

  return(0);
}

int Epetra_CrsGraph::Allocate(const int* numIndicesPerRow, int Inc, bool staticProfile)
{
  if(RowMap().GlobalIndicesInt()) {
    return TAllocate<int>(numIndicesPerRow, Inc, staticProfile);
  }

  if(RowMap().GlobalIndicesLongLong()) {
    TAllocate<int>(numIndicesPerRow, Inc, staticProfile);
    TAllocate<long long>(numIndicesPerRow, Inc, staticProfile);
    return 0;
  }
  
  throw ReportError("Epetra_CrsGraph::Allocate: Internal error.", -1);
}


// private =====================================================================
/*
int Epetra_CrsGraph::ReAllocate() {
  // Reallocate storage that was deleted in OptimizeStorage

  // Since NIPR is in Copy mode, NAIPR will become a Copy as well
  CrsGraphData_->NumAllocatedIndicesPerRow_ = CrsGraphData_->NumIndicesPerRow_;

  CrsGraphData_->StorageOptimized_ = false;

  return(0);
}
*/
//==============================================================================
Epetra_CrsGraph::~Epetra_CrsGraph()
{
  CleanupData();
}

// private =====================================================================
void Epetra_CrsGraph::CleanupData() {
  if(CrsGraphData_ != 0) {
    CrsGraphData_->DecrementReferenceCount();
    if(CrsGraphData_->ReferenceCount() == 0) {
      delete CrsGraphData_;
      CrsGraphData_ = 0;
    }
  }
}

//==============================================================================
template<typename int_type>
int Epetra_CrsGraph::InsertGlobalIndices(int_type Row, int NumIndices, int_type* indices) {
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreGlobal(true);
  int locRow = LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR(InsertIndicesIntoSorted(locRow, NumIndices, indices));

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::InsertGlobalIndices(int Row, int NumIndices, int* indices) {

  if(RowMap().GlobalIndicesInt())
    return InsertGlobalIndices<int>(Row, NumIndices, indices);
  else
    throw ReportError("Epetra_CrsGraph::InsertGlobalIndices int version called for a graph that is not int.", -1);
}

//==============================================================================
int Epetra_CrsGraph::InsertGlobalIndices(long long Row, int NumIndices, long long* indices) {

  if(RowMap().GlobalIndicesLongLong())
    return InsertGlobalIndices<long long>(Row, NumIndices, indices);
  else
    throw ReportError("Epetra_CrsGraph::InsertGlobalIndices long long version called for a graph that is not long long.", -1);
}

//==============================================================================
int Epetra_CrsGraph::InsertMyIndices(int Row, int NumIndices, int* indices) {

  if(IndicesAreGlobal()) {
    EPETRA_CHK_ERR(-2); // Cannot insert local indices into a global graph
  }
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed

  if (CrsGraphData_->HaveColMap_) {
    SetIndicesAreLocal(true);
  }
  else {
     if (!IndicesAreLocal()) {
       EPETRA_CHK_ERR(-4);
     }
  }

  EPETRA_CHK_ERR(InsertIndicesIntoSorted(Row, NumIndices, indices));

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// protected ===================================================================
template<typename int_type>
int Epetra_CrsGraph::InsertIndices(int Row,
				   int NumIndices,
				   int_type* UserIndices)
{
  if (StorageOptimized()) EPETRA_CHK_ERR(-1); // Cannot insert into an optimized graph

  SetSorted(false); // No longer in sorted state.
  CrsGraphData_->NoRedundancies_ = false; // Redundancies possible.
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  int j;
  int ierr = 0;

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-2); // Not in Row range
    
  int& current_numAllocIndices = CrsGraphData_->NumAllocatedIndicesPerRow_[Row];
  int& current_numIndices = CrsGraphData_->NumIndicesPerRow_[Row];

  Epetra_CrsGraphData::IndexData<int_type>& Data = CrsGraphData_->Data<int_type>();

  if(CrsGraphData_->CV_ == View) {
    if(Data.Indices_[Row] != 0) 
      ierr = 2; // This row has been defined already.  Issue warning.
    Data.Indices_[Row] = UserIndices;
    current_numAllocIndices = NumIndices;
    current_numIndices = NumIndices;
  }
  else {
    // if HaveColMap_ is true, UserIndices is copied into a new array,
    // and then modified. The UserIndices pointer is updated to point to this 
    // new array. If HaveColMap_ is false, nothing is done. This way,
    // the same UserIndices pointer can be used later on regardless of whether
    // changes were made.
    if(CrsGraphData_->HaveColMap_) { //only insert indices in col map if defined
      if (CrsGraphData_->NumTempColIndices_ < NumIndices) {
        delete [] Data.TempColIndices_;
        Data.TempColIndices_ = new int_type[NumIndices];
        CrsGraphData_->NumTempColIndices_ = NumIndices;
      }
      int_type * tempIndices = Data.TempColIndices_;
      int loc = 0;
      if(IndicesAreLocal()) {
        for(j = 0; j < NumIndices; ++j)
          if(CrsGraphData_->ColMap_.MyLID(static_cast<int>(UserIndices[j]))) 
            tempIndices[loc++] = UserIndices[j];
      }
      else {
        for(j = 0; j < NumIndices; ++j) {
          const int Index = CrsGraphData_->ColMap_.LID(UserIndices[j]);
          if (Index > -1)
            tempIndices[loc++] = UserIndices[j];
        }

      }
      if(loc != NumIndices) 
        ierr = 2; //Some columns excluded
      NumIndices = loc;
      UserIndices = tempIndices;
    }

    int start = current_numIndices;
    int stop = start + NumIndices;
    if (CrsGraphData_->StaticProfile_) {
      if(stop > current_numAllocIndices)
        EPETRA_CHK_ERR(-2); // Cannot reallocate storage if graph created using StaticProfile
    }
    else {
      if (current_numAllocIndices > 0 && stop > current_numAllocIndices)
        ierr = 3;
      Data.SortedEntries_[Row].entries_.resize(stop, IndexBase() - 1);
      Data.Indices_[Row] = stop>0 ? &Data.SortedEntries_[Row].entries_[0] : NULL;

      current_numAllocIndices =  (int) Data.SortedEntries_[Row].entries_.capacity();    
    }

    current_numIndices = stop;
    int_type* RowIndices = Data.Indices_[Row]+start;
    for(j = 0; j < NumIndices; j++) {
      RowIndices[j] = UserIndices[j];
    }
  }

  if (CrsGraphData_->MaxNumIndices_ < current_numIndices) {
    CrsGraphData_->MaxNumIndices_ = current_numIndices;
  }
  EPETRA_CHK_ERR(ierr);


  if(CrsGraphData_->ReferenceCount() > 1)
    return(1); // return 1 if data is shared
  else
    return(0);
}

// =========================================================================
int Epetra_CrsGraph::InsertIndices(int Row,
				   int NumIndices,
				   int* UserIndices)
{
  if(RowMap().GlobalIndicesTypeValid())
    return InsertIndices<int>(Row, NumIndices, UserIndices);
  else
    throw ReportError("Epetra_CrsGraph::InsertIndices global index type unknown.", -1);
}

// =========================================================================
int Epetra_CrsGraph::InsertIndices(int Row,
				   int NumIndices,
				   long long* UserIndices)
{
  if(RowMap().GlobalIndicesLongLong())
    return InsertIndices<long long>(Row, NumIndices, UserIndices);
  else
    throw ReportError("Epetra_CrsGraph::InsertIndices long long version called for a graph that is not long long.", -1);
}

// =========================================================================
template<typename int_type>
int Epetra_CrsGraph::InsertIndicesIntoSorted(int Row,
              int NumIndices,
              int_type* UserIndices)
{
  // This function is only valid for COPY mode with non-static profile and
  // sorted entries. Otherwise, go to the other function.
  if (!CrsGraphData_->NoRedundancies_ || !CrsGraphData_->Sorted_  || 
      CrsGraphData_->StaticProfile_ || CrsGraphData_->CV_ == View ) 
    return InsertIndices(Row, NumIndices, UserIndices);

  if (StorageOptimized()) EPETRA_CHK_ERR(-1); // Cannot insert into an optimized graph

  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  int ierr = 0;

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-2); // Not in Row range
    
  int& current_numAllocIndices = CrsGraphData_->NumAllocatedIndicesPerRow_[Row];
  int& current_numIndices = CrsGraphData_->NumIndicesPerRow_[Row];

  Epetra_CrsGraphData::IndexData<int_type>& Data = CrsGraphData_->Data<int_type>();

  // if HaveColMap_ is true, UserIndices filters out excluded indices,
  // and then modified. The UserIndices pointer is updated to point to this 
  // new array. If HaveColMap_ is false, nothing is done. This way,
  // the same UserIndices pointer can be used later on regardless of whether
  // changes were made.
  if(CrsGraphData_->HaveColMap_) { //only insert indices in col map if defined
    if (CrsGraphData_->NumTempColIndices_ < NumIndices) {
      delete [] Data.TempColIndices_;
      Data.TempColIndices_ = new int_type[NumIndices];
      CrsGraphData_->NumTempColIndices_ = NumIndices;
    }
    int_type * tempIndices = Data.TempColIndices_;
    int loc = 0;
    if(IndicesAreLocal()) {
      for(int j = 0; j < NumIndices; ++j)
        if(CrsGraphData_->ColMap_.MyLID(static_cast<int>(UserIndices[j])))
          tempIndices[loc++] = UserIndices[j];
    }
    else {
      for(int j = 0; j < NumIndices; ++j) {
        if (CrsGraphData_->ColMap_.MyGID(UserIndices[j])) {
          tempIndices[loc++] = UserIndices[j];
        }
      }
    }
    if(loc != NumIndices) 
      ierr = 2; //Some columns excluded
    NumIndices = loc;
    UserIndices = tempIndices;
  }

  // for non-static profile, directly insert into a list that we always
  // keep sorted.
  Data.SortedEntries_[Row].AddEntries(NumIndices, UserIndices);
  current_numIndices = (int) Data.SortedEntries_[Row].entries_.size();
  current_numAllocIndices = (int) Data.SortedEntries_[Row].entries_.capacity();
  // reset the pointer to the respective data
  Data.Indices_[Row] = current_numIndices > 0 ? &Data.SortedEntries_[Row].entries_[0] : NULL;

  if (CrsGraphData_->MaxNumIndices_ < current_numIndices) {
    CrsGraphData_->MaxNumIndices_ = current_numIndices;
  }
  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1); // return 1 if data is shared
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::InsertIndicesIntoSorted(int Row,
              int NumIndices,
              int* UserIndices)
{
  if(RowMap().GlobalIndicesTypeValid())
    return InsertIndicesIntoSorted<int>(Row, NumIndices, UserIndices);
  else
    throw ReportError("Epetra_CrsGraph::InsertIndicesIntoSorted global index type unknown.", -1);
}

//==============================================================================
int Epetra_CrsGraph::InsertIndicesIntoSorted(int Row,
              int NumIndices,
              long long* UserIndices)
{
  if(RowMap().GlobalIndicesLongLong())
    return InsertIndicesIntoSorted<long long>(Row, NumIndices, UserIndices);
  else
    throw ReportError("Epetra_CrsGraph::InsertIndicesIntoSorted long long version called for a graph that is not long long.", -1);
}

//==============================================================================
template<typename int_type>
int Epetra_CrsGraph::RemoveGlobalIndices(int_type Row, int NumIndices, int_type* indices) {
  int j;
  int k;
  int ierr = 0;
  int Loc;

  if(IndicesAreContiguous() || StorageOptimized()) 
    EPETRA_CHK_ERR(-1); // Indices cannot be individually deleted and newed

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove global indices from a filled graph

  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  int locRow = LRID(Row); // Normalize row range
    
  if(locRow < 0 || locRow >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumCurrentIndices = CrsGraphData_->NumIndicesPerRow_[locRow];

  Epetra_CrsGraphData::IndexData<int_type>& Data = CrsGraphData_->Data<int_type>();

  for(j = 0; j < NumIndices; j++) {
    int_type Index = indices[j];
    if(FindGlobalIndexLoc(locRow,Index,j,Loc)) {
      for(k = Loc+1; k < NumCurrentIndices; k++) 
        Data.Indices_[locRow][k-1] = Data.Indices_[locRow][k];
      NumCurrentIndices--;
      CrsGraphData_->NumIndicesPerRow_[locRow]--;
      if (!CrsGraphData_->StaticProfile_)
        Data.SortedEntries_[locRow].entries_.pop_back();
      else
        Data.Indices_[locRow][NumCurrentIndices-1] = IndexBase() - 1;
    }
  }
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::RemoveGlobalIndices(int Row, int NumIndices, int* indices)
{
  if(RowMap().GlobalIndicesInt())
    return RemoveGlobalIndices<int>(Row, NumIndices, indices);
  else
    throw ReportError("Epetra_CrsGraph::RemoveGlobalIndices int version called for a graph that is not int.", -1);
}

//==============================================================================
int Epetra_CrsGraph::RemoveGlobalIndices(long long Row, int NumIndices, long long* indices)
{
  if(RowMap().GlobalIndicesLongLong())
    return RemoveGlobalIndices<long long>(Row, NumIndices, indices);
  else
    throw ReportError("Epetra_CrsGraph::RemoveGlobalIndices long long version called for a graph that is not long long.", -1);
}

//==============================================================================
int Epetra_CrsGraph::RemoveMyIndices(int Row, int NumIndices, int* indices) {

  if(IndicesAreContiguous() || StorageOptimized()) 
    EPETRA_CHK_ERR(-1); // Indices cannot be individually deleted and newed

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph

  int j;
  int k;
  int ierr = 0;
  int Loc;

  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumCurrentIndices = CrsGraphData_->NumIndicesPerRow_[Row];

  Epetra_CrsGraphData::IndexData<int>& Data = CrsGraphData_->Data<int>();

  for(j = 0; j < NumIndices; j++) {
    int Index = indices[j];
    if(FindMyIndexLoc(Row,Index,j,Loc)) {
      for(k = Loc + 1; k < NumCurrentIndices; k++) 
        Data.Indices_[Row][k-1] = Data.Indices_[Row][k];
      NumCurrentIndices--;
      CrsGraphData_->NumIndicesPerRow_[Row]--;
      if (!CrsGraphData_->StaticProfile_)
        Data.SortedEntries_[Row].entries_.pop_back();
      else
        Data.Indices_[Row][NumCurrentIndices-1] = IndexBase() - 1;
    }
  }
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
template<typename int_type>
int Epetra_CrsGraph::TRemoveGlobalIndices(long long Row) {
  int j;
  int ierr = 0;

  if(IndicesAreContiguous() || StorageOptimized()) 
    EPETRA_CHK_ERR(-1); // Indices cannot be individually deleted and newed

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph
  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  int locRow = LRID(Row); // Normalize row range
    
  if(locRow < 0 || locRow >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  Epetra_CrsGraphData::IndexData<int_type>& Data = CrsGraphData_->Data<int_type>();

  if (CrsGraphData_->StaticProfile_) {
    int NumIndices = CrsGraphData_->NumIndicesPerRow_[locRow];
  
    const int indexBaseMinusOne = IndexBase() - 1;
    for(j = 0; j < NumIndices; j++) 
      Data.Indices_[locRow][j] = indexBaseMinusOne; // Set to invalid 
  }
  else
    Data.SortedEntries_[locRow].entries_.resize(0);
 
  CrsGraphData_->NumIndicesPerRow_[locRow] = 0;


  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

int Epetra_CrsGraph::RemoveGlobalIndices(long long Row)
{
	if(RowMap().GlobalIndicesLongLong())
		return TRemoveGlobalIndices<long long>(Row);

	if(RowMap().GlobalIndicesInt())
		return TRemoveGlobalIndices<int>(Row);

	throw ReportError("Epetra_CrsGraph::RemoveGlobalIndices: Internal error.", -1);
}

//==============================================================================
int Epetra_CrsGraph::RemoveMyIndices(int Row)
{
  int ierr = 0;

  if(IndicesAreContiguous() || StorageOptimized()) 
    EPETRA_CHK_ERR(-1); // Indices cannot be individually deleted and newed

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph

  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  Epetra_CrsGraphData::IndexData<int>& Data = CrsGraphData_->Data<int>();

  if (CrsGraphData_->StaticProfile_) {
    int NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];
    for(int j = 0; j < NumIndices; j++) 
      Data.Indices_[Row][j] = -1; // Set to invalid 
  }
  else
    Data.SortedEntries_[Row].entries_.resize(0);

  CrsGraphData_->NumIndicesPerRow_[Row] = 0;

  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// protected ===================================================================
template<typename int_type>
bool Epetra_CrsGraph::FindGlobalIndexLoc(int LocalRow,
					 int_type Index,
					 int Start,
					 int& Loc) const
{
  int NumIndices = NumMyIndices(LocalRow);
  int_type* locIndices = TIndices<int_type>(LocalRow);

  // If we have transformed the column indices, we must map this global Index to local
  if(CrsGraphData_->IndicesAreLocal_) {
    Index = LCID(Index);
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, locIndices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(locIndices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

bool Epetra_CrsGraph::FindGlobalIndexLoc(int LocalRow,
					 int Index,
					 int Start,
					 int& Loc) const
{
  if(RowMap().GlobalIndicesInt())
	return FindGlobalIndexLoc<int>(LocalRow, Index, Start, Loc);
  else
	throw ReportError("Epetra_CrsGraph::FindGlobalIndexLoc int version called for a graph that is not int.", -1);
}

bool Epetra_CrsGraph::FindGlobalIndexLoc(int LocalRow,
					 long long Index,
					 int Start,
					 int& Loc) const
{
  if(RowMap().GlobalIndicesLongLong())
	return FindGlobalIndexLoc<long long>(LocalRow, Index, Start, Loc);
  else
	throw ReportError("Epetra_CrsGraph::FindGlobalIndexLoc long long version called for a graph that is not long long.", -1);
}

// protected ===================================================================
template<typename int_type>
bool Epetra_CrsGraph::FindGlobalIndexLoc(int NumIndices,
					 const int_type* indices,
					 int_type Index,
					 int Start,
					 int& Loc) const
{
  // If we have transformed the column indices, we must map this global Index to local
  if(CrsGraphData_->IndicesAreLocal_) {
    Index = LCID(Index);
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(indices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

bool Epetra_CrsGraph::FindGlobalIndexLoc(int NumIndices,
					 const int* indices,
					 int Index,
					 int Start,
					 int& Loc) const
{
  if(RowMap().GlobalIndicesInt())
	return FindGlobalIndexLoc<int>(NumIndices, indices, Index, Start, Loc);
  else
	throw ReportError("Epetra_CrsGraph::FindGlobalIndexLoc int version called for a graph that is not int.", -1);
}

bool Epetra_CrsGraph::FindGlobalIndexLoc(int NumIndices,
					 const long long* indices,
					 long long Index,
					 int Start,
					 int& Loc) const
{
  if(RowMap().GlobalIndicesLongLong())
	return FindGlobalIndexLoc<long long>(NumIndices, indices, Index, Start, Loc);
  else
	throw ReportError("Epetra_CrsGraph::FindGlobalIndexLoc long long version called for a graph that is not long long.", -1);
}

// protected ===================================================================
bool Epetra_CrsGraph::FindMyIndexLoc(int LocalRow,
				     int Index,
				     int Start,
				     int& Loc) const
{
  int NumIndices = NumMyIndices(LocalRow);
  int* locIndices = Indices(LocalRow);

  if(!CrsGraphData_->IndicesAreLocal_) {
    throw ReportError("Epetra_CrsGraph::FindMyIndexLoc", -1);// Indices must be local
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, locIndices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(locIndices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

// protected ===================================================================
bool Epetra_CrsGraph::FindMyIndexLoc(int NumIndices,
				     const int* indices,
				     int Index,
				     int Start,
				     int& Loc) const
{
  if(!CrsGraphData_->IndicesAreLocal_) {
    throw ReportError("Epetra_CrsGraph::FindMyIndexLoc", -1);// Indices must be local
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(indices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

//==============================================================================
int Epetra_CrsGraph::FillComplete() {
  EPETRA_CHK_ERR(FillComplete(RowMap(), RowMap()));
  return(0); // uses FillComplete(...)'s shared-ownership test. 
}

//==============================================================================
int Epetra_CrsGraph::FillComplete(const Epetra_BlockMap& domainMap, const Epetra_BlockMap& rangeMap) {
  if(!domainMap.GlobalIndicesMatch(rangeMap))
     throw ReportError("Epetra_CrsGraph::FillComplete: cannot be called with different indices types for domainMap and rangeMap", -1);

  if(!RowMap().GlobalIndicesMatch(domainMap))
    throw ReportError("Epetra_CrsGraph::FillComplete: cannot be called with different indices types for row map and incoming rangeMap", -1);

  CrsGraphData_->DomainMap_ = domainMap;
  CrsGraphData_->RangeMap_ = rangeMap;

  MakeIndicesLocal(domainMap, rangeMap); // Convert global indices to local indices to on each processor
  SortIndices();  // Sort column entries from smallest to largest
  RemoveRedundantIndices(); // Get rid of any redundant index values
  DetermineTriangular(); //determine if lower/upper triangular
  CrsGraphData_->MakeImportExport(); // Build Import or Export objects
  ComputeGlobalConstants(); // Compute constants that require communication
  SetFilled(true);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::TransformToLocal() {
  return(FillComplete());
}

//==============================================================================
int Epetra_CrsGraph::TransformToLocal(const Epetra_BlockMap* domainMap, const Epetra_BlockMap* rangeMap) {
  return(FillComplete(*domainMap, *rangeMap));
}

// private =====================================================================
int Epetra_CrsGraph::ComputeGlobalConstants()
{
  if(GlobalConstantsComputed()) 
    return(0);

  Epetra_LongLongSerialDenseVector tempvec(8); // Temp space

  const int numMyBlockRows = NumMyBlockRows();

  CrsGraphData_->NumMyEntries_ = 0; // Compute Number of Nonzero entries and max
  CrsGraphData_->MaxNumIndices_ = 0;
  {for(int i = 0; i < numMyBlockRows; i++) {
    CrsGraphData_->NumMyEntries_ += CrsGraphData_->NumIndicesPerRow_[i];
    CrsGraphData_->MaxNumIndices_ = EPETRA_MAX(CrsGraphData_->MaxNumIndices_,CrsGraphData_->NumIndicesPerRow_[i]);
  }}
  
  // Case 1:  Constant block size (including blocksize = 1)
  if(RowMap().ConstantElementSize()) {
    tempvec[0] = CrsGraphData_->NumMyEntries_;
    tempvec[1] = CrsGraphData_->NumMyBlockDiagonals_;

    Comm().SumAll(&tempvec[0], &tempvec[2], 2);
    int tmp_MaxNumIndices = CrsGraphData_->MaxNumIndices_;
    Comm().MaxAll(&tmp_MaxNumIndices, &CrsGraphData_->GlobalMaxNumIndices_, 1);
    
    CrsGraphData_->NumGlobalEntries_ = tempvec[2];
    CrsGraphData_->NumGlobalBlockDiagonals_ = tempvec[3];

    int RowElementSize = RowMap().MaxElementSize();
    int ColElementSize = RowElementSize;
    CrsGraphData_->NumGlobalDiagonals_ = tempvec[3] * RowElementSize;
    CrsGraphData_->NumMyNonzeros_ = CrsGraphData_->NumMyEntries_ * RowElementSize * ColElementSize;
    CrsGraphData_->NumGlobalNonzeros_ = CrsGraphData_->NumGlobalEntries_ * RowElementSize * ColElementSize;
    CrsGraphData_->MaxNumNonzeros_ = CrsGraphData_->MaxNumIndices_ * RowElementSize * ColElementSize;
    CrsGraphData_->GlobalMaxNumNonzeros_ = CrsGraphData_->GlobalMaxNumIndices_ * RowElementSize * ColElementSize;
  }

  // Case 2:  Variable block size (more work)
  else {
    CrsGraphData_->NumMyNonzeros_ = 0;  // We will count the number of nonzeros here
    CrsGraphData_->MaxNumNonzeros_ = 0;  // We will determine the max number of nonzeros in any one block row
    int* RowElementSizeList = RowMap().ElementSizeList();
    int* ColElementSizeList = RowElementSizeList;
    if(Importer() != 0) 
      ColElementSizeList = ColMap().ElementSizeList();
    for(int i = 0; i < numMyBlockRows; i++){
      int NumEntries = CrsGraphData_->NumIndicesPerRow_[i];
      int* indices = CrsGraphData_->Data<int>().Indices_[i];
      if(NumEntries > 0) {
	int CurNumNonzeros = 0;
	int RowDim = RowElementSizeList[i];
	for(int j = 0; j < NumEntries; j++) {
	  int ColDim = ColElementSizeList[indices[j]];
	  CurNumNonzeros += RowDim*ColDim;
	  CrsGraphData_->MaxColDim_ = EPETRA_MAX(CrsGraphData_->MaxColDim_, ColDim);
	}
	CrsGraphData_->MaxNumNonzeros_ = EPETRA_MAX(CrsGraphData_->MaxNumNonzeros_, CurNumNonzeros);
	CrsGraphData_->NumMyNonzeros_ += CurNumNonzeros;
      }
    }
    
    // Sum Up all nonzeros
    
    tempvec[0] = CrsGraphData_->NumMyEntries_;
    tempvec[1] = CrsGraphData_->NumMyBlockDiagonals_;
    tempvec[2] = CrsGraphData_->NumMyDiagonals_;
    tempvec[3] = CrsGraphData_->NumMyNonzeros_;
    
    Comm().SumAll(&tempvec[0], &tempvec[4], 4);
    
    CrsGraphData_->NumGlobalEntries_ = tempvec[4];
    CrsGraphData_->NumGlobalBlockDiagonals_ = tempvec[5];
    CrsGraphData_->NumGlobalDiagonals_ = tempvec[6];
    CrsGraphData_->NumGlobalNonzeros_ = tempvec[7];

    tempvec[0] = CrsGraphData_->MaxNumIndices_;
    tempvec[1] = CrsGraphData_->MaxNumNonzeros_;

    Comm().MaxAll(&tempvec[0], &tempvec[2], 2);

    CrsGraphData_->GlobalMaxNumIndices_ = (int) tempvec[2];
    CrsGraphData_->GlobalMaxNumNonzeros_ = (int) tempvec[3];
  }
  
  CrsGraphData_->NumGlobalRows_ = CrsGraphData_->RangeMap_.NumGlobalPoints();
  CrsGraphData_->NumGlobalCols_ = DomainMap().NumGlobalPoints();

  CrsGraphData_->GlobalConstantsComputed_ = true;
  
  return(0);
}

void epetra_shellsort(int* list, int length)
{
  int i, j, j2, temp, istep;
  unsigned step;

  step = length/2;
  while (step > 0)
  {
    for (i=step; i < length; i++)
    {
      istep = step;
      j = i;
      j2 = j-istep;
      temp = list[i];
      if (list[j2] > temp) {
        while ((j >= istep) && (list[j2] > temp))
        {
          list[j] = list[j2];
          j = j2;
          j2 -= istep;
        }
        list[j] = temp;
      }
    }

    if (step == 2)
      step = 1;
    else
      step = (int) (step / 2.2);
  }
}

//==============================================================================
int Epetra_CrsGraph::SortIndices() {
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-1);
  if(Sorted())
    return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort, which is fast if indices are already sorted.

  const int numMyBlockRows = NumMyBlockRows();
  for(int i = 0; i < numMyBlockRows; i++){
    int n = CrsGraphData_->NumIndicesPerRow_[i];
    int* const list = CrsGraphData_->Data<int>().Indices_[i];

    epetra_shellsort(list, n);
//    int m = n/2;
    
//    while(m > 0) {
 //     int max = n - m;
//      for(int j = 0; j < max; j++) {
//        int k = j;
//        while(k>-1) {
//	  if(list[k+m] >= list[k])
//	    break;
//	  int itemp = list[k+m];
//	  list[k+m] = list[k];
//	  list[k] = itemp;
//          k-=m;
//	}
//      }
//      m = m/2;
//    }
  }
  SetSorted(true);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

void epetra_crsgraph_compress_out_duplicates(int len, int* list, int& newlen)
{
  //
  //This function runs the array ('list') checking for
  //duplicate entries. Any duplicates that are found are
  //removed by sliding subsequent data down in the array,
  //over-writing the duplicates. Finally, the new length
  //of the array (i.e., the number of unique entries) is
  //placed in the output parameter 'newlen'. The array is
  //**not** re-allocated.
  //
  //!*!*!*!
  //Important assumption: The array contents are assumed to
  //be sorted before this function is called. If the array
  //contents are not sorted, then the behavior of this
  //function is undefined.
  //!*!*!*!
  //

  if (len < 2) return;

  int* ptr0 = &list[0];
  int* ptr1 = &list[1];

  int* ptr_end = &list[len-1];

  while(*ptr0 != *ptr1 && ptr1 < ptr_end) {
    ++ptr0;
    ++ptr1;
  }

  if (ptr1 < ptr_end) {
    //if ptr1 < ptr_end we've found a duplicate...

    ++ptr0;
    ++ptr1;

    while(*ptr0 == *ptr1 && ptr1 < ptr_end) ++ptr1;

    while(ptr1 < ptr_end) {

      int val = *ptr1++;

      while(val == *ptr1 && ptr1 < ptr_end) {
        ++ptr1;
      }

      *ptr0++ = val;
    }

    if (*(ptr0-1) != *ptr1) *ptr0++ = *ptr1;

    int num_removed = (int)(ptr_end - ptr0 + 1);
    newlen = len - num_removed;
  }
  else {
    if (*ptr0 == *ptr1) newlen = len - 1;
    else newlen = len;
  }
}

//==============================================================================
int Epetra_CrsGraph::RemoveRedundantIndices()
{
  if(!Sorted())
    EPETRA_CHK_ERR(-1);  // Must have sorted index set
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Indices must be local

  // Note:  This function assumes that SortIndices was already called.
  // For each row, remove column indices that are repeated.

  const int numMyBlockRows = NumMyBlockRows();

  if(NoRedundancies() == false) {
    int* numIndicesPerRow = CrsGraphData_->NumIndicesPerRow_.Values();
    for(int i=0; i<numMyBlockRows; ++i) {
      int NumIndices = numIndicesPerRow[i];
      int* col_indices = this->Indices(i);
  
      if(NumIndices > 1) {
        epetra_crsgraph_compress_out_duplicates(NumIndices, col_indices,
                                                numIndicesPerRow[i]);
      }
      // update vector size and address in memory
      if (!CrsGraphData_->StaticProfile_) {
        CrsGraphData_->Data<int>().SortedEntries_[i].entries_.assign(col_indices, col_indices+numIndicesPerRow[i]);
        if (numIndicesPerRow[i] > 0) {
          CrsGraphData_->Data<int>().Indices_[i] = &(CrsGraphData_->Data<int>().SortedEntries_[i].entries_[0]);
        }
        else {
          CrsGraphData_->Data<int>().Indices_[i] = NULL;
        }
      }
    }
  }

  SetNoRedundancies(true);
  return 0;
}

//==============================================================================
int Epetra_CrsGraph::DetermineTriangular()
{
  // determine if graph is upper or lower triangular or has no diagonal
  
  if(!Sorted())
    EPETRA_CHK_ERR(-1);  // Must have sorted index set
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Indices must be local

  CrsGraphData_->NumMyDiagonals_ = 0;
  CrsGraphData_->NumMyBlockDiagonals_ = 0;

  const Epetra_BlockMap& rowMap = RowMap();
  const Epetra_BlockMap& colMap = ColMap();

  const int numMyBlockRows = NumMyBlockRows();

  for(int i = 0; i < numMyBlockRows; i++) {
    int NumIndices = NumMyIndices(i);
    if(NumIndices > 0) {
      long long ig = rowMap.GID(i);
      int* col_indices = this->Indices(i);

      int jl_0 = col_indices[0];
      int jl_n = col_indices[NumIndices-1];

      if(jl_n > i) CrsGraphData_->LowerTriangular_ = false;
      if(jl_0 < i) CrsGraphData_->UpperTriangular_ = false;

      //jl will be the local-index for the diagonal that we
      //want to search for.
      int jl = colMap.LID(ig);

      int insertPoint = -1;
      if (Epetra_Util_binary_search(jl, col_indices, NumIndices, insertPoint)>-1) {
        CrsGraphData_->NumMyBlockDiagonals_++;
        CrsGraphData_->NumMyDiagonals_ += rowMap.ElementSize(i);
      }
    }
  }

  CrsGraphData_->NoDiagonal_ = (CrsGraphData_->NumMyBlockDiagonals_ == 0);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// private =====================================================================
int Epetra_CrsGraph::MakeColMap_int(const Epetra_BlockMap& domainMap,
				const Epetra_BlockMap& rangeMap)
{
  (void)rangeMap;
  int i;
  int j;

  if(CrsGraphData_->HaveColMap_) 
    return(0); // Already have a Column Map

  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-1); // Return error: Indices must be global
  
  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  int numDomainElements = domainMap.NumMyElements();
  bool * LocalGIDs  = 0;
  if (numDomainElements>0) LocalGIDs  = new bool[numDomainElements];
  for (i=0; i<numDomainElements; i++) LocalGIDs[i] = false; // Assume domain GIDs are not local

  // In principle it is good to have RemoteGIDs and RemotGIDList be as long as the number of remote GIDs
  // on this processor, but this would require two passes through the column IDs, so we make it the max of 100
  // and the number of block rows.
  const int numMyBlockRows = NumMyBlockRows();
  int  hashsize = numMyBlockRows; if (hashsize < 100) hashsize = 100;
  //cout << "numMyBlockRows = " << numMyBlockRows << " hashsize = " << hashsize << endl;
  Epetra_HashTable RemoteGIDs(hashsize); 
  Epetra_HashTable RemoteGIDList(hashsize);

  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;
  for(i = 0; i < numMyBlockRows; i++) {
    const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
    int* ColIndices = CrsGraphData_->Data<int>().Indices_[i];
    for(j = 0; j < NumIndices; j++) {
      int GID = ColIndices[j];
      // Check if GID matches a row GID
      int LID = domainMap.LID(GID);
      if(LID != -1) {
	bool alreadyFound = LocalGIDs[LID];
	if (!alreadyFound) {
          LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
          NumLocalColGIDs++;
	}
      }
      else {
	if(RemoteGIDs.Get(GID) == -1) { // This means its a new remote GID
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.Add(NumRemoteColGIDs++, GID);
	}
      }
    }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (domainMap.Comm().NumProc()==1) { 
    
    if (NumRemoteColGIDs!=0) {
      throw ReportError("Some column IDs are not in domainMap.  If matrix is rectangular, you must pass in domainMap to FillComplete",-1); // Sanity test: When one processor,there can be no remoteGIDs
    }
    if (NumLocalColGIDs==numDomainElements) {
      CrsGraphData_->ColMap_ = domainMap;
      CrsGraphData_->HaveColMap_ = true;
      if (LocalGIDs!=0) delete [] LocalGIDs; 
      return(0); 
    }
  }
      
  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int numMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  Epetra_IntSerialDenseVector ColIndices;
  if(numMyBlockCols > 0) 
    ColIndices.Size(numMyBlockCols);

  int* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList.Get(i); 

  int NLists = 1;
  Epetra_IntSerialDenseVector PIDList;
  Epetra_IntSerialDenseVector SizeList;
  int* RemoteSizeList = 0;
  bool DoSizes = !domainMap.ConstantElementSize(); // If not constant element size, then we must exchange
      
  if(NumRemoteColGIDs > 0) 
    PIDList.Size(NumRemoteColGIDs);

  if(DoSizes) {
    if(numMyBlockCols > 0) 
      SizeList.Size(numMyBlockCols);
    RemoteSizeList = SizeList.Values() + NumLocalColGIDs;
    NLists++;
  }
  EPETRA_CHK_ERR(domainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
      
  // Sort External column indices so that all columns coming from a given remote processor are contiguous

  Epetra_Util Util;
  int* SortLists[2]; // this array is allocated on the stack, and so we won't need to delete it.bb
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  Util.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, NLists, SortLists, 0, 0);
  if (CrsGraphData_->SortGhostsAssociatedWithEachProcessor_) {
    // Sort external column indices so that columns from a given remote processor are not only contiguous
    // but also in ascending order. NOTE: I don't know if the number of externals associated
    // with a given remote processor is known at this point ... so I count them here.

    NLists--;
    int StartCurrent, StartNext;
    StartCurrent = 0; StartNext = 1;
    while ( StartNext < NumRemoteColGIDs ) {
      if ((PIDList.Values())[StartNext]==(PIDList.Values())[StartNext-1]) StartNext++;
      else {
        if(DoSizes) SortLists[0] = &(RemoteSizeList[StartCurrent]);
        Util.Sort(true,StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]),0,0,NLists,SortLists, 0, 0);
        StartCurrent = StartNext; StartNext++;
      }
    }
    if(DoSizes) SortLists[0] = &(RemoteSizeList[StartCurrent]);
    Util.Sort(true, StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]), 0, 0, NLists, SortLists, 0, 0);
  }

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise 
  // (2) We step through the GIDs of the domainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if(NumLocalColGIDs == domainMap.NumMyElements()) {
    domainMap.MyGlobalElements(ColIndices.Values()); // Load Global Indices into first numMyBlockCols elements column GID list
    if(DoSizes) 
      domainMap.ElementSizeList(SizeList.Values()); // Load ElementSizeList too
  }
  else {
    int NumMyElements = domainMap.NumMyElements();
    int* MyGlobalElements = domainMap.MyGlobalElements();
    int* ElementSizeList = 0;
    if(DoSizes) 
      ElementSizeList = domainMap.ElementSizeList();
    int NumLocalAgain = 0;
    for(i = 0; i < NumMyElements; i++) {
      if(LocalGIDs[i]) {
	if(DoSizes) 
	  SizeList[NumLocalAgain] = ElementSizeList[i];
	ColIndices[NumLocalAgain++] = MyGlobalElements[i];
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }

  // Done with this array
  if (LocalGIDs!=0) delete [] LocalGIDs; 


  // Make Column map with same element sizes as Domain map

  if(domainMap.MaxElementSize() == 1) { // Simple map
    Epetra_Map temp((int) -1, numMyBlockCols, ColIndices.Values(), domainMap.IndexBase(), domainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  else if(domainMap.ConstantElementSize()) { // Constant Block size map
    Epetra_BlockMap temp((int) -1, numMyBlockCols, ColIndices.Values(), domainMap.MaxElementSize(),domainMap.IndexBase(), domainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  else { // Most general case where block size is variable.
    Epetra_BlockMap temp((int) -1, numMyBlockCols, ColIndices.Values(), SizeList.Values(), domainMap.IndexBase(), domainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  CrsGraphData_->HaveColMap_ = true;

  return(0);
}

//==============================================================================
int Epetra_CrsGraph::MakeColMap_LL(const Epetra_BlockMap& domainMap,
				const Epetra_BlockMap& rangeMap)
{
  (void)rangeMap;
  int i;
  int j;

  if(CrsGraphData_->HaveColMap_) 
    return(0); // Already have a Column Map

  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-1); // Return error: Indices must be global
  
  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  int numDomainElements = domainMap.NumMyElements();
  bool * LocalGIDs  = 0;
  if (numDomainElements>0) LocalGIDs  = new bool[numDomainElements];
  for (i=0; i<numDomainElements; i++) LocalGIDs[i] = false; // Assume domain GIDs are not local

  // In principle it is good to have RemoteGIDs and RemotGIDList be as long as the number of remote GIDs
  // on this processor, but this would require two passes through the column IDs, so we make it the max of 100
  // and the number of block rows.
  const int numMyBlockRows = NumMyBlockRows();
  int  hashsize = numMyBlockRows; if (hashsize < 100) hashsize = 100;
  //cout << "numMyBlockRows = " << numMyBlockRows << " hashsize = " << hashsize << endl;
  Epetra_HashTable RemoteGIDs(hashsize); 
  Epetra_HashTable RemoteGIDList(hashsize);

  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;

  if(IndicesAreLocal())
  {
	  Epetra_CrsGraphData::IndexData<int>& intData = CrsGraphData_->Data<int>();

	  for(i = 0; i < numMyBlockRows; i++) {
		const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
		int* ColIndices = intData.Indices_[i];
		for(j = 0; j < NumIndices; j++) {
		  int GID = ColIndices[j];
		  // Check if GID matches a row GID
		  int LID = domainMap.LID(GID);
		  if(LID != -1) {
		bool alreadyFound = LocalGIDs[LID];
		if (!alreadyFound) {
			  LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
			  NumLocalColGIDs++;
		}
		  }
		  else {
		if(RemoteGIDs.Get(GID) == -1) { // This means its a new remote GID
		  RemoteGIDs.Add(GID, NumRemoteColGIDs);
		  RemoteGIDList.Add(NumRemoteColGIDs++, GID);
		}
		  }
		}
	  }
  }
  else if(IndicesAreGlobal())
  {
	  Epetra_CrsGraphData::IndexData<long long>& LLData = CrsGraphData_->Data<long long>();

	  for(i = 0; i < numMyBlockRows; i++) {
		const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
		long long* ColIndices = LLData.Indices_[i];
		for(j = 0; j < NumIndices; j++) {
		  long long GID = ColIndices[j];
		  // Check if GID matches a row GID
		  int LID = domainMap.LID(GID);
		  if(LID != -1) {
		bool alreadyFound = LocalGIDs[LID];
		if (!alreadyFound) {
			  LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
			  NumLocalColGIDs++;
		}
		  }
		  else {
		if(RemoteGIDs.Get(GID) == -1) { // This means its a new remote GID
		  RemoteGIDs.Add(GID, NumRemoteColGIDs);
		  RemoteGIDList.Add(NumRemoteColGIDs++, GID);
		}
		  }
		}
	  }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (domainMap.Comm().NumProc()==1) { 
    
    if (NumRemoteColGIDs!=0) {
      throw ReportError("Some column IDs are not in domainMap.  If matrix is rectangular, you must pass in domainMap to FillComplete",-1); // Sanity test: When one processor,there can be no remoteGIDs
    }
    if (NumLocalColGIDs==numDomainElements) {
      CrsGraphData_->ColMap_ = domainMap;
      CrsGraphData_->HaveColMap_ = true;
      if (LocalGIDs!=0) delete [] LocalGIDs; 
      return(0); 
    }
  }
      
  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int numMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  Epetra_LongLongSerialDenseVector ColIndices;
  if(numMyBlockCols > 0) 
    ColIndices.Size(numMyBlockCols);

  long long* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList.Get(i); 

  int NLists = 1;
  Epetra_IntSerialDenseVector PIDList;
  Epetra_IntSerialDenseVector SizeList;
  int* RemoteSizeList = 0;
  bool DoSizes = !domainMap.ConstantElementSize(); // If not constant element size, then we must exchange
      
  if(NumRemoteColGIDs > 0) 
    PIDList.Size(NumRemoteColGIDs);

  if(DoSizes) {
    if(numMyBlockCols > 0) 
      SizeList.Size(numMyBlockCols);
    RemoteSizeList = SizeList.Values() + NumLocalColGIDs;
    NLists++;
  }
  EPETRA_CHK_ERR(domainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
      
  // Sort External column indices so that all columns coming from a given remote processor are contiguous

  Epetra_Util Util;
  //int* SortLists[2]; // this array is allocated on the stack, and so we won't need to delete it.bb
  //SortLists[0] = RemoteColIndices;
  //SortLists[1] = RemoteSizeList;
  Util.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, 1, &RemoteSizeList, 1, &RemoteColIndices);
  if (CrsGraphData_->SortGhostsAssociatedWithEachProcessor_) {
    // Sort external column indices so that columns from a given remote processor are not only contiguous
    // but also in ascending order. NOTE: I don't know if the number of externals associated
    // with a given remote processor is known at this point ... so I count them here.

	int* SortLists[1];

    NLists--;
    int StartCurrent, StartNext;
    StartCurrent = 0; StartNext = 1;
    while ( StartNext < NumRemoteColGIDs ) {
      if ((PIDList.Values())[StartNext]==(PIDList.Values())[StartNext-1]) StartNext++;
      else {
        if(DoSizes) SortLists[0] = &(RemoteSizeList[StartCurrent]);
        Util.Sort(true,StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]),0,0,NLists,SortLists, 0, 0);
        StartCurrent = StartNext; StartNext++;
      }
    }
    if(DoSizes) SortLists[0] = &(RemoteSizeList[StartCurrent]);
    Util.Sort(true, StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]), 0, 0, NLists, SortLists, 0, 0);
  }

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise 
  // (2) We step through the GIDs of the domainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if(NumLocalColGIDs == domainMap.NumMyElements()) {
    domainMap.MyGlobalElements(ColIndices.Values()); // Load Global Indices into first numMyBlockCols elements column GID list
    if(DoSizes) 
      domainMap.ElementSizeList(SizeList.Values()); // Load ElementSizeList too
  }
  else {
    int NumMyElements = domainMap.NumMyElements();
    long long* MyGlobalElements = domainMap.MyGlobalElements_LL();
    int* ElementSizeList = 0;
    if(DoSizes) 
      ElementSizeList = domainMap.ElementSizeList();
    int NumLocalAgain = 0;
    for(i = 0; i < NumMyElements; i++) {
      if(LocalGIDs[i]) {
	if(DoSizes) 
	  SizeList[NumLocalAgain] = ElementSizeList[i];
	ColIndices[NumLocalAgain++] = MyGlobalElements[i];
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }

  // Done with this array
  if (LocalGIDs!=0) delete [] LocalGIDs; 


  // Make Column map with same element sizes as Domain map

  if(domainMap.MaxElementSize() == 1) { // Simple map
	Epetra_Map temp((long long) -1, numMyBlockCols, ColIndices.Values(), domainMap.IndexBase(), domainMap.Comm());
	CrsGraphData_->ColMap_ = temp;
  }
  else if(domainMap.ConstantElementSize()) { // Constant Block size map
	Epetra_BlockMap temp((long long) -1, numMyBlockCols, ColIndices.Values(), domainMap.MaxElementSize(),domainMap.IndexBase(), domainMap.Comm());
	CrsGraphData_->ColMap_ = temp;
  }
  else { // Most general case where block size is variable.
	Epetra_BlockMap temp((long long) -1, numMyBlockCols, ColIndices.Values(), SizeList.Values(), domainMap.IndexBase(), domainMap.Comm());
	CrsGraphData_->ColMap_ = temp;
  }
  CrsGraphData_->HaveColMap_ = true;

  return(0);
}

int Epetra_CrsGraph::MakeColMap(const Epetra_BlockMap& domainMap,
				const Epetra_BlockMap& rangeMap)
{
  if(!domainMap.GlobalIndicesMatch(rangeMap))
    throw ReportError("Epetra_CrsGraph::MakeColMap: cannot be called with different indices types for domainMap and rangeMap", -1);
  
  if(!RowMap().GlobalIndicesMatch(domainMap))
    throw ReportError("Epetra_CrsGraph::MakeColMap: cannot be called with different indices types for row map and incoming rangeMap", -1);
  
  if(RowMap().GlobalIndicesInt())
    return MakeColMap_int(domainMap, rangeMap);
  
  if(RowMap().GlobalIndicesLongLong())
    return MakeColMap_LL(domainMap, rangeMap);

  throw ReportError("Epetra_CrsGraph::MakeColMap: Internal error, unable to determine global index type of maps", -1);
}

// protected ===================================================================
int Epetra_CrsGraph::MakeIndicesLocal(const Epetra_BlockMap& domainMap, const Epetra_BlockMap& rangeMap) {
  if(!domainMap.GlobalIndicesMatch(rangeMap))
     throw ReportError("Epetra_CrsGraph::MakeIndicesLocal: cannot be called with different indices types for domainMap and rangeMap", -1);

  if(!RowMap().GlobalIndicesMatch(domainMap))
    throw ReportError("Epetra_CrsGraph::MakeIndicesLocal: cannot be called with different indices types for row map and incoming rangeMap", -1);

  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if(IndicesAreLocal() && IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-1); // Return error: Indices must not be both local and global

  MakeColMap(domainMap, rangeMap); // If user has not prescribed column map, create one from indices
  const Epetra_BlockMap& colmap = ColMap();

  // Store number of local columns
  CrsGraphData_->NumMyCols_ = ColMap().NumMyPoints();
  CrsGraphData_->NumMyBlockCols_ = ColMap().NumMyElements();

  // Transform indices to local index space
  const int numMyBlockRows = NumMyBlockRows();

  if(IndicesAreGlobal()) {
    // Check if ColMap is monotone. If not, the list will get unsorted.
    bool mapMonotone = true;
    {
      long long oldGID = colmap.GID(0);
      for (int i=1; i<colmap.NumMyElements(); ++i) {
        if (oldGID > colmap.GID(i)) {
          mapMonotone = false;
          break;
        }
        oldGID = colmap.GID(i);
      }
    }
    if (Sorted())
      SetSorted(mapMonotone);

	Epetra_CrsGraphData::IndexData<int>& intData = CrsGraphData_->Data<int>();

	if(RowMap().GlobalIndicesInt())
	{
		// now comes the actual transformation
		for(int i = 0; i < numMyBlockRows; i++) {
		  const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
		  int* ColIndices = intData.Indices_[i];
		  for(int j = 0; j < NumIndices; j++) {
			int GID = ColIndices[j];
			int LID = colmap.LID(GID);
			if(LID != -1) 
			  ColIndices[j] = LID;
			else 
			  throw ReportError("Internal error in FillComplete ",-1); 
		  }
		}
	}
	else if(RowMap().GlobalIndicesLongLong())
	{
		Epetra_CrsGraphData::IndexData<long long>& LL_Data = CrsGraphData_->Data<long long>();

		// now comes the actual transformation
		for(int i = 0; i < numMyBlockRows; i++) {
		  const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
		  long long* ColIndices = LL_Data.Indices_[i];
		  int* intColIndices = intData.Indices_[i];
		  for(int j = 0; j < NumIndices; j++) {
			long long GID = ColIndices[j];
			int LID = colmap.LID(GID);
			if(LID != -1) 
			  intColIndices[j] = LID;
			else 
			  throw ReportError("Internal error in FillComplete ",-1); 
		  }
		}

		LL_Data.Deallocate(); // deallocate long long data since indices are local now.
	}

  }
	
  SetIndicesAreLocal(true);
  SetIndicesAreGlobal(false);
	
  if(CrsGraphData_->ReferenceCount() > 1)
    return(1); // return 1 if data is shared
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::OptimizeStorage() {
  int NumIndices;
  const int numMyBlockRows = NumMyBlockRows();

  Epetra_CrsGraphData::IndexData<int>& Data = CrsGraphData_->Data<int>();

  if(StorageOptimized()) 
    return(0); // Have we been here before?
  if (!Filled()) EPETRA_CHK_ERR(-1); // Cannot optimize storage before calling FillComplete()

  bool Contiguous = true; // Assume contiguous is true
  for(int i = 1; i < numMyBlockRows; i++) {
    NumIndices = CrsGraphData_->NumIndicesPerRow_[i-1];
    int NumAllocateIndices = CrsGraphData_->NumAllocatedIndicesPerRow_[i-1];

    // Check if NumIndices is same as NumAllocatedIndices and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if((NumIndices != NumAllocateIndices) || 
       (Data.Indices_[i] != Data.Indices_[i-1] + NumIndices)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumIndices and NumAllocatedIndices
  //        for the last row could be different, but I don't think it matters.


  if((CrsGraphData_->CV_ == View) && !Contiguous) 
    return(3);  // This is user data, it's not contiguous and we can't make it so.

  // This next step constructs the scan sum of the number of indices per row.  Note that the
  // storage associated with NumIndicesPerRow is used by IndexOffset, so we have to be
  // careful with this sum operation

  if(CrsGraphData_->IndexOffset_.Values() != CrsGraphData_->NumIndicesPerRow_.Values())
    CrsGraphData_->IndexOffset_.MakeViewOf(CrsGraphData_->NumIndicesPerRow_);

  int * numIndicesPerRow = CrsGraphData_->NumIndicesPerRow_.Values();
  int curNumIndices = numIndicesPerRow[0];
  numIndicesPerRow[0] = 0;
  for (int i=0; i<numMyBlockRows; ++i) {
    int nextNumIndices = numIndicesPerRow[i+1];
    numIndicesPerRow[i+1] = numIndicesPerRow[i]+curNumIndices;
    curNumIndices = nextNumIndices;
  }

// *******************************
// Data NOT contiguous, make it so
// *******************************
  if(!Contiguous) { // Must pack indices if not already contiguous

    // Allocate one big integer array for all index values
    if (!(StaticProfile())) { // If static profile, All_Indices_ is already allocated, only need to pack data
      int errorcode = Data.All_Indices_.Size(CrsGraphData_->NumMyNonzeros_);
      if(errorcode != 0) throw ReportError("Error with All_Indices_ allocation.", -99);
    }
      // Pack indices into All_Indices_

    int* all_indices = Data.All_Indices_.Values();
    int * indexOffset = CrsGraphData_->IndexOffset_.Values();
    int ** indices = Data.Indices_;
		
    if (!(StaticProfile())) {
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(indexOffset,all_indices,indices)
#endif   
      for(int i = 0; i < numMyBlockRows; i++) {
	int numColIndices = indexOffset[i+1] - indexOffset[i];
        int* ColIndices = indices[i];
        int *newColIndices = all_indices+indexOffset[i];
        for(int j = 0; j < numColIndices; j++) newColIndices[j] = ColIndices[j];
      }
      for(int i = 0; i < numMyBlockRows; i++) {
        if (indices[i]!=0) {
          Data.SortedEntries_[i].entries_.clear();
          indices[i] = 0;
        }
     }
   } // End of non-contiguous non-static section	 
   else {

     for(int i = 0; i < numMyBlockRows; i++) {
       int numColIndices = indexOffset[i+1] - indexOffset[i];
       int* ColIndices = indices[i];
       int *newColIndices = all_indices+indexOffset[i];
       if (ColIndices!=newColIndices) // No need to copy if pointing to same space
       for(int j = 0; j < numColIndices; j++) newColIndices[j] = ColIndices[j];
       indices[i] = 0;
     }
   } // End of non-Contiguous static section
  } // End of non-Contiguous section
  else { // Start of Contiguous section
    // if contiguous, set All_Indices_ from CrsGraphData_->Indices_[0].
    // Execute the assignment block in parallel using the same pattern as SpMV
    // in order to improve page placement
    if (numMyBlockRows > 0 && !(StaticProfile())) {
      const int numMyNonzeros = NumMyNonzeros();
      int errorcode = Data.All_Indices_.Size(numMyNonzeros);
      if(errorcode != 0)  throw ReportError("Error with All_Indices_ allocation.", -99);
      int* new_all_indices = Data.All_Indices_.Values();
      int* old_all_indices = Data.Indices_[0];
      int * indexOffset = CrsGraphData_->IndexOffset_.Values();

#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(indexOffset,old_all_indices,new_all_indices)
#endif
     for(int i = 0; i < numMyBlockRows; i++) {
       int numColIndices = indexOffset[i+1] - indexOffset[i];
       int *oldColIndices = old_all_indices+indexOffset[i];
       int *newColIndices = new_all_indices+indexOffset[i];
       for(int j = 0; j < numColIndices; j++) newColIndices[j] = oldColIndices[j];
     }

//#ifdef EPETRA_HAVE_OMP
//#pragma omp parallel for default(none) shared(all_indices_values,indices_values)
//#endif   
//      for(int ii=0; ii<numMyNonzeros; ++ii) {
//        all_indices_values[ii] = indices_values[ii];
//      }
    }
  }

  // Delete unneeded storage
  CrsGraphData_->NumAllocatedIndicesPerRow_.Resize(0);
  delete [] Data.Indices_; Data.Indices_=0;
  Data.SortedEntries_.clear();

  SetIndicesAreContiguous(true); // Can no longer dynamically add or remove indices
  CrsGraphData_->StorageOptimized_ = true;

  return(0);
}

//==============================================================================
template<typename int_type>
int Epetra_CrsGraph::ExtractGlobalRowCopy(int_type Row, int LenOfIndices, int& NumIndices, int_type* targIndices) const 
{
  int j;

  int locRow = LRID(Row); // Normalize row range

  if(locRow < 0 || locRow >= NumMyBlockRows())
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = NumMyIndices(locRow);
  if(LenOfIndices < NumIndices)
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices

  if(IndicesAreLocal())
  {
    int * srcIndices = TIndices<int>(locRow);
    // static_cast is ok because global indices were created from int values and hence must fit ints
    for(j = 0; j < NumIndices; j++)
      targIndices[j] = static_cast<int_type>(GCID(srcIndices[j]));
  }
  else
  {
    int_type * srcIndices = TIndices<int_type>(locRow);
    for(j = 0; j < NumIndices; j++)
      targIndices[j] = srcIndices[j];
  }
  
  return(0);
}

int Epetra_CrsGraph::ExtractGlobalRowCopy(int Row, int LenOfIndices, int& NumIndices, int* targIndices) const 
{
  if(RowMap().GlobalIndicesInt())
    return ExtractGlobalRowCopy<int>(Row, LenOfIndices, NumIndices, targIndices);
  else
	throw ReportError("Epetra_CrsGraph::ExtractGlobalRowCopy int version called for a graph that is not int.", -1);
}

int Epetra_CrsGraph::ExtractGlobalRowCopy(long long Row, int LenOfIndices, int& NumIndices, long long* targIndices) const 
{
  if(RowMap().GlobalIndicesLongLong())
    return ExtractGlobalRowCopy<long long>(Row, LenOfIndices, NumIndices, targIndices);
  else
	throw ReportError("Epetra_CrsGraph::ExtractGlobalRowCopy long long version called for a graph that is not long long.", -1);
}

//==============================================================================
template<typename int_type>
int Epetra_CrsGraph::ExtractMyRowCopy(int Row, int LenOfIndices, int& NumIndices, int_type* targIndices) const
{
  int j;

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = NumMyIndices(Row);
  if(LenOfIndices < NumIndices) 
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-3); // There are no local indices yet

  int * srcIndices = TIndices<int>(Row);
  for(j = 0; j < NumIndices; j++)
    targIndices[j] = srcIndices[j];
  
  return(0);
}

int Epetra_CrsGraph::ExtractMyRowCopy(int Row, int LenOfIndices, int& NumIndices, int* targIndices) const
{
	if(RowMap().GlobalIndicesTypeValid())
		return ExtractMyRowCopy<int>(Row, LenOfIndices, NumIndices, targIndices);
	else
		throw ReportError("Epetra_CrsGraph::ExtractMyRowCopy graph global index type unknown.", -1);
}

//==============================================================================
int Epetra_CrsGraph::ExtractGlobalRowView(int Row, int& NumIndices, int*& targIndices) const 
{
  if(!RowMap().GlobalIndicesInt())
    throw ReportError("Epetra_CrsGraph::ExtractGlobalRowView int version called for a graph that is not int.", -1);

  int locRow = LRID(Row); // Normalize row range

  if(locRow < 0 || locRow >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // There are no global indices

  NumIndices = NumMyIndices(locRow);

  targIndices = TIndices<int>(locRow);
  
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractGlobalRowView(long long Row, int& NumIndices, long long*& targIndices) const 
{
  if(!RowMap().GlobalIndicesLongLong())
    throw ReportError("Epetra_CrsGraph::ExtractGlobalRowView long long version called for a graph that is not long long.", -1);

  int locRow = LRID(Row); // Normalize row range

  if(locRow < 0 || locRow >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // There are no global indices

  NumIndices = NumMyIndices(locRow);

  targIndices = TIndices<long long>(locRow);
  
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractMyRowView(int Row, int& NumIndices, int*& targIndices) const 
{
  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // There are no local indices

  NumIndices = NumMyIndices(Row);

  targIndices = TIndices<int>(Row);
	
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::NumGlobalIndices(long long Row) const {
  int locRow = LRID(Row);
  if(locRow != -1) 
    return(NumMyIndices(locRow));
  else 
    return(0); // No indices for this row on this processor
}

//==============================================================================
int Epetra_CrsGraph::NumAllocatedGlobalIndices(long long Row) const {
  int locRow = LRID(Row);
  if(locRow != -1) 
    return(NumAllocatedMyIndices(locRow));
  else 
    return(0); // No indices allocated for this row on this processor
}

//==============================================================================
int Epetra_CrsGraph::ReplaceRowMap(const Epetra_BlockMap& newmap)
{
  if (RowMap().PointSameAs(newmap)) {
    Epetra_DistObject::Map_ = newmap;
    CrsGraphData_->RowMap_ = newmap;
    CrsGraphData_->MakeImportExport();
    return(0);
  }

  return(-1);
}

//==============================================================================
int Epetra_CrsGraph::ReplaceColMap(const Epetra_BlockMap& newmap)
{
  if (!HaveColMap() && !IndicesAreLocal() && !IndicesAreGlobal() && newmap.GlobalIndicesMatch(RowMap())) {
    CrsGraphData_->ColMap_            = newmap;
    CrsGraphData_->NumGlobalBlockCols_= newmap.NumGlobalElements();
    CrsGraphData_->NumMyBlockCols_    = newmap.NumMyElements();
    CrsGraphData_->MaxColDim_         = newmap.MaxElementSize();
    CrsGraphData_->GlobalMaxColDim_   = newmap.MaxElementSize();
    CrsGraphData_->NumGlobalCols_     = newmap.NumGlobalPoints();
    CrsGraphData_->NumMyCols_         = newmap.NumMyPoints();
    CrsGraphData_->HaveColMap_        = true;
    return(0);
  }

  if(ColMap().PointSameAs(newmap)) {
    CrsGraphData_->ColMap_ = newmap;
    CrsGraphData_->MakeImportExport();
    return(0);
  }
  
  return(-1);
}

// private =====================================================================
int Epetra_CrsGraph::CheckSizes(const Epetra_SrcDistObject& Source) {
  try {
    const Epetra_CrsGraph& A = dynamic_cast<const Epetra_CrsGraph&>(Source); // downcast Source from SrcDistObject to CrsGraph
    if(!A.GlobalConstantsComputed()) 
      EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  }
  catch(...) {
    return(0); // No error at this point, object could be a RowMatrix
  }
  return(0);
}

  // TODO the functions below need to be changed for long long

// private =====================================================================
int Epetra_CrsGraph::CopyAndPermute(const Epetra_SrcDistObject& Source,
				    int NumSameIDs, 
				    int NumPermuteIDs, 
				    int* PermuteToLIDs,
				    int* PermuteFromLIDs,
                                    const Epetra_OffsetIndex * Indexor)
{ 
  try {
    const Epetra_CrsGraph& A = dynamic_cast<const Epetra_CrsGraph&>(Source);
    EPETRA_CHK_ERR(CopyAndPermuteCrsGraph(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					  PermuteFromLIDs,Indexor));
  }
  catch(...) {
    try {
      const Epetra_RowMatrix& A = dynamic_cast<const Epetra_RowMatrix&>(Source);
      EPETRA_CHK_ERR(CopyAndPermuteRowMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					     PermuteFromLIDs,Indexor));
    }
    catch(...) {
      EPETRA_CHK_ERR(-1); // Incompatible SrcDistObject
    }
  }
  
  return(0);
}

// private =====================================================================
template<typename int_type>
int Epetra_CrsGraph::CopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
					     int NumSameIDs, 
					     int NumPermuteIDs, 
					     int* PermuteToLIDs,
					     int* PermuteFromLIDs,
                                             const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor;
  int i;
  int j;
  int NumIndices;
  int FromRow;
  int_type ToRow;
  int maxNumIndices = A.MaxNumEntries();
  Epetra_IntSerialDenseVector local_indices_vec;
  Epetra_LongLongSerialDenseVector global_indices_vec;
  Epetra_SerialDenseVector Values;

  int* local_indices = 0;
  int_type* global_indices = 0;

  if(maxNumIndices > 0) {
    local_indices_vec.Size(maxNumIndices);
	local_indices = local_indices_vec.Values();

	if(A.RowMatrixRowMap().GlobalIndicesLongLong())
	{
		global_indices_vec.Size(maxNumIndices);
		global_indices = reinterpret_cast<int_type*>(global_indices_vec.Values());
	}
	else
	{
		global_indices = reinterpret_cast<int_type*>(local_indices);
	}

    Values.Size(maxNumIndices); // Must extract values even though we discard them
  }

  const Epetra_Map& rowMap = A.RowMatrixRowMap();
  const Epetra_Map& colMap = A.RowMatrixColMap();
  
  // Do copy first
  for(i = 0; i < NumSameIDs; i++) {
    ToRow = (int) rowMap.GID(i);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, maxNumIndices, NumIndices, Values.Values(), local_indices));
    for(j = 0; j < NumIndices; j++) 
      global_indices[j] = (int_type) colMap.GID(local_indices[j]); // convert to GIDs
    // Place into target graph.  
    int ierr = InsertGlobalIndices(ToRow, NumIndices, global_indices);
    if(ierr < 0) EPETRA_CHK_ERR(ierr);
  }
  
  // Do local permutation next
  for(i = 0; i < NumPermuteIDs; i++) {
    FromRow = PermuteFromLIDs[i];
    ToRow = (int_type) GRID(PermuteToLIDs[i]);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumIndices, NumIndices, Values.Values(), local_indices));
    for(j = 0; j < NumIndices; j++) 
      global_indices[j] = (int_type) colMap.GID(local_indices[j]); // convert to GIDs
    int ierr = InsertGlobalIndices(ToRow, NumIndices, global_indices); // Place into target graph.
    if(ierr < 0) EPETRA_CHK_ERR(ierr);
  }
  
  return(0);
}

int Epetra_CrsGraph::CopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
					     int NumSameIDs, 
					     int NumPermuteIDs, 
					     int* PermuteToLIDs,
					     int* PermuteFromLIDs,
                                             const Epetra_OffsetIndex * Indexor)
{
	if(A.RowMatrixRowMap().GlobalIndicesInt())
		return CopyAndPermuteRowMatrix<int>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);

	if(A.RowMatrixRowMap().GlobalIndicesLongLong())
		return CopyAndPermuteRowMatrix<long long>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);

	throw ReportError("Epetra_CrsGraph::CopyAndPermuteRowMatrix: Unable to determine global index type of map", -1);
}

// private =====================================================================
template<typename int_type>
int Epetra_CrsGraph::CopyAndPermuteCrsGraph(const Epetra_CrsGraph& A,
					    int NumSameIDs, 
					    int NumPermuteIDs, 
					    int* PermuteToLIDs,
					    int* PermuteFromLIDs,
                                            const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor;
  int i;
  int_type Row;
  int NumIndices;
  int_type* indices = 0;
  int_type FromRow, ToRow;
  int maxNumIndices = A.MaxNumIndices();
  Epetra_IntSerialDenseVector int_IndicesVector;
  Epetra_LongLongSerialDenseVector LL_IndicesVector;

  if(maxNumIndices > 0 && A.IndicesAreLocal()) {
	  if(A.RowMap().GlobalIndicesInt())
	  {
        int_IndicesVector.Size(maxNumIndices);
        indices = reinterpret_cast<int_type*>(int_IndicesVector.Values());
	  }
	  else if(A.RowMap().GlobalIndicesLongLong())
	  {
		LL_IndicesVector.Size(maxNumIndices);
	    indices = reinterpret_cast<int_type*>(LL_IndicesVector.Values());
	  }
  }

  // Do copy first
  if(NumSameIDs > 0) {
    if(A.IndicesAreLocal()) {
      for(i = 0; i < NumSameIDs; i++) {
        Row = (int_type) GRID(i);
        EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumIndices, NumIndices, indices));
        // Place into target graph.  
        int ierr = InsertGlobalIndices(Row, NumIndices, indices); 
        if(ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for(i = 0; i < NumSameIDs; i++) {
        Row = (int_type) GRID(i);
        EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumIndices, indices));
        // Place into target graph.  
        int ierr = InsertGlobalIndices(Row, NumIndices, indices); 
        if(ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }	
  }

  // Do local permutation next
  if(NumPermuteIDs > 0) {
    if(A.IndicesAreLocal()) {
      for(i = 0; i < NumPermuteIDs; i++) {
        FromRow = (int_type) A.GRID(PermuteFromLIDs[i]);
        ToRow = (int_type) GRID(PermuteToLIDs[i]);
        EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumIndices, NumIndices, indices));
        // Place into target graph.
        int ierr = InsertGlobalIndices(ToRow, NumIndices, indices); 
        if (ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for(i = 0; i < NumPermuteIDs; i++) {
        FromRow = (int_type) A.GRID(PermuteFromLIDs[i]);
        ToRow = (int_type) GRID(PermuteToLIDs[i]);
        EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumIndices, indices));
        // Place into target graph.
        int ierr = InsertGlobalIndices(ToRow, NumIndices, indices); 
        if (ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
  }	
	
  return(0);
}

int Epetra_CrsGraph::CopyAndPermuteCrsGraph(const Epetra_CrsGraph& A,
					    int NumSameIDs, 
					    int NumPermuteIDs, 
					    int* PermuteToLIDs,
					    int* PermuteFromLIDs,
                                            const Epetra_OffsetIndex * Indexor)
{
	if(A.RowMap().GlobalIndicesInt())
		return CopyAndPermuteCrsGraph<int>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);

	if(A.RowMap().GlobalIndicesLongLong())
		return CopyAndPermuteCrsGraph<long long>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);

	throw ReportError("Epetra_CrsGraph::CopyAndPermuteCrsGraph: Unable to determine global index type of map", -1);
}

// private =====================================================================
int Epetra_CrsGraph::PackAndPrepare(const Epetra_SrcDistObject& Source, 
				    int NumExportIDs, 
				    int* ExportLIDs,
				    int& LenExports, 
				    char*& Exports, 
				    int& SizeOfPacket, 
                                    int * Sizes,
                                    bool& VarSizes,
				    Epetra_Distributor& Distor) 
{
  int globalMaxNumIndices = 0;
  int TotalSendSize = 0;

  VarSizes = true;

  SizeOfPacket = (int)sizeof(int); 

  if(NumExportIDs <= 0) return(0);

  try {
    const Epetra_CrsGraph& A = dynamic_cast<const Epetra_CrsGraph&>(Source);
    globalMaxNumIndices = A.GlobalMaxNumIndices();
    for( int i = 0; i < NumExportIDs; ++i )
    {
      Sizes[i] = (A.NumMyIndices( ExportLIDs[i] ) + 2);
      TotalSendSize += Sizes[i];
    }
  }
  catch(...) {
    try {
      const Epetra_RowMatrix& A = dynamic_cast<const Epetra_RowMatrix&>(Source);
      int maxNumIndices = A.MaxNumEntries();
      A.Comm().MaxAll(&maxNumIndices, &globalMaxNumIndices, 1);
      for( int i = 0; i < NumExportIDs; ++i )
      {
        int NumEntries;
        A.NumMyRowEntries( ExportLIDs[i], NumEntries );
        Sizes[i] = (NumEntries + 2);
        TotalSendSize += Sizes[i];
      }
    }
    catch(...) {
      EPETRA_CHK_ERR(-1); // Bad cast
    }
  }

  CrsGraphData_->ReAllocateAndCast(Exports, LenExports, TotalSendSize*SizeOfPacket);

  try {
    const Epetra_CrsGraph& A = dynamic_cast<const Epetra_CrsGraph&>(Source);
    EPETRA_CHK_ERR(PackAndPrepareCrsGraph(A, NumExportIDs, ExportLIDs, LenExports, Exports,
					  SizeOfPacket, Sizes, VarSizes, Distor));
  }
  catch(...) {
    const Epetra_RowMatrix& A = dynamic_cast<const Epetra_RowMatrix&>(Source);
    EPETRA_CHK_ERR(PackAndPrepareRowMatrix(A, NumExportIDs, ExportLIDs, LenExports, Exports,
					   SizeOfPacket, Sizes, VarSizes, Distor));
  }
  return(0);
}

// private =====================================================================
int Epetra_CrsGraph::PackAndPrepareCrsGraph(const Epetra_CrsGraph& A, 
					    int NumExportIDs, 
					    int* ExportLIDs,
					    int& LenExports, 
					    char*& Exports, 
					    int& SizeOfPacket, 
                                            int* Sizes,
                                            bool& VarSizes,
					    Epetra_Distributor& Distor)
{
  (void)LenExports;
  (void)SizeOfPacket;
  (void)Sizes;
  (void)VarSizes;
  (void)Distor;
  int i;
  int NumIndices;
  int* indices = 0;
  int FromRow;
  int* intptr;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int maxNumIndices = A.MaxNumIndices();
  //if( maxNumIndices ) indices = new int[maxNumIndices];

  intptr = (int*) Exports;
  for(i = 0; i < NumExportIDs; i++) {
    FromRow = A.GRID(ExportLIDs[i]);
    *intptr = FromRow;
    indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumIndices, NumIndices, indices));
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += (NumIndices+2); // Point to next segment
  }

  //if( indices ) delete [] indices;
    
  return(0);
}

// private =====================================================================
int Epetra_CrsGraph::PackAndPrepareRowMatrix(const Epetra_RowMatrix& A, 
					     int NumExportIDs, 
					     int* ExportLIDs,
					     int& LenExports, 
					     char*& Exports, 
					     int& SizeOfPacket, 
                                             int* Sizes,
                                             bool& VarSizes,
					     Epetra_Distributor& Distor)
{
  (void)LenExports;
  (void)SizeOfPacket;
  (void)Sizes;
  (void)VarSizes;
  (void)Distor;
  int i;
  int j;
  int NumIndices;
  int* indices = 0;
  long long FromRow;
  int* intptr;
  Epetra_SerialDenseVector Values;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int maxNumIndices = A.MaxNumEntries();
  if(maxNumIndices > 0) {
    Values.Size(maxNumIndices);
//    indices = new int[maxNumIndices];
  }
  const Epetra_Map& rowMap = A.RowMatrixRowMap();
  const Epetra_Map& colMap = A.RowMatrixColMap();

  intptr = (int*) Exports;
  for(i = 0; i < NumExportIDs; i++) {
    FromRow = rowMap.GID(ExportLIDs[i]);
    *intptr = FromRow;
    indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumIndices, NumIndices, Values.Values(), indices));
    for(j = 0; j < NumIndices; j++) indices[j] = colMap.GID(indices[j]); // convert to GIDs
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += (NumIndices+2); // Point to next segment
  }

//  if( indices ) delete [] indices;
 
  return(0);
}

// private =====================================================================
int Epetra_CrsGraph::UnpackAndCombine(const Epetra_SrcDistObject& Source, 
                                      int NumImportIDs, 
                                      int* ImportLIDs, 
                                      int LenImports, 
                                      char* Imports, 
                                      int& SizeOfPacket, 
                                      Epetra_Distributor& Distor,
                                      Epetra_CombineMode CombineMode,
                                      const Epetra_OffsetIndex * Indexor) 
{
  (void)Source;
  (void)LenImports;
  (void)SizeOfPacket;
  (void)Distor;
  (void)CombineMode;
  (void)Indexor;
  if(NumImportIDs <= 0) 
    return(0);

  int NumIndices;
  int* indices;
  int ToRow;
  int i;
  
  int* intptr;
  // Unpack it...

  // Each segment of Sends will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.

  intptr = (int*) Imports;
    
  for(i = 0; i < NumImportIDs; i++) {
    ToRow = GRID(ImportLIDs[i]);
    assert((intptr[0])==ToRow); // Sanity check
    NumIndices = intptr[1];
    indices = intptr + 2; 
    // Insert indices
    int ierr = InsertGlobalIndices(ToRow, NumIndices, indices);
    if(ierr < 0) 
      EPETRA_CHK_ERR(ierr);
    intptr += (NumIndices+2); // Point to next segment
  }

  //destroy buffers since this operation is usually only done once
  if( LenExports_ ) {
    delete [] Exports_;
    Exports_ = 0;
    LenExports_ = 0;
  }
  if( LenImports_ ) {
    delete [] Imports_;
    Imports_ = 0;
    LenImports_ = 0;
  }
  
  return(0);
}

// protected ===================================================================
bool Epetra_CrsGraph::GlobalConstantsComputed() const {
  int mineComputed = 0;
  int allComputed;
  if(CrsGraphData_->GlobalConstantsComputed_) 
    mineComputed = 1;
  RowMap().Comm().MinAll(&mineComputed, &allComputed, 1); // Find out if any PEs changed constants info
  // If allComputed comes back zero then at least one PE need global constants recomputed.
  return(allComputed==1);
}

// private =====================================================================
void Epetra_CrsGraph::ComputeIndexState() {
  int myIndicesAreLocal = 0;
  int myIndicesAreGlobal = 0;
  if(CrsGraphData_->IndicesAreLocal_) 
    myIndicesAreLocal = 1;
  if(CrsGraphData_->IndicesAreGlobal_) 
    myIndicesAreGlobal = 1;
  int allIndicesAreLocal;
  int allIndicesAreGlobal;
  RowMap().Comm().MaxAll(&myIndicesAreLocal, &allIndicesAreLocal, 1); // Find out if any PEs changed Local Index info
  RowMap().Comm().MaxAll(&myIndicesAreGlobal, &allIndicesAreGlobal, 1); // Find out if any PEs changed Global Index info
  CrsGraphData_->IndicesAreLocal_ = (allIndicesAreLocal==1); // If indices are local on one PE, should be local on all
  CrsGraphData_->IndicesAreGlobal_ = (allIndicesAreGlobal==1);  // If indices are global on one PE should be local on all
}

//==============================================================================
void Epetra_CrsGraph::Print (ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for(int iproc = 0; iproc < NumProc; iproc++) {
    if(MyPID == iproc) {
      if(MyPID == 0) {
	os << "\nNumber of Global Block Rows  = " << NumGlobalBlockRows()      << endl;
	os <<   "Number of Global Block Cols  = " << NumGlobalBlockCols()      << endl;
	os <<   "Number of Global Block Diags = " << NumGlobalBlockDiagonals() << endl;
	os <<   "Number of Global Entries     = " << NumGlobalEntries()        << endl;
	os << "\nNumber of Global Rows        = " << NumGlobalRows()           << endl;
	os <<   "Number of Global Cols        = " << NumGlobalCols()           << endl;
	os <<   "Number of Global Diagonals   = " << NumGlobalDiagonals()      << endl;
	os <<   "Number of Global Nonzeros    = " << NumGlobalNonzeros()       << endl;
	os << "\nGlobal Maximum Block Row Dim = " << GlobalMaxRowDim()         << endl;
	os <<   "Global Maximum Block Col Dim = " << GlobalMaxColDim()         << endl;
	os <<   "Global Maximum Num Indices   = " << GlobalMaxNumIndices()     << endl;
	if(LowerTriangular()) os << " ** Matrix is Lower Triangular **"        << endl;
	if(UpperTriangular()) os << " ** Matrix is Upper Triangular **"        << endl;
	if(NoDiagonal())      os << " ** Matrix has no diagonal     **"        << endl << endl;
      }
      os << "\nNumber of My Block Rows  = " << NumMyBlockRows()      << endl;
      os <<   "Number of My Block Cols  = " << NumMyBlockCols()      << endl;
      os <<   "Number of My Block Diags = " << NumMyBlockDiagonals() << endl;
      os <<   "Number of My Entries     = " << NumMyEntries()        << endl;
      os << "\nNumber of My Rows        = " << NumMyRows()           << endl;
      os <<   "Number of My Cols        = " << NumMyCols()           << endl;
      os <<   "Number of My Diagonals   = " << NumMyDiagonals()      << endl;
      os <<   "Number of My Nonzeros    = " << NumMyNonzeros()       << endl;
      os << "\nMy Maximum Block Row Dim = " << MaxRowDim()           << endl;
      os <<   "My Maximum Block Col Dim = " << MaxColDim()           << endl;
      os <<   "My Maximum Num Indices   = " << MaxNumIndices()       << endl << endl;

      int NumMyBlockRows1 = NumMyBlockRows();
      int MaxNumIndices1 = MaxNumIndices();
      Epetra_IntSerialDenseVector Indices1(MaxNumIndices1);
      int NumIndices1;
      int i;
      int j;
      
      os.width(14);
      os <<  "       Row Index "; os << " ";
      for(j = 0; j < MaxNumIndices(); j++) {   
	os.width(12);
	os <<  "Col Index"; os << "      ";
      }
      os << endl;
      for(i = 0; i < NumMyBlockRows1; i++) {
	long long Row = GRID(i); // Get global row number
	ExtractGlobalRowCopy(Row, MaxNumIndices1, NumIndices1, Indices1.Values());
				
	os.width(14);
	os <<  Row ; os << "    ";	
	for(j = 0; j < NumIndices1 ; j++) {   
	  os.width(12);
	  os <<  Indices1[j]; os << "    ";
	}
	os << endl;
      }      
      os << flush;
    }
    // Do a few global ops to give I/O a chance to complete
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
  }
}

//==============================================================================
Epetra_CrsGraph& Epetra_CrsGraph::operator = (const Epetra_CrsGraph& Source) {
  if ((this == &Source) || (CrsGraphData_ == Source.CrsGraphData_))
    return(*this); // this and Source are same Graph object, or both point to same CrsGraphData object

  CleanupData();
  CrsGraphData_ = Source.CrsGraphData_;
  CrsGraphData_->IncrementReferenceCount();

  return(*this);
}
