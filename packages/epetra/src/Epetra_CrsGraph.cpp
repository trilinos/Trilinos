
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
#include "Epetra_SerialDenseVector.h"
#include "Epetra_OffsetIndex.h"

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int* NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, RowMap))
{
  Allocate(NumIndicesPerRow, 1);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, RowMap))
{
  Allocate(&NumIndicesPerRow, 0);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
				 const Epetra_BlockMap& RowMap, 
				 const Epetra_BlockMap& ColMap, 
				 int* NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, RowMap, ColMap))
{
  Allocate(NumIndicesPerRow, 1);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, 
				 const Epetra_BlockMap& RowMap, 
				 const Epetra_BlockMap& ColMap, 
				 int NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    CrsGraphData_(new Epetra_CrsGraphData(CV, RowMap, ColMap))
{
  Allocate(&NumIndicesPerRow, 0);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph& Graph)
  : Epetra_DistObject(Graph),
    CrsGraphData_(Graph.CrsGraphData_)
{
  CrsGraphData_->IncrementReferenceCount();
}

// private =====================================================================
int Epetra_CrsGraph::Allocate(int* NumIndicesPerRow, int Inc) {
  int i;
  const int numMyBlockRows = CrsGraphData_->NumMyBlockRows_;

  //*** portions specific to ISDVs
  int errorcode = CrsGraphData_->NumIndicesPerRow_.Size(CrsGraphData_->NumMyBlockRows_);
  if(errorcode != 0) 
    throw ReportError("Error with NumIndicesPerRow_ allocation.", -99);
  errorcode = CrsGraphData_->NumAllocatedIndicesPerRow_.Size(CrsGraphData_->NumMyBlockRows_);
  if(errorcode != 0) 
    throw ReportError("Error with NumAllocatedIndicesPerRow_ allocation.", -99);

  if(CrsGraphData_->CV_ == Copy) {
    for(i = 0; i < numMyBlockRows; i++) {
      CrsGraphData_->NumAllocatedIndicesPerRow_[i] = NumIndicesPerRow[i*Inc];
      //CrsGraphData_->NumIndicesPerRow_[i] = 0; // already 0 from sizing
    }
  }
  else {
    //CrsGraphData_->NumAllocatedIndicesPerRow_[i] = 0; // already 0 from sizing
    //CrsGraphData_->NumIndicesPerRow_[i] = 0; // already 0 from sizing
  }
  //***

  CrsGraphData_->MaxNumIndices_ = 0;
  
  // Allocate and initialize entries if we are copying data
  if(CrsGraphData_->CV_ == Copy) {
    for(i = 0; i < numMyBlockRows; i++) {
      const int NumIndices = NumIndicesPerRow[i*Inc]; // Inc is either 0 or 1.  Allows use of same Allocate function
                                                      // for all constructors.
      if(NumIndices > 0) {
	Epetra_IntSerialDenseVector temp(NumIndices);
	CrsGraphData_->Indices_[i] = temp;
      }
      else {
	Epetra_IntSerialDenseVector temp(0);
	CrsGraphData_->Indices_[i] = temp;
      }

      CrsGraphData_->NumAllocatedIndicesPerRow_[i] = NumIndices;
      const int indexBaseMinusOne = IndexBase() - 1;
      int* ColIndices = CrsGraphData_->Indices_[i].Values();
      for(int j = 0; j < NumIndices; j++) {
	ColIndices[j] = indexBaseMinusOne; // Fill column indices with out-of-range values
      }
    }
  }	 
  else { // CV_ == View
    for(i = 0; i < numMyBlockRows; i++) {
      Epetra_IntSerialDenseVector temp(0);
      CrsGraphData_->Indices_[i] = temp;
    }
  }

  SetAllocated(true);
  CrsGraphData_->UpdateSidekick();

  return(0);
}

// private =====================================================================
int Epetra_CrsGraph::ReAllocate() {
  // Reallocate storage that was deleted in OptimizeStorage

  // Since NIPR is in Copy mode, NAIPR will become a Copy as well
  CrsGraphData_->NumAllocatedIndicesPerRow_ = CrsGraphData_->NumIndicesPerRow_;

  CrsGraphData_->StorageOptimized_ = false;

  return(0);
}

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
int Epetra_CrsGraph::InsertGlobalIndices(int Row, int NumIndices, int* Indices) {
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreGlobal(true);
  Row = LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR(InsertIndices(Row, NumIndices, Indices));

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::InsertMyIndices(int Row, int NumIndices, int* Indices) {

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert local indices into a global graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreLocal(true);
  EPETRA_CHK_ERR(InsertIndices(Row, NumIndices, Indices));

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// protected ===================================================================
int Epetra_CrsGraph::InsertIndices(int Row,
				   int NumIndices,
				   int* UserIndices)
{
  SetSorted(false); // No longer in sorted state.
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  int j;
  int ierr = 0;

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  if(CrsGraphData_->CV_ == View) {
    if(CrsGraphData_->Indices_[Row].Values() != 0) 
      ierr = 2; // This row has be defined already.  Issue warning.
    Epetra_IntSerialDenseVector temp(View, UserIndices, NumIndices);
    CrsGraphData_->Indices_[Row] = temp;
    CrsGraphData_->NumAllocatedIndicesPerRow_[Row] = NumIndices;
    CrsGraphData_->NumIndicesPerRow_[Row] = NumIndices;
  }
  else {
    // if HaveColMap_ is true, UserIndices is copied into a new array,
    // and then modified. The UserIndices pointer is updated to point to this 
    // new array. If HaveColMap_ is false, nothing is done. This way,
    // the same UserIndices pointer can be used later on regardless of whether
    // changes were made.
    Epetra_IntSerialDenseVector tempVector;
    if(CrsGraphData_->HaveColMap_) { //only insert indices in col map if defined
      tempVector.Size(NumIndices);
      int loc = 0;
      if(IndicesAreLocal()) {
        for(j = 0; j < NumIndices; ++j)
          if(CrsGraphData_->ColMap_.MyLID(UserIndices[j])) 
	    tempVector[loc++] = UserIndices[j];
      }
      else {
        for(j = 0; j < NumIndices; ++j)
          if(CrsGraphData_->ColMap_.MyGID(UserIndices[j])) 
	    tempVector[loc++] = UserIndices[j];
      }
      if(loc != NumIndices) 
	ierr = 2; //Some columns excluded
      NumIndices = loc;
      UserIndices = tempVector.Values();
    }

    int start = CrsGraphData_->NumIndicesPerRow_[Row];
    int stop = start + NumIndices;
    int NumAllocatedIndices = CrsGraphData_->NumAllocatedIndicesPerRow_[Row];
    if(stop > NumAllocatedIndices) {
      if(NumAllocatedIndices == 0) {
	Epetra_IntSerialDenseVector temp(NumIndices);
	CrsGraphData_->Indices_[Row] = temp;
      }
      else {
	ierr = 1; // Out of room.  Must delete and allocate more space...
	CrsGraphData_->Indices_[Row].Resize(stop);
      }
      CrsGraphData_->NumAllocatedIndicesPerRow_[Row] = stop;
    }
    
    CrsGraphData_->NumIndicesPerRow_[Row] = stop;
    int* RowIndices = CrsGraphData_->Indices_[Row].Values();
    for(j = start; j < stop; j++) {
      RowIndices[j] = UserIndices[j-start];
    }
  }
  CrsGraphData_->MaxNumIndices_ = EPETRA_MAX(CrsGraphData_->MaxNumIndices_, CrsGraphData_->NumIndicesPerRow_[Row]);
  CrsGraphData_->UpdateSidekick(Row);
  EPETRA_CHK_ERR(ierr);


  if(CrsGraphData_->ReferenceCount() > 1)
    return(1); // return 1 if data is shared
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::RemoveGlobalIndices(int Row, int NumIndices, int* Indices) {
  int j;
  int k;
  int ierr = 0;
  int Loc;

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove global indices from a filled graph

  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumCurrentIndices = CrsGraphData_->NumIndicesPerRow_[Row];
  
  for(j = 0; j < NumIndices; j++) {
    int Index = Indices[j];
    if(FindGlobalIndexLoc(Row,Index,j,Loc)) {
      for(k = Loc+1; k < NumCurrentIndices; k++) 
	CrsGraphData_->Indices_[Row][k-1] = CrsGraphData_->Indices_[Row][k];
      NumCurrentIndices--;
      CrsGraphData_->NumIndicesPerRow_[Row]--;
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
int Epetra_CrsGraph::RemoveMyIndices(int Row, int NumIndices, int* Indices) {
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
  
  for(j = 0; j < NumIndices; j++) {
    int Index = Indices[j];
    if(FindMyIndexLoc(Row,Index,j,Loc)) {
      for(k = Loc + 1; k < NumCurrentIndices; k++) 
	CrsGraphData_->Indices_[Row][k-1] = CrsGraphData_->Indices_[Row][k];
      NumCurrentIndices--;
      CrsGraphData_->NumIndicesPerRow_[Row]--;
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
int Epetra_CrsGraph::RemoveGlobalIndices(int Row) {
  int j;
  int ierr = 0;

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph
  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];
  CrsGraphData_->NumIndicesPerRow_[Row] = 0;
  
  const int indexBaseMinusOne = IndexBase() - 1;
  for(j = 0; j < NumIndices; j++) 
    CrsGraphData_->Indices_[Row][j] = indexBaseMinusOne; // Set to invalid 

  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::RemoveMyIndices(int Row)
{
  int ierr = 0;

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph
  if(CrsGraphData_->CV_ == View) 
    EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];
  CrsGraphData_->NumIndicesPerRow_[Row] = 0;
  
  for(int j = 0; j < NumIndices; j++) 
    CrsGraphData_->Indices_[Row][j] = -1; // Set to invalid 

  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// protected ===================================================================
bool Epetra_CrsGraph::FindGlobalIndexLoc(int LocalRow,
					 int Index,
					 int Start,
					 int& Loc) const
{
  int NumIndices = CrsGraphData_->NumIndicesPerRow_[LocalRow];
  int* Indices = CrsGraphData_->Indices_[LocalRow].Values();

  // If we have transformed the column indices, we must map this global Index to local
  if(CrsGraphData_->IndicesAreLocal_) {
    Index = LCID(Index);
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, Indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(Indices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

// protected ===================================================================
bool Epetra_CrsGraph::FindGlobalIndexLoc(int NumIndices,
					 const int* Indices,
					 int Index,
					 int Start,
					 int& Loc) const
{
  // If we have transformed the column indices, we must map this global Index to local
  if(CrsGraphData_->IndicesAreLocal_) {
    Index = LCID(Index);
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, Indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(Indices[j0] == Index) {
	Loc = j0;
	return(true);
      }
      j0++;
    }
  }
  return(false);
}

// protected ===================================================================
bool Epetra_CrsGraph::FindMyIndexLoc(int LocalRow,
				     int Index,
				     int Start,
				     int& Loc) const
{
  int NumIndices = CrsGraphData_->NumIndicesPerRow_[LocalRow];
  int* Indices = CrsGraphData_->Indices_[LocalRow].Values();

  if(!CrsGraphData_->IndicesAreLocal_) {
    throw ReportError("Epetra_CrsGraph::FindMyIndexLoc", -1);// Indices must be local
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, Indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(Indices[j0] == Index) {
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
				     const int* Indices,
				     int Index,
				     int Start,
				     int& Loc) const
{
  if(!CrsGraphData_->IndicesAreLocal_) {
    throw ReportError("Epetra_CrsGraph::FindMyIndexLoc", -1);// Indices must be local
  }

  if (CrsGraphData_->Sorted_) {
    int insertPoint;
    Loc = Epetra_Util_binary_search(Index, Indices, NumIndices, insertPoint);
    return( Loc > -1 );
  }
  else {
    int j, j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
    for(j = 0; j < NumIndices; j++) {
      if(j0 >= NumIndices) 
	j0 = 0; // wrap around
      if(Indices[j0] == Index) {
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
int Epetra_CrsGraph::FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap) {
  CrsGraphData_->DomainMap_ = DomainMap;
  CrsGraphData_->RangeMap_ = RangeMap;

  MakeIndicesLocal(DomainMap, RangeMap); // Convert indices to zero based on each processor
  SortIndices();  // Sort column entries from smallest to largest
  RemoveRedundantIndices(); // Get rid of any redundant index values
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
int Epetra_CrsGraph::TransformToLocal(const Epetra_BlockMap* DomainMap, const Epetra_BlockMap* RangeMap) {
  return(FillComplete(*DomainMap, *RangeMap));
}

// private =====================================================================
int Epetra_CrsGraph::ComputeGlobalConstants() {
  int i;
  int j;

  if(GlobalConstantsComputed()) 
    return(0);

  Epetra_IntSerialDenseVector tempvec(8); // Temp space

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
    Comm().MaxAll(&CrsGraphData_->MaxNumIndices_, &CrsGraphData_->GlobalMaxNumIndices_, 1);
    
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
    for(i = 0; i < numMyBlockRows; i++){
      int NumEntries = CrsGraphData_->NumIndicesPerRow_[i];
      int* Indices = CrsGraphData_->Indices_[i].Values();
      if(NumEntries > 0) {
	int CurNumNonzeros = 0;
	int RowDim = RowElementSizeList[i];
	for(j = 0; j < NumEntries; j++) {
	  int ColDim = ColElementSizeList[Indices[j]];
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

    CrsGraphData_->GlobalMaxNumIndices_ = tempvec[2];
    CrsGraphData_->GlobalMaxNumNonzeros_ = tempvec[3];
  }
  
  CrsGraphData_->NumGlobalRows_ = CrsGraphData_->RangeMap_.NumGlobalPoints();
  CrsGraphData_->NumGlobalCols_ = DomainMap().NumGlobalPoints();

  CrsGraphData_->GlobalConstantsComputed_ = true;
  
  return(0);
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
    int* const list = CrsGraphData_->Indices_[i].Values();
    int m = n/2;
    
    while(m > 0) {
      int max = n - m;
      for(int j = 0; j < max; j++) {
	for(int k = j; k >= 0; k-=m) {
	  if(list[k+m] >= list[k])
	    break;
	  int itemp = list[k+m];
	  list[k+m] = list[k];
	  list[k] = itemp;
	}
      }
      m = m/2;
    }
  }
  SetSorted(true);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::RemoveRedundantIndices() {
  int i;
  int j;
  int k;
  int ig;
  int jg, jl;

  if(NoRedundancies()) 
    return(0);
  if(!Sorted())
    EPETRA_CHK_ERR(-1);  // Must have sorted index set
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Indices must be local

  // For each row, remove column indices that are repeated.
  // Also, determine if graph is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortIndices was already called.
  
  CrsGraphData_->NumMyDiagonals_ = 0;
  CrsGraphData_->NumMyBlockDiagonals_ = 0;
  const int numMyBlockRows = NumMyBlockRows();
  for(i = 0; i < numMyBlockRows; i++) {
    bool diagfound = false;
    int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
    if(NumIndices > 0) {
      int* const Indices = CrsGraphData_->Indices_[i].Values();
      int j0 = 0;
      ig = GRID(i);
      jl = Indices[0];
      jg = GCID(jl);
      if(jl > i) CrsGraphData_->LowerTriangular_ = false;
      if(jl < i) CrsGraphData_->UpperTriangular_ = false;
      if(jg == ig) diagfound = true;
      for(j = 1; j < NumIndices; j++) {
	jl = Indices[j];
	jg = GCID(jl);
	if (jl > i) CrsGraphData_->LowerTriangular_ = false;
	if (jl < i) CrsGraphData_->UpperTriangular_ = false;
	if (jg == ig) diagfound = true;
	if (jl == Indices[j0]) { // Check if index is repeated
	  CrsGraphData_->NumIndicesPerRow_[i]--; // Decrement NumIndices count
	  CrsGraphData_->NumMyNonzeros_ --;
	  for(k = j; k < NumIndices-1; k++) 
	    Indices[k] = Indices[k+1]; // Shift indices
	  NumIndices--;
	  j--;
	}
	else 
	  j0=j; // Redefine comparison index value
      }
      if(diagfound) {
	CrsGraphData_->NumMyBlockDiagonals_++;
	CrsGraphData_->NumMyDiagonals_ += RowMap().ElementSize(i);
      }
    }
  }
	
  CrsGraphData_->NoDiagonal_ = (CrsGraphData_->NumMyBlockDiagonals_ == 0);

  SetNoRedundancies(true);

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

// private =====================================================================
int Epetra_CrsGraph::MakeColMap(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap) {
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

  int hashSize = DomainMap.NumMyElements();
  if (hashSize<5) hashSize = 5; // Make sure there are at least a few elements in this hash table
  Epetra_HashTable LocalGIDs(hashSize);
  Epetra_HashTable RemoteGIDs(hashSize);
  Epetra_HashTable RemoteGIDList(hashSize);

  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;
  const int numMyBlockRows = NumMyBlockRows();
  for(i = 0; i < numMyBlockRows; i++) {
    const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
    int* ColIndices = CrsGraphData_->Indices_[i].Values();
    for(j = 0; j < NumIndices; j++) {
      int GID = ColIndices[j];
      // Check if GID matches a row GID
      if(DomainMap.MyGID(GID)) {
	if(LocalGIDs.Get(GID) == -1) // This means its a new local GID
	  LocalGIDs.Add(GID, NumLocalColGIDs++);
      }
      else {
	if(RemoteGIDs.Get(GID) == -1) { // This means its a new remote GID
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.Add(NumRemoteColGIDs++, GID);
	}
      }
    }
  }

  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int NumMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  Epetra_IntSerialDenseVector ColIndices;
  if(NumMyBlockCols > 0) 
    ColIndices.Size(NumMyBlockCols);

  int* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList.Get(i); 

  int NLists = 1;
  Epetra_IntSerialDenseVector PIDList;
  Epetra_IntSerialDenseVector SizeList;
  int* RemoteSizeList = 0;
  bool DoSizes = !DomainMap.ConstantElementSize(); // If not constant element size, then we must exchange
      
  if(NumRemoteColGIDs > 0) 
    PIDList.Size(NumRemoteColGIDs);

  if(DoSizes) {
    if(NumMyBlockCols > 0) 
      SizeList.Size(NumMyBlockCols);
    RemoteSizeList = SizeList.Values() + NumLocalColGIDs;
    NLists++;
  }
  EPETRA_CHK_ERR(DomainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
      
  // Sort External column indices so that all columns coming from a given remote processor are contiguous

  Epetra_Util Util;
  int* SortLists[2]; // this array is allocated on the stack, and so we won't need to delete it.bb
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  Util.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, NLists, SortLists);

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise 
  // (2) We step through the GIDs of the DomainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if(NumLocalColGIDs == DomainMap.NumMyElements()) {
    DomainMap.MyGlobalElements(ColIndices.Values()); // Load Global Indices into first NumMyBlockCols elements column GID list
    if(DoSizes) 
      DomainMap.ElementSizeList(SizeList.Values()); // Load ElementSizeList too
  }
  else {
    int NumMyElements = DomainMap.NumMyElements();
    int* MyGlobalElements = DomainMap.MyGlobalElements();
    int* ElementSizeList = 0;
    if(DoSizes) 
      ElementSizeList = DomainMap.ElementSizeList();
    int NumLocalAgain = 0;
    for(i = 0; i < NumMyElements; i++) {
      int GID = MyGlobalElements[i];
      if(LocalGIDs.Get(GID) != -1) {
	if(DoSizes) 
	  SizeList[NumLocalAgain] = ElementSizeList[i];
	ColIndices[NumLocalAgain++] = GID;
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }


  // Make Column map with same element sizes as Domain map

  if(DomainMap.MaxElementSize() == 1) { // Simple map
    Epetra_Map temp(-1, NumMyBlockCols, ColIndices.Values(), DomainMap.IndexBase(), DomainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  else if(DomainMap.ConstantElementSize()) { // Constant Block size map
    Epetra_BlockMap temp(-1, NumMyBlockCols, ColIndices.Values(), DomainMap.MaxElementSize(),DomainMap.IndexBase(), DomainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  else { // Most general case where block size is variable.
    Epetra_BlockMap temp(-1, NumMyBlockCols, ColIndices.Values(), SizeList.Values(), DomainMap.IndexBase(), DomainMap.Comm());
    CrsGraphData_->ColMap_ = temp;
  }
  CrsGraphData_->HaveColMap_ = true;

  return(0);
}

// protected ===================================================================
int Epetra_CrsGraph::MakeIndicesLocal(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap) {
  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if(IndicesAreLocal() && IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-1); // Return error: Indices must not be both local and global

  MakeColMap(DomainMap, RangeMap); // If user has not prescribed column map, create one from indices
  const Epetra_BlockMap& colmap = ColMap();

  // Store number of local columns
  CrsGraphData_->NumMyCols_ = ColMap().NumMyPoints();
  CrsGraphData_->NumMyBlockCols_ = ColMap().NumMyElements();
  // Transform indices to local index space

  const int numMyBlockRows = NumMyBlockRows();

  if(IndicesAreGlobal()) {
    for(int i = 0; i < numMyBlockRows; i++) {
      const int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
      int* ColIndices = CrsGraphData_->Indices_[i].Values();
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
	
  SetIndicesAreLocal(true);
  SetIndicesAreGlobal(false);
	
  if(CrsGraphData_->ReferenceCount() > 1)
    return(1); // return 1 if data is shared
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::OptimizeStorage() {
  int i;
  int j;
  int NumIndices;
  const int numMyBlockRows = NumMyBlockRows();

  if(StorageOptimized()) 
    return(0); // Have we been here before?

  bool Contiguous = true; // Assume contiguous is true
  for(i = 1; i < numMyBlockRows; i++) {
    NumIndices = CrsGraphData_->NumIndicesPerRow_[i-1];
    int NumAllocateIndices = CrsGraphData_->NumAllocatedIndicesPerRow_[i-1];

    // Check if NumIndices is same as NumAllocatedIndices and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if((NumIndices != NumAllocateIndices) || 
       (CrsGraphData_->Indices_[i].Values() != CrsGraphData_->Indices_[i-1].Values() + NumIndices)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumIndices and NumAllocatedIndices
  //        for the last row could be different, but I don't think it matters.


  if((CrsGraphData_->CV_ == View) && !Contiguous) 
    return(1);  // This is user data, it's not contiguous and we can't make it so.

  if(CrsGraphData_->NumAllocatedIndicesPerRow_ .Values() != CrsGraphData_->NumIndicesPerRow_.Values())
    CrsGraphData_->NumAllocatedIndicesPerRow_.MakeViewOf(CrsGraphData_->NumIndicesPerRow_); 
  // This space is not needed after construction
  // Once things are contiguous, allocated equals used.

  CrsGraphData_->StorageOptimized_ = true; // We can do it

  if(Contiguous) 
    return(0); // Everything is done.  Return
	
  // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
  CrsGraphData_->NumMyNonzeros_ = 0;
  for(i = 0; i < numMyBlockRows; i++) 
    CrsGraphData_->NumMyNonzeros_ += CrsGraphData_->NumIndicesPerRow_[i];

  // Allocate one big integer array for all index values
  int errorcode = CrsGraphData_->All_Indices_.Size(CrsGraphData_->NumMyNonzeros_);
  if(errorcode != 0) 
    throw ReportError("Error with All_Indices_ allocation.", -99);

  // Set Indices_ to point into All_Indices_

  int* tmp = CrsGraphData_->All_Indices_.Values();
  for(i = 0; i < numMyBlockRows; i++) {
    int NumIndices = CrsGraphData_->NumIndicesPerRow_[i];
    int* ColIndices = CrsGraphData_->Indices_[i].Values();
    for(j = 0; j < NumIndices; j++) 
      tmp[j] = ColIndices[j];
    Epetra_IntSerialDenseVector tempVector(View, tmp, NumIndices);
    CrsGraphData_->Indices_[i] = tempVector;
    tmp += NumIndices; 	// tmp points to the offset in All_Indices_ where Indices_[i] starts.
  }

  SetIndicesAreContiguous(true); // Can no longer dynamically add indices
  CrsGraphData_->UpdateSidekick();

  if(CrsGraphData_->ReferenceCount() > 1)
    return(1);
  else
    return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractGlobalRowCopy(int Row, int LenOfIndices, int& NumIndices, int* Indices) const 
{
  int j;

  Row = LRID(Row); // Normalize row range

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];
  if(LenOfIndices < NumIndices) 
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices

  if(IndicesAreLocal())  
    for(j = 0; j < NumIndices; j++) 
      Indices[j] = GCID(CrsGraphData_->Indices_[Row][j]);
  else 
    for(j = 0; j < NumIndices; j++)
      Indices[j] = CrsGraphData_->Indices_[Row][j];
  
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractMyRowCopy(int Row, int LenOfIndices, int& NumIndices, int* Indices) const 
{
  int j;

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];
  if(LenOfIndices < NumIndices) 
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-3); // There are no local indices yet

  for(j = 0; j < NumIndices; j++)
    Indices[j] = CrsGraphData_->Indices_[Row][j];
  
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractGlobalRowView(int Row, int& NumIndices, int*& Indices) const 
{
  Row = LRID(Row); // Normalize row range

  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // There are no global indices

  NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];

  Indices = CrsGraphData_->Indices_[Row].Values();
  
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::ExtractMyRowView(int Row, int& NumIndices, int*& Indices) const 
{
  if(Row < 0 || Row >= NumMyBlockRows()) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // There are no local indices

  NumIndices = CrsGraphData_->NumIndicesPerRow_[Row];

  Indices = CrsGraphData_->Indices_[Row].Values();
	
  return(0);
}

//==============================================================================
int Epetra_CrsGraph::NumGlobalIndices(int Row) const {
  Row = LRID(Row);
  if(Row != -1) 
    return(CrsGraphData_->NumIndicesPerRow_[Row]);
  else 
    return(0); // No indices for this row on this processor
}

//==============================================================================
int Epetra_CrsGraph::NumAllocatedGlobalIndices(int Row) const {
  Row = LRID(Row);
  if(Row != -1) 
    return(CrsGraphData_->NumAllocatedIndicesPerRow_[Row]);
  else 
    return(0); // No indices allocated for this row on this processor
}

//==============================================================================
int Epetra_CrsGraph::ReplaceRowMap(const Epetra_BlockMap& newmap)
{
  if (RowMap().PointSameAs(newmap)) {
    Epetra_DistObject::Map_ = newmap;
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
int Epetra_CrsGraph::CopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
					     int NumSameIDs, 
					     int NumPermuteIDs, 
					     int* PermuteToLIDs,
					     int* PermuteFromLIDs,
                                             const Epetra_OffsetIndex * Indexor)
{
  int i;
  int j;
  int NumIndices;
  int FromRow;
  int ToRow;
  int MaxNumIndices = A.MaxNumEntries();
  Epetra_IntSerialDenseVector Indices;
  Epetra_SerialDenseVector Values;

  if(MaxNumIndices > 0) {
    Indices.Size(MaxNumIndices);
    Values.Size(MaxNumIndices); // Must extract values even though we discard them
  }

  const Epetra_Map& RowMap = A.RowMatrixRowMap();
  const Epetra_Map& ColMap = A.RowMatrixColMap();
  
  // Do copy first
  for(i = 0; i < NumSameIDs; i++) {
    ToRow = RowMap.GID(i);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, MaxNumIndices, NumIndices, Values.Values(), Indices.Values()));
    for(j = 0; j < NumIndices; j++) 
      Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    // Place into target graph.  
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices.Values());
    if(ierr < 0) EPETRA_CHK_ERR(ierr);
  }
  
  // Do local permutation next
  for(i = 0; i < NumPermuteIDs; i++) {
    FromRow = PermuteFromLIDs[i];
    ToRow = GRID(PermuteToLIDs[i]);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, MaxNumIndices, NumIndices, Values.Values(), Indices.Values()));
    for(j = 0; j < NumIndices; j++) 
      Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices.Values()); // Place into target graph.
    if(ierr < 0) EPETRA_CHK_ERR(ierr);
  }
  
  return(0);
}

// private =====================================================================
int Epetra_CrsGraph::CopyAndPermuteCrsGraph(const Epetra_CrsGraph& A,
					    int NumSameIDs, 
					    int NumPermuteIDs, 
					    int* PermuteToLIDs,
					    int* PermuteFromLIDs,
                                            const Epetra_OffsetIndex * Indexor)
{
  int i;
  int Row;
  int NumIndices;
  int* Indices = 0;
  int FromRow, ToRow;
  int MaxNumIndices = A.MaxNumIndices();
  Epetra_IntSerialDenseVector IndicesVector;

  if(MaxNumIndices > 0 && A.IndicesAreLocal()) {
    IndicesVector.Size(MaxNumIndices);
    Indices = IndicesVector.Values();
  }
  
  // Do copy first
  if(NumSameIDs > 0) {
    if(A.IndicesAreLocal()) {
      for(i = 0; i < NumSameIDs; i++) {
        Row = GRID(i);
        EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Indices));
        // Place into target graph.  
        int ierr = InsertGlobalIndices(Row, NumIndices, Indices); 
        if(ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for(i = 0; i < NumSameIDs; i++) {
        Row = GRID(i);
        EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumIndices, Indices));
        // Place into target graph.  
        int ierr = InsertGlobalIndices(Row, NumIndices, Indices); 
        if(ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }	
  }

  // Do local permutation next
  if(NumPermuteIDs > 0) {
    if(A.IndicesAreLocal()) {
      for(i = 0; i < NumPermuteIDs; i++) {
        FromRow = A.GRID(PermuteFromLIDs[i]);
        ToRow = GRID(PermuteToLIDs[i]);
        EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices));
        // Place into target graph.
        int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices); 
        if (ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for(i = 0; i < NumPermuteIDs; i++) {
        FromRow = A.GRID(PermuteFromLIDs[i]);
        ToRow = GRID(PermuteToLIDs[i]);
        EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumIndices, Indices));
        // Place into target graph.
        int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices); 
        if (ierr < 0) EPETRA_CHK_ERR(ierr); 
      }
    }
  }	
	
  return(0);
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
  int GlobalMaxNumIndices = 0;
  int TotalSendSize = 0;

  VarSizes = true;

  if(NumExportIDs <= 0) 
    return(0);

  try {
    const Epetra_CrsGraph& A = dynamic_cast<const Epetra_CrsGraph&>(Source);
    GlobalMaxNumIndices = A.GlobalMaxNumIndices();
    for( int i = 0; i < NumExportIDs; ++i )
    {
      Sizes[i] = (A.NumMyIndices( ExportLIDs[i] ) + 2);
      TotalSendSize += Sizes[i];
    }
  }
  catch(...) {
    try {
      const Epetra_RowMatrix& A = dynamic_cast<const Epetra_RowMatrix&>(Source);
      int MaxNumIndices = A.MaxNumEntries();
      A.Comm().MaxAll(&MaxNumIndices, &GlobalMaxNumIndices, 1);
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

  SizeOfPacket = sizeof(int); 
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
  int i;
  int NumIndices;
  int* Indices = 0;
  int FromRow;
  int* intptr;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int MaxNumIndices = A.MaxNumIndices();
  //if( MaxNumIndices ) Indices = new int[MaxNumIndices];

  intptr = (int*) Exports;
  for(i = 0; i < NumExportIDs; i++) {
    FromRow = A.GRID(ExportLIDs[i]);
    *intptr = FromRow;
    Indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices));
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += (NumIndices+2); // Point to next segment
  }

  //if( Indices ) delete [] Indices;
    
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
  int i;
  int j;
  int NumIndices;
  int* Indices = 0;
  int FromRow;
  int* intptr;
  Epetra_SerialDenseVector Values;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int MaxNumIndices = A.MaxNumEntries();
  if(MaxNumIndices > 0) {
    Values.Size(MaxNumIndices);
//    Indices = new int[MaxNumIndices];
  }
  const Epetra_Map& RowMap = A.RowMatrixRowMap();
  const Epetra_Map& ColMap = A.RowMatrixColMap();

  intptr = (int*) Exports;
  for(i = 0; i < NumExportIDs; i++) {
    FromRow = RowMap.GID(ExportLIDs[i]);
    *intptr = FromRow;
    Indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], MaxNumIndices, NumIndices, Values.Values(), Indices));
    for(j = 0; j < NumIndices; j++) Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += (NumIndices+2); // Point to next segment
  }

//  if( Indices ) delete [] Indices;
 
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
  if(NumImportIDs <= 0) 
    return(0);

  int NumIndices;
  int* Indices;
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
    Indices = intptr + 2; 
    // Insert indices
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices);
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
	int Row = GRID(i); // Get global row number
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
