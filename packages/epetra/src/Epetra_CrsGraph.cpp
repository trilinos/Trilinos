
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"
#include "Epetra_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_HashTable.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int *NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    ColMap_(0),
    CV_(CV) {
  InitializeDefaults();
  Allocate(NumIndicesPerRow, 1);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    ColMap_(0),
    CV_(CV) {
  InitializeDefaults();
  Allocate(&NumIndicesPerRow, 0);
}
//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
				 const Epetra_BlockMap& ColMap, int *NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    ColMap_(0),
    CV_(CV) {
  InitializeDefaults();
  Allocate(NumIndicesPerRow, 1);
  ColMap_ = new Epetra_BlockMap(ColMap);
}

//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
				 const Epetra_BlockMap& ColMap, int NumIndicesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsGraph"),
    ColMap_(0),
    CV_(CV) {
  InitializeDefaults();
  Allocate(&NumIndicesPerRow, 0);
  ColMap_ = new Epetra_BlockMap(ColMap);
}
//==============================================================================
Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph & Graph) 
  : Epetra_DistObject(Graph),
    DomainMap_(Graph.DomainMap_),
    RangeMap_(Graph.RangeMap_),
    ColMap_(Graph.ColMap_),
    Importer_(Graph.Importer_),
    Exporter_(Graph.Exporter_),
    Filled_(Graph.Filled_),
    Allocated_(false),
    Sorted_(Graph.Sorted_),
    StorageOptimized_(false),
    NoRedundancies_(Graph.NoRedundancies_),
    IndicesAreGlobal_(Graph.IndicesAreGlobal_),
    IndicesAreLocal_(Graph.IndicesAreLocal_),
    IndicesAreContiguous_(false),
    LowerTriangular_(Graph.LowerTriangular_),
    UpperTriangular_(Graph.UpperTriangular_),
    NoDiagonal_(Graph.NoDiagonal_),
    GlobalConstantsComputed_(false),
    IndexBase_(Graph.IndexBase_),
    NumGlobalEntries_(Graph.NumGlobalEntries_),
    NumGlobalBlockRows_(Graph.NumGlobalBlockRows_),
    NumGlobalBlockCols_(Graph.NumGlobalBlockCols_),
    NumGlobalBlockDiagonals_(Graph.NumGlobalBlockDiagonals_),
    NumMyEntries_(Graph.NumMyEntries_),
    NumMyBlockRows_(Graph.NumMyBlockRows_),
    NumMyBlockCols_(Graph.NumMyBlockCols_),
    NumMyBlockDiagonals_(Graph.NumMyBlockDiagonals_),
    MaxRowDim_(Graph.MaxRowDim_),
    MaxColDim_(Graph.MaxColDim_),
    GlobalMaxRowDim_(Graph.GlobalMaxRowDim_),
    GlobalMaxColDim_(Graph.GlobalMaxColDim_),
    MaxNumNonzeros_(Graph.MaxNumNonzeros_),
    GlobalMaxNumNonzeros_(Graph.GlobalMaxNumNonzeros_),
    NumGlobalNonzeros_(Graph.NumGlobalNonzeros_),
    NumGlobalRows_(Graph.NumGlobalRows_),
    NumGlobalCols_(Graph.NumGlobalCols_),
    NumGlobalDiagonals_(Graph.NumGlobalDiagonals_),
    NumMyNonzeros_(Graph.NumMyNonzeros_),
    NumMyRows_(Graph.NumMyRows_),
    NumMyCols_(Graph.NumMyCols_),
    NumMyDiagonals_(Graph.NumMyDiagonals_),
    Indices_(0),
    NumAllocatedIndicesPerRow_(0),
    NumIndicesPerRow_(0),
    MaxNumIndices_(Graph.MaxNumIndices_),
    GlobalMaxNumIndices_(Graph.GlobalMaxNumIndices_),
    All_Indices_(0),
    CV_(Copy)
{
  Allocate(Graph.NumIndicesPerRow(), 1);
  for (int i=0; i<NumMyBlockRows_; i++) {
    NumIndicesPerRow_[i] = NumAllocatedIndicesPerRow_[i];
    for (int j=0; j< NumIndicesPerRow_[i]; j++) Indices_[i][j] = Graph.Indices_[i][j];
  }
  MaxNumIndices_ = Graph.MaxNumIndices();
  if (ColMap_ != 0) ColMap_ = new Epetra_BlockMap(Graph.ColMap());
  if (Importer_ != 0) Importer_ = new Epetra_Import(*Graph.Importer()); // Non-trivial importer, must copy it.

  if (Exporter_ != 0) Exporter_ = new Epetra_Export(*Graph.Exporter()); // Non-trivial exporter, must copy it.
}

//==============================================================================
void Epetra_CrsGraph::InitializeDefaults() { // Initialize all attributes that have trivial default values

  IndexBase_ = RowMap().IndexBase();

  DomainMap_ = 0;
  RangeMap_ = 0;
  Importer_ = 0;
  Exporter_ = 0;

  Filled_ = false;
  Allocated_ = false;
  Sorted_ = false;
  StorageOptimized_ = false;
  NoRedundancies_ = false;
  IndicesAreGlobal_ = false;
  IndicesAreLocal_ = false;
  IndicesAreContiguous_ = false;

  LowerTriangular_ = true;
  UpperTriangular_ = true;
  NoDiagonal_ = true;

  NumGlobalBlockRows_ = RowMap().NumGlobalElements();
  NumGlobalBlockCols_ = NumGlobalBlockRows_;
  NumMyBlockRows_ = RowMap().NumMyElements();
  NumMyBlockCols_ = NumMyBlockRows_;

  NumGlobalRows_ = RowMap().NumGlobalPoints();
  NumGlobalCols_ = NumGlobalRows_;
  NumMyRows_ = RowMap().NumMyPoints();
  NumMyCols_ = NumMyRows_;

  GlobalMaxRowDim_ = RowMap().MaxElementSize();
  MaxRowDim_ = RowMap().MaxElementSize();
  GlobalMaxColDim_ = GlobalMaxRowDim_;
  MaxColDim_ = MaxRowDim_;

  GlobalConstantsComputed_ = false;

  NumGlobalEntries_ = 0;
  NumGlobalBlockDiagonals_ = 0;
  NumMyEntries_ = 0;
  NumMyBlockDiagonals_ = 0;

  NumGlobalNonzeros_ = 0;
  NumGlobalDiagonals_ = 0;
  NumMyNonzeros_ = 0;
  NumMyDiagonals_ = 0;
  
  MaxNumNonzeros_ = 0;
  GlobalMaxNumNonzeros_ = 0;
  MaxNumIndices_ = 0;
  GlobalMaxNumIndices_ = 0;

  All_Indices_ = 0;
}

//==============================================================================
int Epetra_CrsGraph::Allocate(int * NumIndicesPerRow, int Inc ) {

  int i, j;

  // Allocate Index pointers and sizes
  MaxNumIndices_ = 0;
  Indices_ = new int*[NumMyBlockRows_];
  NumIndicesPerRow_ = new int[NumMyBlockRows_];
  NumAllocatedIndicesPerRow_ = new int[NumMyBlockRows_];
  
  // Allocate and initialize entries if we are copying data
  if (CV_==Copy) {
    for (i=0; i<NumMyBlockRows_; i++) {
      const int NumIndices = NumIndicesPerRow[i*Inc]; // Inc is either 0 or 1.  Allows use of same Allocate function
                                                      // for all constructors.
      if (NumIndices > 0) Indices_[i] = new int[NumIndices];
      else Indices_[i] = 0;

      NumAllocatedIndicesPerRow_[i] = NumIndices;
      NumIndicesPerRow_[i] = 0;
      for (j=0; j< NumIndices; j++) Indices_[i][j] = IndexBase_ - 1; // Fill column indices with out-of-range values
    }
  }	 
  else {
    for (i=0; i<NumMyBlockRows_; i++) {
      Indices_[i] = 0;
      NumAllocatedIndicesPerRow_[i] = 0;
      NumIndicesPerRow_[i] = 0;
    }
  }

    SetAllocated(true);
    return(0);
}
//==============================================================================
int Epetra_CrsGraph::ReAllocate() {

  // Reallocate storage that was deleted in OptimizeStorage
  NumAllocatedIndicesPerRow_ = new int[NumMyBlockRows_];
  
  for (int i=0; i<NumMyBlockRows_; i++) NumAllocatedIndicesPerRow_[i] = NumIndicesPerRow_[i];

  StorageOptimized_ = false;
  return(0);
}
//==============================================================================
Epetra_CrsGraph::~Epetra_CrsGraph()
{
  int i;

  if (CV_==Copy) {
    if (All_Indices_!=0) delete [] All_Indices_;
    else for (i=0; i<NumMyBlockRows_; i++) if (NumAllocatedIndicesPerRow_[i] >0) delete  [] Indices_[i];
  }
    
  delete [] Indices_;
  if (!StorageOptimized()) delete [] NumAllocatedIndicesPerRow_;
  delete [] NumIndicesPerRow_;

  if (ColMap_!=0) delete ColMap_;
  if (Importer_!=0) delete Importer_;

  if (Exporter_!=0) delete Exporter_;

  NumMyBlockRows_ = 0;
  
  Filled_ = false;
  Allocated_ = false;
}

//==========================================================================
int Epetra_CrsGraph::InsertGlobalIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreGlobal(true);
  Row = LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR(InsertIndices(Row, NumIndices, Indices));
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::InsertMyIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot insert local indices into a global graph
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreLocal(true);
  EPETRA_CHK_ERR(InsertIndices(Row, NumIndices, Indices));
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::InsertIndices(int Row, int NumIndices, int *Indices) {

  SetSorted(false); // No longer in sorted state.
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  int j;
  int * tmp_Indices;
  int ierr = 0;

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  if (CV_==View) {
    if (Indices_[Row]!=0) ierr = 2; // This row has be defined already.  Issue warning.
    Indices_[Row] = Indices;
    NumAllocatedIndicesPerRow_[Row] = NumIndices;
    NumIndicesPerRow_[Row] = NumIndices;
  }
  else {

    if( ColMap_ ) { //only insert indices in col map if defined
      int * tmpIndices = Indices;
      Indices = new int[NumIndices];
      int loc = 0;
      if( IndicesAreLocal() ) {
        for( j = 0; j < NumIndices; ++j )
          if( ColMap_->MyLID(tmpIndices[j]) ) Indices[loc++] = tmpIndices[j];
      }
      else {
        for( j = 0; j < NumIndices; ++j )
          if( ColMap_->MyGID(tmpIndices[j]) ) Indices[loc++] = tmpIndices[j];
      }
      if( loc != NumIndices ) ierr = 2; //Some columns excluded
      NumIndices = loc;
    }

    int start = NumIndicesPerRow_[Row];
    int stop = start + NumIndices;
    int NumAllocatedIndices = NumAllocatedIndicesPerRow_[Row];
    if (stop > NumAllocatedIndices){
      if (NumAllocatedIndices==0) Indices_[Row] = new int[NumIndices]; // Row was never allocated, so do it
      else {
	ierr = 1; // Out of room.  Must delete and allocate more space...
	tmp_Indices = new int[stop];
	for (j=0; j< start; j++) tmp_Indices[j] = Indices_[Row][j]; // Copy existing entries
	delete [] Indices_[Row]; // Delete old storage
	Indices_[Row] = tmp_Indices; // Set pointer to new storage
      }
      NumAllocatedIndicesPerRow_[Row] = stop;
    }
    
    NumIndicesPerRow_[Row] = stop;
    for (j=start; j<stop; j++) Indices_[Row][j ] = Indices[j-start];

    if( ColMap_ ) delete [] Indices;
  }
  MaxNumIndices_ = EPETRA_MAX(MaxNumIndices_, NumIndicesPerRow_[Row]);
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::RemoveGlobalIndices(int Row, int NumIndices, int *Indices) {


  int j, k;
  int ierr = 0;
  int Loc;

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot remove global indices from a filled graph

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumCurrentIndices = NumIndicesPerRow_[Row];
  
  for (j=0; j<NumIndices; j++) {
    int Index = Indices[j];
    if (FindGlobalIndexLoc(Row,Index,j,Loc)) {
      for (k=Loc+1; k<NumCurrentIndices; k++) Indices_[Row][k-1] = Indices_[Row][k];
      NumCurrentIndices--;
      NumIndicesPerRow_[Row]--;
    }
  }
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::RemoveMyIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph

  int j, k;
  int ierr = 0;
  int Loc;

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumCurrentIndices = NumIndicesPerRow_[Row];
  
  for (j=0; j<NumIndices; j++) {
    int Index = Indices[j];
    if (FindMyIndexLoc(Row,Index,j,Loc)) {
      for (k=Loc+1; k<NumCurrentIndices; k++) Indices_[Row][k-1] = Indices_[Row][k];
      NumCurrentIndices--;
      NumIndicesPerRow_[Row]--;
    }
  }
  SetGlobalConstantsComputed(false); // No longer have valid global constants.

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::RemoveGlobalIndices(int Row) {

  int j;
  int ierr = 0;

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph
  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumIndices = NumIndicesPerRow_[Row];
  NumIndicesPerRow_[Row] = 0;
  
  for (j=0; j<NumIndices; j++) Indices_[Row][j] = IndexBase_ -1; // Set to invalid 

  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::RemoveMyIndices(int Row) {

  int j;
  int ierr = 0;

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot remove from a filled graph
  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  int NumIndices = NumIndicesPerRow_[Row];
  NumIndicesPerRow_[Row] = 0;
  
  for (j=0; j<NumIndices; j++) Indices_[Row][j] = -1; // Set to invalid 

  SetGlobalConstantsComputed(false); // No longer have valid global constants.
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
bool Epetra_CrsGraph::FindGlobalIndexLoc(int LocalRow, int Index, int Start, int & Loc) {
  int j;
  int NumIndices = NumIndicesPerRow_[LocalRow];
  int * Indices = Indices_[LocalRow];

  // If we have transformed the column indices, we must map this global Index to local
  if (IndicesAreLocal()) Index = LCID(Index);

  int j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
  for (j=0; j< NumIndices; j++) {
    if (j0>=NumIndices) j0 = 0; // wrap around

    if (Indices[j0]==Index) {
      Loc = j0;
      return(true);
    }
    j0++;
  }
  return(false);
}

//==========================================================================
bool Epetra_CrsGraph::FindMyIndexLoc(int LocalRow, int Index, int Start, int & Loc) {
  int j;
  int NumIndices = NumIndicesPerRow_[LocalRow];
  int * Indices = Indices_[LocalRow];

  // If we have transformed the column indices, we must map this global Index to local
  if (IndicesAreGlobal()) 
		throw ReportError("Epetra_CrsGraph::FindMyIndexLoc", -1);// Indices must be local

  int j0 = Start; // Start search at index Start (must be >= 0 and < NumIndices)
  for (j=0; j< NumIndices; j++) {
    if (j0>=NumIndices) j0 -= NumIndices; // wrap around
 
    if (Indices[j0]==Index) {
      Loc = j0;
      return(true);
    }
    j0++;
  }
  return(false);
}

//==========================================================================
int Epetra_CrsGraph::TransformToLocal() {
  EPETRA_CHK_ERR(TransformToLocal((Epetra_BlockMap *) (&RowMap()), (Epetra_BlockMap *) (&RowMap())));
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::TransformToLocal(const Epetra_BlockMap *DomainMap, const Epetra_BlockMap *RangeMap) {

  DomainMap_ = DomainMap;
  RangeMap_ = RangeMap;

  MakeIndicesLocal(*DomainMap, *RangeMap); // Convert indices to zero based on each processor
  SortIndices();  // Sort column entries from smallest to largest
  RemoveRedundantIndices(); // Get rid of any redundant index values
  MakeImportExport();  // Build Import or Export objects
  ComputeGlobalConstants(); // Compute constants that require communication
  SetFilled(true);

  return(0);
}

//==========================================================================
int Epetra_CrsGraph::ComputeGlobalConstants() {

  int i, j;

  if (GlobalConstantsComputed()) return(0);

  int * tempvec = new int[8]; // Temp space


  NumMyEntries_ = 0; // Compute Number of Nonzero entries and max
  MaxNumIndices_ = 0;
  {for (int i=0; i< NumMyBlockRows_; i++) {
    NumMyEntries_ += NumIndicesPerRow_[i];
    MaxNumIndices_ = EPETRA_MAX(MaxNumIndices_,NumIndicesPerRow_[i]);
  }}
  

  // Case 1:  Constant block size (including blocksize = 1)
  if (RowMap().ConstantElementSize()) {

    tempvec[0] = NumMyEntries_;
    tempvec[1] = NumMyBlockDiagonals_;

    Comm().SumAll(tempvec, tempvec+2, 2);
    Comm().MaxAll(&MaxNumIndices_, &GlobalMaxNumIndices_, 1);
    
    NumGlobalEntries_ = tempvec[2];
    NumGlobalBlockDiagonals_ = tempvec[3];

    int RowElementSize = RowMap().MaxElementSize();
    int ColElementSize = RowElementSize;
    NumGlobalDiagonals_ = tempvec[3] * RowElementSize;
    NumMyNonzeros_ = NumMyEntries_ * RowElementSize * ColElementSize;
    NumGlobalNonzeros_ = NumGlobalEntries_ * RowElementSize * ColElementSize;
    MaxNumNonzeros_ = MaxNumIndices_ * RowElementSize * ColElementSize;
    GlobalMaxNumNonzeros_ = GlobalMaxNumIndices_ * RowElementSize * ColElementSize;
  }

  // Case 2:  Variable block size (more work)
  else {

    NumMyNonzeros_ = 0;  // We will count the number of nonzeros here
    MaxNumNonzeros_ = 0;  // We will determine the max number of nonzeros in any one block row
    int * RowElementSizeList = RowMap().ElementSizeList();
    int * ColElementSizeList = RowElementSizeList;
    if (Importer()!=0) ColElementSizeList = ColMap().ElementSizeList();
    for (i=0; i<NumMyBlockRows_; i++){
      int NumEntries = NumIndicesPerRow_[i];
      int * Indices = Indices_[i];
      if (NumEntries>0) {
	int CurNumNonzeros = 0;
	int RowDim = RowElementSizeList[i];
	for (j=0; j<NumEntries; j++) {
	  int ColDim = ColElementSizeList[Indices[j]];
	  CurNumNonzeros += RowDim*ColDim;
	  MaxColDim_ = EPETRA_MAX(MaxColDim_, ColDim);
	}
	MaxNumNonzeros_ = EPETRA_MAX(MaxNumNonzeros_, CurNumNonzeros);
	NumMyNonzeros_ += CurNumNonzeros;
      }
    }
    
    // Sum Up all nonzeros

    
    tempvec[0] = NumMyEntries_;
    tempvec[1] = NumMyBlockDiagonals_;
    tempvec[2] = NumMyDiagonals_;
    tempvec[3] = NumMyNonzeros_;
    
    Comm().SumAll(tempvec, tempvec+4, 4);
    
    NumGlobalEntries_ = tempvec[4];
    NumGlobalBlockDiagonals_ = tempvec[5];
    NumGlobalDiagonals_ = tempvec[6];
    NumGlobalNonzeros_ = tempvec[7];

    tempvec[0] = MaxNumIndices_;
    tempvec[1] = MaxNumNonzeros_;

    Comm().MaxAll(tempvec, tempvec+2, 2);

    GlobalMaxNumIndices_ = tempvec[2];
    GlobalMaxNumNonzeros_ = tempvec[3];
  }
  
  NumGlobalRows_ = RangeMap_->NumGlobalPoints();
  NumGlobalCols_ = DomainMap_->NumGlobalPoints();

  GlobalConstantsComputed_ = true;

  delete [] tempvec;
  
  return(0);


}
//==========================================================================
int Epetra_CrsGraph::SortIndices() {

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-1);
  if (Sorted()) return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort, which is fast if indices are already sorted.

  for (int i=0; i<NumMyBlockRows_; i++){

    int n = NumIndicesPerRow_[i];
    int * const list = Indices_[i];
    int m = n/2;
    
    while (m > 0) {
      int max = n - m;
      for (int j=0; j<max; j++)
        {
	  for (int k=j; k>=0; k-=m)
            {
	      if (list[k+m] >= list[k])
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
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::RemoveRedundantIndices() {

  int i, j, k, jj;

  if (NoRedundancies()) return(0);
  if (!Sorted()) EPETRA_CHK_ERR(-1);  // Must have sorted index set
  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Indices must be local

  // For each row, remove column indices that are repeated.
  // Also, determine if graph is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortIndices was already called.
  

  NumMyDiagonals_ = 0;
  NumMyBlockDiagonals_ = 0;
  for (i=0; i<NumMyBlockRows_; i++){
    bool diagfound = false;
    int NumIndices = NumIndicesPerRow_[i];
    if (NumIndices>0) {
      int * const Indices = Indices_[i];
      int j0 = 0;
      jj = Indices[0];
      if (jj > i) LowerTriangular_ = false;
      if (jj < i) UpperTriangular_ = false;
      if (jj ==i) diagfound = true;
      for (j=1; j<NumIndices; j++) {
	jj = Indices[j];
	if (jj > i) LowerTriangular_ = false;
	if (jj < i) UpperTriangular_ = false;
	if (jj ==i) diagfound = true;
	if (jj==Indices[j0]) { // Check if index is repeated
	  NumIndicesPerRow_[i]--; // Decrement NumIndices count
	  NumMyNonzeros_ --;
	  for (k=j; k<NumIndices-1; k++) Indices[k] = Indices[k+1]; // Shift indices
	  NumIndices--;
	  j--;
	}
	else j0=j; // Redefine comparison index value
      }
      if (diagfound) {
	NumMyBlockDiagonals_++;
	NumMyDiagonals_ += RowMap().ElementSize(i);
      }
    }
  }

  NoDiagonal_ = (NumMyBlockDiagonals_==0);

  SetNoRedundancies(true);
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::MakeColMap(const Epetra_BlockMap & DomainMap, const Epetra_BlockMap & RangeMap) {

  int i, j;

  if (ColMap_!=0) return(0); // Already have a Column Map

  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if (IndicesAreLocal()) EPETRA_CHK_ERR(-1); // Return error: Indices must be global
  
  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.

  Epetra_HashTable LocalGIDs(DomainMap.NumMyElements());
  Epetra_HashTable RemoteGIDs(DomainMap.NumMyElements());
  Epetra_HashTable RemoteGIDList(DomainMap.NumMyElements());

  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;
  for (i=0; i<NumMyBlockRows_; i++){
    const int NumIndices = NumIndicesPerRow_[i];
    for (j=0; j<NumIndices; j++) {
      int GID = Indices_[i][j];
      // Check if GID matches a row GID
      if (DomainMap.MyGID(GID)) {
	if (LocalGIDs.Get(GID)==-1) // This means its a new local GID
	  LocalGIDs.Add(GID, NumLocalColGIDs++);
      }
      else {
	if (RemoteGIDs.Get(GID)==-1) { // This means its a new remote GID
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.Add(NumRemoteColGIDs++, GID);
	}
      }
    }
  }

  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int NumMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  int *ColIndices = 0;
  if (NumMyBlockCols>0) ColIndices = new int[NumMyBlockCols];

  int *RemoteColIndices = ColIndices+NumLocalColGIDs; // Points to back end of ColIndices

  for (i = 0; i<NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList.Get(i); 

  int NLists = 1;
  int *PIDList = 0;
  int *SizeList = 0;
  int *RemoteSizeList = 0;
  bool DoSizes = !DomainMap.ConstantElementSize(); // If not constant element size, then we must exchange
      
  if (NumRemoteColGIDs>0) PIDList = new int[NumRemoteColGIDs];

  if (DoSizes) {
    if (NumMyBlockCols>0) SizeList = new int[NumMyBlockCols];
    RemoteSizeList = SizeList+NumLocalColGIDs;
    NLists++;
  }
  EPETRA_CHK_ERR(DomainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList, 0, RemoteSizeList));
      
  // Sort External column indices so that all columns coming from a given remote processor are contiguous
      
  Epetra_Util Util;      
  int **SortLists = new int*[2];
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  Util.Sort(true, NumRemoteColGIDs, PIDList, 0, 0, NLists, SortLists);
  delete [] SortLists;
  if (PIDList!=0) delete [] PIDList;

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise 
  // (2) We step through the GIDs of the DomainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if (NumLocalColGIDs==DomainMap.NumMyElements()) {
    DomainMap.MyGlobalElements(ColIndices); // Load Global Indices into first NumMyBlockCols elements column GID list
    if (DoSizes) DomainMap.ElementSizeList(SizeList); // Load ElementSizeList too
  }
  else {
    int NumMyElements = DomainMap.NumMyElements();
    int * MyGlobalElements = DomainMap.MyGlobalElements();
    int * ElementSizeList = 0;
    if (DoSizes) ElementSizeList = DomainMap.ElementSizeList();
    int NumLocalAgain = 0;
    for (i=0; i<NumMyElements; i++) {
      int GID = MyGlobalElements[i];
      if (LocalGIDs.Get(GID)!=-1) {
	if (DoSizes) SizeList[NumLocalAgain] = ElementSizeList[i];
	ColIndices[NumLocalAgain++] = GID;
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }


  // Make Column map with same element sizes as Domain map

  if (DomainMap.MaxElementSize()==1) // Simple map
    ColMap_ = new Epetra_Map(-1, NumMyBlockCols, ColIndices,
				  DomainMap.IndexBase(), DomainMap.Comm());
  else if (DomainMap.ConstantElementSize()) // Constant Block size map
    ColMap_ = new Epetra_BlockMap(-1, NumMyBlockCols, ColIndices, DomainMap.MaxElementSize(),
				  DomainMap.IndexBase(), DomainMap.Comm());

  // Most general case where block size is variable.
  else
    ColMap_ = new Epetra_BlockMap(-1, NumMyBlockCols, ColIndices, SizeList,
				  DomainMap.IndexBase(), DomainMap.Comm());

  if (NumMyBlockCols>0) delete [] ColIndices; // Delete workspace
  if (SizeList!=0) delete [] SizeList;

  return(0);
}

//==========================================================================
int Epetra_CrsGraph::MakeImportExport() {

 // Create Import object for use by matrix classes.    This is only needed if ColMap and DomainMap are different
  if (!ColMap().SameAs(DomainMap()))
    Importer_ = new Epetra_Import(ColMap(), DomainMap());
  
  // Now see if we need to define an export map.  This is only needed if RowMap and RangeMap are different
  
  if (!RowMap().SameAs(RangeMap()))
    Exporter_ = new Epetra_Export(RowMap(), RangeMap()); // Create Export object. 
   
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::MakeIndicesLocal(const Epetra_BlockMap & DomainMap, const Epetra_BlockMap & RangeMap) {

  ComputeIndexState(); // Update index state by checking IndicesAreLocal/Global on all PEs
  if (IndicesAreLocal() && IndicesAreGlobal()) EPETRA_CHK_ERR(-1); // Return error: Indices must not be both local and global

  MakeColMap(DomainMap, RangeMap); // If user has not prescribed column map, create one from indices
  
  // Transform indices to local index space

  if (IndicesAreGlobal()) {

    for (int i=0; i<NumMyBlockRows_; i++) {
      const int NumIndices = NumIndicesPerRow_[i];
      for (int j=0; j<NumIndices; j++) {
	int GID = Indices_[i][j];
	int LID = LCID(GID);
	if (LID!=-1) Indices_[i][j] = LID;
	else 
	  throw ReportError("Internal error in TransformToLocal ",-1); 
	
      }
    }
  }
  else { // Indices are already local

    // Do a sanity check on column indices.  They must all be in the range 0 to NumMyBlockCols_
    // Indices will be sorted so we only need to check the last one

    if (!Sorted()) SortIndices();  // Must have sorted index set

    for (int i=0; i<NumMyBlockRows_; i++) {
      int NumIndices = NumIndicesPerRow_[i];
      if (NumIndices>0) 
	if (Indices_[i][NumIndices-1] >=NumMyBlockCols_) EPETRA_CHK_ERR(-1);
    }
  }

      
  // Store number of local columns
  NumMyCols_ = ColMap().NumMyPoints();
  NumMyBlockCols_ = ColMap().NumMyElements();

  SetIndicesAreLocal(true);
  SetIndicesAreGlobal(false);

  return(0);
}

//==========================================================================
int Epetra_CrsGraph::OptimizeStorage() {

  int i, j;
  int NumIndices;

  if (StorageOptimized()) return(0); // Have we been here before?

  bool Contiguous = true; // Assume contiguous is true
  for (i=1; i<NumMyBlockRows_; i++){
    NumIndices = NumIndicesPerRow_[i-1];
    int NumAllocateIndices = NumAllocatedIndicesPerRow_[i-1];

    // Check if NumIndices is same as NumAllocatedIndices and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if ((NumIndices!=NumAllocateIndices) || (Indices_[i]!=Indices_[i-1]+NumIndices)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumIndices and NumAllocatedIndices
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) return(1);  // This is user data, it's not contiguous and we can't make it so.

  if (NumAllocatedIndicesPerRow_ != NumIndicesPerRow_)
    delete [] NumAllocatedIndicesPerRow_; // This space is not needed after construction

  NumAllocatedIndicesPerRow_ = NumIndicesPerRow_; // Once things are contiguous, allocated equals used.

  StorageOptimized_ = true; // We can do it

  if (Contiguous) return(0); // Everything is done.  Return

 // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
  NumMyNonzeros_ = 0;
  for (i=0; i< NumMyBlockRows_; i++) NumMyNonzeros_ += NumIndicesPerRow_[i];

  // Allocate one big integer array for all index values
  All_Indices_ = new int[NumMyNonzeros_];
  
  // Set Indices_ to point into All_Indices_
  
  int * tmp = All_Indices_;
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumIndices = NumIndicesPerRow_[i];
    for (j=0; j<NumIndices; j++) tmp[j] = Indices_[i][j];
    if (Indices_[i] !=0) delete [] Indices_[i];
    Indices_[i] = tmp;
    tmp += NumIndices;
  }
  
  SetIndicesAreContiguous(true); // Can no longer dynamically add indices

    return(0);
}
//==========================================================================
int Epetra_CrsGraph::ExtractGlobalRowCopy(int Row, int LenOfIndices, int & NumIndices, int * Indices) const 
{
  int j;

  Row = LRID(Row); // Normalize row range

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = NumIndicesPerRow_[Row];
  if (LenOfIndices < NumIndices) EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices

  if (IndicesAreLocal())  for (j=0; j<NumIndices; j++) Indices[j] = GCID(Indices_[Row][j]);
  else for(j=0; j<NumIndices; j++)Indices[j] = Indices_[Row][j];
  
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::ExtractMyRowCopy(int Row, int LenOfIndices, int & NumIndices, int * Indices) const 
{
  int j;

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  NumIndices = NumIndicesPerRow_[Row];
  if (LenOfIndices < NumIndices) EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumIndices


  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-3); // There are no local indices yet

  for(j=0; j<NumIndices; j++)Indices[j] = Indices_[Row][j];
  
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::ExtractGlobalRowView(int Row, int & NumIndices, int *& Indices) const 
{
  
  Row = LRID(Row); // Normalize row range

  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // There are no global indices

  NumIndices = NumIndicesPerRow_[Row];

  Indices = Indices_[Row];
  
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::ExtractMyRowView(int Row, int & NumIndices, int *& Indices) const 
{
  
  if (Row < 0 || Row >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // There are no local indices

  NumIndices = NumIndicesPerRow_[Row];

  Indices = Indices_[Row];
  
  return(0);
}

//==========================================================================
int Epetra_CrsGraph::NumGlobalIndices(int Row) const {
  Row = LRID(Row);
  if (Row!=-1) return(NumIndicesPerRow_[Row]);
  else return(0); // No indices for this row on this processor
}
//==========================================================================
int Epetra_CrsGraph::NumAllocatedGlobalIndices(int Row) const {
  Row = LRID(Row);
  if (Row!=-1) return(NumAllocatedIndicesPerRow_[Row]);
  else return(0); // No indices allocated for this row on this processor
}
//=========================================================================
int Epetra_CrsGraph::CheckSizes(const Epetra_SrcDistObject & Source) {
  try {
    const Epetra_CrsGraph & A = dynamic_cast<const Epetra_CrsGraph &>(Source);
    if (!A.GlobalConstantsComputed()) EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  }
  catch (...) {
    return(0); // No error at this point, object could be a RowMatrix
  }
  return(0);
}
//=========================================================================
int Epetra_CrsGraph::CopyAndPermute(const Epetra_SrcDistObject & Source,
					 int NumSameIDs, 
					 int NumPermuteIDs, int * PermuteToLIDs,
					 int *PermuteFromLIDs) {
 
  try {
    const Epetra_CrsGraph & A = dynamic_cast<const Epetra_CrsGraph &>(Source);
    EPETRA_CHK_ERR(CopyAndPermuteCrsGraph(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					  PermuteFromLIDs));
  }
  catch (...) {
    try {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
      EPETRA_CHK_ERR(CopyAndPermuteRowMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					  PermuteFromLIDs));
    }
    catch (...) {
      EPETRA_CHK_ERR(-1); // Incompatible SrcDistObject
    }
  }
  
  return(0);
}

//=========================================================================
int Epetra_CrsGraph::CopyAndPermuteCrsGraph(const Epetra_CrsGraph & A,
					    int NumSameIDs, 
					    int NumPermuteIDs, int * PermuteToLIDs,
					    int *PermuteFromLIDs) {
 
  int i;
  
  int Row, NumIndices;
  int * Indices = 0;
  int FromRow, ToRow;

  int MaxNumIndices = A.MaxNumIndices();

  if( MaxNumIndices>0 && A.IndicesAreLocal()) Indices = new int[MaxNumIndices];
  
  // Do copy first
  if (NumSameIDs>0) {
    if (A.IndicesAreLocal()) {
      for (i=0; i<NumSameIDs; i++) {
	Row = GRID(i);
	EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Indices));
	// Place into target graph.  
	int ierr = InsertGlobalIndices(Row, NumIndices, Indices); 
	if (ierr<0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for (i=0; i<NumSameIDs; i++) {
	Row = GRID(i);
	EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumIndices, Indices)); // Set pointer	
	// Place into target graph.
	  int ierr = InsertGlobalIndices(Row, NumIndices, Indices); 
	  if (ierr<0) EPETRA_CHK_ERR(ierr); 
      }
    }	
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (A.IndicesAreLocal()) {
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = A.GRID(PermuteFromLIDs[i]);
	ToRow = GRID(PermuteToLIDs[i]);
	EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices));
	// Place into target graph.
	int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices); 
	if (ierr<0) EPETRA_CHK_ERR(ierr); 
      }
    }
    else { // A.IndiceAreGlobal()
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = A.GRID(PermuteFromLIDs[i]);
	ToRow = GRID(PermuteToLIDs[i]);
	EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumIndices, Indices)); // Set pointer
	// Place into target graph.
	int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices); 
 	if (ierr<0) EPETRA_CHK_ERR(ierr); 
     }
    }
  }	

  if (MaxNumIndices>0 && A.IndicesAreLocal()) delete [] Indices;
    
  return(0);
}

//=========================================================================
int Epetra_CrsGraph::CopyAndPermuteRowMatrix(const Epetra_RowMatrix & A,
					     int NumSameIDs, 
					     int NumPermuteIDs, int * PermuteToLIDs,
					     int *PermuteFromLIDs) {
 
  int i, j;
  int NumIndices;
  int * Indices = 0;
  double * Values = 0;
  int FromRow, ToRow;
  
  int MaxNumIndices = A.MaxNumEntries();

  if (MaxNumIndices>0) {
    Indices = new int[MaxNumIndices];
    Values = new double[MaxNumIndices]; // Must extract values even though we discard them
  }

  const Epetra_Map & RowMap = A.RowMatrixRowMap();
  const Epetra_Map & ColMap = A.RowMatrixColMap();
  
  // Do copy first
  for (i=0; i<NumSameIDs; i++) {
    ToRow = RowMap.GID(i);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, MaxNumIndices, NumIndices, Values, Indices));
    for (j=0; j<NumIndices; j++) Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    // Place into target graph.  
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices);
    if (ierr<0) EPETRA_CHK_ERR(ierr);
  }
  
  // Do local permutation next
  for (i=0; i<NumPermuteIDs; i++) {
    FromRow = PermuteFromLIDs[i];
    ToRow = GRID(PermuteToLIDs[i]);
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, MaxNumIndices, NumIndices, Values, Indices));
    for (j=0; j<NumIndices; j++) Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices); // Place into target graph.
    if (ierr<0) EPETRA_CHK_ERR(ierr);
  }
  
  if(MaxNumIndices>0) {
    delete [] Indices;
    delete [] Values;
  }
  
  return(0);
}

//=========================================================================
int Epetra_CrsGraph::PackAndPrepare(const Epetra_SrcDistObject & Source, 
				     int NumExportIDs, int * ExportLIDs,
				     int Nsend, int Nrecv,
				     int & LenExports, char * & Exports, int & LenImports, 
				     char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor) {
  
  int GlobalMaxNumIndices = 0;

  try {
    const Epetra_CrsGraph & A = dynamic_cast<const Epetra_CrsGraph &>(Source);
    GlobalMaxNumIndices = A.GlobalMaxNumIndices();
  }
  catch (...) {
    try {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
      int MaxNumIndices = A.MaxNumEntries();
      A.Comm().MaxAll(&MaxNumIndices, &GlobalMaxNumIndices, 1);
    }
    catch (...) {
    EPETRA_CHK_ERR(-1); // Bad cast
    }
  }

  int * IntExports = 0;
  int * IntImports = 0;
  int IntPacketSize = GlobalMaxNumIndices + 2;
  SizeOfPacket = IntPacketSize * sizeof(int); 

  if (IntPacketSize*Nsend>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = IntPacketSize*Nsend;
    IntExports = new int[LenExports];
    Exports = (char *) IntExports;
  }

  if (IntPacketSize*Nrecv>LenImports) {
    if (LenImports>0) delete [] Imports;
    LenImports = IntPacketSize*Nrecv;
    IntImports = new int[LenImports];
    Imports = (char *) IntImports;
  }
  if (NumExportIDs<=0) return(0);

  try {
    const Epetra_CrsGraph & A = dynamic_cast<const Epetra_CrsGraph &>(Source);
    EPETRA_CHK_ERR(PackAndPrepareCrsGraph(A, NumExportIDs, ExportLIDs, Nsend, Nrecv, LenExports, Exports,
			   LenImports, Imports, SizeOfPacket, Distor, IntPacketSize));
  }
  catch (...) {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
    EPETRA_CHK_ERR(PackAndPrepareRowMatrix(A, NumExportIDs, ExportLIDs, Nsend, Nrecv, LenExports, Exports,
			   LenImports, Imports, SizeOfPacket, Distor, IntPacketSize));
  }
return(0);
}
//=========================================================================
int Epetra_CrsGraph::PackAndPrepareCrsGraph(const Epetra_CrsGraph & A, 
				     int NumExportIDs, int * ExportLIDs,
				     int Nsend, int Nrecv,
				     int & LenExports, char * & Exports, int & LenImports, 
				     char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor, int IntPacketSize) {
  
  int i;
  
  int NumIndices;
  int * Indices;
  int FromRow;
  int * intptr;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int MaxNumIndices = A.MaxNumIndices();
  intptr = (int *) Exports;
  for (i=0; i<NumExportIDs; i++) {
    FromRow = A.GRID(ExportLIDs[i]);
    *intptr = FromRow;
    Indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices));
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += IntPacketSize; // Point to next segment
  }
    
return(0);
}
//=========================================================================
int Epetra_CrsGraph::PackAndPrepareRowMatrix(const Epetra_RowMatrix & A, 
				     int NumExportIDs, int * ExportLIDs,
				     int Nsend, int Nrecv,
				     int & LenExports, char * & Exports, int & LenImports, 
				     char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor, int IntPacketSize) {
  

  int i, j;
  
  int NumIndices;
  int * Indices;
  int FromRow;
  int * intptr;
  double * Values = 0;
  
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
  //   sized segments for current communication routines.
  int MaxNumIndices = A.MaxNumEntries();
  if (MaxNumIndices>0) Values = new double[MaxNumIndices];
  const Epetra_Map & RowMap = A.RowMatrixRowMap();
  const Epetra_Map & ColMap = A.RowMatrixColMap();

  intptr = (int *) Exports;
  for (i=0; i<NumExportIDs; i++) {
    FromRow = RowMap.GID(ExportLIDs[i]);
    *intptr = FromRow;
    Indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], MaxNumIndices, NumIndices, Values, Indices));
    for (j=0; j<NumIndices; j++) Indices[j] = ColMap.GID(Indices[j]); // convert to GIDs
    intptr[1] = NumIndices; // Load second slot of segment
    intptr += IntPacketSize; // Point to next segment
  }
    
  if (Values!=0) delete [] Values;
return(0);
}
//=========================================================================
int Epetra_CrsGraph::UnpackAndCombine(const Epetra_SrcDistObject & Source, 
				       int NumImportIDs, int * ImportLIDs, 
				       char * Imports, int & SizeOfPacket, 
				      Epetra_Distributor & Distor, Epetra_CombineMode CombineMode) {



  if (NumImportIDs<=0) return(0);

  int GlobalMaxNumIndices = 0;

  try {
    const Epetra_CrsGraph & A = dynamic_cast<const Epetra_CrsGraph &>(Source);
    GlobalMaxNumIndices = A.GlobalMaxNumIndices();
  }
  catch (...) {
    try {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
      int MaxNumIndices = A.MaxNumEntries();
      A.Comm().MaxAll(&MaxNumIndices, &GlobalMaxNumIndices, 1);
    }
    catch (...) {
    EPETRA_CHK_ERR(-1); // Bad cast
    }
  }

  int IntPacketSize = GlobalMaxNumIndices + 2;



  int NumIndices;
  int * Indices;
  int ToRow;
  int i;
  
  int * intptr;
  // Unpack it...

  // Each segment of Sends will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source graph
  // next int:  NumIndices, Number of indices in row.
  // next NumIndices: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be 
  //  wasted but we need fixed sized segments for current communication routines.

  intptr = (int *) Imports;
    
  for (i=0; i<NumImportIDs; i++) {
    ToRow = GRID(ImportLIDs[i]);
    assert((intptr[0])==ToRow); // Sanity check
    NumIndices = intptr[1];
    Indices = intptr + 2; 
    // Insert indices
    int ierr = InsertGlobalIndices(ToRow, NumIndices, Indices);
    if (ierr<0) EPETRA_CHK_ERR(ierr);
    intptr += IntPacketSize; // Point to next segment
  }
  
  return(0);
}
//=========================================================================

bool Epetra_CrsGraph::GlobalConstantsComputed() const {

  int mineComputed = 0;
  int allComputed;
  if (GlobalConstantsComputed_) mineComputed = 1;
  RowMap().Comm().MinAll(&mineComputed, &allComputed, 1); // Find out if any PEs changed constants info
  // If allComputed comes back zero then at least one PE need global constants recomputed.
  return(allComputed==1);
}
//=========================================================================

void Epetra_CrsGraph::ComputeIndexState() {

  int myIndicesAreLocal = 0;
  int myIndicesAreGlobal = 0;
  if (IndicesAreLocal_) myIndicesAreLocal = 1;
  if (IndicesAreGlobal_) myIndicesAreGlobal = 1;
  int allIndicesAreLocal;
  int allIndicesAreGlobal;
  RowMap().Comm().MaxAll(&myIndicesAreLocal, &allIndicesAreLocal, 1); // Find out if any PEs changed Local Index info
  RowMap().Comm().MaxAll(&myIndicesAreGlobal, &allIndicesAreGlobal, 1); // Find out if any PEs changed Global Index info
  IndicesAreLocal_ = (allIndicesAreLocal==1); // If indices are local on one PE, should be local on all
  IndicesAreGlobal_ = (allIndicesAreGlobal==1);  // If indices are global on one PE should be local on all
  return;
}
//=========================================================================

void Epetra_CrsGraph::Print (ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (MyPID==0) {
	os <<  "\nNumber of Global Block Rows  = "; os << NumGlobalBlockRows(); os << endl;
	os <<    "Number of Global Block Cols  = "; os << NumGlobalBlockCols(); os << endl;
	os <<    "Number of Global Block Diags = "; os << NumGlobalBlockDiagonals(); os << endl;
	os <<    "Number of Global Entries     = "; os << NumGlobalEntries(); os << endl;
	os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros(); os << endl;
	os <<  "\nGlobal Maximum Block Row Dim = "; os << GlobalMaxRowDim(); os << endl;
	os <<    "Global Maximum Block Col Dim = "; os << GlobalMaxColDim(); os << endl;
	os <<    "Global Maximum Num Indices   = "; os << GlobalMaxNumIndices(); os << endl;
	if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Block Rows  = "; os << NumMyBlockRows(); os << endl;
      os <<    "Number of My Block Cols  = "; os << NumMyBlockCols(); os << endl;
      os <<    "Number of My Block Diags = "; os << NumMyBlockDiagonals(); os << endl;
      os <<    "Number of My Entries     = "; os << NumMyEntries(); os << endl;
      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl;
      os <<  "\nMy Maximum Block Row Dim = "; os << MaxRowDim(); os << endl;
      os <<    "My Maximum Block Col Dim = "; os << MaxColDim(); os << endl;
      os <<    "My Maximum Num Indices   = "; os << MaxNumIndices(); os << endl; os << endl;

      int NumMyBlockRows1 = NumMyBlockRows();
      int MaxNumIndices1 = MaxNumIndices();
      int * Indices1  = new int[MaxNumIndices1];
      int NumIndices1;
      int i, j;
      
      os.width(14);
      os <<  "       Row Index "; os << " ";
      for (j = 0; j < MaxNumIndices(); j++) {   
	os.width(12);
	os <<  "Col Index"; os << "      ";
      }
      os << endl;
      for (i=0; i<NumMyBlockRows1; i++) {
	int Row = GRID(i); // Get global row number
	ExtractGlobalRowCopy(Row, MaxNumIndices1, NumIndices1, Indices1);
	
	os.width(14);
	os <<  Row ; os << "    ";	
	for (j = 0; j < NumIndices1 ; j++)
	  {   
	    os.width(12);
	    os <<  Indices1[j]; os << "    ";
	  }
	os << endl;
      }

      delete [] Indices1;
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
  }

  return;
}
