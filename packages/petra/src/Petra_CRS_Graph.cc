#include "Petra_CRS_Graph.h"

//==============================================================================
Petra_CRS_Graph::Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int *NumIndicesPerRow) 
  : RowMap_(RowMap),
    ColMap_(RowMap),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate(NumIndicesPerRow, 1);
}

//==============================================================================
Petra_CRS_Graph::Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int NumIndicesPerRow) 
  : RowMap_(RowMap),
    ColMap_(RowMap),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate(&NumIndicesPerRow, 0);
}


//==============================================================================
Petra_CRS_Graph::Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int *NumIndicesPerRow) 
  : RowMap_(RowMap),
    ColMap_(ColMap),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate(NumIndicesPerRow, 1);
}

//==============================================================================
Petra_CRS_Graph::Petra_CRS_Graph(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int NumIndicesPerRow) 
  : RowMap_(RowMap),
    ColMap_(ColMap),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate(&NumIndicesPerRow, 0);
}

//==============================================================================
Petra_CRS_Graph::Petra_CRS_Graph(const Petra_CRS_Graph & Graph) 
  : RowMap_(Graph.RowMap_),
    ColMap_(Graph.ColMap_),
    DomainMap_(Graph.DomainMap_),
    RangeMap_(Graph.RangeMap_),
    ImportMap_(Graph.ImportMap_),
    Importer_(Graph.Importer_),
    ExportMap_(Graph.ExportMap_),
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
    All_Indices_(0),
    LenImports_(0),
    LenExports_(0),
    Imports_(0),
    Exports_(0),
    CV_(Copy)
{
  int ierr = Allocate(Graph.NumIndicesPerRow(), 1);
  for (int i=0; i<NumMyBlockRows_; i++) {
    NumIndicesPerRow_[i] = NumAllocatedIndicesPerRow_[i];
    for (int j=0; j< NumIndicesPerRow_[i]; j++) Indices_[i][j] = Graph.Indices_[i][j];
  }
  MaxNumIndices_ = Graph.MaxNumIndices();
  GlobalMaxNumIndices_ = Graph.GlobalMaxNumIndices();
  if (ImportMap_ != 0) ImportMap_ = new Petra_BlockMap(Graph.ImportMap()); // Non-trivial import column map, must copy it.
  if (Importer_ != 0) Importer_ = new Petra_Import(*Graph.Importer()); // Non-trivial importer, must copy it.

  if (ExportMap_ != 0) ExportMap_ = new Petra_BlockMap(Graph.ExportMap()); // Non-trivial export row map, must copy it.
  if (Exporter_ != 0) Exporter_ = new Petra_Export(*Graph.Exporter()); // Non-trivial exporter, must copy it.
}

//==============================================================================
void Petra_CRS_Graph::InitializeDefaults() { // Initialize all attributes that have trivial default values

  IndexBase_ = RowMap().IndexBase();

  DomainMap_ = 0;
  RangeMap_ = 0;
  ImportMap_ = 0;
  Importer_ = 0;
  ExportMap_ = 0;
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
  NumGlobalBlockCols_ = ColMap().NumGlobalElements();
  NumMyBlockRows_ = RowMap().NumMyElements();
  NumMyBlockCols_ = ColMap().NumMyElements();

  NumGlobalRows_ = RowMap().NumGlobalEquations();
  NumGlobalCols_ = ColMap().NumGlobalEquations();
  NumMyRows_ = RowMap().NumMyEquations();
  NumMyCols_ = ColMap().NumMyEquations();

  GlobalMaxRowDim_ = RowMap().MaxElementSize();
  MaxRowDim_ = RowMap().MaxElementSize();
  GlobalMaxColDim_ = ColMap().MaxElementSize();
  MaxColDim_ = ColMap().MaxElementSize();

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
  LenImports_ = 0;
  LenExports_ = 0;
  Imports_ = 0;
  Exports_ = 0;
}

//==============================================================================
int Petra_CRS_Graph::Allocate(int * NumIndicesPerRow, int Inc ) {

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
#ifdef PETRA_LEVELSCHEDULING
  NumThreads_= 0;
  ThreadStartRows_= 0;
  NumLevels_= 0;
  LevelIndices_= 0;
  LevelOrder_= 0;
#endif

    SetAllocated(true);
    return(0);
}
//==============================================================================
int Petra_CRS_Graph::ReAllocate() {

  // Reallocate storage that was deleted in OptimizeStorage
  NumAllocatedIndicesPerRow_ = new int[NumMyBlockRows_];
  
  for (int i=0; i<NumMyBlockRows_; i++) NumAllocatedIndicesPerRow_[i] = NumIndicesPerRow_[i];

  StorageOptimized_ = false;
  return(0);
}
//==============================================================================
Petra_CRS_Graph::~Petra_CRS_Graph()
{
  int i;

  if (CV_==Copy) {
    if (All_Indices_!=0) delete [] All_Indices_;
    else for (i=0; i<NumMyBlockRows_; i++) if (NumAllocatedIndicesPerRow_[i] >0) delete  [] Indices_[i];
  }
    
  delete [] Indices_;
  if (!StorageOptimized()) delete [] NumAllocatedIndicesPerRow_;
  delete [] NumIndicesPerRow_;

  if (ImportMap_!=0) delete ImportMap_;
  if (Importer_!=0) delete Importer_;

  if (ExportMap_!=0 && ExportMap_ != RangeMap_) delete ExportMap_;
  if (Exporter_!=0) delete Exporter_;

  NumMyBlockRows_ = 0;
  
  Filled_ = false;
  Allocated_ = false;
}

//==========================================================================
int Petra_CRS_Graph::InsertGlobalIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreLocal()) return(-2); // Cannot insert global values into local graph
  if (IndicesAreContiguous()) return(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreGlobal(true);
  Row = LRID(Row); // Find local row number for this global row index

  return(InsertIndices(Row, NumIndices, Indices));
}

//==========================================================================
int Petra_CRS_Graph::InsertMyIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreGlobal()) return(-2); // Cannot insert local indices into a global graph
  if (IndicesAreContiguous()) return(-3); // Indices cannot be individually deleted and newed
  SetIndicesAreLocal(true);
  return(InsertIndices(Row, NumIndices, Indices));
}

//==========================================================================
int Petra_CRS_Graph::InsertIndices(int Row, int NumIndices, int *Indices) {

  SetSorted(false); // No longer in sorted state.

  int j;
  int * tmp_Indices;
  int ierr = 0;

  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range
    
  if (CV_==View) {
    if (Indices_[Row]!=0) ierr = 2; // This row has be defined already.  Issue warning.
    Indices_[Row] = Indices;
    NumAllocatedIndicesPerRow_[Row] = NumIndices;
    NumIndicesPerRow_[Row] = NumIndices;
  }
  else {
    
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
  }
  MaxNumIndices_ = maxfn(MaxNumIndices_, NumIndicesPerRow_[Row]);
  return(ierr);
}

//==========================================================================
int Petra_CRS_Graph::RemoveGlobalIndices(int Row, int NumIndices, int *Indices) {


  int j, k;
  int ierr = 0;
  int Loc;

  if (IndicesAreLocal()) return(-2); // Cannot remove global indices from a filled graph

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range
    
  int NumCurrentIndices = NumIndicesPerRow_[Row];
  
  for (j=0; j<NumIndices; j++) {
    int Index = Indices[j];
    if (FindGlobalIndexLoc(Row,Index,j,Loc)) {
      for (k=Loc+1; k<NumCurrentIndices; k++) Indices_[Row][k-1] = Indices_[Row][k];
      NumCurrentIndices--;
      NumIndicesPerRow_[Row]--;
    }
  }

  return(ierr);
}

//==========================================================================
int Petra_CRS_Graph::RemoveMyIndices(int Row, int NumIndices, int *Indices) {

  if (IndicesAreGlobal()) return(-2); // Cannot insert global values into filled graph

  int j, k;
  int ierr = 0;
  int Loc;

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range
    
  int NumCurrentIndices = NumIndicesPerRow_[Row];
  
  for (j=0; j<NumIndices; j++) {
    int Index = Indices[j];
    if (FindMyIndexLoc(Row,Index,j,Loc)) {
      for (k=Loc+1; k<NumCurrentIndices; k++) Indices_[Row][k-1] = Indices_[Row][k];
      NumCurrentIndices--;
      NumIndicesPerRow_[Row]--;
    }
  }

  return(ierr);
}

//==========================================================================
int Petra_CRS_Graph::RemoveGlobalIndices(int Row) {

  int j;
  int ierr = 0;

  if (IndicesAreLocal()) return(-2); // Cannot remove from a filled graph
  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

  Row = LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range
    
  int NumIndices = NumIndicesPerRow_[Row];
  NumIndicesPerRow_[Row] = 0;
  
  for (j=0; j<NumIndices; j++) Indices_[Row][j] = IndexBase_ -1; // Set to invalid 

  return(ierr);
}

//==========================================================================
int Petra_CRS_Graph::RemoveMyIndices(int Row) {

  int j;
  int ierr = 0;

  if (IndicesAreGlobal()) return(-2); // Cannot remove from a filled graph
  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.

  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range
    
  int NumIndices = NumIndicesPerRow_[Row];
  NumIndicesPerRow_[Row] = 0;
  
  for (j=0; j<NumIndices; j++) Indices_[Row][j] = -1; // Set to invalid 

  return(ierr);
}

//==========================================================================
bool Petra_CRS_Graph::FindGlobalIndexLoc(int LocalRow, int Index, int Start, int & Loc) {
  int j;
  int NumIndices = NumIndicesPerRow_[LocalRow];
  int * Indices = Indices_[LocalRow];

  // If we have transformed the column indices, we must map this global Index to local
  if (IndicesAreLocal()) Index = LCID(Index);

  for (j=0; j< NumIndices; j++) {
    int j0 = (j+Start)%NumIndices; // Start search at index Start
    if (Indices[j0]==Index) {
      Loc = j0;
      return(true);
    }
  }
  return(false);
}

//==========================================================================
bool Petra_CRS_Graph::FindMyIndexLoc(int LocalRow, int Index, int Start, int & Loc) {
  int j;
  int NumIndices = NumIndicesPerRow_[LocalRow];

  // If we have transformed the column indices, we must map this global Index to local
  if (IndicesAreGlobal()) return(-2); // Indices must be local
  for (j=0; j< NumIndices; j++) {
    int j0 = (j+Start)%NumIndices; // Start search at index Start
    if (Indices_[LocalRow][j0]==Index) {
      Loc = j0;
      return(true);
    }
  }
  return(false);
}

//==========================================================================
int Petra_CRS_Graph::TransformToLocal() {
  return(TransformToLocal((Petra_BlockMap *) (&RowMap_), (Petra_BlockMap *) (&ColMap_)));
}

//==========================================================================
int Petra_CRS_Graph::TransformToLocal(Petra_BlockMap *DomainMap, Petra_BlockMap *RangeMap) {

  DomainMap_ = DomainMap;
  RangeMap_ = RangeMap;
  
  MakeIndicesLocal(*DomainMap_, *RangeMap_); // Convert indices to zero based on each processor

  SortIndices();  // Sort column entries from smallest to largest

  RemoveRedundantIndices(); // Get rid of any redundant index values

  ComputeGlobalConstants(); // Compute constants that require communication

  SetFilled(true);

  return(0);
}

//==========================================================================
int Petra_CRS_Graph::ComputeGlobalConstants() {

  int i, j, k;

  if (GlobalConstantsComputed_) return(0);

  int * tempvec = new int[6]; // Temp space


  NumMyEntries_ = 0; // Compute Number of Nonzero entries and max
  MaxNumIndices_ = 0;
  for (int i=0; i< NumMyBlockRows_; i++) {
    NumMyEntries_ += NumIndicesPerRow_[i];
    MaxNumIndices_ = maxfn(MaxNumIndices_,NumIndicesPerRow_[i]);
  }
  

  // Case 1:  Constant block size (including blocksize = 1)
  if (RowMap().ConstantElementSize() && ColMap().ConstantElementSize()) {

    tempvec[0] = NumMyEntries_;
    tempvec[1] = NumMyBlockDiagonals_;

    Comm().SumAll(tempvec, tempvec+2, 2);
    Comm().MaxAll(&MaxNumIndices_, &GlobalMaxNumIndices_, 1);
    
    NumGlobalEntries_ = tempvec[2];
    NumGlobalBlockDiagonals_ = tempvec[3];

    int RowElementSize = RowMap().MaxElementSize();
    int ColElementSize = ColMap().MaxElementSize();
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
    if (Importer()!=0) ColElementSizeList = ImportMap().ElementSizeList();
    for (i=0; i<NumMyBlockRows_; i++){
      int NumEntries = NumIndicesPerRow_[i];
      int * Indices = Indices_[i];
      if (NumEntries>0) {
	int CurNumNonzeros = 0;
	int RowDim = RowElementSizeList[i];
	for (j=0; j<NumEntries; j++) {
	  int ColDim = ColElementSizeList[Indices[j]];
	  CurNumNonzeros += RowDim*ColDim;
	  MaxColDim_ = maxfn(MaxColDim_, ColDim);
	}
	MaxNumNonzeros_ = maxfn(MaxNumNonzeros_, CurNumNonzeros);
	NumMyNonzeros_ += CurNumNonzeros;
      }
    }
    
    // Sum Up all nonzeros

    
    tempvec[0] = NumMyEntries_;
    tempvec[1] = NumMyBlockDiagonals_;
    tempvec[2] = NumMyNonzeros_;
    
    Comm().SumAll(tempvec, tempvec+3, 3);
    
    NumGlobalEntries_ = tempvec[3];
    NumGlobalBlockDiagonals_ = tempvec[4];
    NumGlobalNonzeros_ = tempvec[5];

    tempvec[0] = MaxNumIndices_;
    tempvec[1] = MaxNumNonzeros_;

    Comm().MaxAll(tempvec, tempvec+2, 2);

    GlobalMaxNumIndices_ = tempvec[2];
    GlobalMaxNumNonzeros_ = tempvec[3];
  }
  
  GlobalConstantsComputed_ = true;

  delete [] tempvec;
  
  return(0);


}
//==========================================================================
int Petra_CRS_Graph::SortIndices() {

  if (!IndicesAreLocal()) return(-1);
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
int Petra_CRS_Graph::RemoveRedundantIndices() {

  int i, j, k, jj;

  if (NoRedundancies()) return(0);
  if (!Sorted()) return(-1);  // Must have sorted index set
  if (IndicesAreGlobal()) return(-2); // Indices must be local

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
int Petra_CRS_Graph::MakeIndicesLocal(const Petra_BlockMap & DomainMap, const Petra_BlockMap & RangeMap) {

  int i, j, k;

  if (!IndicesAreLocal()) {

    DomainMap_ = (Petra_BlockMap *) &DomainMap;
    RangeMap_ = (Petra_BlockMap *) &RangeMap;

    // For each row, check if column indices are owned or not.
    // If owned, transform global index to local index.  
    // If not owned, add to ImportMap for later use

    NumMyBlockCols_ = DomainMap_->NumMyElements(); // Redefine NumMyBlockCols_ to number of local domain map elements
    NumMyCols_ = DomainMap_->NumMyEquations(); // Redefine NumMyCols_ to number of local domain map equations
    
    int IncBlockCols = maxfn(minfn(NumMyBlockCols_/4,100),10);
    int MaxBlockCols = 0;
    int *ColIndices = 0;
    if (DomainMap_->DistributedGlobal()) {
      MaxBlockCols = NumMyBlockCols_;
      ColIndices = new int[MaxBlockCols];
    }
    int *new_ColIndices = 0;

    int NewNumMyBlockCols = NumMyBlockCols_;
    for (i=0; i<NumMyBlockRows_; i++){
      const int NumIndices = NumIndicesPerRow_[i];
      for (j=0; j<NumIndices; j++) {
	int GID = Indices_[i][j];
	int LID = LCID(GID);
	if (LID==-1) {
	  bool ExtNew = true;
	  for (k=NumMyBlockCols_; k<NewNumMyBlockCols; k++) {
	    if (ColIndices[k]==GID) {
	      ExtNew = false;
	      break;
	    }
	  }
	  if (ExtNew) {
	    if (NewNumMyBlockCols >= MaxBlockCols) { // Need to expand...
	      MaxBlockCols = maxfn(NumMyBlockCols_+IncBlockCols, MaxBlockCols+IncBlockCols); // Increment column space 
	      new_ColIndices = new int[MaxBlockCols];
	      for (k=NumMyBlockCols_; k<NewNumMyBlockCols; k++) new_ColIndices[k] = ColIndices[k];
	      if (ColIndices!=0) delete [] ColIndices;
	      ColIndices = new_ColIndices;
	    }
	    ColIndices[NewNumMyBlockCols] = GID;
	    /*
         if (Comm().MyPID()==0) {
            cout << " i, j, NumIndices = " << i << "  " << j << "  " << NumIndices << endl;
            cout << "New Column = " << NewNumMyBlockCols << " with GID = " << GID << endl;
         }
	    */
	    NewNumMyBlockCols++;
	  }
	}
      }
    }


    // Create ImportMap.  This map will be used to facilitate communication in matrix classes
    
    if (DomainMap_->DistributedGlobal()) {
      
      // Find processors that own the off-processor GIDs
      int NumRemote = NewNumMyBlockCols - NumMyBlockCols_;
      int *RemoteColIndices = ColIndices+NumMyBlockCols_;
      int NLists = 1;
      int *PIDList = 0;
      int *SizeList = 0;
      int *RemoteSizeList = 0;
      bool DoSizes = !DomainMap_->ConstantElementSize(); // If not constant element size, then we must exchange
      
      if (NumRemote>0) PIDList = new int[NumRemote];

      if (DoSizes) {
	if (NewNumMyBlockCols>0) SizeList = new int[NewNumMyBlockCols];
	RemoteSizeList = SizeList+NumMyBlockCols_;
	NLists++;
      }
      DomainMap_->RemoteIDList(NumRemote, RemoteColIndices, PIDList, 0, RemoteSizeList);
      
      // Sort External column indices so that all columns coming from a given remote processor are contiguous
      
      Petra_Util Util;      
      int **SortLists = new int*[2];
      SortLists[0] = RemoteColIndices;
      SortLists[1] = RemoteSizeList;
      Util.Sort(true, NumRemote, PIDList, 0, 0, NLists, SortLists);
      delete [] SortLists;

      if (NumRemote>0) delete []PIDList;
      
      DomainMap_->MyGlobalElements(ColIndices); // Load Global Indices into first NumMyBlockCols_ elements of import column map
      if (DoSizes)DomainMap_->ElementSizeList(SizeList); // Load ElementSizeList into first NumMyBlockCols_ elements of import size list

      NumMyBlockCols_ = NewNumMyBlockCols; // Redefine NumMyBlockCols_ based on local columns plus number of columns needing imported elements

      // Make Import map with same element sizes as Domain map
      if (DomainMap_->ConstantElementSize()) // Constant Block size map
	ImportMap_ = new Petra_BlockMap(-1, NewNumMyBlockCols, ColIndices, DomainMap_->MaxElementSize(),
					DomainMap_->IndexBase(), DomainMap_->Comm());

      // Most general case where block size is variable.
      else
	ImportMap_ = new Petra_BlockMap(-1, NewNumMyBlockCols, ColIndices, SizeList,
					DomainMap_->IndexBase(), DomainMap_->Comm());

      Importer_ = new Petra_Import(ImportMap(), *DomainMap_); // Create Import object for use by matrix classes.   
    
      delete [] ColIndices; // Delete workspace
      if (DoSizes && NewNumMyBlockCols>0) delete [] SizeList;
      
      // Recompute number of local columns
      NumMyCols_ = ImportMap_->NumMyEquations();
    }

    // Now see if we need to define an export map.  This is only needed if RowMap and RangeMap are different

    if (!RowMap_.SameAs(*RangeMap_)) {
      Exporter_ = new Petra_Export(RowMap(), *RangeMap_); // Create Export object. 
      ExportMap_ = RangeMap_;
    }
   
    SetIndicesAreLocal(true);
    SetIndicesAreGlobal(false);

    // Transform indices to local index space

    for (i=0; i<NumMyBlockRows_; i++) {
      const int NumIndices = NumIndicesPerRow_[i];
      for (j=0; j<NumIndices; j++) {
	int GID = Indices_[i][j];
	int LID = LCID(GID);
	if (LID!=-1) Indices_[i][j] = LID;
	else {
	  cout << "Internal error in TransformToLocal " << endl; 
	  abort();
	}
	
      }
    }
  }
  else { // Indices are already local

    // Do a sanity check on column indices.  They must all be in the range 0 to NumMyBlockCols_
    // Indices will be sorted so we only need to check the last one

    if (!Sorted()) SortIndices();  // Must have sorted index set

    for (i=0; i<NumMyBlockRows_; i++) {
      int NumIndices = NumIndicesPerRow_[i];
      if (NumIndices>0) 
	if (Indices_[i][NumIndices-1] >=NumMyBlockCols_) return(-1);
    }
  }


  return(0);
}

//==========================================================================
int Petra_CRS_Graph::OptimizeStorage() {

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


  if ((CV_==View) && !Contiguous) return(-1);  // This is user data, it's not contiguous and we can't make it so.

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
int Petra_CRS_Graph::ExtractGlobalRowCopy(int Row, int LenOfIndices, int & NumIndices, int * Indices) const 
{
  int j;

  Row = LRID(Row); // Normalize row range


  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range

  NumIndices = NumIndicesPerRow_[Row];
  if (LenOfIndices < NumIndices) return(-2); // Not enough space for copy. Needed size is passed back in NumIndices


  if (IndicesAreLocal())  for (j=0; j<NumIndices; j++) Indices[j] = GCID(Indices_[Row][j]);
  else for(j=0; j<NumIndices; j++)Indices[j] = Indices_[Row][j];
  
  return(0);
}

//==========================================================================
int Petra_CRS_Graph::ExtractMyRowCopy(int Row, int LenOfIndices, int & NumIndices, int * Indices) const 
{
  int j;

  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range

  NumIndices = NumIndicesPerRow_[Row];
  if (LenOfIndices < NumIndices) return(-2); // Not enough space for copy. Needed size is passed back in NumIndices


  if (IndicesAreGlobal()) return(-3); // There are no local indices yet

  for(j=0; j<NumIndices; j++)Indices[j] = Indices_[Row][j];
  
  return(0);
}

//==========================================================================
int Petra_CRS_Graph::ExtractGlobalRowView(int Row, int & NumIndices, int *& Indices) const 
{
  
  Row = LRID(Row); // Normalize row range

  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range

  if (IndicesAreLocal()) return(-2); // There are no global indices

  NumIndices = NumIndicesPerRow_[Row];

  Indices = Indices_[Row];
  
  return(0);
}

//==========================================================================
int Petra_CRS_Graph::ExtractMyRowView(int Row, int & NumIndices, int *& Indices) const 
{
  
  if (Row < 0 || Row >= NumMyBlockRows_) return(-1); // Not in Row range

  if (IndicesAreGlobal()) return(-2); // There are no local indices

  NumIndices = NumIndicesPerRow_[Row];

  Indices = Indices_[Row];
  
  return(0);
}

//==========================================================================
int Petra_CRS_Graph::NumGlobalIndices(int Row) const {
  Row = LRID(Row);
  if (Row!=-1) return(NumIndicesPerRow_[Row]);
  else return(0); // No indices for this row on this processor
}
//==========================================================================
int Petra_CRS_Graph::NumAllocatedGlobalIndices(int Row) const {
  Row = LRID(Row);
  if (Row!=-1) return(NumAllocatedIndicesPerRow_[Row]);
  else return(0); // No indices allocated for this row on this processor
}
//==========================================================================
int Petra_CRS_Graph::LCID( int GCID) const {

  int Index = DomainMap_->LID(GCID); 
  if (Index!=-1) return(Index);
  if (ImportMap_==0) return(-1);
  return(ImportMap_->LID(GCID)); // Check col map 

} 
//==========================================================================
int Petra_CRS_Graph::GCID( int LCID) const {

  int Index = DomainMap_->GID(LCID); // Check row map first
  if (Index!=IndexBase_-1) return(Index);
  if (ImportMap_==0) return(IndexBase_-1);
  return(ImportMap_->GID(LCID)); // Check col map 

} 
//=========================================================================
int Petra_CRS_Graph::Import(const Petra_CRS_Graph& SourceGraph, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode) {

  if (!RowMap().SameAs(Importer.TargetMap())) return(-2);
  if (!SourceGraph.RowMap().SameAs(Importer.SourceMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *RemoteLIDs = Importer.RemoteLIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int SizeOfPacket = SourceGraph.GlobalMaxNumIndices() + 2;
  int Nsend = SizeOfPacket * Importer.NumSend();
  int Nrecv = SizeOfPacket * Importer.NumRecv();

  return(DoTransfer(SourceGraph, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenExports_, Exports_, LenImports_, Imports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_CRS_Graph::Export(const Petra_CRS_Graph& SourceGraph, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode) {

  if (!RowMap().SameAs(Exporter.TargetMap())) return(-2);
  if (!SourceGraph.RowMap().SameAs(Exporter.SourceMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumRemoteIDs();
  int NumExportIDs = Exporter.NumExportIDs();
  int *ExportLIDs = Exporter.ExportLIDs();
  int *RemoteLIDs = Exporter.RemoteLIDs();
  int *PermuteToLIDs = Exporter.PermuteToLIDs();
  int *PermuteFromLIDs = Exporter.PermuteFromLIDs();
  int SizeOfPacket = SourceGraph.GlobalMaxNumIndices() + 2;
  int Nsend = SizeOfPacket * Exporter.NumSend();
  int Nrecv = SizeOfPacket * Exporter.NumRecv();

  return(DoTransfer(SourceGraph, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenExports_, Exports_, LenImports_, Imports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_CRS_Graph::Import(const Petra_CRS_Graph& SourceGraph, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Exporter.SourceMap())) return(-2);
  if (!SourceGraph.RowMap().SameAs(Exporter.TargetMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumExportIDs();
  int NumExportIDs = Exporter.NumRemoteIDs();
  int *ExportLIDs = Exporter.RemoteLIDs();
  int *RemoteLIDs = Exporter.ExportLIDs();
  int *PermuteToLIDs = Exporter.PermuteFromLIDs();
  int *PermuteFromLIDs = Exporter.PermuteToLIDs();
  int SizeOfPacket = SourceGraph.GlobalMaxNumIndices() + 2;
  int Nsend = SizeOfPacket * Exporter.NumRecv();
  int Nrecv = SizeOfPacket * Exporter.NumSend();

  return(DoTransfer(SourceGraph, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenImports_, Imports_, LenExports_, Exports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_CRS_Graph::Export(const Petra_CRS_Graph& SourceGraph, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode) {


  if (!RowMap().SameAs(Importer.SourceMap())) return(-2);
  if (!SourceGraph.RowMap().SameAs(Importer.TargetMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumExportIDs();
  int NumExportIDs = Importer.NumRemoteIDs();
  int *ExportLIDs = Importer.RemoteLIDs();
  int *RemoteLIDs = Importer.ExportLIDs();
  int *PermuteToLIDs = Importer.PermuteFromLIDs();
  int *PermuteFromLIDs = Importer.PermuteToLIDs();
  int SizeOfPacket = SourceGraph.GlobalMaxNumIndices() + 2;
  int Nsend = SizeOfPacket * Importer.NumRecv();
  int Nrecv = SizeOfPacket * Importer.NumSend();

  return(DoTransfer(SourceGraph, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, Nsend, Nrecv, SizeOfPacket,
		    LenImports_, Imports_, LenExports_, Exports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_CRS_Graph::DoTransfer(const Petra_CRS_Graph& SourceGraph, 
				     Petra_CombineMode CombineMode,
				     int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, 
				     int NumExportIDs, 
				     int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
				     int * ExportLIDs,
				     int Nsend, int Nrecv, int SizeOfPacket,
				     int & LenExports, int * & Exports,
				     int & LenImports, int * & Imports,
#ifdef PETRA_MPI
				     GSComm_Plan & Plan, 
#endif
				     bool DoReverse){
  int ierr = 0;
  ierr = 
    CopyAndPermute(*this, SourceGraph, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs);

  if (ierr!=0) return(ierr);
  if (CombineMode==Zero) return(0); // All done if CombineMode only involves copying and permuting


  // Need to send global row ID, NumIndices for the row, then the indices.
  if (Nsend>LenExports) {
    if (LenExports>0) {
      delete [] Exports;
    }
    Exports = new int[Nsend];
    LenExports = Nsend;
  }

  ierr =  Pack(SourceGraph, NumExportIDs, ExportLIDs, Exports);

  if (ierr!=0) return(ierr);
  
  if (Nrecv>LenImports) {
    if (LenImports>0) {
      delete Imports;
    }
    Imports = new int[Nrecv];
    LenImports = Nrecv;
  }
  
#ifdef PETRA_MPI
  if (RowMap().DistributedGlobal()) {
  int msgtag = 32765;
  bool GSComm_OK;
  GSComm_Comm GSComm;

  if (DoReverse)
  // Do the exchange of remote index data
    GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				  reinterpret_cast<char *> (Exports), 
				  SizeOfPacket * sizeof( int ),
				  reinterpret_cast<char *> (Imports) );
  else
    GSComm_OK = GSComm.Do( Plan, msgtag, 
			   reinterpret_cast<char *> (Exports), 
			   SizeOfPacket * sizeof( int ),
			   reinterpret_cast<char *> (Imports) );

  if (!GSComm_OK) return(-1);

  ierr = 
    UnpackAndCombine( *this, SizeOfPacket, NumRemoteIDs, RemoteLIDs, Imports, CombineMode);
  }
#endif

  return(ierr);
  
}
 
//=========================================================================
int Petra_CRS_Graph::CopyAndPermute(Petra_CRS_Graph & Target, 
					 const Petra_CRS_Graph & Source,
					 int NumSameIDs, 
					 int NumPermuteIDs, int * PermuteToLIDs,
					 int *PermuteFromLIDs) {
  
  int i;
  
  int Row, NumIndices;
  int * Indices;
  int FromRow, ToRow;
  
  // Do copy first
  if (NumSameIDs>0) {
    if (Source.IndicesAreLocal()) {
      int MaxNumIndices = Source.MaxNumIndices();
      Indices = new int[MaxNumIndices];  // Need some temporary space
      
      for (i=0; i<NumSameIDs; i++) {
	Row = Target.GRID(i);
	assert(Source.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Indices)==0);
	// Place into target graph.  
	assert(Target.InsertGlobalIndices(Row, NumIndices, Indices)==0); 
      }
      delete [] Indices;
    }
    else { // Source.IndiceAreGlobal()
      for (i=0; i<NumSameIDs; i++) {
	Row = Target.GRID(i);
	assert(Source.ExtractGlobalRowView(Row, NumIndices, Indices)==0); // Set pointer	
	// Place into target graph.
	  assert(Target.InsertGlobalIndices(Row, NumIndices, Indices)==0); 
      }
    }	
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (Source.IndicesAreLocal()) {
      int MaxNumIndices = Source.MaxNumIndices();
      Indices = new int[MaxNumIndices];  // Need some temporary space
      
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = Source.GRID(PermuteFromLIDs[i]);
	ToRow = Target.GRID(PermuteToLIDs[i]);
	assert(Source.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices)==0);
	// Place into target graph.
	assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0); 
      }
      delete [] Indices;
    }
    else { // Source.IndiceAreGlobal()
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = Source.GRID(PermuteFromLIDs[i]);
	ToRow = Target.GRID(PermuteToLIDs[i]);
	assert(Source.ExtractGlobalRowView(FromRow, NumIndices, Indices)==0); // Set pointer
	// Place into target graph.
	assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0); 
      }
    }
  }	
    
  return(0);
}

//=========================================================================
int Petra_CRS_Graph::Pack(const Petra_CRS_Graph & Source,
			       int NumSendIDs, int * SendLIDs, int * Sends) {
  
  int i;
  
  int NumIndices;
  int * Indices;
  int FromRow;
  int * intptr;
  
  if (NumSendIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source graph
    // next int:  NumIndices, Number of indices in row.
    // next NumIndices: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
    //   sized segments for current communication routines.

    intptr = Sends;
    int GlobalMaxNumIndices = Source.GlobalMaxNumIndices();
    int Inc = GlobalMaxNumIndices + 2;
    for (i=0; i<NumSendIDs; i++) {
      FromRow = Source.GRID(SendLIDs[i]);
      *intptr = FromRow;
      Indices = intptr + 2;
      assert(Source.ExtractGlobalRowCopy(FromRow, GlobalMaxNumIndices, NumIndices, Indices)==0);
      intptr[1] = NumIndices; // Load second slot of segment
      intptr += Inc; // Point to next segment
    }
    // assert(ptr-Sends==Nsend); // Sanity check on send count

  }
    
  return(0);
}

//=========================================================================
int Petra_CRS_Graph::UnpackAndCombine(Petra_CRS_Graph & Target, int SizeOfPacket,
					   int NumRecvIDs, int * RecvLIDs, 
					   int * Recvs, Petra_CombineMode CombineMode) {
  int NumIndices;
  int * Indices;
  int ToRow;
  int i;
  
  int * intptr;
  // Unpack it...

  if (NumRecvIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source graph
    // next int:  NumIndices, Number of indices in row.
    // next NumIndices: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be 
    //  wasted but we need fixed sized segments for current communication routines.

    intptr = Recvs;
    
    int Inc = SizeOfPacket;
    for (i=0; i<NumRecvIDs; i++) {
      ToRow = Target.GRID(RecvLIDs[i]);
      assert((intptr[0])==ToRow); // Sanity check
      NumIndices = intptr[1];
      Indices = intptr + 2; 
      // Insert indices
      assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0);
      intptr += Inc; // Point to next segment
    }
    // assert(ptr-Recvs==Nrecv); // Sanity check on receive count
  }
  
  return(0);
}
/*
//=========================================================================
int Petra_CRS_Graph::Import(const Petra_CRS_Graph& SourceGraph, const Petra_Import & Importer) {

  int ierr = 0;
  int Nsend, Nrecv;

  if (!RowMap_.SameAs(Importer.TargetMap())) return(-2);
  if (!SourceGraph.RowMap().SameAs(Importer.SourceMap())) return(-3);
  
  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int * ExportLIDs = Importer.ExportLIDs();
  int *RemoteLIDs = Importer.RemoteLIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();

  // Need to send global row ID, NumIndices for the row, then the indices.
  int SizeOfPacket = SourceGraph.GlobalMaxNumIndices() + 2;
  Nsend = SizeOfPacket * Importer.NumSend();
  if (Nsend>LenExports_) {
    if (LenExports_>0) delete [] Exports_;
    Exports_ = new int[Nsend];
    LenExports_ = Nsend;
  }

  ierr = 
    CopyPermuteAndPack(*this, SourceGraph, 
		       NumSameIDs, 
		       NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, 
		       NumExportIDs, ExportLIDs, Exports_);

  if (ierr!=0) return(ierr);
  
  Nrecv = SizeOfPacket * Importer.NumRecv();
  if (Nrecv>LenImports_) {
    if (LenImports_>0) delete [] Imports_;
    Imports_ = new int[Nrecv];
    LenImports_ = Nrecv;
  }
  
#ifdef PETRA_MPI
  int msgtag = 32765;

  GSComm_Comm * GSComm = new GSComm_Comm();


  // Do the exchange of remote data
  bool GSComm_OK = GSComm->Do( Importer.GSPlan(), msgtag, 
			reinterpret_cast<char *> (Exports_), 
			SizeOfPacket * sizeof( int ),
			reinterpret_cast<char *> (Imports_) );

  delete GSComm;

  if (!GSComm_OK) return(-1);
#endif

  ierr = 
    UnpackAndAdd( *this, SizeOfPacket, NumRemoteIDs, RemoteLIDs, Imports_);

  return(ierr);

}
 
//=========================================================================
int Petra_CRS_Graph::CopyPermuteAndPack(Petra_CRS_Graph & Target, const Petra_CRS_Graph & Source,
					int NumSameIDs, 
					int NumPermuteIDs, int * PermuteToLIDs, int *PermuteFromLIDs,
					int NumSendIDs, int * SendLIDs, int * Sends) {
  
  int i;
  
  int Row, NumIndices;
  int * Indices;
  int FromRow, ToRow;
  int * ptr;
  
  // Do copy first
  if (NumSameIDs>0) {
    if (Target.IndicesAreLocal()) {
      if (!Source.IndicesAreLocal()) return(-1); // Cannot import a local graph from a global graph
      for (i=0; i<NumSameIDs; i++) {
	Row = i;
	assert(Source.ExtractMyRowView(Row, NumIndices, Indices)==0); // Set pointers
	assert(Target.InsertMyIndices(Row, NumIndices, Indices)==0); // Place into target graph (copy or view depending on 
	                                                      // Petra_DataAccess argument used in target construction).
      }
    }
    else { // Target.IndicesAreGlobal()
      if (Source.IndicesAreLocal()) {
	int MaxNumIndices = Source.MaxNumIndices();
	Indices = new int[MaxNumIndices];  // Need some temporary space
					    
	for (i=0; i<NumSameIDs; i++) {
	  Row = Target.GRID(i);
	  assert(Source.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Indices)==0); // Set pointers
	  assert(Target.InsertGlobalIndices(Row, NumIndices, Indices)==0); // Place into target graph (copy or view depending on 
	                                                            // Petra_DataAccess argument used in target construction).
	}
	delete [] Indices;
      }
      else { // Source.IndiceAreGlobal()
	for (i=0; i<NumSameIDs; i++) {
	  Row = Target.GRID(i);
	  assert(Source.ExtractGlobalRowView(Row, NumIndices, Indices)==0); // Set pointers
	  assert(Target.InsertGlobalIndices(Row, NumIndices, Indices)==0); // Place into target graph (copy or view depending on 
	                                                            // Petra_DataAccess argument used in target construction).
	}
      }
    }
  }	


  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (Target.IndicesAreLocal()) {
      if (!Source.IndicesAreLocal()) return(-1); // Cannot import a local graph from a global graph
      for (i=0; i<NumPermuteIDs; i++) {
	FromRow = PermuteFromLIDs[i];
	ToRow = PermuteToLIDs[i];
	assert(Source.ExtractMyRowView(FromRow, NumIndices, Indices)==0); // Set pointers
	assert(Target.InsertMyIndices(ToRow, NumIndices, Indices)==0); // Place into target graph (copy or view depending on 
	                                                      // Petra_DataAccess argument used in target construction).
      }
    }
    else { // Target.IndicesAreGlobal()
      if (Source.IndicesAreLocal()) {
	int MaxNumIndices = Source.MaxNumIndices();
	Indices = new int[MaxNumIndices];  // Need some temporary space
					    
	for (i=0; i<NumSameIDs; i++) {
	  FromRow = Source.GRID(PermuteFromLIDs[i]);
	  ToRow = Target.GRID(PermuteToLIDs[i]);
	  assert(Source.ExtractGlobalRowCopy(FromRow, MaxNumIndices, NumIndices, Indices)==0); // Set pointers
	  assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0); // Place into target graph copy or view depending on 
	                                                            // Petra_DataAccess argument used in target construction).
	}
	delete [] Indices;
      }
      else { // Source.IndiceAreGlobal()
	for (i=0; i<NumSameIDs; i++) {
	  FromRow = Source.GRID(PermuteFromLIDs[i]);
	  ToRow = Target.GRID(PermuteToLIDs[i]);
	  assert(Source.ExtractGlobalRowView(FromRow, NumIndices, Indices)==0); // Set pointers
	  assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0); // Place into target graph (copy or view depending on 
	                                                            // Petra_DataAccess argument used in target construction).
	}
      }
    }
  }	
  
  if (NumSendIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source graph
    // next int:  NumIndices, Number of indices in row.
    // next NumIndices: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
    //   sized segments for current communication routines.

    ptr = Sends;
    int GlobalMaxNumIndices = Source.GlobalMaxNumIndices();
    int Inc = GlobalMaxNumIndices + 2;
    if (Target.IndicesAreLocal()) {
      if (!Source.IndicesAreLocal()) return(-1); // Cannot import a local graph from a global graph
      for (i=0; i<NumSendIDs; i++) {
	FromRow = SendLIDs[i];
	*ptr = Source.GRID(FromRow);
	Indices = ptr + 2; 
	assert(Source.ExtractMyRowCopy(FromRow, GlobalMaxNumIndices, NumIndices, Indices)==0); // Copy Indices
	ptr[1] = NumIndices; // Load second slot of segment
	ptr += Inc; // Point to next segment
      }
    }
    else { // IndicesAreGlobal
      for (i=0; i<NumSendIDs; i++) {
	FromRow = Source.GRID(SendLIDs[i]);
	*ptr = FromRow;
	Indices = ptr + 2;
	assert(Source.ExtractGlobalRowCopy(FromRow, GlobalMaxNumIndices, NumIndices, Indices)==0); // Set pointers
	ptr[1] = NumIndices; // Load second slot of segment
	ptr += Inc; // Point to next segment
      }
    }
    // assert(ptr-Sends==Nsend); // Sanity check on send count

  }
    
  return(0);
}

//=========================================================================
int Petra_CRS_Graph::UnpackAndAdd(Petra_CRS_Graph & Target, int SizeOfPacket,
				  int NumRecvIDs, int * RecvLIDs, 
				  int * Recvs) {
  int Row, NumIndices;
  int * Indices;
  int ToRow;
  int i;
  
  int * ptr;
  // Unpack it...

  if (NumRecvIDs>0) {

    // Each segment of Sends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source graph
    // next int:  NumIndices, Number of indices in row.
    // next NumIndices: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumIndices - NumIndices ints) will be wasted but we need fixed
    //   sized segments for current communication routines.

    ptr = Recvs;

    int Inc = SizeOfPacket;
    if (Target.IndicesAreLocal()) {
      for (i=0; i<NumRecvIDs; i++) {
	ToRow = RecvLIDs[i];
	assert(ptr[0]==Target.GRID(ToRow)); // Sanity check
	NumIndices = ptr[1];
	Indices = ptr + 2; 
	assert(Target.InsertMyIndices(ToRow, NumIndices, Indices)==0); // Put in values
	ptr += Inc; // Point to next segment
      }
    }
    else { // IndicesAreGlobal
      for (i=0; i<NumRecvIDs; i++) {
	ToRow = Target.GRID(RecvLIDs[i]);
	assert(ptr[0]==ToRow); // Sanity check
	NumIndices = ptr[1];
	Indices = ptr + 2; 
	assert(Target.InsertGlobalIndices(ToRow, NumIndices, Indices)==0); // Set pointers
	ptr += Inc; // Point to next segment
      }
    }
    // assert(ptr-Recvs==Nrecv); // Sanity check on receive count
  }
  
  return(0);
}
*/
#ifdef PETRA_LEVELSCHEDULING
//==============================================================================
int Petra_CRS_Graph::ComputeLevels(int NumThreads) {

/* -------------------------------------------------------------

 Description:
    Compute a level schedule of a sparse triangular matrix.
    See the description of parameter iopt and the remarks section.

 Remarks:
 1. Forward level scheduling depth is defined as follows:

               |  0,                                    if L(i,j) = 0 for all j < i
    depth(i) = |
               |   max { depth(i), L(i,j) is nonzero }, otherwise
                  j < i

 2. Backward level scheduling depth is defined as follows:

               |  0,                                   if L(j,i) = 0 for all j > i
    depth(i) = |
               |  max { depth(i), L(j,i) is nonzero }, otherwise
                 j > i

------------------------------------------------------------- */


  // Quick test for upper/lower triangular property
  if (!LowerTriangular() && !UpperTriangular()) return(-1);

  // Quick test for valid number of threads
  if (NumThreads<1) return(-2);


  // Local variables

  int i, j, depth;


  // Initialize a few things

  int n = NumMyBlockRows();
  if (LevelIndices_==0) LevelIndices_ = new int[n+1];
  if (LevelOrder_==0) LevelOrder_ = new int[n];
  int * head = new int[n];
  int * link = new int[n];
  for (i = 0; i <n; i++) {
    LevelIndices_[i] = 0;
    head[i] = 0;
  }
  LevelIndices_[n] = 0;

  if (LowerTriangular()) {

    // ---------------------------------------------------------------
    // Forward level scheduling based on the row structure of L
    // ---------------------------------------------------------------

    head[0] = 0;
    LevelOrder_[0] = 0;

    // ---------------------------------------------------------
    // Temporary storage usage:
    // - LevelOrder_(i) for the level numbers (depth) for node i
    // - LevelIndices_(k+1) to keep the number of nodes at each level k
    // ---------------------------------------------------------
    int * ColIndices = 0;
    for (i = 0; i < n; i++) {
      depth = -1;
      int NumIndices = NumIndicesPerRow()[i];
      if (NumIndices>0) {
	ColIndices = Indices()[i];
	for (j=0 ; j <NumIndices; j++)
	  depth= maxfn(depth,LevelOrder_[ColIndices[j]]);
      }
      depth++;
      LevelOrder_[i] = depth;
      LevelIndices_[depth+1]++;
      link[i] = head[depth];
      head[depth] = i;
    }

  } else if(UpperTriangular()) {

    // ---------------------------------------------------------------
    // Forward level scheduling based on the row structure of L(t) or
    // U starting at the bottom.  (pointer, index) must contain the
    // structure of the diagonal entries.
    // ---------------------------------------------------------------

    head[0] = n-1;
    LevelOrder_[n-1] = 0;

    // ---------------------------------------------------------
    // Temporary storage usage:
    // - LevelOrder_(i) for the level numbers (depth) for node i
    // - LevelIndices_(k+1) to keep the number of nodes at each level k
    // ---------------------------------------------------------

    int * ColIndices = 0;
    for (i = n-1; i >= 0; i--) {
      depth = -1;
      int NumIndices = NumIndicesPerRow()[i];
      if (NumIndices>0) {
	ColIndices = Indices()[i];
	for (j=0; j< NumIndices; j++)
	  depth= maxfn(depth,LevelOrder_[ColIndices[j]]);
      }
      depth++;
      LevelOrder_[i] = depth;
      LevelIndices_[depth + 1]++;
      link[i] = head[depth];
      head[depth] = i;
		  
    }
  }
  
  // ------------------------------------
  // From length of each level, calculate
  // the pointers to each level in LevelOrder_
  // ------------------------------------
  for (i = 1; i <n+1; i++) {
    LevelIndices_[i] += LevelIndices_[i-1];
    if (LevelIndices_[i] == n) {
      NumLevels_ = i;
      break;
    }
  }
  
  // -------------------------------------------------------------
  // Determine LevelOrder_, the ordering of the nodes by level numbers.
  // This can be done in parallel since we already have LevelIndices_.
  // -------------------------------------------------------------
  for (i=0; i<NumLevels_; i++) {
    int i2 = LevelIndices_[i+1] - LevelIndices_[i];
    listfill(i2, head[i], link, &LevelOrder_[LevelIndices_[i]]);
  }
  
  
  delete [] head;
  delete [] link;

  // Now determine starting row for each thread.  Try to split work so each
  // thread gets about the same number of nonzeros at each level.

  // At the end of these loops, ThreadStartRows_[i][j] = the first row that thread j
  // will compute with at level i.

  NumThreads_ = NumThreads;
  ThreadStartRows_ = new int*[NumLevels_];

  bool debug = true;
  if (debug) cout << " Number of Levels = "<< NumLevels_ << endl;
  

    if (debug) cout << " Number of nonzeros per level" << endl;
  for (i=0; i<NumLevels_; i++) {
    ThreadStartRows_[i] = new int[NumThreads_+1];
    int levelstart = LevelIndices_[i];
    int nextlevelstart = LevelIndices_[i+1];
    int nnzlevel = 0;
    for (j=levelstart; j<nextlevelstart; j++) 
      nnzlevel += NumIndicesPerRow()[LevelOrder_[j]]+1;
    if (debug) cout << i <<"  "<< nnzlevel << endl;

    int avg_nnzlevel = maxfn(nnzlevel/NumThreads_, 1);
    ThreadStartRows_[i][0] = levelstart;
    
    int curthread = 0;
    int curnnz = 0;
    for (j=levelstart; j<nextlevelstart; j++) {
      //  if (debug) cout << " Row = "<< LevelOrder_[j] << endl;
      int nnzrow = NumIndicesPerRow()[LevelOrder_[j]];
      curnnz += nnzrow;
      if (curnnz>=avg_nnzlevel) {
	if (curthread==NumThreads_) break;
	curthread++;
	ThreadStartRows_[i][curthread] = j+1;
	curnnz = 0;
      }
    }
    // Set end pointer (and any remaining thread start rows for pathological cases)
    curthread++;
    for (j=curthread; j<=NumThreads; j++) ThreadStartRows_[i][j] = nextlevelstart;
  }
    
	
	
  delete [] LevelIndices_; // For now we do not need this info since the ThreadStartRows_ array serves the same purpose
  
  return(0);
} 
//------------------------------------------------------------
void Petra_CRS_Graph::listfill(int n, int head, int * link, int * vector) {
  
  
  for (int i=0; i<n; i++) {
    vector[n-i-1] = head;
    head = link[head];
  }
}
#endif

//==========================================================================
// Non-member functions

ostream& operator << (ostream& os, const Petra_CRS_Graph& A)
{
  int MyPID = A.RowMap().Comm().MyPID();
  int NumProc = A.RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      if (MyPID==0) {
	os <<  "\nNumber of Global Block Rows  = "; os << A.NumGlobalBlockRows(); os << endl;
	os <<    "Number of Global Block Cols  = "; os << A.NumGlobalBlockCols(); os << endl;
	os <<    "Number of Global Block Diags = "; os << A.NumGlobalBlockDiagonals(); os << endl;
	os <<    "Number of Global Entries     = "; os << A.NumGlobalEntries(); os << endl;
	os <<  "\nNumber of Global Rows        = "; os << A.NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << A.NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << A.NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << A.NumGlobalNonzeros(); os << endl;
	os <<  "\nGlobal Maximum Block Row Dim = "; os << A.GlobalMaxRowDim(); os << endl;
	os <<    "Global Maximum Block Col Dim = "; os << A.GlobalMaxColDim(); os << endl;
	os <<    "Global Maximum Num Indices   = "; os << A.GlobalMaxNumIndices(); os << endl;
	if (A.LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (A.UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (A.NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Block Rows  = "; os << A.NumMyBlockRows(); os << endl;
      os <<    "Number of My Block Cols  = "; os << A.NumMyBlockCols(); os << endl;
      os <<    "Number of My Block Diags = "; os << A.NumMyBlockDiagonals(); os << endl;
      os <<    "Number of My Entries     = "; os << A.NumMyEntries(); os << endl;
      os <<  "\nNumber of My Rows        = "; os << A.NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << A.NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << A.NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << A.NumMyNonzeros(); os << endl;
      os <<  "\nMy Maximum Block Row Dim = "; os << A.MaxRowDim(); os << endl;
      os <<    "My Maximum Block Col Dim = "; os << A.MaxColDim(); os << endl;
      os <<    "My Maximum Num Indices   = "; os << A.MaxNumIndices(); os << endl; os << endl;

      int NumMyBlockRows = A.NumMyBlockRows();
      int MaxNumIndices = A.MaxNumIndices();
      int * Indices  = new int[MaxNumIndices];
      int NumIndices;
      int i, j;
      
      os.width(14);
      os <<  "       Row Index "; os << " ";
      for (j = 0; j < A.MaxNumIndices(); j++) {   
	os.width(12);
	os <<  "Col Index"; os << "      ";
      }
      os << endl;
      for (i=0; i<NumMyBlockRows; i++) {
	int Row = A.GRID(i); // Get global row number
	A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Indices);
	
	os.width(14);
	os <<  Row ; os << "    ";	
	for (j = 0; j < NumIndices ; j++)
	  {   
	    os.width(12);
	    os <<  Indices[j]; os << "    ";
	  }
	os << endl;
      }

      delete [] Indices;
      
      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
  }

  return os;
}
