#include "Petra_RDP_VBR_Matrix.h"

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int *NumBlockEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, NumBlockEntriesPerRow);
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int NumBlockEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, NumBlockEntriesPerRow);
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, 
					   const Petra_BlockMap& ColMap, int *NumBlockEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  cerr << "Not Implemented." << endl; abort();
  InitializeDefaults();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, ColMap, NumBlockEntriesPerRow);
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, 
					   const Petra_BlockMap& ColMap, int NumBlockEntriesPerRow) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  cerr << "Not Implemented." << endl; abort();
  InitializeDefaults();
  Graph_ = new Petra_CRS_Graph(CV, RowMap, ColMap, NumBlockEntriesPerRow);
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_CRS_Graph & Graph) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_((Petra_CRS_Graph*) &Graph),
    Allocated_(false),
    StaticGraph_(true),
    NumMyBlockRows_(Graph.RowMap().NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  int ierr = Allocate();
}

//==============================================================================
Petra_RDP_VBR_Matrix::Petra_RDP_VBR_Matrix(const Petra_RDP_VBR_Matrix & Source) 
  : Petra_Flops(),
    Petra_BLAS(),
    Graph_(Source.Graph_),
    Allocated_(Source.Allocated_),
    StaticGraph_(Source.StaticGraph_),
    NumMyBlockRows_(Source.NumMyBlockRows_),
    CV_(Copy) {
  InitializeDefaults();
  if (!Source.StaticGraph()) Graph_ = new Petra_CRS_Graph(Source.Graph());
  int ierr = Allocate();

  int i, j;
  
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumBlockEntries = NumBlockEntriesPerRow_[i];
    int RowDim = ElementSizeList_[i];
    for (j=0; j < NumBlockEntries; j++) {
      int LDA = Source.LDAs_[i][j];
      int ColDim = Source.ColDims_[i][j];
      ColDims_[i][j] = ColDim;
      LDAs_[i][j] = LDA;
      Values_[i][j] = new double[LDA*ColDim];
      CopyMat(Source.Values_[i][j], LDA, RowDim, ColDim, Values_[i][j], LDA, false);
    }
  }

}

//==============================================================================
void Petra_RDP_VBR_Matrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  Values_ = 0;
  ColDims_ = 0;
  LDAs_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  ImportVector_ = 0;
  LenIntImports_ = 0;
  LenIntExports_ = 0;
  LenImports_ = 0;
  LenExports_ = 0;
  Imports_ = 0;
  Exports_ = 0;
  IntImports_ = 0;
  IntExports_ = 0;

  NumBlockEntriesPerRow_  = 0;
  NumAllocatedBlockEntriesPerRow_ = 0;
  Indices_ = 0;
  ElementSizeList_ = 0;
  FirstElementEntryList_ = 0;
  
  // State variables needed for constructing matrix entry-by-entry

  TempRowDims_ = 0;
  TempColDims_ = 0;
  TempLDAs_ = 0;
  TempValues_ = 0;
  LenTemps_ = 0;
  CurBlockRow_ = 0;
  CurNumBlockEntries_ = 0;
  CurBlockIndices_ = 0;
  CurEntry_ = -1; // Set to -1 to allow a simple sanity check when submitting entries
  CurIndicesAreLocal_ = false;
  CurSubmitMode_ = Insert;
  
  // State variables needed for extracting entries
  CurExtractBlockRow_ = 0;
  CurExtractEntry_ = -1; // Set to -1 to allow a simple sanity check when extracting entries
  CurExtractNumBlockEntries_ = 0;
  CurExtractIndicesAreLocal_ = false;
  CurExtractView_ = false;
  CurRowDim_ = 0;

  // State variable for extracting block diagonal entries
  CurBlockDiag_ = -1; // Set to -1 to allow a simple sanity check when extracting entries

  return;
}

//==============================================================================
int Petra_RDP_VBR_Matrix::Allocate() {

  int i, j;
  
  // Set direct access pointers to graph info (needed for speed)
  NumBlockEntriesPerRow_ = Graph_->NumIndicesPerRow();
  NumAllocatedBlockEntriesPerRow_ = Graph_->NumAllocatedIndicesPerRow();
  Indices_ = Graph_->Indices();

  ElementSizeList_ = RowMap().ElementSizeList();

  FirstElementEntryList_ = RowMap().FirstElementEntryList();
  

  // Allocate Values array
  Values_ = new double**[NumMyBlockRows_];
  ColDims_ = new int *[NumMyBlockRows_];
  LDAs_ = new int *[NumMyBlockRows_];
  // Allocate and initialize entries
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumAllocatedBlockEntries = NumAllocatedBlockEntriesPerRow_[i];
    
    if (NumAllocatedBlockEntries > 0) {
      Values_[i] = new double*[NumAllocatedBlockEntries];
      ColDims_[i] = new int[NumAllocatedBlockEntries];
      LDAs_[i] = new int[NumAllocatedBlockEntries];
      for (j=0; j < NumAllocatedBlockEntries; j++) {
	Values_[i][j] = 0;
	ColDims_[i][j] = 0;
	LDAs_[i] [j] = 0;
      }
    }
    else {
      Values_[i] = 0;
      ColDims_[i] = 0;
      LDAs_[i] = 0;
    }
  }
  SetAllocated(true);
    return(0);
}
//==============================================================================
Petra_RDP_VBR_Matrix::~Petra_RDP_VBR_Matrix(){

  int i;

  for (i=0; i<NumMyBlockRows_; i++) {
    int NumAllocatedBlockEntries = NumAllocatedBlockEntriesPerRow_[i];
    
    if (NumAllocatedBlockEntries >0) {
      if (All_Values_!=0) delete [] All_Values_[i];
      else if (CV_==Copy)
	for (int j=0; j < NumAllocatedBlockEntries; j++) 
	  if (Values_[i][j]!=0) {
	    delete [] Values_[i][j];
	    Values_[i][j] = 0;
	  }
      
      delete Values_[i];
      delete [] ColDims_[i];
      delete [] LDAs_[i];
    }
  }

  if (All_Values_!=0)   delete [] All_Values_;

  if (Values_!=0)       delete [] Values_;
  if (ColDims_!=0)      delete [] ColDims_;
  if (LDAs_!=0)         delete [] LDAs_;


  if (ImportVector_!=0) delete ImportVector_;
  if (Imports_!=0)      delete Imports_;
  if (Exports_!=0)      delete Exports_;  


  NumMyBlockRows_ = 0;

  if (LenTemps_>0) {
    delete [] TempRowDims_;
    delete [] TempColDims_;
    delete [] TempLDAs_;
    delete [] TempValues_;
  }
  
  InitializeDefaults(); // Reset all basic pointers to zero
  Allocated_ = false;

  if (!StaticGraph())   delete Graph_; // We created the graph, so must delete it.
  
}

//==============================================================================
int Petra_RDP_VBR_Matrix::PutScalar(double Scalar) 
{
  for (int i=0; i<NumMyBlockRows_; i++) {
    int NumBlockEntries = NumBlockEntriesPerRow_[i];
    int RowDim = ElementSizeList_[i];
    for (int j=0; j< NumBlockEntries; j++) {
      int LDA = LDAs_[i][j];
      int ColDim = ColDims_[i][j];
      for (int col=0; col < ColDim; col++) {
	double * Values = Values_[i][j]+col*LDA;
	for (int row=0; row < RowDim; row++)
	  *Values++ = Scalar;
      }
    }
  }
  return(0);
}
//==========================================================================
int Petra_RDP_VBR_Matrix::BeginInsertGlobalValues(int BlockRow, int NumBlockEntries, int * BlockIndices) {

  if (IndicesAreLocal()) return(-2); // Cannot insert global values into local graph
  Graph_->SetIndicesAreGlobal(true);
  BlockRow = LRID(BlockRow); // Find local row number for this global row index
  
  bool IndicesAreLocal = false;
  return(BeginInsertValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));

}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginInsertMyValues(int  BlockRow, int NumBlockEntries, int * BlockIndices) {

  if (IndicesAreGlobal()) return(-2); // Cannot insert global values into filled graph
  Graph_->SetIndicesAreLocal(true);
  bool IndicesAreLocal = true;
  return(BeginInsertValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));

}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginInsertValues(int BlockRow, int NumBlockEntries, 
					    int * BlockIndices, bool IndicesAreLocal) {

  if (StaticGraph()) return(-2); // If the matrix graph is fully constructed, we cannot insert new values

  int j;
  int ierr = 0;

  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) return(-1); // Not in BlockRow range    
  if (CV_==View && Values_[BlockRow]!=0) ierr = 2; // This row has be defined already. Issue warning.    
  if (IndicesAreContiguous()) return(-3); // Indices cannot be individually deleted and new

  // Set up pointers, make sure enough temp space for this round of submits

  Petra_CombineMode SubmitMode = Insert;

  return(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
}
//==========================================================================
int Petra_RDP_VBR_Matrix::BeginReplaceGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

   BlockRow = LRID(BlockRow); // Normalize row range
   bool IndicesAreLocal = false;
  return(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginReplaceMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  bool IndicesAreLocal = true;
  return(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginReplaceValues(int BlockRow, int NumBlockEntries, 
					     int *BlockIndices, bool IndicesAreLocal) {

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.
  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) return(-1); // Not in BlockRow range

  Petra_CombineMode SubmitMode = Replace;
  return(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginSumIntoGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

   BlockRow = LRID(BlockRow); // Normalize row range
   bool IndicesAreLocal = false;
  return(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginSumIntoMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  bool IndicesAreLocal = true;
  return(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginSumIntoValues(int BlockRow, int NumBlockEntries, 
					     int *BlockIndices, bool IndicesAreLocal) {

  if (CV_==View) return(-3); // This is a view only.  Cannot remove entries.
  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) return(-1); // Not in BlockRow range

  Petra_CombineMode SubmitMode = Add;
  return(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::SetupForSubmits(int BlockRow, int NumBlockEntries, int * BlockIndices, 
					  bool IndicesAreLocal, Petra_CombineMode SubmitMode) {

  if (NumBlockEntries>LenTemps_) {
    if (LenTemps_>0) {
      delete [] TempRowDims_;
      delete [] TempColDims_;
      delete [] TempLDAs_;
      delete [] TempValues_;
    }
    TempRowDims_ = new int[NumBlockEntries];
    TempColDims_ = new int[NumBlockEntries];
    TempLDAs_ = new int[NumBlockEntries];
    TempValues_ = new double*[NumBlockEntries];
    LenTemps_ = NumBlockEntries;
  }

  CurBlockRow_ = BlockRow;
  CurNumBlockEntries_ = NumBlockEntries;
  CurBlockIndices_ = BlockIndices;
  CurEntry_ = 0;
  CurIndicesAreLocal_ = IndicesAreLocal;
  CurSubmitMode_ = SubmitMode;
    
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::SubmitBlockEntry(double *Values, int LDA, int NumRows, int NumCols) {

  if (CurEntry_==-1) return(-1); // This means that a Begin routine was not called
  if (CurEntry_>=CurNumBlockEntries_) return(-4); // Exceeded the number of entries that can be submitted

  // Fill up temp space with entry

  TempColDims_[CurEntry_] = NumCols;
  TempRowDims_[CurEntry_] = NumRows;
  TempLDAs_[CurEntry_] = LDA;
  TempValues_[CurEntry_] = Values;
  CurEntry_++;
  
  return(0);
}
//==========================================================================
int Petra_RDP_VBR_Matrix::EndSubmitEntries() {

  if (CurEntry_!=CurNumBlockEntries_) return(-6); // Did not submit right number of entries

  if (CurSubmitMode_==Insert) return(EndInsertValues());
  else return(EndReplaceSumIntoValues());
}
//==========================================================================
int Petra_RDP_VBR_Matrix::EndReplaceSumIntoValues() {

  int j;
  int ierr = 0;
  int Loc;

  int RowDim = ElementSizeList_[CurBlockRow_];

  bool SumInto = (CurSubmitMode_==Add);

  if (CurIndicesAreLocal_) {
    for (j=0; j<CurNumBlockEntries_; j++) {
      int BlockIndex = CurBlockIndices_[j];
      if (Graph_->FindMyIndexLoc(CurBlockRow_,BlockIndex,j,Loc)) 
	CopyMat(TempValues_[j], TempLDAs_[j], RowDim, TempColDims_[j], 
		Values_[CurBlockRow_][Loc], LDAs_[CurBlockRow_][Loc], SumInto);
      else return(-2); // Value not found
    }
  }
  else {
    for (j=0; j<CurNumBlockEntries_; j++) {
      int BlockIndex = CurBlockIndices_[j];
      if (Graph_->FindGlobalIndexLoc(CurBlockRow_,BlockIndex,j,Loc)) 
	CopyMat(TempValues_[j], TempLDAs_[j], RowDim, TempColDims_[j], 
		Values_[CurBlockRow_][Loc], LDAs_[CurBlockRow_][Loc], SumInto);
      else return(-2); // Value not found
    }
  }
  return(ierr);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::EndInsertValues() {

  int ierr = 0;
  int j;
  int start = NumBlockEntriesPerRow_[CurBlockRow_];
  int stop = start + CurNumBlockEntries_;
  int NumAllocatedEntries = NumAllocatedBlockEntriesPerRow_[CurBlockRow_];
  if (stop > NumAllocatedEntries){
    if (NumAllocatedEntries==0) { // BlockRow was never allocated, so do it
      Values_[CurBlockRow_] = new double*[CurNumBlockEntries_];
      ColDims_[CurBlockRow_] = new int[CurNumBlockEntries_];
      LDAs_[CurBlockRow_] = new int[CurNumBlockEntries_];
    }
    else {
      ierr = 1; // Out of room.  Must delete and allocate more space...
      double ** tmp_Values = new double*[stop];
      int * tmp_ColDims = new int[stop];
      int * tmp_LDAs = new int[stop];
      for (j=0; j< start; j++) {
	tmp_Values[j] = Values_[CurBlockRow_][j]; // Copy existing entries
	tmp_ColDims[j] = ColDims_[CurBlockRow_][j];
	tmp_LDAs[j] = LDAs_[CurBlockRow_][j];
      }
      delete [] Values_[CurBlockRow_]; // Delete old storage
      delete [] ColDims_[CurBlockRow_];
      delete [] LDAs_[CurBlockRow_];
      
      Values_[CurBlockRow_] = tmp_Values; // Set pointer to new storage
      ColDims_[CurBlockRow_] = tmp_ColDims;
      LDAs_[CurBlockRow_] = tmp_LDAs;
    }
  }
  if (CV_==View) {
    for (j=start; j<stop; j++) {
      Values_[CurBlockRow_][j] = TempValues_[j-start];
      ColDims_[CurBlockRow_][j] = TempColDims_[j-start];
      LDAs_[CurBlockRow_][j] = TempLDAs_[j-start];
    }
  }
  else { // Copy not view
    int RowDim = ElementSizeList_[CurBlockRow_];
    for (j=start; j<stop; j++) {
      int ColDim =  TempColDims_[j-start];
      int LDA = TempLDAs_[j-start];
      Values_[CurBlockRow_][j] = new double[RowDim*ColDim];
      CopyMat(TempValues_[j-start], LDA, RowDim, ColDim, Values_[CurBlockRow_][j], RowDim, false);
      ColDims_[CurBlockRow_][j] = ColDim;
      LDAs_[CurBlockRow_][j] = RowDim;
    }
  }
    
  Graph_->InsertIndices(CurBlockRow_, CurNumBlockEntries_, CurBlockIndices_); // Update graph
  
  return(ierr);
  
}

//=============================================================================
int Petra_RDP_VBR_Matrix::CopyMat(double * A, int LDA, int NumRows, int NumCols, 
					double * B, int LDB, bool SumInto) const {

  int i, j;
  double * ptr1 = B;
  double * ptr2;

  if (LDB<NumRows) return(-1); // Stride of B is not large enough

  if (SumInto) { // Add to existing values
    for (j=0; j<NumCols; j++) {
      ptr1 = B + j*LDB;
      ptr2 = A + j*LDA;
      for (i=0; i<NumRows; i++) *ptr1++ += *ptr2++;
    }
  }
  else {  // Replace values
    for (j=0; j<NumCols; j++) {
      ptr1 = B + j*LDB;
      ptr2 = A + j*LDA;
      for (i=0; i<NumRows; i++) *ptr1++ = *ptr2++;
    }
  }
  return(0);
}
//==========================================================================
int Petra_RDP_VBR_Matrix::TransformToLocal() {
  return(TransformToLocal((Petra_BlockMap *) (&RowMap()), (Petra_BlockMap *) (&ColMap())));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::TransformToLocal(Petra_BlockMap *DomainMap, Petra_BlockMap *RangeMap) {
  
  if (!StaticGraph()) Graph_->MakeIndicesLocal(*DomainMap, *RangeMap);
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) Graph_->TransformToLocal(DomainMap, RangeMap);

  // NumMyCols_ = Graph_->NumMyCols(); // Redefine based on local number of cols

  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::SortEntries() {

  if (!IndicesAreLocal()) return(-1);
  if (Sorted()) return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for (int i=0; i<NumMyBlockRows_; i++){

    double ** Values = Values_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
    int * ColDims = ColDims_[i];
    int * LDAs = LDAs_[i];
    int n = NumEntries;
    int m = n/2;
    
    while (m > 0) {
      int max = n - m;
      for (int j=0; j<max; j++)
        {
	  for (int k=j; k>=0; k-=m)
            {
	      if (Indices[k+m] >= Indices[k])
		break;
	      double *dtemp = Values[k+m];
	      Values[k+m] = Values[k];
	      Values[k] = dtemp;

	      int itemp = Indices[k+m];
	      Indices[k+m] = Indices[k];
	      Indices[k] = itemp;

	      itemp = ColDims[k+m];
	      ColDims[k+m] = ColDims[k];
	      ColDims[k] = itemp;

	      itemp = LDAs[k+m];
	      LDAs[k+m] = LDAs[k];
	      LDAs[k] = itemp;
            }
        }
      m = m/2;
    }
  }
  Graph_->SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::MergeRedundantEntries() {

  int i, j, k;

  if (NoRedundancies()) return(0);
  if (!Sorted()) return(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortEntries was already called.

  bool SumInto = true;
  for (i=0; i<NumMyBlockRows_; i++){
    int NumEntries = NumBlockEntriesPerRow_[i];
    if (NumEntries>0) {
      double ** const Values = Values_[i];
      int * const Indices = Indices_[i];
      int * const LDAs = LDAs_[i];
      int * const ColDims = ColDims_[i];
      int RowDim = ElementSizeList_[i];
      int j0 = 0;
      int jj0 = Indices[j0];
      for (j=1; j<NumEntries; j++) {
	int jj = Indices[j];
	int ColDim = ColDims[j];
	if (jj==jj0) {// Check if index is repeated
	  // Values[j0] += Values[j];
	  CopyMat(Values[j], LDAs[j], RowDim, ColDim, Values[j0], LDAs[j0], SumInto);
	  delete [] Values[j]; // Get rid of space
	  for (k=j; k<NumEntries-1; k++) Values[k] = Values[k+1]; // Sum up values
	  Values[NumEntries-1] = 0; // No space associated with this pointer anymore
	  NumEntries--;
	}
	else {
	  j0=j; // Redefine comparison index value
	  jj0=Indices[j0];
	}
      }
    }
  }
    
  return(Graph_->RemoveRedundantIndices()); // Remove redundant indices and then return


}

//==========================================================================
int Petra_RDP_VBR_Matrix::OptimizeStorage() {

  /* Work on later...
  int i, j;

  // The purpose of this routine is to make the block entries in each row contiguous in memory
  // so that a single call to GEMV or GEMM call be used to compute an entire block row.

  if (StorageOptimized()) return(0); // Have we been here before?

  bool Contiguous = true; // Assume contiguous is true
  for (i=1; i<NumMyBlockRows_; i++){
    int NumEntries = NumBlockEntriesPerRow_[i];
    int NumAllocatedEntries = NumAllocatedBlockEntriesPerRow_[i];
      
    // Check if NumEntries is same as NumAllocatedEntries and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if ((NumEntries!=NumAllocatedEntries) || (Values_[i]!=Values_[i-1]+NumEntries)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumEntries and NumAllocatedEntries
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) return(-1);  // This is user data, it's not contiguous and we can't make it so.

  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  if (ierr) return(ierr);

  if (Contiguous) return(0); // Everything is done.  Return

  // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
  int NumMyNonzeros = Graph_->NumMyNonzeros();

  // Allocate one big integer array for all index values
  All_Values_ = new double[NumMyNonzeros];
  
  // Set Entries_ to point into All_Entries_
  
  double * tmp = All_Values_;
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumEntries = NumEntriesPerBlockRow_[i];
    for (j=0; j<NumEntries; j++) tmp[j] = Values_[i][j];
    if (Values_[i] !=0) delete [] Values_[i];
    Values_[i] = tmp;
    tmp += NumEntries;
  }
  */
  return(0);
}
//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractGlobalRowCopy(int GlobalRow, int Length, 
                   int & NumEntries, double *Values, int * Indices) const {

  cout << "Must implement..." << endl;
  return(0);
}
//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractGlobalBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
							    int & RowDim, int & NumBlockEntries, 
							    int * BlockIndices, int * & ColDims, 
							    int * & LDAs, double ** & Values) const {

  bool IndicesAreLocal = false;
  return(ExtractBlockRowPointers(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
				     ColDims, LDAs, Values, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractMyBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
							int & RowDim, int & NumBlockEntries, 
							int * BlockIndices, int * & ColDims, 
							int * & LDAs, double ** & Values) const {

  bool IndicesAreLocal = true;
  return(ExtractBlockRowPointers(BlockRow,MaxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
				     ColDims, LDAs, Values, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
						      int & RowDim, int & NumBlockEntries, 
						      int * BlockIndices, int * & ColDims, 
						      int * & LDAs, double ** & Values, 
						      bool IndicesAreLocal) const {
  int ierr = 0;
  if (!IndicesAreLocal) {
    ierr = Graph_->ExtractGlobalRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);
    BlockRow = LRID(BlockRow);
  }
  else
    ierr = Graph_->ExtractMyRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);

  if (ierr) return(ierr);

  RowDim = ElementSizeList_[BlockRow];

  ColDims = ColDims_[BlockRow];
  LDAs = LDAs_[BlockRow];
  Values = Values_[BlockRow];


  return(ierr);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractGlobalBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
							 int & RowDim, int & NumBlockEntries, 
							 int * BlockIndices, int * ColDims) const {

  bool IndicesAreLocal = false;
  return(BeginExtractBlockRowCopy(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractMyBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
						     int & RowDim, int & NumBlockEntries, 
						     int * BlockIndices, int * ColDims) const {

  bool IndicesAreLocal = true;
  return(BeginExtractBlockRowCopy(BlockRow,MaxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
						   int & RowDim, int & NumBlockEntries, 
						   int * BlockIndices, int * ColDims, 
						   bool IndicesAreLocal) const  {
  int ierr = 0;
  if (!IndicesAreLocal)
    ierr = Graph_->ExtractGlobalRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);
  if (ierr) return(ierr);

  bool ExtractView = false;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, IndicesAreLocal);
  if (ierr) return(ierr);

  return(ExtractBlockDimsCopy(NumBlockEntries, ColDims));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::SetupForExtracts(int BlockRow, int & RowDim, int NumBlockEntries, bool ExtractView, 
						   bool IndicesAreLocal) const {

  if (!IndicesAreLocal) BlockRow = LRID(BlockRow); // Normalize row range
  CurExtractBlockRow_ = BlockRow;
  CurExtractEntry_ = 0;
  CurExtractNumBlockEntries_ = NumBlockEntries;
  CurExtractIndicesAreLocal_ = IndicesAreLocal;
  CurExtractView_ = ExtractView;
  CurRowDim_ = ElementSizeList_[CurBlockRow_];
  RowDim = CurRowDim_;
    
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractBlockDimsCopy(int NumBlockEntries, int * ColDims) const {

  int * CurColDims = ColDims_[CurExtractBlockRow_];

  for (int i=0; i<NumBlockEntries; i++) ColDims[i] = CurColDims[i];
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const
{
  if (CurExtractEntry_==-1) return(-1); // No BeginCopy routine was called
  int CurColDim = ColDims_[CurExtractBlockRow_][CurExtractEntry_];
  if (LDA*CurColDim>SizeOfValues) return(-2);  // Not enough space

  double * CurValues = Values_[CurExtractBlockRow_][CurExtractEntry_];
  int CurLDA = LDAs_[CurExtractBlockRow_][CurExtractEntry_];

  CurExtractEntry_++; // Increment Entry Pointer

  return(CopyMat(CurValues, CurLDA, CurRowDim_, CurColDim, Values, LDA, SumInto));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
							 int * & BlockIndices, int * & ColDims, int * & LDAs) const
{

  bool IndicesAreLocal = false;
  return(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						     int * & BlockIndices, int * & ColDims, int * & LDAs)  const
{

  bool IndicesAreLocal = true;
  return(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::BeginExtractBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						   int * & BlockIndices, int * & ColDims, 
						   int * & LDAs, bool IndicesAreLocal) const
{
  int ierr = 0;
  if (!IndicesAreLocal)
    ierr = Graph_->ExtractGlobalRowView(BlockRow, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowView(BlockRow,  NumBlockEntries, BlockIndices);
  if (ierr) return(ierr);

  bool ExtractView = true;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, IndicesAreLocal);
  if (ierr) return(ierr);

  return(ExtractBlockDimsView(NumBlockEntries, ColDims, LDAs));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractBlockDimsView(int NumBlockEntries, int * & ColDims, int * & LDAs) const {

  ColDims = ColDims_[CurExtractBlockRow_];
  LDAs = LDAs_[CurExtractBlockRow_];
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractEntryView(double * & Values) const
{
  Values = Values_[CurExtractBlockRow_][CurExtractEntry_];

  CurExtractEntry_++; // Increment Entry Pointer
  return(0);
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						    int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const
{

  Values = Values_[LRID(BlockRow)]; // Pointer to Array of pointers for this row's block entries
  bool IndicesAreLocal = false;
  return(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
}

//==========================================================================
int Petra_RDP_VBR_Matrix::ExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const
{

  Values = Values_[BlockRow]; // Pointer to Array of pointers for this row's block entries
  bool IndicesAreLocal = true;
  return(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
}

//==============================================================================
int Petra_RDP_VBR_Matrix::ExtractDiagonalCopy(Petra_RDP_Vector & Diagonal) const {
	
  if (!Filled()) return(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) return(-2); // Maps must be the same
  double * diagptr = Diagonal.Values();
  int Base = IndexBase();
  for(int i=0; i<NumMyBlockRows_; i++){
    int BlockRow = i + Base;
    int RowDim = ElementSizeList_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
    for (int j=0; j<NumEntries; j++) {
      int BlockCol = Indices[j];
      if (BlockRow==BlockCol) {
	CopyMatDiag(Values_[i][j], LDAs_[i][j], RowDim, ColDims_[i][j], 
		    diagptr+FirstElementEntryList_[i]);
	break;
      }
    }
  }
  return(0);
}
//==============================================================================
int Petra_RDP_VBR_Matrix::BeginExtractBlockDiagonalCopy(int MaxNumBlockDiagonalEntries, 
							int & NumBlockDiagonalEntries, int * RowColDims ) const{
	
  if (!Filled()) return(-1); // Can't get diagonal unless matrix is filled
  CurBlockDiag_ = 0; // Initialize pointer
  NumBlockDiagonalEntries = NumMyBlockRows_;
  if (NumBlockDiagonalEntries>MaxNumBlockDiagonalEntries) return(-1);
  return(RowMap().ElementSizeList(RowColDims));
}
//==============================================================================
int Petra_RDP_VBR_Matrix::ExtractBlockDiagonalEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto ) const {
	
  if (CurBlockDiag_==-1) return(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow =  + IndexBase();
  int RowDim = ElementSizeList_[i];
  int NumEntries = NumBlockEntriesPerRow_[i];
  int * Indices = Indices_[i];
  for (int j=0; j<NumEntries; j++) {
    int Col = Indices[j];
    if (BlockRow==Col) {
      int ColDim = ColDims_[i][j];
      if (LDA*ColDim>SizeOfValues) return(-2); // Not enough room in Values
      CopyMat(Values_[i][j], LDAs_[i][j], RowDim, ColDim, Values, LDA, SumInto);
      break;
    }
  }
  CurBlockDiag_++; // Increment counter
  return(0);
}
//==============================================================================
int Petra_RDP_VBR_Matrix::BeginExtractBlockDiagonalView(int & NumBlockDiagonalEntries, 
							int * & RowColDims ) const {
	
  if (!Filled()) return(-1); // Can't get diagonal unless matrix is filled
  CurBlockDiag_ = 0; // Initialize pointer
  NumBlockDiagonalEntries = NumMyBlockRows_;
  RowColDims = ElementSizeList_;
  return(0);
}
//==============================================================================
int Petra_RDP_VBR_Matrix::ExtractBlockDiagonalEntryView(double * & Values, int & LDA) const {
	
  if (CurBlockDiag_==-1) return(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow =  + IndexBase();
  int RowDim = ElementSizeList_[i];
  int NumEntries = NumBlockEntriesPerRow_[i];
  int * Indices = Indices_[i];
  for (int j=0; j<NumEntries; j++) {
    int Col = Indices[j];
    if (BlockRow==Col) {
      Values = Values_[i][j];
      LDA = LDAs_[i][j];
      break;
    }
  }
  CurBlockDiag_++; // Increment counter
  return(0);
}
//=============================================================================
int Petra_RDP_VBR_Matrix::CopyMatDiag(double * A, int LDA, int NumRows, int NumCols, 
					double * Diagonal) const {

  int i, j;
  double * ptr1 = Diagonal;
  double * ptr2;
  int ndiags = minfn(NumRows,NumCols);

  for (i=0; i<ndiags; i++) {
    ptr2 = A + i*LDA+i;
    *ptr1++ = *ptr2;
  }
  return(0);
}
//=============================================================================
int Petra_RDP_VBR_Matrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {
  if (!Filled()) return(-1); // Can't row unless matrix is filled
  if (!IndicesAreLocal()) return(-2);

  int ierr = 0;
  int BlockRow, BlockOffset;
  ierr = RowMap().FindLocalBlockID(MyRow, BlockRow, BlockOffset);  if (ierr!=0) return(ierr);

  int RowDim, NumBlockEntries;
  int * BlockIndices, * ColDims, * LDAs;
  double ** ValBlocks;
  ierr = ExtractMyBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, ColDims, LDAs, ValBlocks);
  if (ierr!=0) return(ierr);

  int * ColFirstElementEntryList = FirstElementEntryList_;
  if (Importer()!=0) ColFirstElementEntryList = ImportMap().FirstElementEntryList();
  NumEntries = 0;
  for (int i=0; i<NumBlockEntries; i++) {
    int ColDim = ColDims[i];
    NumEntries += ColDim;
    if (NumEntries>Length) return(-3); // Not enough space
    double * A = ValBlocks[i] + BlockOffset; // Point to first element in row
    int LDA = LDAs[i];
    int Index = ColFirstElementEntryList[BlockIndices[i]];
    for (int j=0; j < ColDim; j++) {
      *Values++ = *A;
      A += LDA;
      *Indices++ = Index++;
    }
  }
  return(0);
      
}
//=============================================================================
int Petra_RDP_VBR_Matrix::Multiply1(bool TransA, const Petra_RDP_Vector& x, Petra_RDP_Vector& y) const {
//
// This function forms the product y = A * x or y = A' * x
//

  if (!Filled()) return (-1); // Matrix must be filled.
  
  int i, j;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * FirstElementEntry = FirstElementEntryList_;
  int * ElementSize = ElementSizeList_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;

  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstElementEntryList = FirstElementEntryList_;



  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),1); // Create import vector if needed
      ImportVector_->Import(x, *Importer(), Insert);
      xp = (double*)ImportVector_->Values();
      ColElementSizeList = ImportMap().ElementSizeList(); // The Import map will always have an existing ElementSizeList
      ColFirstElementEntryList = ImportMap().FirstElementEntryList(); // Import map will always have an existing ...
    }
    

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*)ExportVector_->Values();
    }
    
    // Do actual computation
    int NumMyRows_ = NumMyRows();
    for (i=0; i< NumMyRows_; i++) yp[i] = 0.0;  // Initialize y
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      double * cury = yp + *FirstElementEntry++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
	//sum += BlockRowValues[j] * xp[BlockRowIndices[j]];
	double * A = BlockRowValues[j];
	int LDA = BlockRowLDAs[j];
	int Index = BlockRowIndices[j];
	double * curx = xp + ColFirstElementEntryList[Index];
	int ColDim = ColElementSizeList[Index];
	GEMV('N', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);
	
      }
      if (Exporter()!=0) y.Export(*ExportVector_, *Exporter(), Add); // Fill y with Values from export vector
    }
  }
  
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),1); // Create Export vector if needed
      ExportVector_->Import(x, *Exporter(), Insert);
      xp = (double*)ExportVector_->Values();
    }
    
    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),1); // Create import vector if needed
      yp = (double*)ImportVector_->Values();
      ColElementSizeList = ImportMap().ElementSizeList(); // The Import map will always have an existing ElementSizeList
      ColFirstElementEntryList = ImportMap().FirstElementEntryList(); // Import map will always have an existing ...
    }
    
    // Do actual computation
    int NumMyCols_ = NumMyCols();
    for (i=0; i < NumMyCols_; i++) yp[i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      double * curx = xp + *FirstElementEntry++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
	//yp[BlockRowIndices[j]] += BlockRowValues[j] * xp[i];
	double * A = BlockRowValues[j];
	int LDA = BlockRowLDAs[j];
	int Index = BlockRowIndices[j];
	double * cury = yp + ColFirstElementEntryList[Index];
	int ColDim = ColElementSizeList[Index];
	GEMV('T', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);
	
      }
    }
    if (Importer()!=0) y.Export(*ImportVector_, *Importer(), Add); // Fill y with Values from export vector
  }
  
  UpdateFlops(2*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_VBR_Matrix::Multiply(bool TransA, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  
  if (!Filled()) return (-1); // Matrix must be filled.
  
  int i, j, k;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstElementEntryList = FirstElementEntryList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstElementEntryList = FirstElementEntryList_;

   
  int NumVectors = X.NumVectors();
  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  
  if (!TransA) {
    
    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()<NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
       // Create import vector if needed
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),NumVectors);

      ImportVector_->Import(X, *Importer(), Insert);
      Xp = (double**)ImportVector_->Pointers();
      ColElementSizeList = ImportMap().ElementSizeList();
      ColFirstElementEntryList = ImportMap().FirstElementEntryList();
    }
    
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
       // Create Export vector if needed
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),NumVectors);

      ExportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ExportVector_->Pointers();
      RowElementSizeList = ImportMap().ElementSizeList();
      RowFirstElementEntryList = ImportMap().FirstElementEntryList();
    }
    else
      Y.PutScalar(0.0); // Zero y values
    
    // Do actual computation
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      int  yoff = *RowFirstElementEntryList++;
      int RowDim = *RowElementSizeList++;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
		       ColFirstElementEntryList, ColElementSizeList, 
		       1.0, BlockRowValues, BlockRowLDAs, Xp, 1.0, Yp, NumVectors);
    }
    if (Exporter()!=0) Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
       // Create Export vector if needed
      if (ExportVector_==0) ExportVector_ = new Petra_RDP_MultiVector(RowMap(),NumVectors);

      ExportVector_->Import(X, *Exporter(), Insert);
      Xp = (double**)ExportVector_->Pointers();
      ColElementSizeList = ExportMap().ElementSizeList();
      ColFirstElementEntryList = ExportMap().FirstElementEntryList();
    }
  
    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()<NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      // Create import vector if needed
      if (ImportVector_==0) ImportVector_ = new Petra_RDP_MultiVector(ImportMap(),NumVectors);

      ImportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ImportVector_->Pointers();
      RowElementSizeList = ImportMap().ElementSizeList();
      RowFirstElementEntryList = ImportMap().FirstElementEntryList();
    }
    else
      Y.PutScalar(0.0); // Zero y values

    // Do actual computation
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      int  xoff = *ColFirstElementEntryList++;
      int RowDim = *ColElementSizeList++;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, xoff, 
		       RowFirstElementEntryList, RowElementSizeList, 
		       1.0, BlockRowValues, BlockRowLDAs, Xp, 1.0, Yp, NumVectors);
    }

    if (Importer()!=0) Y.Export(*ImportVector_, *Importer(), Add); // Fill Y with Values from export vector
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}
//=============================================================================
void Petra_RDP_VBR_Matrix::BlockRowMultiply(bool TransA, int RowDim, int NumEntries, 
		      int * BlockIndices, int RowOff,
		      int * FirstElementEntryList, int * ElementSizeList,
		      double Alpha, double ** As, int * LDAs, 
		      double ** X, double Beta, double ** Y, int NumVectors) const {
  int j, k;

  if (!TransA) {
    for (j=0; j < NumEntries; j++) {
      double * A = As[j];
      int LDA = LDAs[j];
      int BlockIndex = BlockIndices[j];
      int xoff = FirstElementEntryList[BlockIndex];
      int ColDim = ElementSizeList[BlockIndex];
      for (k=0; k<NumVectors; k++) {
	double * curx = X[k] + xoff;
	double * cury = Y[k] + RowOff;
	GEMV('N', RowDim, ColDim, Alpha, A, LDA, curx, Beta, cury);
      }
    }
  }
  else {
    for (j=0; j < NumEntries; j++) {
      double * A = As[j];
      int LDA = LDAs[j];
      int BlockIndex = BlockIndices[j];
      int yoff = FirstElementEntryList[BlockIndex];
      int ColDim = ElementSizeList[BlockIndex];
      for (k=0; k<NumVectors; k++) {
	double * curx = X[k] + RowOff;
	double * cury = Y[k] + yoff;
	GEMV('T', RowDim, ColDim, Alpha, A, LDA, curx, Beta, cury);
      }
    }
  }

    return;
}
//=============================================================================
int Petra_RDP_VBR_Matrix::Solve(bool Upper, bool TransA, bool UnitDiagonal, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const {
  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //

  if (!Filled()) return (-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) return (-2);
  if ((!Upper) && (!LowerTriangular())) return (-3);
  if (!NoDiagonal()) return (-4); // We must use UnitDiagonal

  int i, j, j0, k;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * FirstElementEntry = FirstElementEntryList_;
  int * ElementSize = ElementSizeList_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;

  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstElementEntryList = FirstElementEntryList_;

  // If upper, point to last row
  if (Upper) {
    NumBlockEntriesPerRow += NumMyBlockRows_-1;
    FirstElementEntry += NumMyBlockRows_-1;
    ElementSize += NumMyBlockRows_-1;
    LDAs += NumMyBlockRows_-1;
    Indices += NumMyBlockRows_-1;
    Values += NumMyBlockRows_-1;
  }

  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();

  if (X.Pointers() != Y.Pointers()) Y = X; // Copy X into Y if they are not the same multivector

  bool Case1 = (((!TransA) && Upper) ||  (TransA && !Upper)); 
  // Case 2 = (((TransA) && Upper) || (!TransA) && Lower);
  if (Case1) {
      for (i=0; i < NumMyBlockRows_; i++) {
	int      NumEntries = *NumBlockEntriesPerRow--;
	int *    BlockRowIndices = *Indices--;
	double ** BlockRowValues  = *Values--;
	int *    BlockRowLDAs = *LDAs--;
	int  yoff = *FirstElementEntry--;
	int RowDim = *ElementSize--;
	BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
			 ColFirstElementEntryList, ColElementSizeList, 
			 1.0, BlockRowValues, BlockRowLDAs, Yp, -1.0, Yp, NumVectors);
      }
    }
    else {
      for (i=0; i < NumMyBlockRows_; i++) {
	int      NumEntries = *NumBlockEntriesPerRow++;
	int *    BlockRowIndices = *Indices++;
	double ** BlockRowValues  = *Values++;
	int *    BlockRowLDAs = *LDAs++;
	int  yoff = *FirstElementEntry++;
	int RowDim = *ElementSize++;
	BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
			 ColFirstElementEntryList, ColElementSizeList, 
			 1.0, BlockRowValues, BlockRowLDAs, Yp, -1.0, Yp, NumVectors);
      }
    }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Petra_RDP_VBR_Matrix::InvRowSums(Petra_RDP_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//
  return(InverseSums(true, x));
}

//=============================================================================
int Petra_RDP_VBR_Matrix::InvColSums(Petra_RDP_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//
  return(InverseSums(false, x));
}
//=============================================================================
int Petra_RDP_VBR_Matrix::InverseSums(bool DoRows, Petra_RDP_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!Filled()) return (-1); // Matrix must be filled.
  if (DoRows) {
    if (!Graph().RangeMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the range of A
  }
  else {
    if (!Graph().DomainMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the domain of A
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstElementEntryList = FirstElementEntryList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstElementEntryList = FirstElementEntryList_;
  if (Importer()!=0) {
    ColElementSizeList = ImportMap().ElementSizeList();
    ColFirstElementEntryList = ImportMap().FirstElementEntryList();
  }

  x.PutScalar(0.0); // Zero out result vector

  double * xp = (double*)x.Values();

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  Petra_RDP_Vector * x_tmp = 0;
  if (!DoRows) {
    if (Importer()!=0) {
      x_tmp = new Petra_RDP_Vector(ImportMap()); // Create import vector if needed
      xp = (double*)x_tmp->Values();
    }
  }

  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    double ** BlockRowValues  = *Values++;
    int *    BlockRowLDAs = *LDAs++;
    int xoff = *RowFirstElementEntryList++;
    int RowDim = *RowElementSizeList++;
    if (DoRows) {
      for (int ii=0; ii < NumEntries; ii++) {
	double * x = xp+xoff;
	double * A = BlockRowValues[ii];
	int LDA = BlockRowLDAs[ii];
	int BlockIndex = BlockRowIndices[ii];
	int ColDim = ColElementSizeList[BlockIndex];
	for (int j=0; j<ColDim; j++) {
	  double * curEntry = A + j*LDA;
	  for (int k=0; k<RowDim; k++)
	    x[k] += fabs(*curEntry++);
	}
      }
    }
    else {
      for (int ii=0; ii < NumEntries; ii++) {
	double * A = BlockRowValues[ii];
	int LDA = BlockRowLDAs[ii];
	int BlockIndex = BlockRowIndices[ii];
	int off = ColFirstElementEntryList[BlockIndex];
	int ColDim = ColElementSizeList[BlockIndex];
	double * curx = xp+off;
	for (int j=0; j<ColDim; j++) {
	  double * curEntry = A + j*LDA;
	  for (int k=0; k<RowDim; k++)
	    curx[j] += fabs(*curEntry++);
	}
      }
    }
  }

  if (!DoRows) {
    if (Importer()!=0){
      x.Export(*x_tmp, *Importer(), Add); // Fill x with Values from import vector
      delete x_tmp;
      xp = (double*) x.Values();
    }
  }
  int NumMyRows_ = NumMyRows();
  for (int i=0; i < NumMyRows_; i++) {
    double scale = xp[i];
    if (scale<Petra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero row/col sum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Petra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());

  return(ierr);
}
//=============================================================================
int Petra_RDP_VBR_Matrix::LeftScale(const Petra_RDP_Vector& x) {
//
// Multiply the ith row of A by x[i].
//
  return(Scale(true, x));
}

//=============================================================================
int Petra_RDP_VBR_Matrix::RightScale(const Petra_RDP_Vector& x) {
//
// Multiply the jth column of A by x[j].
//
  return(Scale (false, x));
}
//=============================================================================
int Petra_RDP_VBR_Matrix::Scale(bool DoRows, const Petra_RDP_Vector& x) {

  if (!Filled()) return (-1); // Matrix must be filled.
  if (DoRows) {
    if (!Graph().RangeMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the range of A
  }
  else {
    if (!Graph().DomainMap().SameAs(x.Map())) return(-2); // x must have the same distribution as the domain of A
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstElementEntryList = FirstElementEntryList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstElementEntryList = FirstElementEntryList_;
  if (Importer()!=0) {
    ColElementSizeList = ImportMap().ElementSizeList();
    ColFirstElementEntryList = ImportMap().FirstElementEntryList();
  }

  double * xp = (double*)x.Values();

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  Petra_RDP_Vector * x_tmp = 0;
  if (!DoRows) {
    if (Importer()!=0) {
      x_tmp = new Petra_RDP_Vector(ImportMap()); // Create import vector if needed
      x_tmp->Import(x,*Importer(), Insert); // x_tmp will have all the values we need
      xp = (double*)x_tmp->Values();
    }
  }

  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    double ** BlockRowValues  = *Values++;
    int *    BlockRowLDAs = *LDAs++;
    int xoff = *RowFirstElementEntryList++;
    int RowDim = *RowElementSizeList++;
    if (DoRows) {
      for (int ii=0; ii < NumEntries; ii++) {
	double * x = xp+xoff;
	double * A = BlockRowValues[ii];
	int LDA = BlockRowLDAs[ii];
	int BlockIndex = BlockRowIndices[ii];
	int ColDim = ColElementSizeList[BlockIndex];
	for (int j=0; j<ColDim; j++) {
	  double * curEntry = A + j*LDA;
	  for (int k=0; k<RowDim; k++)
	    *curEntry++ *= x[k];
	}
      }
    }
    else {
      for (int ii=0; ii < NumEntries; ii++) {
	double * A = BlockRowValues[ii];
	int LDA = BlockRowLDAs[ii];
	int BlockIndex = BlockRowIndices[ii];
	int off = ColFirstElementEntryList[BlockIndex];
	int ColDim = ColElementSizeList[BlockIndex];
	double * curx = xp+off;
	for (int j=0; j<ColDim; j++) {
	  double * curEntry = A + j*LDA;
	  for (int k=0; k<RowDim; k++)
	    *curEntry++ *= curx[j];
	}
      }
    }
  }

  if (x_tmp!=0) delete x_tmp;
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());

  return(ierr);
}
//=============================================================================
double Petra_RDP_VBR_Matrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  if (!Filled()) return (-1); // Matrix must be filled.

  int MaxRowDim_ = MaxRowDim();
  double * tempv = new double[MaxRowDim_];

  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * ElementSize = ElementSizeList_;
  int ** ColDims = ColDims_;
  int ** LDAs = LDAs_;
  double *** Values = Values_;

  double Local_NormInf = 0.0;
  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++ ;
    int RowDim = *ElementSize++;
    int * BlockRowColDims = *ColDims++;
    int * BlockRowLDAs = *LDAs++;
    double ** BlockRowValues  = *Values++;
    BlockRowNormInf(RowDim, NumEntries, BlockRowColDims, 
		    BlockRowLDAs, BlockRowValues, tempv);
    for (int j=0; j < RowDim; j++) Local_NormInf = maxfn(Local_NormInf, tempv[j]);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  delete [] tempv;
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
void Petra_RDP_VBR_Matrix::BlockRowNormInf(int RowDim, int NumEntries, 
					   int * ColDims, int * LDAs, double ** As, 
					   double * Y) const {
  int i, j, k;
  for (k=0; k<RowDim; k++) Y[k] = 0.0;

  for (i=0; i < NumEntries; i++) {
    double * A = As[i];
    int LDA = LDAs[i];
    int ColDim = ColDims[i];
    for (j=0; j<ColDim; j++) {
      for (k=0; k<RowDim; k++) Y[k] += fabs(A[k]);
      A += LDA;
    }
  }
  return;
}
//=============================================================================
double Petra_RDP_VBR_Matrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);
  int * ColFirstElementEntryList = FirstElementEntryList_;
  if (Importer()!=0) {
    ColFirstElementEntryList = ImportMap().FirstElementEntryList();
  }

  Petra_RDP_Vector * x = new Petra_RDP_Vector(RowMap()); // Need temp vector for column sums
  
  double * xp = (double*)x->Values();
  Petra_RDP_MultiVector * x_tmp = 0;

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Petra_RDP_Vector(ImportMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }

  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * ElementSize = ElementSizeList_;
  int ** ColDims = ColDims_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;

  for (int i=0; i < NumMyBlockRows_; i++) {
    int NumEntries = *NumBlockEntriesPerRow++;
    int RowDim = *ElementSize++;
    int * BlockRowColDims = *ColDims++;
    int * BlockRowLDAs = *LDAs++;
    int *    BlockRowIndices = *Indices++;
    double ** BlockRowValues  = *Values++;
    BlockRowNormOne(RowDim, NumEntries, BlockRowIndices, BlockRowColDims, 
		    BlockRowLDAs, BlockRowValues,  ColFirstElementEntryList, xp);
  }
  if (Importer()!=0) x->Export(*x_tmp, *Importer(), Add); // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=============================================================================
void Petra_RDP_VBR_Matrix::BlockRowNormOne(int RowDim, int NumEntries, int * BlockRowIndices,
					   int * ColDims, int * LDAs, double ** As, 
					   int * ColFirstElementEntryList, double * x) const {
  int i, j, k;

  for (i=0; i < NumEntries; i++) {
    double * A = As[i];
    int LDA = LDAs[i];
    int ColDim = ColDims[i];
    double * curx = x + ColFirstElementEntryList[BlockRowIndices[i]];
    for (j=0; j<ColDim; j++) {
      for (k=0; k<RowDim; k++) curx[j] += fabs(A[k]);
      A += LDA;
    }
  }
  return;
}
//=========================================================================
int Petra_RDP_VBR_Matrix::Import(const Petra_RDP_VBR_Matrix& SourceMatrix, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode) {

  int ierr = 0;

  if (!RowMap().SameAs(Importer.TargetMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Importer.SourceMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumRemoteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *RemoteLIDs = Importer.RemoteLIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumNonzeros();
  int IntSizeOfPacket = 2*SourceMatrix.GlobalMaxNumBlockEntries() + 3;
  int Nsend = SizeOfPacket * Importer.NumSend();
  int Nrecv = SizeOfPacket * Importer.NumRecv();
  int IntNsend = IntSizeOfPacket * Importer.NumSend();
  int IntNrecv = IntSizeOfPacket * Importer.NumRecv();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, 
		    Nsend, Nrecv, SizeOfPacket, IntNsend, IntNrecv, IntSizeOfPacket,
		    LenExports_, Exports_, LenIntExports_, IntExports_, 
		    LenImports_, Imports_, LenIntImports_, IntImports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_RDP_VBR_Matrix::Export(const Petra_RDP_VBR_Matrix& SourceMatrix, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode){

  int ierr = 0;

  if (!RowMap().SameAs(Exporter.TargetMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Exporter.SourceMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumRemoteIDs();
  int NumExportIDs = Exporter.NumExportIDs();
  int *ExportLIDs = Exporter.ExportLIDs();
  int *RemoteLIDs = Exporter.RemoteLIDs();
  int *PermuteToLIDs = Exporter.PermuteToLIDs();
  int *PermuteFromLIDs = Exporter.PermuteFromLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumNonzeros();
  int IntSizeOfPacket = 2*SourceMatrix.GlobalMaxNumBlockEntries() + 3;
  int Nsend = SizeOfPacket * Exporter.NumSend();
  int Nrecv = SizeOfPacket * Exporter.NumRecv();
  int IntNsend = IntSizeOfPacket * Exporter.NumSend();
  int IntNrecv = IntSizeOfPacket * Exporter.NumRecv();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs,
		    Nsend, Nrecv, SizeOfPacket, IntNsend, IntNrecv, IntSizeOfPacket,
		    LenExports_, Exports_, LenIntExports_, IntExports_, 
		    LenImports_, Imports_, LenIntImports_, IntImports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    false));
}
//=========================================================================
int Petra_RDP_VBR_Matrix::Import(const Petra_RDP_VBR_Matrix& SourceMatrix, 
				 const Petra_Export & Exporter, Petra_CombineMode CombineMode) {

  int ierr = 0;

  if (!RowMap().SameAs(Exporter.SourceMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Exporter.TargetMap())) return(-3);

  int NumSameIDs = Exporter.NumSameIDs();
  int NumPermuteIDs = Exporter.NumPermuteIDs();
  int NumRemoteIDs = Exporter.NumExportIDs();
  int NumExportIDs = Exporter.NumRemoteIDs();
  int *ExportLIDs = Exporter.RemoteLIDs();
  int *RemoteLIDs = Exporter.ExportLIDs();
  int *PermuteToLIDs = Exporter.PermuteFromLIDs();
  int *PermuteFromLIDs = Exporter.PermuteToLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumNonzeros();
  int IntSizeOfPacket = 2*SourceMatrix.GlobalMaxNumBlockEntries() + 3;
  int Nsend = SizeOfPacket * Exporter.NumRecv();
  int Nrecv = SizeOfPacket * Exporter.NumSend();
  int IntNsend = IntSizeOfPacket * Exporter.NumRecv();
  int IntNrecv = IntSizeOfPacket * Exporter.NumSend();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs,
		    Nsend, Nrecv, SizeOfPacket, IntNsend, IntNrecv, IntSizeOfPacket,
		    LenImports_, Imports_, LenIntImports_, IntImports_, 
		    LenExports_, Exports_, LenIntExports_, IntExports_,
#ifdef PETRA_MPI
		    Exporter.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_RDP_VBR_Matrix::Export(const Petra_RDP_VBR_Matrix& SourceMatrix, 
				 const Petra_Import & Importer, Petra_CombineMode CombineMode){

  int ierr = 0;

  if (!RowMap().SameAs(Importer.SourceMap())) return(-2);
  if (!SourceMatrix.RowMap().SameAs(Importer.TargetMap())) return(-3);

  int NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumRemoteIDs = Importer.NumExportIDs();
  int NumExportIDs = Importer.NumRemoteIDs();
  int *ExportLIDs = Importer.RemoteLIDs();
  int *RemoteLIDs = Importer.ExportLIDs();
  int *PermuteToLIDs = Importer.PermuteFromLIDs();
  int *PermuteFromLIDs = Importer.PermuteToLIDs();
  int SizeOfPacket = SourceMatrix.GlobalMaxNumNonzeros();
  int IntSizeOfPacket = 2*SourceMatrix.GlobalMaxNumBlockEntries() + 3;
  int Nsend = SizeOfPacket * Importer.NumRecv();
  int Nrecv = SizeOfPacket * Importer.NumSend();
  int IntNsend = IntSizeOfPacket * Importer.NumRecv();
  int IntNrecv = IntSizeOfPacket * Importer.NumSend();

  return(DoTransfer(SourceMatrix, CombineMode, NumSameIDs, NumPermuteIDs, NumRemoteIDs, NumExportIDs,
		    PermuteToLIDs, PermuteFromLIDs, RemoteLIDs, ExportLIDs, 
		    Nsend, Nrecv, SizeOfPacket, IntNsend, IntNrecv, IntSizeOfPacket,
		    LenImports_, Imports_, LenIntImports_, IntImports_, 
		    LenExports_, Exports_, LenIntExports_, IntExports_,
#ifdef PETRA_MPI
		    Importer.GSPlan(), 
#endif
		    true));
}
//=========================================================================
int Petra_RDP_VBR_Matrix::DoTransfer(const Petra_RDP_VBR_Matrix& SourceMatrix, 
				     Petra_CombineMode CombineMode,
				     int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, 
				     int NumExportIDs, 
				     int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
				     int * ExportLIDs,
				     int Nsend, int Nrecv, int SizeOfPacket, 
				     int IntNsend, int IntNrecv, int IntSizeOfPacket,
				     int & LenExports, double * & Exports, int & LenIntExports, int * & IntExports,
				     int & LenImports, double * & Imports, int & LenIntImports, int * & IntImports,
#ifdef PETRA_MPI
				     GSComm_Plan & Plan, 
#endif
				     bool DoReverse){
  int ierr = 0;
  ierr = 
    CopyAndPermute(*this, SourceMatrix, 
		   NumSameIDs, 
		   NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs);

  if (ierr!=0) return(ierr);
  if (CombineMode==Zero) return(0); // All done if CombineMode only involves copying and permuting


  // Need to send global row ID, NumIndices for the row, then the indices.
  if (Nsend>LenExports) {
    if (LenExports>0) delete [] Exports;
    Exports = new double[Nsend];
    LenExports = Nsend;
  }
  // Need to send global row ID, NumIndices for the row, then the indices.
  if (IntNsend>LenIntExports) {
    if (LenIntExports>0) delete [] IntExports;
    IntExports = new int[IntNsend];
    LenIntExports = IntNsend;
  }

  ierr =  Pack(SourceMatrix, SizeOfPacket, IntSizeOfPacket, NumExportIDs, ExportLIDs, Exports, IntExports);

  if (ierr!=0) return(ierr);
  
  if (Nrecv>LenImports) {
    if (LenImports>0) delete Imports;
    Imports = new double[Nrecv];
    LenImports = Nrecv;
  }
  
  if (IntNrecv>LenIntImports) {
    if (LenIntImports>0) delete IntImports;
    IntImports = new int[IntNrecv];
    LenIntImports = IntNrecv;
  }
  
#ifdef PETRA_MPI
  if (RowMap().DistributedGlobal()) {
    int msgtag = 32765;
    bool GSComm_OK;
    GSComm_Comm GSComm;

    if (DoReverse)
      // Do the exchange of remote values data
      GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				    reinterpret_cast<char *> (Exports), 
				    SizeOfPacket * sizeof( double ),
				    reinterpret_cast<char *> (Imports) );
    else
      GSComm_OK = GSComm.Do( Plan, msgtag, 
			     reinterpret_cast<char *> (Exports), 
			     SizeOfPacket * sizeof( double ),
			     reinterpret_cast<char *> (Imports) );

    if (!GSComm_OK) return(-1);

    if (DoReverse)
      // Do the exchange of remote index data
      GSComm_OK = GSComm.DoReverse( Plan, msgtag, 
				    reinterpret_cast<char *> (IntExports), 
				    SizeOfPacket * sizeof( int ),
				    reinterpret_cast<char *> (IntImports) );
    else
      GSComm_OK = GSComm.Do( Plan, msgtag, 
			     reinterpret_cast<char *> (IntExports), 
			     SizeOfPacket * sizeof( int ),
			     reinterpret_cast<char *> (IntImports) );
  
    if (!GSComm_OK) return(-1);

    ierr = 
      UnpackAndCombine( *this, SizeOfPacket, IntSizeOfPacket, NumRemoteIDs, RemoteLIDs,
			Imports, IntImports, CombineMode);
  }
#endif

  return(ierr);
  
}
 
//=========================================================================
int Petra_RDP_VBR_Matrix::CopyAndPermute(Petra_RDP_VBR_Matrix & Target, 
					 const Petra_RDP_VBR_Matrix & Source,
					 int NumSameIDs, 
					 int NumPermuteIDs, int * PermuteToLIDs,
					 int *PermuteFromLIDs){
  
  int i, j;
  
  int BlockRow, NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims, * LDAs;
  double ** Values;
  int FromBlockRow, ToBlockRow;
  
  // Do copy first
  if (NumSameIDs>0) {
      int MaxNumBlockEntries = Source.MaxNumBlockEntries();
      BlockIndices = new int[MaxNumBlockEntries];  // Need some temporary space
      
      
      for (i=0; i<NumSameIDs; i++) {
	BlockRow = Target.GRID(i);
	assert(Source.ExtractGlobalBlockRowPointers(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, 
					        BlockIndices, ColDims, LDAs, Values)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.BeginReplaceGlobalValues(BlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  assert(Target.BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices)==0); 
	// Insert block entries one-at-a-time
	for (j=0; j<NumBlockEntries; j++) Target.SubmitBlockEntry(Values[j], LDAs[j], RowDim, ColDims[j]);
	Target.EndSubmitEntries(); // Complete this block row
      }
      delete [] BlockIndices;
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
      int MaxNumBlockEntries = Source.MaxNumBlockEntries();
      BlockIndices = new int[MaxNumBlockEntries];  // Need some temporary space
      
      for (i=0; i<NumPermuteIDs; i++) {
	FromBlockRow = Source.GRID(PermuteFromLIDs[i]);
	ToBlockRow = Target.GRID(PermuteToLIDs[i]);
	assert(Source.ExtractGlobalBlockRowPointers(FromBlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, 
					        BlockIndices, ColDims, LDAs, Values)==0); // Set pointers
	// Place into target matrix.  Depends on Petra_DataAccess copy/view and static/dynamic graph.
	if (Target.StaticGraph())
	  assert(Target.BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  assert(Target.BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0); 
	// Insert block entries one-at-a-time
	for (j=0; j<NumBlockEntries; j++) Target.SubmitBlockEntry(Values[j], LDAs[j], RowDim, ColDims[j]);
	Target.EndSubmitEntries(); // Complete this block row
      }
      delete [] BlockIndices;
    }
    
  return(0);
}

//=========================================================================
int Petra_RDP_VBR_Matrix::Pack(const Petra_RDP_VBR_Matrix & Source, int SizeOfPacket, int IntSizeOfPacket,
			       int NumSendIDs, int * SendLIDs, double * Sends, 
			       int * IntSends){
  
  int i, j;
  
  int BlockRow, NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * Values;
  int FromBlockRow, ToBlockRow;
  double * valptr;
  int * intptr;
  
  if (NumSendIDs>0) {

    // Each segment of IntSends will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source matrix
    // next int:  RowDim of Block Row
    // next int: NumBlockEntries, Number of indices in row.
    // next NumBlockEntries: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumNonzeros - NumBlockEntries ints) will be wasted but we need fixed
    //   sized segments for current communication routines.

    // Each segment of Sends will be filled with values.

    valptr = Sends;
    intptr = IntSends;
    int GlobalMaxNumBlockEntries = Source.GlobalMaxNumBlockEntries();
    bool NoSumInto = false;
    for (i=0; i<NumSendIDs; i++) {
      FromBlockRow = Source.GRID(SendLIDs[i]);
      BlockIndices = intptr + 3;
      ColDims = BlockIndices + GlobalMaxNumBlockEntries;
      assert(Source.BeginExtractGlobalBlockRowCopy(FromBlockRow, GlobalMaxNumBlockEntries, RowDim,
						   NumBlockEntries, BlockIndices, ColDims)==0);
      // Now extract each block entry into send buffer
      Values = valptr;
      for (j=0; j<NumBlockEntries; j++) {
	int SizeOfValues = RowDim*ColDims[j];
	Source.ExtractEntryCopy(SizeOfValues, Values, RowDim, NoSumInto);
	Values += SizeOfValues;
      }
      // Fill first three slots of intptr with info
      intptr[0] = FromBlockRow;
      intptr[1] = RowDim;
      intptr[2] = NumBlockEntries;
      valptr += SizeOfPacket; // Point to next segment
      intptr += IntSizeOfPacket; // Point to next segment
    }
    // assert(ptr-Sends==Nsend); // Sanity check on send count

  }
    
  return(0);
}

//=========================================================================
int Petra_RDP_VBR_Matrix::UnpackAndCombine(Petra_RDP_VBR_Matrix & Target, int SizeOfPacket,int IntSizeOfPacket,
					   int NumRecvIDs, int * RecvLIDs, 
					   double * Recvs, int * IntRecvs, Petra_CombineMode CombineMode){
  int BlockRow, NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * Values;
  int ToBlockRow;
  int i, j;
  
  double * valptr;
  int * intptr;
  // Unpack it...

  if (NumRecvIDs>0) {

    // Each segment of IntRecvs will be filled by a packed row of information for each row as follows:
    // 1st int: GRID of row where GRID is the global row ID for the source matrix
    // next int:  NumBlockEntries, Number of indices in row.
    // next int:  RowDim of Block Row
    // next NumBlockEntries: The actual indices for the row.
    // Any remaining space (of length GlobalMaxNumNonzeros - NumBlockEntries ints) will be 
    //  wasted but we need fixed sized segments for current communication routines.

    valptr = Recvs;
    intptr = IntRecvs;
    int GlobalMaxNumBlockEntries = Target.GlobalMaxNumBlockEntries();
    
    for (i=0; i<NumRecvIDs; i++) {
      ToBlockRow = Target.GRID(RecvLIDs[i]);
      assert((intptr[0])==ToBlockRow); // Sanity check
      RowDim = Target.RowMap().ElementSize(RecvLIDs[i]);
      assert((intptr[1])==RowDim); // Sanity check
      NumBlockEntries = intptr[2];
      BlockIndices = intptr + 3; 
      ColDims = BlockIndices + GlobalMaxNumBlockEntries;
      if (CombineMode==Add) {
	if (Target.StaticGraph())
	  // Replace any current values
	  assert(Target.BeginSumIntoGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  // Insert values
	  assert(Target.BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
      }
      else if (CombineMode==Insert) {
	if (Target.StaticGraph())
	  // Replace any current values
	  assert(Target.BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  // Insert values
	  assert(Target.BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
      }
      // Now extract each block entry into send buffer
      Values = valptr;
      for (j=0; j<NumBlockEntries; j++) {
	int LDA = RowDim;
	int ColDim = ColDims[j];
	Target.SubmitBlockEntry(Values, LDA, RowDim, ColDim);
	Values += (LDA*ColDim);
      }
      valptr += SizeOfPacket; // Point to next segment
      intptr += IntSizeOfPacket; // Point to next segment
    }
    // assert(ptr-Recvs==Nrecv); // Sanity check on receive count
  }
  
  return(0);
}

// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_VBR_Matrix& A) {
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
	os <<    "Number of Global Blk Entries = "; os << A.NumGlobalBlockEntries(); os << endl;
	os <<    "Global Max Num Block Entries = "; os << A.GlobalMaxNumBlockEntries(); os << endl;
	os <<  "\nNumber of Global Rows        = "; os << A.NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << A.NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << A.NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << A.NumGlobalNonzeros(); os << endl;
	os <<    "Global Maximum Num Entries   = "; os << A.GlobalMaxNumNonzeros(); os << endl;
	if (A.LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (A.UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (A.NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Block Rows  = "; os << A.NumMyBlockRows(); os << endl;
      os <<    "Number of My Block Cols  = "; os << A.NumMyBlockCols(); os << endl;
      os <<    "Number of My Block Diags = "; os << A.NumMyBlockDiagonals(); os << endl;
      os <<    "Number of My Blk Entries = "; os << A.NumMyBlockEntries(); os << endl;
      os <<    "My Max Num Block Entries = "; os << A.MaxNumBlockEntries(); os << endl;
      os <<  "\nNumber of My Rows        = "; os << A.NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << A.NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << A.NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << A.NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << A.MaxNumBlockEntries(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    A.Comm().Barrier();
    A.Comm().Barrier();
    A.Comm().Barrier();
  }

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      long olda = os.setf(ios::right,ios::adjustfield);
      long oldf = os.setf(ios::scientific,ios::floatfield);
      int oldp = os.precision(12);
      int NumMyRows = A.NumMyRows();
      int MaxNumIndices = A.MaxNumBlockEntries();
      int * Indices  = new int[MaxNumIndices];
      double * Values  = new double[MaxNumIndices];
      int NumIndices;
      int i, j;

      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (i=0; i<NumMyRows; i++) {
	int Row = A.GRID(i); // Get global row number
	A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);
	
	for (j = 0; j < NumIndices ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  Indices[j]; os << "    ";
	  os.width(20);
	  os <<  Values[j]; os << "    ";
	  os << endl;
	}
      }

      delete [] Indices;
      delete [] Values;
      
      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    A.Comm().Barrier();
    A.Comm().Barrier();
    A.Comm().Barrier();
  }

  return os;
}
