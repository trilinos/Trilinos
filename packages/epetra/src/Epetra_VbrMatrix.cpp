
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

#include "Epetra_VbrMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "Epetra_SerialDenseMatrix.h"

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int *NumBlockEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, RowMap, NumBlockEntriesPerRow);
  assert(Allocate()==0);
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumBlockEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, RowMap, NumBlockEntriesPerRow);
  assert(Allocate()==0);
}
//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
				   const Epetra_BlockMap& ColMap, int *NumBlockEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, RowMap, ColMap, NumBlockEntriesPerRow);
  assert(Allocate()==0);
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
				   const Epetra_BlockMap& ColMap, int NumBlockEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyBlockRows_(RowMap.NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, RowMap, ColMap, NumBlockEntriesPerRow);
  assert(Allocate()==0);
}
//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph & Graph) 
  : Epetra_DistObject(Graph.RowMap(), "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_((Epetra_CrsGraph*) &Graph),
    Allocated_(false),
    StaticGraph_(true),
    NumMyBlockRows_(Graph.RowMap().NumMyElements()),
    CV_(CV)
{
  InitializeDefaults();
  assert(Allocate()==0);
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(const Epetra_VbrMatrix & Source) 
  : Epetra_DistObject(Source),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(Source.Graph_),
    Allocated_(Source.Allocated_),
    StaticGraph_(Source.StaticGraph_),
    UseTranspose_(Source.UseTranspose_),
    NumMyBlockRows_(Source.NumMyBlockRows_),
    CV_(Copy) {
  InitializeDefaults();
  if (!Source.StaticGraph()) Graph_ = new Epetra_CrsGraph(Source.Graph());
  assert(Allocate()==0);

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
void Epetra_VbrMatrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  UseTranspose_ = false;
  Values_ = 0;
  ColDims_ = 0;
  LDAs_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  ImportVector_ = 0;

  NumBlockEntriesPerRow_  = 0;
  NumAllocatedBlockEntriesPerRow_ = 0;
  Indices_ = 0;
  ElementSizeList_ = 0;
  FirstPointInElementList_ = 0;
  
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

	// Atributes that support the Epetra_RowMatrix and Epetra_Operator interfaces
	RowMatrixRowMap_ = 0;
	RowMatrixColMap_ = 0;
	RowMatrixImporter_ = 0;
	OperatorDomainMap_ = 0;
	OperatorRangeMap_ = 0;
	HavePointObjects_ = false;

	OperatorX_ = 0;
	OperatorY_ = 0;
  return;
}

//==============================================================================
int Epetra_VbrMatrix::Allocate() {

  int i, j;
  
  // Set direct access pointers to graph info (needed for speed)
  NumBlockEntriesPerRow_ = Graph_->NumIndicesPerRow();
  NumAllocatedBlockEntriesPerRow_ = Graph_->NumAllocatedIndicesPerRow();
  Indices_ = Graph_->Indices();

  ElementSizeList_ = RowMap().ElementSizeList();

  FirstPointInElementList_ = RowMap().FirstPointInElementList();
  

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
Epetra_VbrMatrix::~Epetra_VbrMatrix(){

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
      
      delete [] Values_[i];
      delete [] ColDims_[i];
      delete [] LDAs_[i];
    }
  }

  if (All_Values_!=0)   delete [] All_Values_;

  if (Values_!=0)       delete [] Values_;
  if (ColDims_!=0)      delete [] ColDims_;
  if (LDAs_!=0)         delete [] LDAs_;


  if (ImportVector_!=0) delete ImportVector_;


  NumMyBlockRows_ = 0;

  if (LenTemps_>0) {
    delete [] TempRowDims_;
    delete [] TempColDims_;
    delete [] TempLDAs_;
    delete [] TempValues_;
  }

	// Delete any objects related to supporting the RowMatrix and Operator interfaces
	if (HavePointObjects_) {
		if (!RowMatrixColMap().SameAs(RowMatrixRowMap())) delete RowMatrixColMap_;
		if (!OperatorDomainMap().SameAs(RowMatrixRowMap())) delete OperatorDomainMap_;
		if (!OperatorRangeMap().SameAs(RowMatrixRowMap())) delete OperatorRangeMap_;
		delete RowMatrixRowMap_;
		delete RowMatrixImporter_;
		HavePointObjects_ = false;
	}
		
	if (OperatorX_!=0) {
		delete OperatorX_;
		delete OperatorY_;
	}
  
  InitializeDefaults(); // Reset all basic pointers to zero
  Allocated_ = false;

  if (!StaticGraph())   delete Graph_; // We created the graph, so must delete it.
  
}

//==============================================================================
int Epetra_VbrMatrix::PutScalar(double ScalarConstant) 
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
	  *Values++ = ScalarConstant;
      }
    }
  }
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::BeginInsertGlobalValues(int BlockRow, int NumBlockEntries, int * BlockIndices) {

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  Graph_->SetIndicesAreGlobal(true);
  BlockRow = LRID(BlockRow); // Find local row number for this global row index
  
  bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginInsertValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginInsertMyValues(int  BlockRow, int NumBlockEntries, int * BlockIndices) {

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  Graph_->SetIndicesAreLocal(true);
  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginInsertValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);

}

//==========================================================================
int Epetra_VbrMatrix::BeginInsertValues(int BlockRow, int NumBlockEntries, 
					    int * BlockIndices, bool IndicesAreLocal) {

  if (StaticGraph()) EPETRA_CHK_ERR(-2); // If the matrix graph is fully constructed, we cannot insert new values

  int ierr = 0;

  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in BlockRow range    
  if (CV_==View && Values_[BlockRow]!=0) ierr = 2; // This row has be defined already. Issue warning.    
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new

  // Set up pointers, make sure enough temp space for this round of submits

  Epetra_CombineMode SubmitMode = Insert;

  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
  EPETRA_CHK_ERR(ierr);
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::BeginReplaceGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

   BlockRow = LRID(BlockRow); // Normalize row range
   bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginReplaceMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginReplaceValues(int BlockRow, int NumBlockEntries, 
					     int *BlockIndices, bool IndicesAreLocal) {

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.
  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in BlockRow range

  Epetra_CombineMode SubmitMode = Zero; // This is a misuse of Zero mode, fix it later
  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

   BlockRow = LRID(BlockRow); // Normalize row range
   bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoValues(int BlockRow, int NumBlockEntries, 
					     int *BlockIndices, bool IndicesAreLocal) {

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.
  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in BlockRow range

  Epetra_CombineMode SubmitMode = Add;
  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, IndicesAreLocal, SubmitMode));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SetupForSubmits(int BlockRow, int NumBlockEntries, int * BlockIndices, 
					  bool IndicesAreLocal, Epetra_CombineMode SubmitMode) {

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
int Epetra_VbrMatrix::SubmitBlockEntry(double *Values, int LDA, int NumRows, int NumCols) {

  if (CurEntry_==-1) EPETRA_CHK_ERR(-1); // This means that a Begin routine was not called
  if (CurEntry_>=CurNumBlockEntries_) EPETRA_CHK_ERR(-4); // Exceeded the number of entries that can be submitted

  // Fill up temp space with entry

  TempColDims_[CurEntry_] = NumCols;
  TempRowDims_[CurEntry_] = NumRows;
  TempLDAs_[CurEntry_] = LDA;
  TempValues_[CurEntry_] = Values;
  CurEntry_++;
  
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::EndSubmitEntries() {

  if (CurEntry_!=CurNumBlockEntries_) EPETRA_CHK_ERR(-6); // Did not submit right number of entries

  if (CurSubmitMode_==Insert) {
    EPETRA_CHK_ERR(EndInsertValues());
  }
  else {
    EPETRA_CHK_ERR(EndReplaceSumIntoValues());
  }
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::EndReplaceSumIntoValues() {

  int j;
  int ierr = 0;
  int Loc;

  int RowDim = ElementSizeList_[CurBlockRow_];

  bool SumInto = (CurSubmitMode_==Add);

  if (CurIndicesAreLocal_) {
    for (j=0; j<CurNumBlockEntries_; j++) {
      int BlockIndex = CurBlockIndices_[j];
      if (Graph_->FindMyIndexLoc(CurBlockRow_,BlockIndex,j,Loc)) {
	if (Values_[CurBlockRow_][Loc]==0) {
	  LDAs_[CurBlockRow_][Loc] = TempLDAs_[j];
	  ColDims_[CurBlockRow_][Loc] = TempColDims_[j];
	  int LDA = LDAs_[CurBlockRow_][Loc];
	  int NumCols = TempColDims_[j];
	  Values_[CurBlockRow_][Loc] = new double[LDA*NumCols];
	  for (int i=0; i<LDA*NumCols; i++) Values_[CurBlockRow_][Loc][i] = 0.0;
	}
	CopyMat(TempValues_[j], TempLDAs_[j], RowDim, TempColDims_[j], 
		Values_[CurBlockRow_][Loc], LDAs_[CurBlockRow_][Loc], SumInto);
      }
      else {EPETRA_CHK_ERR(-2)}; // Value not found
    }
  }
  else {
    for (j=0; j<CurNumBlockEntries_; j++) {
      int BlockIndex = CurBlockIndices_[j];
      if (Graph_->FindGlobalIndexLoc(CurBlockRow_,BlockIndex,j,Loc)) {
	if (Values_[CurBlockRow_][Loc]==0) {
	  LDAs_[CurBlockRow_][Loc] = TempLDAs_[j];
	  ColDims_[CurBlockRow_][Loc] = TempColDims_[j];
	  int LDA = LDAs_[CurBlockRow_][Loc];
	  int NumCols = TempColDims_[j];
	  Values_[CurBlockRow_][Loc] = new double[LDA*NumCols];
	  for (int i=0; i<LDA*NumCols; i++) Values_[CurBlockRow_][Loc][i] = 0.0;
	}
	CopyMat(TempValues_[j], TempLDAs_[j], RowDim, TempColDims_[j], 
		Values_[CurBlockRow_][Loc], LDAs_[CurBlockRow_][Loc], SumInto);
      }
      else {EPETRA_CHK_ERR(-2)}; // Value not found
    }
  }
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::EndInsertValues() {

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
    
  EPETRA_CHK_ERR(Graph_->InsertIndices(CurBlockRow_, CurNumBlockEntries_, CurBlockIndices_)); // Update graph
  
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::CopyMat(double * A, int LDA, int NumRows, int NumCols, 
					double * B, int LDB, bool SumInto) const {

  int i, j;
  double * ptr1 = B;
  double * ptr2;

  if (LDB<NumRows) EPETRA_CHK_ERR(-1); // Stride of B is not large enough

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
int Epetra_VbrMatrix::TransformToLocal() {
  EPETRA_CHK_ERR(TransformToLocal((Epetra_BlockMap *) (&RowMap()), (Epetra_BlockMap *) (&RowMap())));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::TransformToLocal(const Epetra_BlockMap *DomainMap, const Epetra_BlockMap *RangeMap) {
  
  if (!StaticGraph()) EPETRA_CHK_ERR(Graph_->MakeIndicesLocal(*DomainMap, *RangeMap));
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) EPETRA_CHK_ERR(Graph_->TransformToLocal(DomainMap, RangeMap));

  // NumMyCols_ = Graph_->NumMyCols(); // Redefine based on local number of cols

  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SortEntries() {

  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-1);
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
int Epetra_VbrMatrix::MergeRedundantEntries() {

  int i, j, k;

  if (NoRedundancies()) return(0);
  if (!Sorted()) EPETRA_CHK_ERR(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortEntries was already called.

  bool SumInto = true;
  for (i=0; i<NumMyBlockRows_; i++){
    int NumEntries = NumBlockEntriesPerRow_[i];
    if (NumEntries>1) {
      double ** const Values = Values_[i];
      int * const Indices = Indices_[i];
      int * const LDAs = LDAs_[i];
      int * const ColDims = ColDims_[i];
      int RowDim = ElementSizeList_[i];
      int curEntry =0;
      double* curValues = Values[0];
      int curLDA = LDAs[0];
      int curColDim = ColDims[0];
      for (int k=1; k<NumEntries; k++) {
	if (Indices[k]==Indices[k-1]) {
	  CopyMat(Values[k], LDAs[k], RowDim, ColDims[k],
		  curValues, curLDA, SumInto);
	}
	else {
	  CopyMat(curValues, curLDA, RowDim, curColDim,
		  Values[curEntry], LDAs[curEntry], false);
	  curEntry++;
	  curValues = Values[k];
	  curLDA = LDAs[k];
	  curColDim = ColDims[k];
	}
      }
      CopyMat(curValues, curLDA, RowDim, curColDim,
	      Values[curEntry], LDAs[curEntry], false);
    }
  }
    
  EPETRA_CHK_ERR(Graph_->RemoveRedundantIndices()); // Remove redundant indices and then return
  return(0);

}

//==========================================================================
int Epetra_VbrMatrix::OptimizeStorage() {

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


  if ((CV_==View) && !Contiguous) EPETRA_CHK_ERR(-1);  // This is user data, it's not contiguous and we can't make it so.

  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  if (ierr) EPETRA_CHK_ERR(ierr);

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
int Epetra_VbrMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, 
                   int & NumEntries, double *Values, int * Indices) const {

  cout << "Must implement..." << endl;
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::ExtractGlobalBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
							    int & RowDim, int & NumBlockEntries, 
							    int * BlockIndices, int * & ColDims, 
							    int * & LDAs, double ** & Values) const {

  bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(ExtractBlockRowPointers(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
				     ColDims, LDAs, Values, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractMyBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
							int & RowDim, int & NumBlockEntries, 
							int * BlockIndices, int * & ColDims, 
							int * & LDAs, double ** & Values) const {

  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(ExtractBlockRowPointers(BlockRow,MaxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
				     ColDims, LDAs, Values, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
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

  if (ierr) EPETRA_CHK_ERR(ierr);

  RowDim = ElementSizeList_[BlockRow];

  ColDims = ColDims_[BlockRow];
  LDAs = LDAs_[BlockRow];
  Values = Values_[BlockRow];


  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractGlobalBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
							 int & RowDim, int & NumBlockEntries, 
							 int * BlockIndices, int * ColDims) const {

  bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowCopy(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractMyBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
						     int & RowDim, int & NumBlockEntries, 
						     int * BlockIndices, int * ColDims) const {

  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowCopy(BlockRow,MaxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
						   int & RowDim, int & NumBlockEntries, 
						   int * BlockIndices, int * ColDims, 
						   bool IndicesAreLocal) const  {
  int ierr = 0;
  if (!IndicesAreLocal)
    ierr = Graph_->ExtractGlobalRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowCopy(BlockRow, MaxNumBlockEntries, NumBlockEntries, BlockIndices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  bool ExtractView = false;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, IndicesAreLocal);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractBlockDimsCopy(NumBlockEntries, ColDims));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SetupForExtracts(int BlockRow, int & RowDim, int NumBlockEntries, bool ExtractView, 
						   bool IndicesAreLocal) const {

  if (!IndicesAreLocal) BlockRow = LRID(BlockRow); // Normalize row range
  CurExtractBlockRow_ = BlockRow;
  CurExtractEntry_ = 0;
  CurExtractNumBlockEntries_ = NumBlockEntries;
  CurExtractIndicesAreLocal_ = IndicesAreLocal;
  CurExtractView_ = ExtractView;
  CurRowDim_ = ElementSizeList_[CurExtractBlockRow_];
  RowDim = CurRowDim_;
    
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractBlockDimsCopy(int NumBlockEntries, int * ColDims) const {

  int * CurColDims = ColDims_[CurExtractBlockRow_];

  for (int i=0; i<NumBlockEntries; i++) ColDims[i] = CurColDims[i];
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const
{
  if (CurExtractEntry_==-1) EPETRA_CHK_ERR(-1); // No BeginCopy routine was called
  int CurColDim = ColDims_[CurExtractBlockRow_][CurExtractEntry_];
  if (LDA*CurColDim>SizeOfValues) EPETRA_CHK_ERR(-2);  // Not enough space

  double * CurValues = Values_[CurExtractBlockRow_][CurExtractEntry_];
  int CurLDA = LDAs_[CurExtractBlockRow_][CurExtractEntry_];

  CurExtractEntry_++; // Increment Entry Pointer

  EPETRA_CHK_ERR(CopyMat(CurValues, CurLDA, CurRowDim_, CurColDim, Values, LDA, SumInto));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
							 int * & BlockIndices, int * & ColDims, int * & LDAs) const
{

  bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						     int * & BlockIndices, int * & ColDims, int * & LDAs)  const
{

  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						   int * & BlockIndices, int * & ColDims, 
						   int * & LDAs, bool IndicesAreLocal) const
{
  int ierr = 0;
  if (!IndicesAreLocal)
    ierr = Graph_->ExtractGlobalRowView(BlockRow, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowView(BlockRow,  NumBlockEntries, BlockIndices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  bool ExtractView = true;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, IndicesAreLocal);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractBlockDimsView(NumBlockEntries, ColDims, LDAs));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractBlockDimsView(int NumBlockEntries, int * & ColDims, int * & LDAs) const {

  ColDims = ColDims_[CurExtractBlockRow_];
  LDAs = LDAs_[CurExtractBlockRow_];
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractEntryView(double * & Values) const
{
  Values = Values_[CurExtractBlockRow_][CurExtractEntry_];

  CurExtractEntry_++; // Increment Entry Pointer
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						    int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const
{

  Values = Values_[LRID(BlockRow)]; // Pointer to Array of pointers for this row's block entries
  bool IndicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
						int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const
{

  Values = Values_[BlockRow]; // Pointer to Array of pointers for this row's block entries
  bool IndicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
				  ColDims, LDAs, IndicesAreLocal));
  return(0);
}

//==============================================================================
int Epetra_VbrMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
	
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) EPETRA_CHK_ERR(-2); // Maps must be the same
  double * diagptr = Diagonal.Values();
  for(int i=0; i<NumMyBlockRows_; i++){
    int BlockRow = i;
    int RowDim = ElementSizeList_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
    for (int j=0; j<NumEntries; j++) {
      int BlockCol = Indices[j];
      if (BlockRow==BlockCol) {
	CopyMatDiag(Values_[i][j], LDAs_[i][j], RowDim, ColDims_[i][j], 
		    diagptr+FirstPointInElementList_[i]);
	break;
      }
    }
  }
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::BeginExtractBlockDiagonalCopy(int MaxNumBlockDiagonalEntries, 
							int & NumBlockDiagonalEntries, int * RowColDims ) const{
	
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  CurBlockDiag_ = 0; // Initialize pointer
  NumBlockDiagonalEntries = NumMyBlockRows_;
  if (NumBlockDiagonalEntries>MaxNumBlockDiagonalEntries) EPETRA_CHK_ERR(-2);
  EPETRA_CHK_ERR(RowMap().ElementSizeList(RowColDims));
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::ExtractBlockDiagonalEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto ) const {
	
  if (CurBlockDiag_==-1) EPETRA_CHK_ERR(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow = i;
  int RowDim = ElementSizeList_[i];
  int NumEntries = NumBlockEntriesPerRow_[i];
  int * Indices = Indices_[i];
  for (int j=0; j<NumEntries; j++) {
    int Col = Indices[j];
    if (BlockRow==Col) {
      int ColDim = ColDims_[i][j];
      if (LDA*ColDim>SizeOfValues) EPETRA_CHK_ERR(-2); // Not enough room in Values
      CopyMat(Values_[i][j], LDAs_[i][j], RowDim, ColDim, Values, LDA, SumInto);
      break;
    }
  }
  CurBlockDiag_++; // Increment counter
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::BeginExtractBlockDiagonalView(int & NumBlockDiagonalEntries, 
							int * & RowColDims ) const {
	
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  CurBlockDiag_ = 0; // Initialize pointer
  NumBlockDiagonalEntries = NumMyBlockRows_;
  RowColDims = ElementSizeList_;
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::ExtractBlockDiagonalEntryView(double * & Values, int & LDA) const {
	
  if (CurBlockDiag_==-1) EPETRA_CHK_ERR(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow = i;
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
int Epetra_VbrMatrix::CopyMatDiag(double * A, int LDA, int NumRows, int NumCols, 
					double * Diagonal) const {

  int i;
  double * ptr1 = Diagonal;
  double * ptr2;
  int ndiags = EPETRA_MIN(NumRows,NumCols);

  for (i=0; i<ndiags; i++) {
    ptr2 = A + i*LDA+i;
    *ptr1++ = *ptr2;
  }
  return(0);
}
//=============================================================================
int Epetra_VbrMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {

  int BlockRow, BlockOffset;
  int ierr = RowMap().FindLocalElementID(MyRow, BlockRow, BlockOffset);  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int NumBlockEntries = NumMyBlockEntries(BlockRow);
  NumEntries = 0;
  for (int i=0; i<NumBlockEntries; i++) NumEntries += ColDims_[BlockRow][i];
  return(0);  
}
//=============================================================================
int Epetra_VbrMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't row unless matrix is filled
  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-2);

  int ierr = 0;
  int BlockRow, BlockOffset;
  ierr = RowMap().FindLocalElementID(MyRow, BlockRow, BlockOffset);  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int RowDim, NumBlockEntries;
  int * BlockIndices, * ColDims, * LDAs;
  double ** ValBlocks;
  ierr = ExtractMyBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, ColDims, LDAs, ValBlocks);
  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int * ColFirstPointInElementList = FirstPointInElementList_;
  if (Importer()!=0) ColFirstPointInElementList = ColMap().FirstPointInElementList();
  NumEntries = 0;
  for (int i=0; i<NumBlockEntries; i++) {
    int ColDim = ColDims[i];
    NumEntries += ColDim;
    if (NumEntries>Length) EPETRA_CHK_ERR(-3); // Not enough space
    double * A = ValBlocks[i] + BlockOffset; // Point to first element in row
    int LDA = LDAs[i];
    int Index = ColFirstPointInElementList[BlockIndices[i]];
    for (int j=0; j < ColDim; j++) {
      *Values++ = *A;
      A += LDA;
      *Indices++ = Index++;
    }
  }
  return(0);
      
}
//=============================================================================
int Epetra_VbrMatrix::Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {
//
// This function forms the product y = A * x or y = A' * x
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  
  int i, j;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * FirstPointInElement = FirstPointInElementList_;
  int * ElementSize = ElementSizeList_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;

  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;



  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*)ImportVector_->Values();
      ColElementSizeList = ColMap().ElementSizeList(); // The Import map will always have an existing ElementSizeList
      ColFirstPointInElementList = ColMap().FirstPointInElementList(); // Import map will always have an existing ...
    }
    

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
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
      double * cury = yp + *FirstPointInElement++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
				//sum += BlockRowValues[j] * xp[BlockRowIndices[j]];
				double * A = BlockRowValues[j];
				int LDA = BlockRowLDAs[j];
				int Index = BlockRowIndices[j];
				double * curx = xp + ColFirstPointInElementList[Index];
				int ColDim = ColElementSizeList[Index];
				GEMV('N', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);
				
      }
      if (Exporter()!=0) {
				EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
			}
		}
  }
  
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*)ExportVector_->Values();
    }
    
    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      yp = (double*)ImportVector_->Values();
      ColElementSizeList = ColMap().ElementSizeList(); // The Import map will always have an existing ElementSizeList
      ColFirstPointInElementList = ColMap().FirstPointInElementList(); // Import map will always have an existing ...
    }
    
    // Do actual computation
    int NumMyCols_ = NumMyCols();
    for (i=0; i < NumMyCols_; i++) yp[i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      double * curx = xp + *FirstPointInElement++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
	//yp[BlockRowIndices[j]] += BlockRowValues[j] * xp[i];
	double * A = BlockRowValues[j];
	int LDA = BlockRowLDAs[j];
	int Index = BlockRowIndices[j];
	double * cury = yp + ColFirstPointInElementList[Index];
	int ColDim = ColElementSizeList[Index];
	GEMV('T', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);
	
      }
    }
    if (Importer()!=0) {
			EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
		}
	}
  
  UpdateFlops(2*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  
  int i;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstPointInElementList = FirstPointInElementList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;

   
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
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors);

      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
      ColElementSizeList = ColMap().ElementSizeList();
      ColFirstPointInElementList = ColMap().FirstPointInElementList();
    }
    
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
       // Create Export vector if needed
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors);

      ExportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ExportVector_->Pointers();
      RowElementSizeList = ColMap().ElementSizeList();
      RowFirstPointInElementList = ColMap().FirstPointInElementList();
    }
    else
      Y.PutScalar(0.0); // Zero y values
    
    // Do actual computation
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      int  yoff = *RowFirstPointInElementList++;
      int RowDim = *RowElementSizeList++;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
		       ColFirstPointInElementList, ColElementSizeList, 
		       1.0, BlockRowValues, BlockRowLDAs, Xp, 1.0, Yp, NumVectors);
    }
    if (Exporter()!=0) {
			EPETRA_CHK_ERR(Y.Export(*ExportVector_, *Exporter(), Add)); // Fill Y with Values from export vector
		}
	}
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()<NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
       // Create Export vector if needed
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors);

      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
      ColElementSizeList = RowMap().ElementSizeList();
      ColFirstPointInElementList = RowMap().FirstPointInElementList();
    }
  
    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()<NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      // Create import vector if needed
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors);

      ImportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ImportVector_->Pointers();
      RowElementSizeList = ColMap().ElementSizeList();
      RowFirstPointInElementList = ColMap().FirstPointInElementList();
    }
    else
      Y.PutScalar(0.0); // Zero y values

    // Do actual computation
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      double ** BlockRowValues  = *Values++;
      int *    BlockRowLDAs = *LDAs++;
      int  xoff = *ColFirstPointInElementList++;
      int RowDim = *ColElementSizeList++;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, xoff, 
		       RowFirstPointInElementList, RowElementSizeList, 
		       1.0, BlockRowValues, BlockRowLDAs, Xp, 1.0, Yp, NumVectors);
    }

    if (Importer()!=0) {
			EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
		}
	}

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowMultiply(bool TransA, int RowDim, int NumEntries, 
		      int * BlockIndices, int RowOff,
		      int * FirstPointInElementList, int * ElementSizeList,
		      double Alpha, double ** As, int * LDAs, 
		      double ** X, double Beta, double ** Y, int NumVectors) const {
  int j, k;

  if (!TransA) {
    for (j=0; j < NumEntries; j++) {
      double * A = As[j];
      int LDA = LDAs[j];
      int BlockIndex = BlockIndices[j];
      int xoff = FirstPointInElementList[BlockIndex];
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
      int yoff = FirstPointInElementList[BlockIndex];
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
int Epetra_VbrMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) EPETRA_CHK_ERR(-3);
  if (!NoDiagonal()) EPETRA_CHK_ERR(-4); // We must use UnitDiagonal

  int i;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * FirstPointInElement = FirstPointInElementList_;
  int * ElementSize = ElementSizeList_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;

  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;

  // If upper, point to last row
  if (Upper) {
    NumBlockEntriesPerRow += NumMyBlockRows_-1;
    FirstPointInElement += NumMyBlockRows_-1;
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
	int  yoff = *FirstPointInElement--;
	int RowDim = *ElementSize--;
	BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
			 ColFirstPointInElementList, ColElementSizeList, 
			 1.0, BlockRowValues, BlockRowLDAs, Yp, -1.0, Yp, NumVectors);
      }
    }
    else {
      for (i=0; i < NumMyBlockRows_; i++) {
	int      NumEntries = *NumBlockEntriesPerRow++;
	int *    BlockRowIndices = *Indices++;
	double ** BlockRowValues  = *Values++;
	int *    BlockRowLDAs = *LDAs++;
	int  yoff = *FirstPointInElement++;
	int RowDim = *ElementSize++;
	BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
			 ColFirstPointInElementList, ColElementSizeList, 
			 1.0, BlockRowValues, BlockRowLDAs, Yp, -1.0, Yp, NumVectors);
      }
    }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::InvRowSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//
  EPETRA_CHK_ERR(InverseSums(true, x));
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::InvColSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//
  EPETRA_CHK_ERR(InverseSums(false, x));
  return(0);
}
//=============================================================================
int Epetra_VbrMatrix::InverseSums(bool DoRows, Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (DoRows) {
    if (!Graph().RangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A
  }
  else {
    if (!Graph().DomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstPointInElementList = FirstPointInElementList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;
  if (Importer()!=0) {
    ColElementSizeList = ColMap().ElementSizeList();
    ColFirstPointInElementList = ColMap().FirstPointInElementList();
  }

  x.PutScalar(0.0); // Zero out result vector

  double * xp = (double*)x.Values();

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  Epetra_Vector * x_tmp = 0;
  if (!DoRows) {
    if (Importer()!=0) {
      x_tmp = new Epetra_Vector(ColMap()); // Create import vector if needed
      xp = (double*)x_tmp->Values();
    }
  }

  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    double ** BlockRowValues  = *Values++;
    int *    BlockRowLDAs = *LDAs++;
    int xoff = *RowFirstPointInElementList++;
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
	int off = ColFirstPointInElementList[BlockIndex];
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
      EPETRA_CHK_ERR(x.Export(*x_tmp, *Importer(), Add)); // Fill x with Values from import vector
      delete x_tmp;
      xp = (double*) x.Values();
    }
  }
  int NumMyRows_ = NumMyRows();
  {for (int i=0; i < NumMyRows_; i++) {
    double scale = xp[i];
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero row/col sum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }}
  UpdateFlops(NumGlobalNonzeros());

  EPETRA_CHK_ERR(ierr);
  return(0);
}
//=============================================================================
int Epetra_VbrMatrix::LeftScale(const Epetra_Vector& x) {
//
// Multiply the ith row of A by x[i].
//
  EPETRA_CHK_ERR(Scale(true, x));
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::RightScale(const Epetra_Vector& x) {
//
// Multiply the jth column of A by x[j].
//
  EPETRA_CHK_ERR(Scale (false, x));
  return(0);
}
//=============================================================================
int Epetra_VbrMatrix::Scale(bool DoRows, const Epetra_Vector& x) {

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (DoRows) {
    if (!Graph().RangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A
  }
  else {
    if (!Graph().DomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** LDAs = LDAs_;
  int ** Indices = Indices_;
  double *** Values = Values_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstPointInElementList = FirstPointInElementList_;
  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;
  if (Importer()!=0) {
    ColElementSizeList = ColMap().ElementSizeList();
    ColFirstPointInElementList = ColMap().FirstPointInElementList();
  }

  double * xp = (double*)x.Values();

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  Epetra_Vector * x_tmp = 0;
  if (!DoRows) {
    if (Importer()!=0) {
      x_tmp = new Epetra_Vector(ColMap()); // Create import vector if needed
      EPETRA_CHK_ERR(x_tmp->Import(x,*Importer(), Insert)); // x_tmp will have all the values we need
      xp = (double*)x_tmp->Values();
    }
  }

  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    double ** BlockRowValues  = *Values++;
    int *    BlockRowLDAs = *LDAs++;
    int xoff = *RowFirstPointInElementList++;
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
	int off = ColFirstPointInElementList[BlockIndex];
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

  EPETRA_CHK_ERR(ierr);
  return(0);
}
//=============================================================================
double Epetra_VbrMatrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

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
    for (int j=0; j < RowDim; j++) Local_NormInf = EPETRA_MAX(Local_NormInf, tempv[j]);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  delete [] tempv;
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowNormInf(int RowDim, int NumEntries, 
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
double Epetra_VbrMatrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);
  int * ColFirstPointInElementList = FirstPointInElementList_;
  if (Importer()!=0) {
    ColFirstPointInElementList = ColMap().FirstPointInElementList();
  }

  Epetra_Vector * x = new Epetra_Vector(RowMap()); // Need temp vector for column sums
  
  double * xp = (double*)x->Values();
  Epetra_MultiVector * x_tmp = 0;

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create temporary import vector if needed
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
		    BlockRowLDAs, BlockRowValues,  ColFirstPointInElementList, xp);
  }
  if (Importer()!=0) {
		EPETRA_CHK_ERR(x->Export(*x_tmp, *Importer(), Add));
	} // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowNormOne(int RowDim, int NumEntries, int * BlockRowIndices,
					   int * ColDims, int * LDAs, double ** As, 
					   int * ColFirstPointInElementList, double * x) const {
  int i, j, k;

  for (i=0; i < NumEntries; i++) {
    double * A = As[i];
    int LDA = LDAs[i];
    int ColDim = ColDims[i];
    double * curx = x + ColFirstPointInElementList[BlockRowIndices[i]];
    for (j=0; j<ColDim; j++) {
      for (k=0; k<RowDim; k++) curx[j] += fabs(A[k]);
      A += LDA;
    }
  }
  return;
}
//=========================================================================
int Epetra_VbrMatrix::CheckSizes(const Epetra_DistObject & Source) {
  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  if (!A.Graph().GlobalConstantsComputed()) EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::CopyAndPermute(const Epetra_DistObject & Source,
				     int NumSameIDs, 
				     int NumPermuteIDs, int * PermuteToLIDs,
				     int *PermuteFromLIDs){
  
  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  int i, j;
  
  int BlockRow, NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims, * LDAs;
  double ** Values;
  int FromBlockRow, ToBlockRow;
  
  // Do copy first
  if (NumSameIDs>0) {
      int MaxNumBlockEntries = A.MaxNumBlockEntries();
      BlockIndices = new int[MaxNumBlockEntries];  // Need some temporary space
      
      
      for (i=0; i<NumSameIDs; i++) {
	BlockRow = GRID(i);
	assert(A.ExtractGlobalBlockRowPointers(BlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, 
					        BlockIndices, ColDims, LDAs, Values)==0); // Set pointers
	// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
	if (StaticGraph() || IndicesAreLocal())
	  assert(BeginReplaceGlobalValues(BlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  assert(BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices)==0); 
	// Insert block entries one-at-a-time
	for (j=0; j<NumBlockEntries; j++) SubmitBlockEntry(Values[j], LDAs[j], RowDim, ColDims[j]);
	EndSubmitEntries(); // Complete this block row
      }
      delete [] BlockIndices;
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
      int MaxNumBlockEntries = A.MaxNumBlockEntries();
      BlockIndices = new int[MaxNumBlockEntries];  // Need some temporary space
      
      for (i=0; i<NumPermuteIDs; i++) {
	FromBlockRow = A.GRID(PermuteFromLIDs[i]);
	ToBlockRow = GRID(PermuteToLIDs[i]);
	assert(A.ExtractGlobalBlockRowPointers(FromBlockRow, MaxNumBlockEntries, RowDim, NumBlockEntries, 
					        BlockIndices, ColDims, LDAs, Values)==0); // Set pointers
	// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
	if (StaticGraph() || IndicesAreLocal())
	  assert(BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
	else
	  assert(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0); 
	// Insert block entries one-at-a-time
	for (j=0; j<NumBlockEntries; j++) SubmitBlockEntry(Values[j], LDAs[j], RowDim, ColDims[j]);
	EndSubmitEntries(); // Complete this block row
      }
      delete [] BlockIndices;
    }
    
  return(0);
}

//=========================================================================
int Epetra_VbrMatrix::PackAndPrepare(const Epetra_DistObject & Source,int NumExportIDs, int * ExportLIDs,
				      int Nsend, int Nrecv,
				      int & LenExports, char * & Exports, int & LenImports, 
				      char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor){
  

  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);

  double * DoubleExports = 0;
  double * DoubleImports = 0;
  int GlobalMaxNumNonzeros = A.GlobalMaxNumNonzeros();
  int GlobalMaxNumBlockEntries = A.GlobalMaxNumBlockEntries();
  // Will have GlobalMaxNumEntries doubles, GlobalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = GlobalMaxNumNonzeros +  
                        (((2*GlobalMaxNumBlockEntries+3)+sizeof(int)-1)*sizeof(int))/sizeof(double);
  SizeOfPacket = DoublePacketSize * sizeof(double); 




  if (DoublePacketSize*Nsend>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = DoublePacketSize*Nsend;
    DoubleExports = new double[LenExports];
    Exports = (char *) DoubleExports;
  }

  if (DoublePacketSize*Nrecv>LenImports) {
    if (LenImports>0) delete [] Imports;
    LenImports = DoublePacketSize*Nrecv;
    DoubleImports = new double[LenImports];
    Imports = (char *) DoubleImports;
  }

  if (NumExportIDs<=0) return(0); // All done if nothing to pack

  int i, j;
  
  int NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * Values;
  int FromBlockRow;
  double * valptr, * dintptr;
  int * intptr;
  
  // Each segment of IntExports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  RowDim of Block Row
  // next int: NumBlockEntries, Number of indices in row.
  // next NumBlockEntries: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumNonzeros - NumBlockEntries ints) will be wasted but we need fixed
  //   sized segments for current communication routines.

  // Each segment of Exports will be filled with values.

  valptr = (double *) Exports;
  dintptr = valptr + GlobalMaxNumNonzeros;
  intptr = (int *) dintptr;
  bool NoSumInto = false;
  for (i=0; i<NumExportIDs; i++) {
    FromBlockRow = A.GRID(ExportLIDs[i]);
    BlockIndices = intptr + 3;
    ColDims = BlockIndices + GlobalMaxNumBlockEntries;
    assert(A.BeginExtractGlobalBlockRowCopy(FromBlockRow, GlobalMaxNumBlockEntries, RowDim,
						 NumBlockEntries, BlockIndices, ColDims)==0);
    // Now extract each block entry into send buffer
    Values = valptr;
    for (j=0; j<NumBlockEntries; j++) {
      int SizeOfValues = RowDim*ColDims[j];
      A.ExtractEntryCopy(SizeOfValues, Values, RowDim, NoSumInto);
      Values += SizeOfValues;
    }
    // Fill first three slots of intptr with info
    intptr[0] = FromBlockRow;
    intptr[1] = RowDim;
    intptr[2] = NumBlockEntries;
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + GlobalMaxNumNonzeros;
    intptr = (int *) dintptr;
  }
    
  return(0);
}

//=========================================================================
int Epetra_VbrMatrix::UnpackAndCombine(const Epetra_DistObject & Source, 
				       int NumImportIDs, int * ImportLIDs, 
				       char * Imports, int & SizeOfPacket, 
				       Epetra_Distributor & Distor, 
				       Epetra_CombineMode CombineMode){
  if (NumImportIDs<=0) return(0);

  if(   CombineMode != Add
     && CombineMode != Zero
     && CombineMode != Insert )
    EPETRA_CHK_ERR(-1); // CombineMode not supported, default to mode Zero


  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  int NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * Values;
  int ToBlockRow;
  int i, j;
  
  double * valptr, *dintptr;
  int * intptr;
  int GlobalMaxNumNonzeros = A.GlobalMaxNumNonzeros();
  int GlobalMaxNumBlockEntries = A.GlobalMaxNumBlockEntries();
  // Will have GlobalMaxNumEntries doubles, GlobalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = GlobalMaxNumNonzeros +  
                        (((2*GlobalMaxNumBlockEntries+3)+sizeof(int)-1)*sizeof(int))/sizeof(double);
  // Unpack it...


  // Each segment of IntImports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumBlockEntries, Number of indices in row.
  // next int:  RowDim of Block Row
  // next NumBlockEntries: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumNonzeros - NumBlockEntries ints) will be 
  //  wasted but we need fixed sized segments for current communication routines.

  valptr = (double *) Imports;
  dintptr = valptr + GlobalMaxNumNonzeros;
  intptr = (int *) dintptr;
    
  for (i=0; i<NumImportIDs; i++) {
    ToBlockRow = GRID(ImportLIDs[i]);
    assert((intptr[0])==ToBlockRow); // Sanity check
    RowDim = RowMap().ElementSize(ImportLIDs[i]);
    assert((intptr[1])==RowDim); // Sanity check
    NumBlockEntries = intptr[2];
    BlockIndices = intptr + 3; 
    ColDims = BlockIndices + GlobalMaxNumBlockEntries;
    if (CombineMode==Add) {
      if (StaticGraph() || IndicesAreLocal())
	// Replace any current values
	assert(BeginSumIntoGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
      else
	// Insert values
	assert(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
    }
    else if (CombineMode==Insert) {
      if (StaticGraph() || IndicesAreLocal())
	// Replace any current values
	assert(BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
      else
	// Insert values
	assert(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices)==0);
    }
    // Now extract each block entry into send buffer
    Values = valptr;
    for (j=0; j<NumBlockEntries; j++) {
      int LDA = RowDim;
      int ColDim = ColDims[j];
      SubmitBlockEntry(Values, LDA, RowDim, ColDim);
      Values += (LDA*ColDim);
    }
    EndSubmitEntries(); // Done with this block row
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + GlobalMaxNumNonzeros;
    intptr = (int *) dintptr;
  }
  
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::GeneratePointObjects() const {

	if (HavePointObjects_)  return(0); // Already done

	// Generate a point-wise compatible row map
	EPETRA_CHK_ERR(BlockMap2PointMap(RowMap(), RowMatrixRowMap_));

	// For each of the next maps, first check if the corresponding block map
	// is the same as the block row map for this matrix.  If so, then we simply
	// copy the pointer.  Otherwise we create a new point map.

	// This check can save storage and time since it will often be the case that the
	// domain, range and row block maps will be the same.  Also, in the serial case,
	// the column block map will also often be the same as the row block map.

	if (ColMap().SameAs(RowMap())) {
		RowMatrixColMap_ = RowMatrixRowMap_;
	}
	else {
		EPETRA_CHK_ERR(BlockMap2PointMap(ColMap(), RowMatrixColMap_));
	}

	if (DomainMap().SameAs(RowMap())) {
		OperatorDomainMap_ = RowMatrixRowMap_;
	}
	else {
		EPETRA_CHK_ERR(BlockMap2PointMap(DomainMap(), OperatorDomainMap_));
	}
	if (RangeMap().SameAs(RowMap())) {
		OperatorRangeMap_ = RowMatrixRowMap_;
	}
	else {
		EPETRA_CHK_ERR(BlockMap2PointMap(RangeMap(), OperatorRangeMap_));
	}

	// Finally generate Importer that will migrate needed domain elements to the column map
	// layout.
	RowMatrixImporter_ = new Epetra_Import(*RowMatrixColMap_, *OperatorDomainMap_);

	HavePointObjects_ = true;
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::BlockMap2PointMap(const Epetra_BlockMap & BlockMap, Epetra_Map * & PointMap) const {
	// Generate an Epetra_Map that has the same number and distribution of points
	// as the input Epetra_BlockMap object.  The global IDs for the output PointMap
	// are computed by using the MaxElementSize of the BlockMap.  For variable block
	// sizes this will create gaps in the GID space, but that is OK for Epetra_Maps.

	int MaxElementSize = BlockMap.MaxElementSize();
	int PtNumMyElements = BlockMap.NumMyPoints();
	int * PtMyGlobalElements = 0;
	if (PtNumMyElements>0) PtMyGlobalElements = new int[PtNumMyElements];

	int NumMyElements = BlockMap.NumMyElements();

	int curID = 0;
	for (int i=0; i<NumMyElements; i++) {
		int StartID = BlockMap.GID(i)*MaxElementSize;
		int ElementSize = BlockMap.ElementSize(i);
		for (int j=0; j<ElementSize; j++) PtMyGlobalElements[curID++] = StartID+j;
	}
	assert(curID==PtNumMyElements); // Sanity test

	PointMap = new Epetra_Map(-1, PtNumMyElements, PtMyGlobalElements, BlockMap.IndexBase(), BlockMap.Comm());

	if (PtNumMyElements>0) delete [] PtMyGlobalElements;

	if (!BlockMap.PointSameAs(*PointMap)) {EPETRA_CHK_ERR(-1);} // Maps not compatible
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::UpdateOperatorXY(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

	if (OperatorX_!=0)
		if (OperatorX_->NumVectors()!=X.NumVectors()) {delete OperatorX_; OperatorX_ = 0; delete OperatorY_; OperatorY_=0;}
	if (OperatorX_==0) {
		if (!X.Map().PointSameAs(DomainMap())) EPETRA_CHK_ERR(-1); // X not point-wise compatible with the block domain map
		if (!Y.Map().PointSameAs(RangeMap())) EPETRA_CHK_ERR(-2); // Y not point-wise compatible with the block col map
		OperatorX_ = new Epetra_MultiVector(View, DomainMap(), X.Pointers(), X.NumVectors());
		OperatorY_ = new Epetra_MultiVector(View, RangeMap(), Y.Pointers(), Y.NumVectors());
	}
	else {
		EPETRA_CHK_ERR(OperatorX_->ResetView(X.Pointers()));
		EPETRA_CHK_ERR(OperatorY_->ResetView(Y.Pointers()));
	}
	return(0);
}
//=========================================================================
int Epetra_VbrMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

		EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Epetra_VbrMatrix::Multiply(Epetra_VbrMatrix::UseTranspose(), *OperatorX_, *OperatorY_));
		return(0);
}
//=========================================================================
int Epetra_VbrMatrix::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

		EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Solve(UpperTriangular(), Epetra_VbrMatrix::UseTranspose(), NoDiagonal(), *OperatorX_, *OperatorY_));
		return(0);
}
//=========================================================================
void Epetra_VbrMatrix::Print(ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (MyPID==0) {
	os <<  "\nNumber of Global Block Rows  = "; os << NumGlobalBlockRows(); os << endl;
	os <<    "Number of Global Block Cols  = "; os << NumGlobalBlockCols(); os << endl;
	os <<    "Number of Global Block Diags = "; os << NumGlobalBlockDiagonals(); os << endl;
	os <<    "Number of Global Blk Entries = "; os << NumGlobalBlockEntries(); os << endl;
	os <<    "Global Max Num Block Entries = "; os << GlobalMaxNumBlockEntries(); os << endl;
	os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros(); os << endl;
	os <<    "Global Maximum Num Entries   = "; os << GlobalMaxNumNonzeros(); os << endl;
	if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Block Rows  = "; os << NumMyBlockRows(); os << endl;
      os <<    "Number of My Block Cols  = "; os << NumMyBlockCols(); os << endl;
      os <<    "Number of My Block Diags = "; os << NumMyBlockDiagonals(); os << endl;
      os <<    "Number of My Blk Entries = "; os << NumMyBlockEntries(); os << endl;
      os <<    "My Max Num Block Entries = "; os << MaxNumBlockEntries(); os << endl;
      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << MaxNumBlockEntries(); os << endl; os << endl;

      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }

  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumBlockRows1 = NumMyBlockRows();
      int MaxNumBlockEntries1 = MaxNumBlockEntries();
      int * BlockIndices1  = new int[MaxNumBlockEntries1];
      double ** Values1;
      int RowDim1, NumBlockEntries1, * ColDims1, * LDAs1;
      int i, j;

      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Block Row Index ";
	os.width(10);
	os <<  "   Block Col Index \n";
	os.width(20);
	os <<  "   Values     ";
	os << endl;
      }
      for (i=0; i<NumBlockRows1; i++) {
	int BlockRow1 = GRID(i); // Get global row number
	ExtractGlobalBlockRowPointers(BlockRow1, MaxNumBlockEntries1, RowDim1, 
				      NumBlockEntries1, BlockIndices1, ColDims1,
				      LDAs1, Values1);
	
	for (j = 0; j < NumBlockEntries1 ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  BlockRow1 ; os << "    ";	
	  os.width(10);
	  os <<  BlockIndices1[j]; os << "    " << endl;
	  os.width(20);
	  Epetra_SerialDenseMatrix entry(View, Values1[j], LDAs1[j], RowDim1, ColDims1[j]);
	  os << entry; os << "    ";
	  os << endl;
	}
      }

      delete [] BlockIndices1;
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }}

  return;
}
