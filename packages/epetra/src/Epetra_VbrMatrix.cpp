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

#include "Epetra_ConfigDefs.h"
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES // FIXME
// FIXME long long : whole file

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& rowMap, int *NumBlockEntriesPerRow) 
  : Epetra_DistObject(rowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    NumMyBlockRows_(rowMap.NumMyElements()),
    CV_(CV),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, rowMap, NumBlockEntriesPerRow);
  int err = Allocate();
  assert( err == 0 );
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& rowMap, int NumBlockEntriesPerRow) 
  : Epetra_DistObject(rowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    NumMyBlockRows_(rowMap.NumMyElements()),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, rowMap, NumBlockEntriesPerRow);
  int err = Allocate();
  assert( err == 0 );
}
//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& rowMap, 
           const Epetra_BlockMap& colMap, int *NumBlockEntriesPerRow) 
  : Epetra_DistObject(rowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    NumMyBlockRows_(rowMap.NumMyElements()),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, rowMap, colMap, NumBlockEntriesPerRow);
  int err = Allocate();
  assert( err == 0 );
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_BlockMap& rowMap, 
           const Epetra_BlockMap& colMap, int NumBlockEntriesPerRow) 
  : Epetra_DistObject(rowMap, "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    NumMyBlockRows_(rowMap.NumMyElements()),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Graph_ = new Epetra_CrsGraph(CV, rowMap, colMap, NumBlockEntriesPerRow);
  int err = Allocate();
  assert( err == 0 );
}
//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph & graph) 
  : Epetra_DistObject(graph.RowMap(), "Epetra::VbrMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(new Epetra_CrsGraph(graph)),
    Allocated_(false),
    StaticGraph_(true),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    NumMyBlockRows_(graph.RowMap().NumMyElements()),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  constructedWithFilledGraph_ = graph.Filled();
  InitializeDefaults();
  int err = Allocate();
  assert(err==0);
}

//==============================================================================
Epetra_VbrMatrix::Epetra_VbrMatrix(const Epetra_VbrMatrix & Source) 
  : Epetra_DistObject(Source),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(new Epetra_CrsGraph(Source.Graph())),
    Allocated_(Source.Allocated_),
    StaticGraph_(true),
    UseTranspose_(Source.UseTranspose_),
    constructedWithFilledGraph_(Source.constructedWithFilledGraph_),
    matrixFillCompleteCalled_(Source.matrixFillCompleteCalled_),
    NumMyBlockRows_(0),
    CV_(Copy),
    HavePointObjects_(false),
    squareFillCompleteCalled_(false)
 {

  InitializeDefaults();
  operator=(Source);
}

//==============================================================================
Epetra_VbrMatrix& Epetra_VbrMatrix::operator=(const Epetra_VbrMatrix& src)
{
  if (this == &src) {
    return(*this);
  }

  DeleteMemory();

  Allocated_ = src.Allocated_;
  StaticGraph_ = src.StaticGraph_;
  UseTranspose_ = src.UseTranspose_;
  NumMyBlockRows_ = src.NumMyBlockRows_;
  CV_ = src.CV_;

  InitializeDefaults();

  //our Graph_ member was delete in DeleteMemory() so we don't need to
  //delete it now before re-creating it.
  Graph_ = new Epetra_CrsGraph(src.Graph());

  int err = Allocate();
  assert( err == 0 );

  int i, j;

#if 0  
  bool ConstantShape = true;
  const int NOTSETYET = -13 ; 
  int MyLda = NOTSETYET ; 
  int MyColDim = NOTSETYET ; 
  int MyRowDim = NOTSETYET ; 
  //  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  //  if (ierr) EPETRA_CHK_ERR(ierr);
  //  if ( ierr != 0 ) ConstantShape = false ;   // Do we really need this?
  assert( ConstantShape ) ; 
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumBlockEntries = NumBlockEntriesPerRow_[i];
    for (j=0; j < NumBlockEntries; j++) {
      Epetra_SerialDenseMatrix* ThisBlock = src.Entries_[i][j];
      if ( MyLda == NOTSETYET ) {
  MyLda = ThisBlock->LDA() ; 
        MyColDim = ThisBlock->ColDim() ; 
        MyRowDim = ThisBlock->RowDim() ; 
      } else {
  if ( MyLda !=  ThisBlock->LDA() ) ConstantShape = false ; 
  if ( MyRowDim !=  ThisBlock->RowDim() ) ConstantShape = false ; 
  if ( MyColDim !=  ThisBlock->ColDim() ) ConstantShape = false ; 
      }
    }
  }  
  if ( false &&  ConstantShape ) {

    int NumMyNonzeros = src.Graph_->NumMyNonzeros();

    All_Values_ = new double[NumMyNonzeros];
    All_Values_Orig_ = All_Values_ ;
    for (i=0; i<NumMyBlockRows_; i++) {
      double* Values_ThisBlockRow = All_Values_ ; 
      int NumBlockEntries = NumBlockEntriesPerRow_[i];
      for (j=0; j < NumBlockEntries; j++) {
  double* Values_ThisBlockEntry = All_Values_ ; 
  Epetra_SerialDenseMatrix* M_SDM = src.Entries_[i][j] ;
  for ( int kk = 0 ; kk < MyColDim ; kk++ ) { 
    for ( int ll = 0 ; ll < MyRowDim ; ll++ ) { 
      *All_Values_ = (*M_SDM)(ll,kk); 
      All_Values_++ ; 
    }
  }
  Entries_[i][j] = new Epetra_SerialDenseMatrix(View, Values_ThisBlockEntry,
                  MyLda, MyRowDim, MyColDim );
      }
    }
  } else { 
#endif    
    for (i=0; i<NumMyBlockRows_; i++) {
      int NumBlockEntries = NumBlockEntriesPerRow_[i];
      for (j=0; j < NumBlockEntries; j++) {
        //if src.Entries_[i][j] != 0, set Entries_[i][j] to be
        //a copy of it.
        Entries_[i][j] = src.Entries_[i][j] != 0 ?
          new Epetra_SerialDenseMatrix(*(src.Entries_[i][j])) : 0;
      }
#if 0
    }
#endif
  }

  if ( src.StorageOptimized() ) this->OptimizeStorage() ;
  
  return( *this );
}

//==============================================================================
void Epetra_VbrMatrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  UseTranspose_ = false;
  Entries_ = 0;
  All_Values_ = 0;
  All_Values_Orig_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  NormFrob_ = -1.0;
  ImportVector_ = 0;

  NumBlockEntriesPerRow_  = 0;
  NumAllocatedBlockEntriesPerRow_ = 0;
  Indices_ = 0;
  ElementSizeList_ = 0;
  FirstPointInElementList_ = 0;
  
  // State variables needed for constructing matrix entry-by-entry

  TempRowDims_ = 0;
  TempEntries_ = 0;
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

  // Attributes that support the Epetra_RowMatrix and Epetra_Operator interfaces
  RowMatrixRowMap_ = 0;
  RowMatrixColMap_ = 0;
  RowMatrixImporter_ = 0;
  OperatorDomainMap_ = 0;
  OperatorRangeMap_ = 0;
  HavePointObjects_ = false;
  StorageOptimized_ = false ; 
  
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
  

  // Allocate Entries array
  Entries_ = new Epetra_SerialDenseMatrix**[NumMyBlockRows_];
  // Allocate and initialize entries
  for (i=0; i<NumMyBlockRows_; i++) {
    int NumAllocatedBlockEntries = NumAllocatedBlockEntriesPerRow_[i];
    
    if (NumAllocatedBlockEntries > 0) {
      Entries_[i] = new Epetra_SerialDenseMatrix*[NumAllocatedBlockEntries];
      for (j=0; j < NumAllocatedBlockEntries; j++) {
        Entries_[i][j] = 0;
      }
    }
    else {
      Entries_[i] = 0;
    }
  }
  SetAllocated(true);
  return(0);
}
//==============================================================================
Epetra_VbrMatrix::~Epetra_VbrMatrix()
{
  DeleteMemory();
}

//==============================================================================
void Epetra_VbrMatrix::DeleteMemory()
{
  int i;

  for (i=0; i<NumMyBlockRows_; i++) {
    int NumAllocatedBlockEntries = NumAllocatedBlockEntriesPerRow_[i];

    if (NumAllocatedBlockEntries >0) {

      for (int j=0; j < NumAllocatedBlockEntries; j++) {
        if (Entries_[i][j]!=0) {
          delete Entries_[i][j];
        }
      }
      delete [] Entries_[i];
    }
  }

  if (All_Values_Orig_!=0)   delete [] All_Values_Orig_;
  All_Values_Orig_ = 0;

  if (Entries_!=0)       delete [] Entries_;
  Entries_ = 0;

  if (ImportVector_!=0) delete ImportVector_;
  ImportVector_ = 0;

  NumMyBlockRows_ = 0;

  if (LenTemps_>0) {
    delete [] TempRowDims_;
    delete [] TempEntries_;
  }

  // Delete any objects related to supporting the RowMatrix and Operator interfaces
  if (HavePointObjects_) {
#if 1
    if ( RowMatrixColMap_ != RowMatrixRowMap_ )        delete RowMatrixColMap_;
    if ( OperatorDomainMap_ != RowMatrixRowMap_ )        delete OperatorDomainMap_;
    if ( OperatorRangeMap_ != RowMatrixRowMap_ )        delete OperatorRangeMap_;
#else
    //  this can fail, see bug # 1114
    if (!RowMatrixColMap().SameAs(RowMatrixRowMap())) delete RowMatrixColMap_;
    if (!OperatorDomainMap().SameAs(RowMatrixRowMap())) delete OperatorDomainMap_;
    if (!OperatorRangeMap().SameAs(RowMatrixRowMap())) delete OperatorRangeMap_;
#endif
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

  delete Graph_; // We created the graph, so must delete it.
  Graph_ = 0;
}

//==============================================================================
int Epetra_VbrMatrix::PutScalar(double ScalarConstant) 
{
  if (!Allocated_) return(0);

  for (int i=0; i<NumMyBlockRows_; i++) {
    int NumBlockEntries = NumBlockEntriesPerRow_[i];
    int RowDim = ElementSizeList_[i];
    for (int j=0; j< NumBlockEntries; j++) {
      int LDA = Entries_[i][j]->LDA();
      int ColDim = Entries_[i][j]->N();
      for (int col=0; col < ColDim; col++) {
  double * Entries = Entries_[i][j]->A()+col*LDA;
  for (int row=0; row < RowDim; row++)
    *Entries++ = ScalarConstant;
      }
    }
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;
  return(0);
}

//==============================================================================
int Epetra_VbrMatrix::Scale(double ScalarConstant) 
{
  for (int i=0; i<NumMyBlockRows_; i++) {
    int NumBlockEntries = NumBlockEntriesPerRow_[i];
    int RowDim = ElementSizeList_[i];
    for (int j=0; j< NumBlockEntries; j++) {
      int LDA = Entries_[i][j]->LDA();
      int ColDim = Entries_[i][j]->N();
      for (int col=0; col < ColDim; col++) {
  double * Entries = Entries_[i][j]->A()+col*LDA;
  for (int row=0; row < RowDim; row++)
    *Entries++ *= ScalarConstant;
      }
    }
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::BeginInsertGlobalValues(int BlockRow, int NumBlockEntries, int * BlockIndices) {

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  Graph_->SetIndicesAreGlobal(true);
  int LocalBlockRow = LRID(BlockRow); // Find local row number for this global row index
  
  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginInsertValues(LocalBlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginInsertMyValues(int  BlockRow, int NumBlockEntries,
            int * BlockIndices)
{
  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot insert local values until MakeIndicesLocal() is called either directly or through FillComplete()
  Graph_->SetIndicesAreLocal(true);
  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginInsertValues(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);

}

//==========================================================================
int Epetra_VbrMatrix::BeginInsertValues(int BlockRow, int NumBlockEntries, 
          int * BlockIndices, bool indicesAreLocal)
{
  if (StaticGraph()) EPETRA_CHK_ERR(-2); // If the matrix graph is fully constructed, we cannot insert new values

  int ierr = 0;

  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) {
    EPETRA_CHK_ERR(-1); // Not in BlockRow range    
  }
  if (CV_==View && Entries_[BlockRow]!=0) ierr = 2; // This row has be defined already. Issue warning.    
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new

  // Set up pointers, make sure enough temp space for this round of submits

  Epetra_CombineMode SubmitMode = Insert;

  EPETRA_CHK_ERR(ierr);
  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal, SubmitMode));
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::BeginReplaceGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  BlockRow = LRID(BlockRow); // Normalize row range
  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginReplaceMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  if (!Graph_->IndicesAreLocal()) EPETRA_CHK_ERR(-1);

  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginReplaceValues(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginReplaceValues(int BlockRow,
           int NumBlockEntries, 
           int *BlockIndices,
           bool indicesAreLocal)
{
  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in BlockRow range

  Epetra_CombineMode SubmitMode = Zero; // This is a misuse of Zero mode, fix it later
  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal, SubmitMode));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  BlockRow = LRID(BlockRow); // Normalize row range
  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices) {

  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginSumIntoValues(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginSumIntoValues(int BlockRow, int NumBlockEntries, 
           int *BlockIndices, bool indicesAreLocal) {

  if (BlockRow < 0 || BlockRow >= NumMyBlockRows_) EPETRA_CHK_ERR(-1); // Not in BlockRow range

  Epetra_CombineMode SubmitMode = Add;
  EPETRA_CHK_ERR(SetupForSubmits(BlockRow, NumBlockEntries, BlockIndices, indicesAreLocal, SubmitMode));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SetupForSubmits(int BlockRow, int NumBlockEntries,
              int * BlockIndices, 
              bool indicesAreLocal,
              Epetra_CombineMode SubmitMode) {

  if (NumBlockEntries>LenTemps_) {
    if (LenTemps_>0) {
      delete [] TempRowDims_;
      delete [] TempEntries_;
    }
    TempRowDims_ = new int[NumBlockEntries];
    TempEntries_ = new Epetra_SerialDenseMatrix*[NumBlockEntries];
    LenTemps_ = NumBlockEntries;
  }

  CurBlockRow_ = BlockRow;
  CurNumBlockEntries_ = NumBlockEntries;
  CurBlockIndices_ = BlockIndices;
  CurEntry_ = 0;
  CurIndicesAreLocal_ = indicesAreLocal;
  CurSubmitMode_ = SubmitMode;
    
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::DirectSubmitBlockEntry(int GlobalBlockRow, int GlobalBlockCol,
                                             const double *values, int LDA,
                                             int NumRows, int NumCols, bool sum_into)
{
  int ierr = 0;
  int LocalBlockRow = LRID(GlobalBlockRow); // Normalize row range
  int Loc;
  bool found = Graph_->FindGlobalIndexLoc(LocalBlockRow, GlobalBlockCol,0,Loc);
  if (found) {
    if (Entries_[LocalBlockRow][Loc]==0) {
      Entries_[LocalBlockRow][Loc] =
        new Epetra_SerialDenseMatrix(Copy, const_cast<double*>(values), LDA, NumRows, NumCols, false);
    }
    else {
      Epetra_SerialDenseMatrix* target = Entries_[LocalBlockRow][Loc];
 
      target->CopyMat(values, LDA, NumRows, NumCols,
                      target->A_, target->LDA_, sum_into);
    }
  }
  else {
    ierr = -2;
  }
  EPETRA_CHK_ERR(ierr);
  return ierr;
}

//==========================================================================
int Epetra_VbrMatrix::SubmitBlockEntry(double *values, int LDA,
               int NumRows, int NumCols)
{
  if (CurEntry_==-1) EPETRA_CHK_ERR(-1); // This means that a Begin routine was not called
  if (CurEntry_>=CurNumBlockEntries_) EPETRA_CHK_ERR(-4); // Exceeded the number of entries that can be submitted

  // Fill up temp space with entry

  TempRowDims_[CurEntry_] = NumRows;
  TempEntries_[CurEntry_] = new Epetra_SerialDenseMatrix(CV_, values, LDA,
               NumRows, NumCols, false);
  CurEntry_++;

  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SubmitBlockEntry( Epetra_SerialDenseMatrix &Mat)
{
  return SubmitBlockEntry( Mat.A(), Mat.LDA(), Mat.M(), Mat.N() );
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
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;
 return(0);
}

//==========================================================================
int Epetra_VbrMatrix::EndReplaceSumIntoValues()
{
  int j;
  int ierr = 0;
  int Loc;

  int RowDim = ElementSizeList_[CurBlockRow_];

  bool SumInto = (CurSubmitMode_==Add);

  for (j=0; j<CurNumBlockEntries_; j++) {
    int BlockIndex = CurBlockIndices_[j];

    bool code = false;
    if (CurIndicesAreLocal_) {
      code =Graph_->FindMyIndexLoc(CurBlockRow_,BlockIndex,j,Loc);
    }
    else {
      code =Graph_->FindGlobalIndexLoc(CurBlockRow_,BlockIndex,j,Loc);
    }

    bool need_to_delete_temp_entry = true;

    if (code) {
      Epetra_SerialDenseMatrix* src = TempEntries_[j];

      int ColDim = src->N();

      if (Entries_[CurBlockRow_][Loc]==0) {
        Entries_[CurBlockRow_][Loc] = src;
        need_to_delete_temp_entry = false;
      }
      else {
        Epetra_SerialDenseMatrix* target = Entries_[CurBlockRow_][Loc];

        target->CopyMat(src->A_, src->LDA_, RowDim, ColDim,
                        target->A_, target->LDA_, SumInto);
      }
    }
    else {
      ierr=2; // Block Discarded, Not Found
    }

    if (need_to_delete_temp_entry) {
      delete TempEntries_[j];
    }
  }

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::EndInsertValues()
{
  int ierr = 0;
  int j;

  int NumValidBlockIndices = CurNumBlockEntries_;
  int * ValidBlockIndices = new int[CurNumBlockEntries_];
  for ( j=0; j < CurNumBlockEntries_; ++j ) ValidBlockIndices[j] = j;
    
  if( Graph_->HaveColMap() ) { //test and discard indices not in ColMap
    NumValidBlockIndices = 0;
    const Epetra_BlockMap& map = Graph_->ColMap();

    for ( j = 0; j < CurNumBlockEntries_; ++j ) {
      bool myID = CurIndicesAreLocal_ ?
        map.MyLID(CurBlockIndices_[j]) : map.MyGID(CurBlockIndices_[j]);

      if( myID ) {
        ValidBlockIndices[ NumValidBlockIndices++ ] = j;
      }
      else ierr=2; // Discarding a Block not found in ColMap
    }
  }

  int oldNumBlocks = NumBlockEntriesPerRow_[CurBlockRow_];
  int oldNumAllocBlocks = NumAllocatedBlockEntriesPerRow_[CurBlockRow_];

  // Update graph
  ierr = Graph_->InsertIndices(CurBlockRow_, CurNumBlockEntries_, CurBlockIndices_);

  int newNumAllocBlocks = NumAllocatedBlockEntriesPerRow_[CurBlockRow_];

  if (newNumAllocBlocks > oldNumAllocBlocks) {
    ierr = 3;//positive warning code indicates we expanded allocated memory
    Epetra_SerialDenseMatrix** tmp_Entries = new Epetra_SerialDenseMatrix*[newNumAllocBlocks];
    for(j=0; j<oldNumBlocks; ++j) tmp_Entries[j] = Entries_[CurBlockRow_][j];
    for(j=oldNumBlocks; j<newNumAllocBlocks; ++j) tmp_Entries[j] = NULL;
    delete [] Entries_[CurBlockRow_];
    Entries_[CurBlockRow_] = tmp_Entries;
  }

  int start = oldNumBlocks;
  int stop = start + NumValidBlockIndices;
  int NumAllocatedEntries = NumAllocatedBlockEntriesPerRow_[CurBlockRow_];

  if (stop <= NumAllocatedEntries) {
    for (j=start; j<stop; j++) {
      Epetra_SerialDenseMatrix& mat =
        *(TempEntries_[ValidBlockIndices[j-start]]);
  
      Entries_[CurBlockRow_][j] = new Epetra_SerialDenseMatrix(CV_, mat.A(),
          mat.LDA(),
          mat.M(),
          mat.N());
    }
  }
  else {
    ierr = -4;
  }


  delete [] ValidBlockIndices;

  for (j=0; j<CurNumBlockEntries_; ++j) {
    delete TempEntries_[j];
  }

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::CopyMat(double * A, int LDA, int NumRows, int NumCols, 
            double * B, int LDB, bool SumInto) const
{
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
int Epetra_VbrMatrix::FillComplete() {
  squareFillCompleteCalled_ = true;
  EPETRA_CHK_ERR(FillComplete(RowMap(), RowMap()));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::FillComplete(const Epetra_BlockMap& domain_map,
           const Epetra_BlockMap& range_map)
{
  int returnValue = 0;

  if (Graph_->Filled()) {
    if (!constructedWithFilledGraph_ && !matrixFillCompleteCalled_) {
      returnValue = 2;
    }
  }

 if(!StaticGraph()) {
    EPETRA_CHK_ERR(Graph_->MakeIndicesLocal(domain_map, range_map));
  }

  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values

  if(!StaticGraph()) {
    EPETRA_CHK_ERR(Graph_->FillComplete(domain_map, range_map));
  }

  // NumMyCols_ = Graph_->NumMyCols(); // Redefine based on local number of cols

  matrixFillCompleteCalled_ = true;

  if (squareFillCompleteCalled_) {
    if (DomainMap().NumGlobalElements64() != RangeMap().NumGlobalElements64()) {
      returnValue = 3;
    }
    squareFillCompleteCalled_ = false;
    EPETRA_CHK_ERR(returnValue);
  }

  return(returnValue);
}

//==========================================================================
int Epetra_VbrMatrix::TransformToLocal() {
  EPETRA_CHK_ERR(FillComplete());
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::TransformToLocal(const Epetra_BlockMap* domainMap, const Epetra_BlockMap* rangeMap) {
  EPETRA_CHK_ERR(FillComplete(*domainMap, *rangeMap));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SortEntries() {

  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-1);
  if (Sorted()) return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for (int i=0; i<NumMyBlockRows_; i++){

    Epetra_SerialDenseMatrix ** Entries = Entries_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
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
        Epetra_SerialDenseMatrix *dtemp = Entries[k+m];
        Entries[k+m] = Entries[k];
        Entries[k] = dtemp;

        int itemp = Indices[k+m];
        Indices[k+m] = Indices[k];
        Indices[k] = itemp;
            }
        }
      m = m/2;
    }
  }
  Graph_->SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::MergeRedundantEntries()
{

  if (NoRedundancies()) return(0);
  if (!Sorted()) EPETRA_CHK_ERR(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal
  // Note:  This function assumes that SortEntries was already called.

  bool SumInto = true;
  for (int i=0; i<NumMyBlockRows_; i++){
    int NumEntries = NumBlockEntriesPerRow_[i];
    if (NumEntries>1) {
      Epetra_SerialDenseMatrix ** const Entries = Entries_[i];
      int * const Indices = Indices_[i];
      int RowDim = ElementSizeList_[i];
      int curEntry =0;
      Epetra_SerialDenseMatrix* curBlkEntry = Entries[0];
      for (int k=1; k<NumEntries; k++) {
  if (Indices[k]==Indices[k-1]) {
    CopyMat(Entries[k]->A(), Entries[k]->LDA(), RowDim, Entries[k]->N(),
      curBlkEntry->A(), curBlkEntry->LDA(), SumInto);
  }
  else {
    CopyMat(curBlkEntry->A(), curBlkEntry->LDA(), RowDim, curBlkEntry->N(),
      Entries[curEntry]->A(), Entries[curEntry]->LDA(), false);
    curEntry++;
    curBlkEntry = Entries[k];
  }
      }
      CopyMat(curBlkEntry->A(), curBlkEntry->LDA(), RowDim, curBlkEntry->N(),
        Entries[curEntry]->A(), Entries[curEntry]->LDA(), false);
    }
  }
    
  EPETRA_CHK_ERR(Graph_->RemoveRedundantIndices()); // Remove redundant indices and then return
  return(0);

}

//==========================================================================
int Epetra_VbrMatrix::OptimizeStorage() {

  if (StorageOptimized()) return(0); // Have we been here before?

  // cout << __FILE__ << " " << __LINE__ << " " << Graph_->NumMyNonzeros() << endl ; 

  bool ConstantShape = true;
  int i,j;
  const int NOTSETYET = -13 ; 
  int MyLda = NOTSETYET ; 
  int MyColDim = NOTSETYET ; 
  int MyRowDim = NOTSETYET ; 
  //
  //  I don't know why the underlying graph musthave optimized storage, but 
  //  the test for optimized storage depends upon the graph hainvg optimized 
  //  storage and that is enough for me.  
  //
  //  Ideally, we would like to have optimized storage for the graph.  But, 
  //  for now, this is commented out until bug 1097 is fixed
  //
  //  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  //  if (ierr) EPETRA_CHK_ERR(ierr);
  //  if ( ierr != 0 ) ConstantShape = false ; 
  if( ConstantShape ) {
    for (i=0; i<NumMyBlockRows_; i++) {
      int NumBlockEntries = NumBlockEntriesPerRow_[i];
      for (j=0; j < NumBlockEntries; j++) {
  Epetra_SerialDenseMatrix* ThisBlock = Entries_[i][j];
  if ( MyLda == NOTSETYET ) {
    MyLda = ThisBlock->LDA() ; 
    MyColDim = ThisBlock->ColDim() ; 
    MyRowDim = ThisBlock->RowDim() ; 
  } else {
    if ( MyLda !=  ThisBlock->LDA() ) ConstantShape = false ; 
    if ( MyRowDim !=  ThisBlock->RowDim() ) ConstantShape = false ; 
    if ( MyColDim !=  ThisBlock->ColDim() ) ConstantShape = false ; 
  }
      }
    }    
  }  
  // cout << __FILE__ << " " << __LINE__ << " " << Graph_->NumMyNonzeros() << endl ; 

  if ( ConstantShape ) {

    int numMyNonzeros = Graph_->NumMyEntries();
    int coef_len = MyColDim*MyRowDim*numMyNonzeros;
    All_Values_ = new double[coef_len];
    All_Values_Orig_ = All_Values_ ;
    for (i=0; i<NumMyBlockRows_; i++) {
      int NumBlockEntries = NumBlockEntriesPerRow_[i];
      for (j=0; j < NumBlockEntries; j++) {
  double* Values_ThisBlockEntry = All_Values_ ; 
  Epetra_SerialDenseMatrix* M_SDM = Entries_[i][j] ;
  for ( int kk = 0 ; kk < MyColDim ; kk++ ) { 
    for ( int ll = 0 ; ll < MyRowDim ; ll++ ) { 
      *All_Values_ = (*M_SDM)(ll,kk) ;
      All_Values_++ ; 
    }
  }
  //  Entries_[i][j] = new Epetra_SerialDenseMatrix(*(src.Entries_[i][j]));
  delete Entries_[i][j]; 
  Entries_[i][j] = new Epetra_SerialDenseMatrix(View, Values_ThisBlockEntry,
                  MyLda, MyRowDim, MyColDim );
      }
    }
    StorageOptimized_ = true ; 
  }

  // cout << __FILE__ << " " << __LINE__ << " " << Graph_->NumMyNonzeros() << endl ; 


  /* Work on later...
     int i, j;

     // The purpose of this routine is to make the block entries in each row contiguous in memory
     // so that a single call to GEMV or GEMM call be used to compute an entire block row.


     bool Contiguous = true; // Assume contiguous is true
     for (i=1; i<NumMyBlockRows_; i++){
     int NumEntries = NumBlockEntriesPerRow_[i-1];
     int NumAllocatedEntries = NumAllocatedBlockEntriesPerRow_[i-1];
      
     // Check if NumEntries is same as NumAllocatedEntries and 
     // check if end of beginning of current row starts immediately after end of previous row.
     if ((NumEntries!=NumAllocatedEntries) || (Entries_[i]!=Entries_[i-1]+NumEntries)) {
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
     int numMyNonzeros = Graph_->NumMyNonzeros();

     // Allocate one big integer array for all index values
     All_Values_ = new double[numMyNonzeros];
  
     // Set Entries_ to point into All_Entries_
  
     double * tmp = All_Values_;
     for (i=0; i<NumMyBlockRows_; i++) {
     int NumEntries = NumEntriesPerBlockRow_[i];
     for (j=0; j<NumEntries; j++) tmp[j] = Entries_[i][j];
     if (Entries_[i] !=0) delete [] Entries_[i];
     Entries_[i] = tmp;
     tmp += NumEntries;
     }
  */
  // cout << __FILE__ << " " << __LINE__ << " " << Graph_->NumMyNonzeros() << endl ; 
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::ExtractGlobalRowCopy(int GlobalRow,
             int Length, 
             int & NumEntries,
             double *values,
             int * Indices) const
{
  (void)GlobalRow;
  (void)Length;
  (void)NumEntries;
  (void)values;
  (void)Indices;
  cout << "Must implement..." << endl;
  return(0);
}
//==========================================================================
int Epetra_VbrMatrix::ExtractGlobalBlockRowPointers(int BlockRow, int maxNumBlockEntries, 
                int & RowDim, int & NumBlockEntries, 
                int * BlockIndices,
                Epetra_SerialDenseMatrix** & Entries) const {

  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(ExtractBlockRowPointers(BlockRow, maxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
           Entries, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractMyBlockRowPointers(int BlockRow, int maxNumBlockEntries, 
            int & RowDim, int & NumBlockEntries, 
            int * BlockIndices,
            Epetra_SerialDenseMatrix** & Entries) const {

  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(ExtractBlockRowPointers(BlockRow,maxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
           Entries, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractBlockRowPointers(int BlockRow, int maxNumBlockEntries, 
                int & RowDim, int & NumBlockEntries, 
                int * BlockIndices,
                Epetra_SerialDenseMatrix** & Entries, 
                bool indicesAreLocal) const {
  int ierr = 0;
  if (!indicesAreLocal) {
    ierr = Graph_->ExtractGlobalRowCopy(BlockRow, maxNumBlockEntries,
          NumBlockEntries, BlockIndices);
    BlockRow = LRID(BlockRow);
  }
  else {
    ierr = Graph_->ExtractMyRowCopy(BlockRow, maxNumBlockEntries,
            NumBlockEntries, BlockIndices);
  }
  if (ierr) EPETRA_CHK_ERR(ierr);

  RowDim = ElementSizeList_[BlockRow];

  Entries = Entries_[BlockRow];


  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractGlobalBlockRowCopy(int BlockRow, int maxNumBlockEntries, 
                 int & RowDim, int & NumBlockEntries, 
                 int * BlockIndices, int * ColDims) const {

  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowCopy(BlockRow, maxNumBlockEntries, RowDim, NumBlockEntries, BlockIndices, 
            ColDims, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractMyBlockRowCopy(int BlockRow, int maxNumBlockEntries, 
             int & RowDim, int & NumBlockEntries, 
             int * BlockIndices, int * ColDims) const {

  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowCopy(BlockRow,maxNumBlockEntries , RowDim, NumBlockEntries, BlockIndices, 
            ColDims, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractBlockRowCopy(int BlockRow, int maxNumBlockEntries, 
                 int & RowDim, int & NumBlockEntries, 
                 int * BlockIndices, int * ColDims, 
                 bool indicesAreLocal) const  {
  int ierr = 0;
  if (!indicesAreLocal)
    ierr = Graph_->ExtractGlobalRowCopy(BlockRow, maxNumBlockEntries, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowCopy(BlockRow, maxNumBlockEntries, NumBlockEntries, BlockIndices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  bool ExtractView = false;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, indicesAreLocal);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractBlockDimsCopy(NumBlockEntries, ColDims));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::SetupForExtracts(int BlockRow, int & RowDim, int NumBlockEntries, bool ExtractView, 
               bool indicesAreLocal) const {

  if (!indicesAreLocal) BlockRow = LRID(BlockRow); // Normalize row range
  CurExtractBlockRow_ = BlockRow;
  CurExtractEntry_ = 0;
  CurExtractNumBlockEntries_ = NumBlockEntries;
  CurExtractIndicesAreLocal_ = indicesAreLocal;
  CurExtractView_ = ExtractView;
  CurRowDim_ = ElementSizeList_[CurExtractBlockRow_];
  RowDim = CurRowDim_;
    
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractBlockDimsCopy(int NumBlockEntries, int * ColDims) const {

  for (int i=0; i<NumBlockEntries; i++) {
    ColDims[i] = Entries_[CurExtractBlockRow_][i]->N();
  }
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractEntryCopy(int SizeOfValues,
               double * values,
               int LDA,
               bool SumInto) const
{
  (void)SumInto;
  if (CurExtractEntry_==-1) EPETRA_CHK_ERR(-1); // No BeginCopy routine was called
  int CurColDim = Entries_[CurExtractBlockRow_][CurExtractEntry_]->N();
  if (LDA*CurColDim>SizeOfValues) EPETRA_CHK_ERR(-2);  // Not enough space

  Epetra_SerialDenseMatrix* CurEntries = Entries_[CurExtractBlockRow_][CurExtractEntry_];
  int CurRowDim = CurEntries->M();
  int CurLDA = CurEntries->LDA();

  CurExtractEntry_++; // Increment Entry Pointer

  double* vals = CurEntries->A();
  if (LDA==CurRowDim && CurLDA==CurRowDim) // Columns of the entry are contiguous, can treat as single array
    for (int ii=0; ii<CurRowDim*CurColDim; ++ii)
      values[ii] = vals[ii];
  else {
    double * CurTargetCol = values;
    double * CurSourceCol = vals;

    for (int jj=0; jj<CurColDim; jj++) {
      for (int ii=0; ii<CurRowDim; ++ii) 
  CurTargetCol[ii] = CurSourceCol[ii];
      CurTargetCol += LDA;
      CurSourceCol += CurLDA;
    }
  }

  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractGlobalBlockRowView(int BlockRow,
                 int & RowDim, int & NumBlockEntries, 
                 int * & BlockIndices) const
{

  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries,
            BlockIndices, indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
             int * & BlockIndices)  const
{

  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
            indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::BeginExtractBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
                 int * & BlockIndices,
                 bool indicesAreLocal) const
{
  int ierr = 0;
  if (!indicesAreLocal)
    ierr = Graph_->ExtractGlobalRowView(BlockRow, NumBlockEntries, BlockIndices);
  else
    ierr = Graph_->ExtractMyRowView(BlockRow,  NumBlockEntries, BlockIndices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  bool ExtractView = true;
  ierr = SetupForExtracts(BlockRow, RowDim, NumBlockEntries, ExtractView, indicesAreLocal);
  if (ierr) EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractEntryView(Epetra_SerialDenseMatrix* & entry) const
{
  entry = Entries_[CurExtractBlockRow_][CurExtractEntry_];

  CurExtractEntry_++; // Increment Entry Pointer
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
            int * & BlockIndices,
            Epetra_SerialDenseMatrix** & Entries) const
{

  Entries = Entries_[LRID(BlockRow)]; // Pointer to Array of pointers for this row's block entries
  bool indicesAreLocal = false;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries,
            BlockIndices, 
            indicesAreLocal));
  return(0);
}

//==========================================================================
int Epetra_VbrMatrix::ExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
              int * & BlockIndices,
              Epetra_SerialDenseMatrix** & Entries) const
{

  Entries = Entries_[BlockRow]; // Pointer to Array of pointers for this row's block entries
  bool indicesAreLocal = true;
  EPETRA_CHK_ERR(BeginExtractBlockRowView(BlockRow, RowDim, NumBlockEntries, BlockIndices, 
            indicesAreLocal));
  return(0);
}

//==============================================================================
int Epetra_VbrMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
  
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) EPETRA_CHK_ERR(-2); // Maps must be the same
  double * diagptr = Diagonal.Values();
  for (int i=0; i<NumMyBlockRows_; i++){
    int BlockRow = GRID64(i); // FIXME long long
    int RowDim = ElementSizeList_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
    for (int j=0; j<NumEntries; j++) {
      int BlockCol = GCID64(Indices[j]);// FIXME long long
      if (BlockRow==BlockCol) {
  CopyMatDiag(Entries_[i][j]->A(), Entries_[i][j]->LDA(), RowDim,
        Entries_[i][j]->N(), diagptr+FirstPointInElementList_[i]);
  break;
      }
    }
  }
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::ReplaceDiagonalValues(const Epetra_Vector & Diagonal) {
  
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) EPETRA_CHK_ERR(-2); // Maps must be the same
  int ierr = 0;
  double * diagptr = (double *) Diagonal.Values(); // Dangerous but being lazy
  for (int i=0; i<NumMyBlockRows_; i++){
    int BlockRow = GRID64(i);// FIXME long long
    int RowDim = ElementSizeList_[i];
    int NumEntries = NumBlockEntriesPerRow_[i];
    int * Indices = Indices_[i];
    bool DiagMissing = true;
    for (int j=0; j<NumEntries; j++) {
      int BlockCol = GCID64(Indices[j]);// FIXME long long
      if (BlockRow==BlockCol) {
  ReplaceMatDiag(Entries_[i][j]->A(), Entries_[i][j]->LDA(), RowDim, Entries_[i][j]->N(), 
           diagptr+FirstPointInElementList_[i]);
  DiagMissing = false;
  break;
      }
    }
    if (DiagMissing) ierr = 1; // flag a warning error
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;
  EPETRA_CHK_ERR(ierr);
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::BeginExtractBlockDiagonalCopy(int maxNumBlockDiagonalEntries, 
                int & NumBlockDiagonalEntries, int * RowColDims ) const{
  
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  CurBlockDiag_ = 0; // Initialize pointer
  NumBlockDiagonalEntries = NumMyBlockRows_;
  if (NumBlockDiagonalEntries>maxNumBlockDiagonalEntries) EPETRA_CHK_ERR(-2);
  EPETRA_CHK_ERR(RowMap().ElementSizeList(RowColDims));
  return(0);
}
//==============================================================================
int Epetra_VbrMatrix::ExtractBlockDiagonalEntryCopy(int SizeOfValues, double * values, int LDA, bool SumInto ) const {
  
  if (CurBlockDiag_==-1) EPETRA_CHK_ERR(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow = i;
  int RowDim = ElementSizeList_[i];
  int NumEntries = NumBlockEntriesPerRow_[i];
  int * Indices = Indices_[i];
  for (int j=0; j<NumEntries; j++) {
    int Col = Indices[j];
    if (BlockRow==Col) {
      int ColDim = Entries_[i][j]->N();
      if (LDA*ColDim>SizeOfValues) EPETRA_CHK_ERR(-2); // Not enough room in values
      CopyMat(Entries_[i][j]->A(), Entries_[i][j]->LDA(), RowDim, ColDim, values, LDA, SumInto);
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
int Epetra_VbrMatrix::ExtractBlockDiagonalEntryView(double * & values, int & LDA) const {
  
  if (CurBlockDiag_==-1) EPETRA_CHK_ERR(-1); // BeginExtractBlockDiagonalCopy was not called
  int i = CurBlockDiag_;
  int BlockRow = i;
  int NumEntries = NumBlockEntriesPerRow_[i];
  int * Indices = Indices_[i];
  for (int j=0; j<NumEntries; j++) {
    int Col = Indices[j];
    if (BlockRow==Col) {
      values = Entries_[i][j]->A();
      LDA = Entries_[i][j]->LDA();
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
int Epetra_VbrMatrix::ReplaceMatDiag(double * A, int LDA, int NumRows, int NumCols, 
             double * Diagonal) {

  int i;
  double * ptr1 = Diagonal;
  double * ptr2;
  int ndiags = EPETRA_MIN(NumRows,NumCols);

  for (i=0; i<ndiags; i++) {
    ptr2 = A + i*LDA+i;
    *ptr2 = *ptr1++;
  }
  return(0);
}
//=============================================================================
int Epetra_VbrMatrix::MaxNumEntries() const {

  int outval = 0;

  for (int i=0; i<NumMyBlockRows_; i++){
    int NumBlockEntries = NumMyBlockEntries(i);
    int NumEntries = 0;
    for (int j=0; j<NumBlockEntries; j++) NumEntries += Entries_[i][j]->N();
    outval = EPETRA_MAX(outval,NumEntries);
  }
  return(outval);
}
//=============================================================================
int Epetra_VbrMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {

  int BlockRow, BlockOffset;
  int ierr = RowMap().FindLocalElementID(MyRow, BlockRow, BlockOffset);  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int NumBlockEntries = NumMyBlockEntries(BlockRow);
  NumEntries = 0;
  for (int i=0; i<NumBlockEntries; i++) NumEntries += Entries_[BlockRow][i]->N();
  return(0);  
}
//=============================================================================
int Epetra_VbrMatrix::ExtractMyRowCopy(int MyRow,
               int Length,
               int & NumEntries,
               double *values,
               int * Indices) const
{
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't extract row unless matrix is filled
  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-2);

  int ierr = 0;
  int BlockRow, BlockOffset;
  ierr = RowMap().FindLocalElementID(MyRow, BlockRow, BlockOffset);
  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int RowDim, NumBlockEntries;
  int * BlockIndices;
  Epetra_SerialDenseMatrix ** ValBlocks;
  ierr = ExtractMyBlockRowView(BlockRow, RowDim, NumBlockEntries,
             BlockIndices, ValBlocks);
  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  int * ColFirstPointInElementList = FirstPointInElementList_;
  if (Importer()!=0) {
    ColFirstPointInElementList = ColMap().FirstPointInElementList();
  }
  NumEntries = 0;
  for (int i=0; i<NumBlockEntries; i++) {
    int ColDim = ValBlocks[i]->N();
    NumEntries += ColDim;
    if (NumEntries>Length) EPETRA_CHK_ERR(-3); // Not enough space
    int LDA = ValBlocks[i]->LDA();
    double * A = ValBlocks[i]->A()+BlockOffset; // Point to first element in row
    int Index = ColFirstPointInElementList[BlockIndices[i]];
    for (int j=0; j < ColDim; j++) {
      *values++ = *A;
      A += LDA;
      *Indices++ = Index++;
    }
  }

  return(0);      
}
//=============================================================================
int Epetra_VbrMatrix::Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const
{
  //
  // This function forms the product y = A * x or y = A' * x
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  
  int i, j;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * FirstPointInElement = FirstPointInElementList_;
  int * ElementSize = ElementSizeList_;
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;

  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  int * ColElementSizeList = ColMap().ElementSizeList();
  int * ColFirstPointInElementList = ColMap().FirstPointInElementList();

  bool x_and_y_same = false;
  Epetra_Vector* ytemp = 0;
  if (xp == yp) {
    x_and_y_same = true;
    ytemp = new Epetra_Vector(y.Map());
    yp = (double*)ytemp->Values();
  }

  Epetra_MultiVector& yref = x_and_y_same ? *ytemp : y;
  //
  //Now yref is a reference to ytemp if x_and_y_same == true, otherwise it is a reference
  //to y.
  //

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
  if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*)ImportVector_->Values();
      ColElementSizeList = ColMap().ElementSizeList(); // The Import map will always have an existing ElementSizeList
      ColFirstPointInElementList = ColMap().FirstPointInElementList(); // Import map will always have an existing ...
    }


    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
  if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*)ExportVector_->Values();
    }
    
    // Do actual computation
    int NumMyRows_ = NumMyRows();
    for (i=0; i< NumMyRows_; i++) yp[i] = 0.0;  // Initialize y
    
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
      double * cury = yp + *FirstPointInElement++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
  //sum += BlockRowValues[j] * xp[BlockRowIndices[j]];
  double * A = BlockRowValues[j]->A();
  int LDA = BlockRowValues[j]->LDA();
  int Index = BlockRowIndices[j];
  double * curx = xp + ColFirstPointInElementList[Index];
  int ColDim = ColElementSizeList[Index];
  GEMV('N', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);      
      }
    }
    if (Exporter()!=0) {
      yref.PutScalar(0.0);
      EPETRA_CHK_ERR(yref.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(yref.Reduce());
  }
  
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
  if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*)ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
  if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
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
      Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
      double * curx = xp + *FirstPointInElement++;
      int      RowDim = *ElementSize++;
      for (j=0; j < NumEntries; j++) {
  //yp[BlockRowIndices[j]] += BlockRowValues[j] * xp[i];
  double * A = BlockRowValues[j]->A();
  int LDA = BlockRowValues[j]->LDA();
  int Index = BlockRowIndices[j];
  double * cury = yp + ColFirstPointInElementList[Index];
  int ColDim = ColElementSizeList[Index];
  GEMV('T', RowDim, ColDim, 1.0, A, LDA, curx, 1.0, cury);
  
      }
    }
    if (Importer()!=0) {
      yref.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(yref.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(yref.Reduce());
  }

  if (x_and_y_same == true) {
    y = *ytemp;
    delete ytemp;
  }
  
  UpdateFlops(2*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_VbrMatrix::DoMultiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  
  int i;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;
  
  int * RowElementSizeList = ElementSizeList_;
  int * RowFirstPointInElementList = FirstPointInElementList_;
  int * ColElementSizeList = ColMap().ElementSizeList();
  int * ColFirstPointInElementList = ColMap().FirstPointInElementList();

   
  int NumVectors = X.NumVectors();
  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  bool x_and_y_same = false;
  Epetra_MultiVector* Ytemp = 0;
  if (*Xp == *Yp) {
    x_and_y_same = true;
    Ytemp = new Epetra_MultiVector(Y.Map(), NumVectors);
    Yp = (double**)Ytemp->Pointers();
  }

  Epetra_MultiVector& Yref = x_and_y_same ? *Ytemp : Y;
  //
  //Now Yref is a reference to Ytemp if x_and_y_same == true, otherwise it is a reference
  //to Y.
  //

  if (!TransA) {
    
    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
  if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
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
  if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      // Create Export vector if needed
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors);

      ExportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ExportVector_->Pointers();
      RowElementSizeList = ColMap().ElementSizeList();
      RowFirstPointInElementList = ColMap().FirstPointInElementList();
    }
    else {
      Yref.PutScalar(0.0); // Zero y values
    }

    // Do actual computation
    if ( All_Values_ != 0 ) { // Contiguous Storage 
      if ( ! TransA && *RowElementSizeList <= 4 ) { 

  int RowDim = *RowElementSizeList ;

  Epetra_SerialDenseMatrix* Asub = **Entries;
  double *A = Asub->A_ ;

  if ( NumVectors == 1 ) { 

    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;

      double * xptr = Xp[0];
      double y0 = 0.0;
      double y1 = 0.0;
      double y2 = 0.0;
      double y3 = 0.0;
      switch(RowDim) {
      case 1:
        for (int j=0; j < NumEntries; ++j) {
    int xoff = ColFirstPointInElementList[*BlockRowIndices++];
    
    y0 += A[0]*xptr[xoff];
    A += 1; 
        }
        break;
        
      case 2:
        for (int j=0; j < NumEntries; ++j) {
    int xoff = ColFirstPointInElementList[*BlockRowIndices++];
    y0 += A[0]*xptr[xoff] + A[2]*xptr[xoff+1];
    y1 += A[1]*xptr[xoff] + A[3]*xptr[xoff+1];
    A += 4 ;
        }
        break;
        
      case 3:
        for (int j=0; j < NumEntries; ++j) {
    int xoff = ColFirstPointInElementList[*BlockRowIndices++];
    y0 += A[0]*xptr[xoff+0] + A[3]*xptr[xoff+1] + A[6]*xptr[xoff+2];
    y1 += A[1]*xptr[xoff+0] + A[4]*xptr[xoff+1] + A[7]*xptr[xoff+2];
    y2 += A[2]*xptr[xoff+0] + A[5]*xptr[xoff+1] + A[8]*xptr[xoff+2];
    A += 9 ;
        }
        break;

      case 4:
        for (int j=0; j < NumEntries; ++j) {
    int xoff = ColFirstPointInElementList[*BlockRowIndices++];
    y0 += A[0]*xptr[xoff+0] + A[4]*xptr[xoff+1] + A[8]*xptr[xoff+2] +A[12]*xptr[xoff+3];
    y1 += A[1]*xptr[xoff+0] + A[5]*xptr[xoff+1] + A[9]*xptr[xoff+2] +A[13]*xptr[xoff+3];
    y2 += A[2]*xptr[xoff+0] + A[6]*xptr[xoff+1] + A[10]*xptr[xoff+2] +A[14]*xptr[xoff+3];
    y3 += A[3]*xptr[xoff+0] + A[7]*xptr[xoff+1] + A[11]*xptr[xoff+2] +A[15]*xptr[xoff+3];
    A += 16 ;
        }
        break;
      default:
        assert(false);
      }
      int  yoff = *RowFirstPointInElementList++;
      switch(RowDim) {
      case 4:
        *(Yp[0]+yoff+3) = y3;
      case 3:
        *(Yp[0]+yoff+2) = y2;
      case 2:
        *(Yp[0]+yoff+1) = y1;
      case 1:
        *(Yp[0]+yoff) = y0;
      }
    }
  } else { // NumVectors != 1 
    double *OrigA = A ; 
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int  yoff = *RowFirstPointInElementList++;
      int *    BRI = *Indices++;
    
      for (int k=0; k<NumVectors; k++) {
        int *    BlockRowIndices = BRI;
        A = OrigA ; 
        double * xptr = Xp[k];
        double y0 = 0.0;
        double y1 = 0.0;
        double y2 = 0.0;
        double y3 = 0.0;
        switch(RowDim) {
        case 1:
    for (int j=0; j < NumEntries; ++j) {
      int xoff = ColFirstPointInElementList[*BlockRowIndices++];
      
      y0 += A[0]*xptr[xoff];
      A += 1; 
    }
    break;
    
        case 2:
    for (int j=0; j < NumEntries; ++j) {
      int xoff = ColFirstPointInElementList[*BlockRowIndices++];
      y0 += A[0]*xptr[xoff] + A[2]*xptr[xoff+1];
      y1 += A[1]*xptr[xoff] + A[3]*xptr[xoff+1];
      A += 4 ;
    }
    break;
    
        case 3:
    for (int j=0; j < NumEntries; ++j) {
      int xoff = ColFirstPointInElementList[*BlockRowIndices++];
      y0 += A[0]*xptr[xoff+0] + A[3]*xptr[xoff+1] + A[6]*xptr[xoff+2];
      y1 += A[1]*xptr[xoff+0] + A[4]*xptr[xoff+1] + A[7]*xptr[xoff+2];
      y2 += A[2]*xptr[xoff+0] + A[5]*xptr[xoff+1] + A[8]*xptr[xoff+2];
      A += 9 ;
    }
    break;
      case 4:
        for (int j=0; j < NumEntries; ++j) {
    int xoff = ColFirstPointInElementList[*BlockRowIndices++];
    y0 += A[0]*xptr[xoff+0] + A[4]*xptr[xoff+1] + A[8]*xptr[xoff+2] +A[12]*xptr[xoff+3];
    y1 += A[1]*xptr[xoff+0] + A[5]*xptr[xoff+1] + A[9]*xptr[xoff+2] +A[13]*xptr[xoff+3];
    y2 += A[2]*xptr[xoff+0] + A[6]*xptr[xoff+1] + A[10]*xptr[xoff+2] +A[14]*xptr[xoff+3];
    y3 += A[3]*xptr[xoff+0] + A[7]*xptr[xoff+1] + A[11]*xptr[xoff+2] +A[15]*xptr[xoff+3];
    A += 16 ;
        }
        break;
        default:
    assert(false);
        }
        switch(RowDim) {
        case 4:
    *(Yp[k]+yoff+3) = y3;
        case 3:
    *(Yp[k]+yoff+2) = y2;
        case 2:
    *(Yp[k]+yoff+1) = y1;
        case 1:
    *(Yp[k]+yoff) = y0;
        }
      }
      OrigA = A ; 
    }
  } // end else NumVectors != 1 

      } else {  
  
  for (i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
    int  yoff = *RowFirstPointInElementList++;
    int RowDim = *RowElementSizeList++;
    
    FastBlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
             ColFirstPointInElementList, ColElementSizeList, 
             BlockRowValues, Xp, Yp, NumVectors);
  }
      }
    } else { 
      for (i=0; i < NumMyBlockRows_; i++) {
  int      NumEntries = *NumBlockEntriesPerRow++;
  int *    BlockRowIndices = *Indices++;
  Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
  int  yoff = *RowFirstPointInElementList++;
  int RowDim = *RowElementSizeList++;
  BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
       ColFirstPointInElementList, ColElementSizeList, 
       BlockRowValues, Xp, Yp, NumVectors);
      }
    }

    if (Exporter()!=0) {
      Yref.PutScalar(0.0);
      EPETRA_CHK_ERR(Yref.Export(*ExportVector_, *Exporter(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Yref.Reduce());
  }
  else { // Transpose operation
    
    
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
  if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
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
  if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      // Create import vector if needed
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors);

      ImportVector_->PutScalar(0.0); // Zero y values
      Yp = (double**)ImportVector_->Pointers();
      RowElementSizeList = ColMap().ElementSizeList();
      RowFirstPointInElementList = ColMap().FirstPointInElementList();
    }
    else
      Yref.PutScalar(0.0); // Zero y values

    // Do actual computation

    if ( All_Values_ != 0 ) { // Contiguous Storage 
      for (i=0; i < NumMyBlockRows_; i++) {
  int      NumEntries = *NumBlockEntriesPerRow++;
  int *    BlockRowIndices = *Indices++;
  Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
  int  xoff = *RowFirstPointInElementList++;
  int RowDim = *RowElementSizeList++;
  FastBlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, xoff, 
       ColFirstPointInElementList, ColElementSizeList, 
       BlockRowValues, Xp, Yp, NumVectors);
      }
    } else {
      for (i=0; i < NumMyBlockRows_; i++) {
  int      NumEntries = *NumBlockEntriesPerRow++;
  int *    BlockRowIndices = *Indices++;
  Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
  int  xoff = *RowFirstPointInElementList++;
  int RowDim = *RowElementSizeList++;
  BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, xoff, 
       ColFirstPointInElementList, ColElementSizeList, 
       BlockRowValues, Xp, Yp, NumVectors);
      }
    }

    if (Importer()!=0) {
      Yref.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Yref.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Yref.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());

  if (x_and_y_same == true) {
    Y = *Ytemp;
    delete Ytemp;
  }

  return(0);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowMultiply(bool TransA,
          int RowDim,
          int NumEntries, 
          int * BlockIndices,
          int RowOff,
          int * FirstPointInElementList,
          int * ElementSizeList,
          double Alpha,
          Epetra_SerialDenseMatrix** As,
          double ** X,
          double Beta,
          double ** Y,
          int NumVectors) const
{
  //This overloading of BlockRowMultiply is the same as the one below, except
  //that this one accepts Alpha and Beta arguments. This BlockRowMultiply is
  //called from within the 'solve' methods.
  int j, k;
  if (!TransA) {
    for (j=0; j < NumEntries; j++) {
      Epetra_SerialDenseMatrix* Asub = As[j];
      double * A = Asub->A();
      int LDA = Asub->LDA();
      int BlockIndex = BlockIndices[j];
      int xoff = FirstPointInElementList[BlockIndex];
      int ColDim = ElementSizeList[BlockIndex];

      for (k=0; k<NumVectors; k++) {
  double * curx = X[k] + xoff;
  double * cury = Y[k] + RowOff;

  GEMV('N', RowDim, ColDim, Alpha, A, LDA, curx, Beta, cury);
      }//end for (k
    }//end for (j
  }
  else { //TransA == true
    for (j=0; j < NumEntries; j++) {
      double * A = As[j]->A();
      int LDA = As[j]->LDA();
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
void Epetra_VbrMatrix::FastBlockRowMultiply(bool TransA,
  int RowDim,
  int NumEntries, 
  int * BlockRowIndices,
  int yoff,
  int * ColFirstPointInElementList,
  int * ColElementSizeList,
  Epetra_SerialDenseMatrix** BlockRowValues,
  double ** Xp,
  double ** Yp,
  int NumVectors) const
{
  int j, k;
  if (!TransA) {
    for (k=0; k<NumVectors; k++) {
      double * y = Yp[k] + yoff;
      double * xptr = Xp[k];

      Epetra_SerialDenseMatrix* Asub = BlockRowValues[0];
      double *OrigA = Asub->A_;
      int LDA = Asub->LDA_;
      int OrigBlockIndex = BlockRowIndices[0];
      int ColDim = ColElementSizeList[OrigBlockIndex];
 
      assert( RowDim == ColDim ) ; 
      assert( RowDim == LDA ) ; 
      
      switch(RowDim) {
#if 0
      case 1:
  double y0 = y[0];
  double y1 = y[0];
  double *A = OrigA ;
  for (j=0; j < NumEntries; ++j) {
    
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
    
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    y[0] += A[0]*x[0];
  }//end for (j
  break;
  
      case 2:
  double y0 = y[0];
  double y1 = y[0];
  double *A = OrigA ;
  for (j=0; j < NumEntries; ++j) {
    
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
    
    y0 += A[0]*xptr[xoff] + A[2]*xptr[xoff+1];
    y1 += A[1]*xptr[xoff] + A[3]*xptr[xoff+1];
    A += ColDim * LDA  ;
  }
  y[0] = y0;
  y[1] = y1;
  break;
  
      case 3:
  for (j=0; j < NumEntries; ++j) {
      
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
      
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    y[0] += A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
    y[1] += A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
    y[2] += A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
  }
  break;

      case 4:
  for (j=0; j < NumEntries; ++j) {
      
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
      
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    y[0] += A[0]*x[0] + A[4]*x[1] + A[8]*x[2] + A[12]*x[3];
    y[1] += A[1]*x[0] + A[5]*x[1] + A[9]*x[2] + A[13]*x[3];
    y[2] += A[2]*x[0] + A[6]*x[1] + A[10]*x[2] + A[14]*x[3];
    y[3] += A[3]*x[0] + A[7]*x[1] + A[11]*x[2] + A[15]*x[3];
  }
  break;
#endif
      case 5:
  for (j=0; j < NumEntries; ++j) {
      
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
      
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    y[0] += A[0]*x[0] + A[5]*x[1] + A[10]*x[2] + A[15]*x[3] + A[20]*x[4];
    y[1] += A[1]*x[0] + A[6]*x[1] + A[11]*x[2] + A[16]*x[3] + A[21]*x[4];
    y[2] += A[2]*x[0] + A[7]*x[1] + A[12]*x[2] + A[17]*x[3] + A[22]*x[4];
    y[3] += A[3]*x[0] + A[8]*x[1] + A[13]*x[2] + A[18]*x[3] + A[23]*x[4];
    y[4] += A[4]*x[0] + A[9]*x[1] + A[14]*x[2] + A[19]*x[3] + A[24]*x[4];
  }
  break;

      case 6:
  for (j=0; j < NumEntries; ++j) {
      
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
      
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    y[0] += A[0]*x[0] + A[6]*x[1] + A[12]*x[2] + A[18]*x[3] + A[24]*x[4]
      + A[30]*x[5];
    y[1] += A[1]*x[0] + A[7]*x[1] + A[13]*x[2] + A[19]*x[3] + A[25]*x[4]
      + A[31]*x[5];
    y[2] += A[2]*x[0] + A[8]*x[1] + A[14]*x[2] + A[20]*x[3] + A[26]*x[4]
      + A[32]*x[5];
    y[3] += A[3]*x[0] + A[9]*x[1] + A[15]*x[2] + A[21]*x[3] + A[27]*x[4]
      + A[33]*x[5];
    y[4] += A[4]*x[0] + A[10]*x[1] + A[16]*x[2] + A[22]*x[3] + A[28]*x[4]
      + A[34]*x[5];
    y[5] += A[5]*x[0] + A[11]*x[1] + A[17]*x[2] + A[23]*x[3] + A[29]*x[4]
      + A[35]*x[5];
  }
  break;

      default:
  for (j=0; j < NumEntries; ++j) {
      
    int BlockIndex = BlockRowIndices[j];
    int xoff = ColFirstPointInElementList[BlockIndex];
      
    double *A = OrigA + j * ColDim * LDA  ;
    double *x = xptr + xoff;
    GEMV('N', RowDim, ColDim, 1.0, A, LDA, x, 1.0, y);
  }
      }//end switch

    }//end for (k
  }
  else { //TransA == true
    for (j=0; j < NumEntries; j++) {
      double * A = BlockRowValues[j]->A();
      int LDA = BlockRowValues[j]->LDA();
      int BlockIndex = BlockRowIndices[j];
      int Yyoff = ColFirstPointInElementList[BlockIndex];
      int ColDim = ColElementSizeList[BlockIndex];
      for (k=0; k<NumVectors; k++) {
  double * x = Xp[k] + yoff;
  double * y = Yp[k] + Yyoff;
  GEMV('T', RowDim, ColDim, 1.0, A, LDA, x, 1.0, y);
      }
    }
  }

}
void Epetra_VbrMatrix::BlockRowMultiply(bool TransA,
          int RowDim,
          int NumEntries, 
          int * BlockIndices,
          int RowOff,
          int * FirstPointInElementList,
          int * ElementSizeList,
          Epetra_SerialDenseMatrix** As,
          double ** X,
          double ** Y,
          int NumVectors) const
{
  //This overloading of BlockRowMultiply is the same as the one above, except
  //that this one doesn't accept the Alpha and Beta arguments (it assumes that
  //they are both 1.0) and contains some inlined unrolled code for certain
  //cases (certain block-sizes) rather than calling GEMV every time. This
  //BlockRowMultiply is called from within the 'Multiply' methods.
  //Note: Scott Hutchinson's Aztec method 'dvbr_sparax_basic' was consulted in
  //the optimizing of this method.

  int j, k;
  if (!TransA) {
    for (k=0; k<NumVectors; k++) {
      double * y = Y[k] + RowOff;
      double * xptr = X[k];

      for (j=0; j < NumEntries; ++j) {
  Epetra_SerialDenseMatrix* Asub = As[j];
  double * A = Asub->A_;
  int LDA = Asub->LDA_;
  int BlockIndex = BlockIndices[j];
  int xoff = FirstPointInElementList[BlockIndex];
  int ColDim = ElementSizeList[BlockIndex];

  double * x = xptr + xoff;

  //Call GEMV if sub-block is non-square or if LDA != RowDim.
  if (LDA != RowDim || ColDim != RowDim) {
    GEMV('N', RowDim, ColDim, 1.0, A, LDA, x, 1.0, y);
    continue;
  }

  //It is a big performance win to use inlined, unrolled code for small
  //common block sizes rather than calling GEMV.

  switch(RowDim) {
  case 1:
    y[0] += A[0]*x[0];
    break;

  case 2:
    y[0] += A[0]*x[0] + A[2]*x[1];
    y[1] += A[1]*x[0] + A[3]*x[1];
    break;

  case 3:
    y[0] += A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
    y[1] += A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
    y[2] += A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
    break;

  case 4:
    y[0] += A[0]*x[0] + A[4]*x[1] + A[8]*x[2] + A[12]*x[3];
    y[1] += A[1]*x[0] + A[5]*x[1] + A[9]*x[2] + A[13]*x[3];
    y[2] += A[2]*x[0] + A[6]*x[1] + A[10]*x[2] + A[14]*x[3];
    y[3] += A[3]*x[0] + A[7]*x[1] + A[11]*x[2] + A[15]*x[3];
    break;

  case 5:
    y[0] += A[0]*x[0] + A[5]*x[1] + A[10]*x[2] + A[15]*x[3] + A[20]*x[4];
    y[1] += A[1]*x[0] + A[6]*x[1] + A[11]*x[2] + A[16]*x[3] + A[21]*x[4];
    y[2] += A[2]*x[0] + A[7]*x[1] + A[12]*x[2] + A[17]*x[3] + A[22]*x[4];
    y[3] += A[3]*x[0] + A[8]*x[1] + A[13]*x[2] + A[18]*x[3] + A[23]*x[4];
    y[4] += A[4]*x[0] + A[9]*x[1] + A[14]*x[2] + A[19]*x[3] + A[24]*x[4];
    break;

  case 6:
    y[0] += A[0]*x[0] + A[6]*x[1] + A[12]*x[2] + A[18]*x[3] + A[24]*x[4]
      + A[30]*x[5];
    y[1] += A[1]*x[0] + A[7]*x[1] + A[13]*x[2] + A[19]*x[3] + A[25]*x[4]
      + A[31]*x[5];
    y[2] += A[2]*x[0] + A[8]*x[1] + A[14]*x[2] + A[20]*x[3] + A[26]*x[4]
      + A[32]*x[5];
    y[3] += A[3]*x[0] + A[9]*x[1] + A[15]*x[2] + A[21]*x[3] + A[27]*x[4]
      + A[33]*x[5];
    y[4] += A[4]*x[0] + A[10]*x[1] + A[16]*x[2] + A[22]*x[3] + A[28]*x[4]
      + A[34]*x[5];
    y[5] += A[5]*x[0] + A[11]*x[1] + A[17]*x[2] + A[23]*x[3] + A[29]*x[4]
      + A[35]*x[5];
    break;

  default:
    GEMV('N', RowDim, ColDim, 1.0, A, LDA, x, 1.0, y);
  }//end switch
      }//end for (k
    }//end for (j
  }
  else { //TransA == true
    for (j=0; j < NumEntries; j++) {
      double * A = As[j]->A();
      int LDA = As[j]->LDA();
      int BlockIndex = BlockIndices[j];
      int yoff = FirstPointInElementList[BlockIndex];
      int ColDim = ElementSizeList[BlockIndex];
      for (k=0; k<NumVectors; k++) {
  double * x = X[k] + RowOff;
  double * y = Y[k] + yoff;
  GEMV('T', RowDim, ColDim, 1.0, A, LDA, x, 1.0, y);
      }
    }
  }

  return;
}
//=============================================================================
int Epetra_VbrMatrix::DoSolve(bool Upper, bool TransA,
          bool UnitDiagonal,
          const Epetra_MultiVector& X,
          Epetra_MultiVector& Y) const
{
  (void)UnitDiagonal;
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
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;

  int * ColElementSizeList = ElementSizeList_;
  int * ColFirstPointInElementList = FirstPointInElementList_;

  // If upper, point to last row
  if (Upper) {
    NumBlockEntriesPerRow += NumMyBlockRows_-1;
    FirstPointInElement += NumMyBlockRows_-1;
    ElementSize += NumMyBlockRows_-1;
    Indices += NumMyBlockRows_-1;
    Entries += NumMyBlockRows_-1;
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
      Epetra_SerialDenseMatrix** BlockRowValues  = *Entries--;
      int  yoff = *FirstPointInElement--;
      int RowDim = *ElementSize--;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
           ColFirstPointInElementList, ColElementSizeList, 
           1.0, BlockRowValues, Yp, -1.0, Yp, NumVectors);
    }
  }
  else {
    for (i=0; i < NumMyBlockRows_; i++) {
      int      NumEntries = *NumBlockEntriesPerRow++;
      int *    BlockRowIndices = *Indices++;
      Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
      int  yoff = *FirstPointInElement++;
      int RowDim = *ElementSize++;
      BlockRowMultiply(TransA, RowDim, NumEntries, BlockRowIndices, yoff, 
           ColFirstPointInElementList, ColElementSizeList, 
           1.0, BlockRowValues, Yp, -1.0, Yp, NumVectors);
    }
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
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
  bool hasOperatorMap = false;
  if (DoRows) {
    if ( !Graph().RangeMap().SameAs(x.Map()) ) {
      hasOperatorMap = OperatorRangeMap().SameAs(x.Map());
      if( !hasOperatorMap )
        EPETRA_CHK_ERR(-2);
    }
  }
  else {
    if ( !Graph().DomainMap().SameAs(x.Map()) ) {
      hasOperatorMap = OperatorDomainMap().SameAs(x.Map());
      if( !hasOperatorMap )
        EPETRA_CHK_ERR(-2);
    }
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;
  
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
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
    int xoff = *RowFirstPointInElementList++;
    int RowDim = *RowElementSizeList++;
    if (DoRows) {
      for (int ii=0; ii < NumEntries; ii++) {
  double * xptr = xp+xoff;
  double * A = BlockRowValues[ii]->A();
  int LDA = BlockRowValues[ii]->LDA();
  int BlockIndex = BlockRowIndices[ii];
  int ColDim = ColElementSizeList[BlockIndex];
  for (int j=0; j<ColDim; j++) {
    double * curEntry = A + j*LDA;
    for (int k=0; k<RowDim; k++)
      xptr[k] += std::abs(*curEntry++);
  }
      }
    }
    else {
      for (int ii=0; ii < NumEntries; ii++) {
  double * A = BlockRowValues[ii]->A();
  int LDA = BlockRowValues[ii]->LDA();
  int BlockIndex = BlockRowIndices[ii];
  int off = ColFirstPointInElementList[BlockIndex];
  int ColDim = ColElementSizeList[BlockIndex];
  double * curx = xp+off;
  for (int j=0; j<ColDim; j++) {
    double * curEntry = A + j*LDA;
    for (int k=0; k<RowDim; k++)
      curx[j] += std::abs(*curEntry++);
  }
      }
    }
  }

  if (!DoRows) {
    if (Importer()!=0){
      Epetra_Vector  *x_blocked = 0;
      if(hasOperatorMap)
        x_blocked = new Epetra_Vector( ::View, DomainMap(), &x[0] );
      else
        x_blocked = &x;
      x_blocked->PutScalar(0.0);
      EPETRA_CHK_ERR(x_blocked->Export(*x_tmp, *Importer(), Add)); // Fill x with Values from import vector
      if(hasOperatorMap)
        delete x_blocked;
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
  UpdateFlops(NumGlobalNonzeros64());

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
  bool hasOperatorMap = false;
  if (DoRows) {
    if ( !Graph().RangeMap().SameAs(x.Map()) ) {
      hasOperatorMap = OperatorRangeMap().SameAs(x.Map());
      if( !hasOperatorMap )
        EPETRA_CHK_ERR(-2);
    }
  }
  else {
    if ( !Graph().DomainMap().SameAs(x.Map()) ) {
      hasOperatorMap = OperatorDomainMap().SameAs(x.Map());
      if( !hasOperatorMap )
        EPETRA_CHK_ERR(-2);
    }
  }
  int ierr = 0;
  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_; 
  
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
      Epetra_Vector  *x_blocked = 0;
      if(hasOperatorMap)
        x_blocked = new Epetra_Vector( ::View, Graph().DomainMap(), (double *) &x[0] );
      else
        x_blocked = (Epetra_Vector *) &x;
      x_tmp = new Epetra_Vector(ColMap()); // Create import vector if needed
      EPETRA_CHK_ERR(x_tmp->Import(*x_blocked,*Importer(), Insert)); // x_tmp will have all the values we need
      xp = (double*)x_tmp->Values();
    }
  }

  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++;
    int *    BlockRowIndices = *Indices++;
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
    int xoff = *RowFirstPointInElementList++;
    int RowDim = *RowElementSizeList++;
    if (DoRows) {
      for (int ii=0; ii < NumEntries; ii++) {
  double * xptr = xp+xoff;
  double * A = BlockRowValues[ii]->A();
  int LDA = BlockRowValues[ii]->LDA();
  int BlockIndex = BlockRowIndices[ii];
  int ColDim = ColElementSizeList[BlockIndex];
  for (int j=0; j<ColDim; j++) {
    double * curEntry = A + j*LDA;
    for (int k=0; k<RowDim; k++)
      *curEntry++ *= xptr[k];
  }
      }
    }
    else {
      for (int ii=0; ii < NumEntries; ii++) {
  double * A = BlockRowValues[ii]->A();
  int LDA = BlockRowValues[ii]->LDA();
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
  NormFrob_ = -1.0;
  UpdateFlops(NumGlobalNonzeros64());

  EPETRA_CHK_ERR(ierr);
  return(0);
}
//=============================================================================
double Epetra_VbrMatrix::NormInf() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormInf is called.
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormInf_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
  //  if( NormInf_ >= 0.0) 
    return(NormInf_);
#endif

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int MaxRowDim_ = MaxRowDim();
  double * tempv = new double[MaxRowDim_];

  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * ElementSize = ElementSizeList_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;

  double Local_NormInf = 0.0;
  for (int i=0; i < NumMyBlockRows_; i++) {
    int      NumEntries = *NumBlockEntriesPerRow++ ;
    int RowDim = *ElementSize++;
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
    BlockRowNormInf(RowDim, NumEntries, 
        BlockRowValues, tempv);
    for (int j=0; j < RowDim; j++) Local_NormInf = EPETRA_MAX(Local_NormInf, tempv[j]);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  delete [] tempv;
  UpdateFlops(NumGlobalNonzeros64());
  return(NormInf_);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowNormInf(int RowDim, int NumEntries,
               Epetra_SerialDenseMatrix** As, 
               double * Y) const
{
  int i, j, k;
  for (k=0; k<RowDim; k++) Y[k] = 0.0;

  for (i=0; i < NumEntries; i++) {
    double * A = As[i]->A();
    int LDA = As[i]->LDA();
    int ColDim = As[i]->N();
    for (j=0; j<ColDim; j++) {
      for (k=0; k<RowDim; k++) Y[k] += std::abs(A[k]);
      A += LDA;
    }
  }
  return;
}
//=============================================================================
double Epetra_VbrMatrix::NormOne() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormOne is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormOne_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormOne_);
#endif

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
  int ** Indices = Indices_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;

  for (int i=0; i < NumMyBlockRows_; i++) {
    int NumEntries = *NumBlockEntriesPerRow++;
    int RowDim = *ElementSize++;
    int *    BlockRowIndices = *Indices++;
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;
    BlockRowNormOne(RowDim, NumEntries, BlockRowIndices,
        BlockRowValues,  ColFirstPointInElementList, xp);
  }
  if (Importer()!=0) {
    x->PutScalar(0.0);
    EPETRA_CHK_ERR(x->Export(*x_tmp, *Importer(), Add));
  } // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros64());
  return(NormOne_);
}
//=============================================================================
double Epetra_VbrMatrix::NormFrobenius() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormOne is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormFrob_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormFrob_);
#endif

  int * NumBlockEntriesPerRow = NumBlockEntriesPerRow_;
  int * ElementSize = ElementSizeList_;
  Epetra_SerialDenseMatrix*** Entries = Entries_;

  double local_sum = 0.0;

  for (int i=0; i < NumMyBlockRows_; i++) {
    int NumEntries = *NumBlockEntriesPerRow++;
    int RowDim = *ElementSize++;
    Epetra_SerialDenseMatrix** BlockRowValues  = *Entries++;

    for(int ii=0; ii<NumEntries; ++ii) {
      double* A = BlockRowValues[ii]->A();
      int LDA = BlockRowValues[ii]->LDA();
      int colDim = BlockRowValues[ii]->N();
      for(int j=0; j<colDim; ++j) {
        for(int k=0; k<RowDim; ++k) {
          local_sum += A[k]*A[k];
        }
        A += LDA;
      }
    }

//The following commented-out code performs the calculation by running
//all the way across each point-entry row before going to the next.
//This is the order in which the CrsMatrix frobenius norm is calculated,
//but for a VbrMatrix it is faster to run the block-entries in
//column-major order (which is how they're stored), as the above code
//does.
//But if the same matrix is stored in both Vbr and Crs form, and you wish
//to compare the frobenius norms, it is necessary to run through the
//coefficients in the same order to avoid round-off differences.
//
//    for(int k=0; k<RowDim; ++k) {
//      for(int ii=0; ii<NumEntries; ++ii) {
//        double* A = BlockRowValues[ii]->A();
//        int LDA = BlockRowValues[ii]->LDA();
//        int colDim = BlockRowValues[ii]->N();
//        for(int j=0; j<colDim; ++j) {
//          double Ak = A[k+j*LDA];
//          local_sum += Ak*Ak;
//        }
//      }
//    }

  }

  double global_sum = 0.0;
  Comm().SumAll(&local_sum, &global_sum, 1);

  NormFrob_ = std::sqrt(global_sum);

  UpdateFlops(NumGlobalNonzeros64());

  return(NormFrob_);
}
//=============================================================================
void Epetra_VbrMatrix::BlockRowNormOne(int RowDim, int NumEntries, int * BlockRowIndices,
               Epetra_SerialDenseMatrix** As, 
               int * ColFirstPointInElementList, double * x) const {
  int i, j, k;

  for (i=0; i < NumEntries; i++) {
    double * A = As[i]->A();
    int LDA = As[i]->LDA();
    int ColDim = As[i]->N();
    double * curx = x + ColFirstPointInElementList[BlockRowIndices[i]];
    for (j=0; j<ColDim; j++) {
      for (k=0; k<RowDim; k++) curx[j] += std::abs(A[k]);
      A += LDA;
    }
  }
  return;
}
//=========================================================================
int Epetra_VbrMatrix::CheckSizes(const Epetra_SrcDistObject & Source) {
  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  if (!A.Graph().GlobalConstantsComputed()) EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::CopyAndPermute(const Epetra_SrcDistObject & Source,
                                     int NumSameIDs, 
                                     int NumPermuteIDs,
                                     int * PermuteToLIDs,
                                     int *PermuteFromLIDs,
                                     const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor; 
  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  int i, j;
  
  int BlockRow, NumBlockEntries;
  int * BlockIndices;
  int RowDim;
  Epetra_SerialDenseMatrix ** Entries;
  int FromBlockRow, ToBlockRow;
  
  // Do copy first
  if (NumSameIDs>0) {
    int maxNumBlockEntries = A.MaxNumBlockEntries();
    BlockIndices = new int[maxNumBlockEntries];  // Need some temporary space
 

    for (i=0; i<NumSameIDs; i++) {
      BlockRow = GRID64(i);// FIXME long long
      EPETRA_CHK_ERR( A.ExtractGlobalBlockRowPointers(BlockRow,
                                                      maxNumBlockEntries,
                  RowDim, NumBlockEntries, 
                  BlockIndices, Entries)); // Set pointers
      // Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
      if (StaticGraph() || IndicesAreLocal()) {
  EPETRA_CHK_ERR(BeginReplaceGlobalValues(BlockRow, NumBlockEntries,
                                                BlockIndices));
      }
      else {
  EPETRA_CHK_ERR(BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices));
      }
      // Insert block entries one-at-a-time
      for (j=0; j<NumBlockEntries; j++) SubmitBlockEntry(Entries[j]->A(),
               Entries[j]->LDA(),
               RowDim, Entries[j]->N());
      EndSubmitEntries(); // Complete this block row
    }
    delete [] BlockIndices;
  }

  // Do local permutation next
  if (NumPermuteIDs>0) {
    int maxNumBlockEntries = A.MaxNumBlockEntries();
    BlockIndices = new int[maxNumBlockEntries];  // Need some temporary space
      
    for (i=0; i<NumPermuteIDs; i++) {
      FromBlockRow = A.GRID64(PermuteFromLIDs[i]);// FIXME long long
      ToBlockRow = GRID64(PermuteToLIDs[i]);// FIXME long long
      EPETRA_CHK_ERR(A.ExtractGlobalBlockRowPointers(FromBlockRow, maxNumBlockEntries, RowDim, NumBlockEntries, 
               BlockIndices, Entries)); // Set pointers
      // Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
      if (StaticGraph() || IndicesAreLocal()) {
  EPETRA_CHK_ERR(BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
      else {
  EPETRA_CHK_ERR(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
      // Insert block entries one-at-a-time
      for (j=0; j<NumBlockEntries; j++) SubmitBlockEntry(Entries[j]->A(),
               Entries[j]->LDA(), RowDim, Entries[j]->N());
      EndSubmitEntries(); // Complete this block row
    }
    delete [] BlockIndices;
  }
 
  return(0);
}

//=========================================================================
int Epetra_VbrMatrix::PackAndPrepare(const Epetra_SrcDistObject & Source,
                                     int NumExportIDs,
                                     int * ExportLIDs,
                                     int & LenExports,
                                     char * & Exports,
                                     int & SizeOfPacket,
                                     int * Sizes,
                                     bool & VarSizes,
                                     Epetra_Distributor & Distor)
{
  (void)LenExports;
  (void)Sizes;
  (void)VarSizes;
  (void)Distor; 
  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);

  double * DoubleExports = 0;
  int globalMaxNumNonzeros = A.GlobalMaxNumNonzeros();
  int globalMaxNumBlockEntries = A.GlobalMaxNumBlockEntries();
  // Will have globalMaxNumEntries doubles, globalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = globalMaxNumNonzeros +  
    (((2*globalMaxNumBlockEntries+3)+(int)sizeof(int)-1)*(int)sizeof(int))/(int)sizeof(double);
  SizeOfPacket = DoublePacketSize * (int)sizeof(double); 

  if (DoublePacketSize*NumExportIDs>LenExports_) {
    if (LenExports_>0) delete [] Exports_;
    LenExports_ = DoublePacketSize*NumExportIDs;
    DoubleExports = new double[LenExports_];
    Exports_ = (char *) DoubleExports;
  }

  if (NumExportIDs<=0) return(0); // All done if nothing to pack

  int i, j;
  
  int NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * Entries;
  int FromBlockRow;
  double * valptr, * dintptr;
  int * intptr;
  
  // Each segment of IntExports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  RowDim of Block Row
  // next int: NumBlockEntries, Number of indices in row.
  // next NumBlockEntries: The actual indices for the row.
  // Any remaining space (of length globalMaxNumNonzeros - NumBlockEntries ints) will be wasted but we need fixed
  //   sized segments for current communication routines.

  // Each segment of Exports will be filled with values.

  valptr = (double *) Exports;
  dintptr = valptr + globalMaxNumNonzeros;
  intptr = (int *) dintptr;
  bool NoSumInto = false;
  for (i=0; i<NumExportIDs; i++) {
    FromBlockRow = A.GRID64(ExportLIDs[i]);// FIXME long long
    BlockIndices = intptr + 3;
    ColDims = BlockIndices + globalMaxNumBlockEntries;
    EPETRA_CHK_ERR(A.BeginExtractGlobalBlockRowCopy(FromBlockRow, globalMaxNumBlockEntries, RowDim,
              NumBlockEntries, BlockIndices, ColDims));
    // Now extract each block entry into send buffer
    Entries = valptr;
    for (j=0; j<NumBlockEntries; j++) {
      int SizeOfValues = RowDim*ColDims[j];
      A.ExtractEntryCopy(SizeOfValues, Entries, RowDim, NoSumInto);
      Entries += SizeOfValues;
    }
    // Fill first three slots of intptr with info
    intptr[0] = FromBlockRow;
    intptr[1] = RowDim;
    intptr[2] = NumBlockEntries;
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + globalMaxNumNonzeros;
    intptr = (int *) dintptr;
  }
    
  return(0);
}

//=========================================================================
int Epetra_VbrMatrix::UnpackAndCombine(const Epetra_SrcDistObject & Source, 
                                       int NumImportIDs,
                                       int * ImportLIDs, 
                                       int LenImports,
                                       char * Imports,
                                       int & SizeOfPacket, 
                                       Epetra_Distributor & Distor, 
                                       Epetra_CombineMode CombineMode,
                                       const Epetra_OffsetIndex * Indexor)
{
  (void)LenImports;
  (void)SizeOfPacket;
  (void)Distor;
  (void)Indexor;
  if (NumImportIDs<=0) return(0);

  if(   CombineMode != Add
  && CombineMode != Zero
  && CombineMode != Insert )
    EPETRA_CHK_ERR(-1); // CombineMode not supported, default to mode Zero

  const Epetra_VbrMatrix & A = dynamic_cast<const Epetra_VbrMatrix &>(Source);
  int NumBlockEntries;
  int * BlockIndices;
  int RowDim, * ColDims;
  double * values;
  int ToBlockRow;
  int i, j;
  
  double * valptr, *dintptr;
  int * intptr;
  int globalMaxNumNonzeros = A.GlobalMaxNumNonzeros();
  int globalMaxNumBlockEntries = A.GlobalMaxNumBlockEntries();
  // Will have globalMaxNumEntries doubles, globalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = globalMaxNumNonzeros +  
    (((2*globalMaxNumBlockEntries+3)+(int)sizeof(int)-1)*(int)sizeof(int))/(int)sizeof(double);
  // Unpack it...


  // Each segment of IntImports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumBlockEntries, Number of indices in row.
  // next int:  RowDim of Block Row
  // next NumBlockEntries: The actual indices for the row.
  // Any remaining space (of length globalMaxNumNonzeros - NumBlockEntries ints) will be 
  //  wasted but we need fixed sized segments for current communication routines.

  valptr = (double *) Imports;
  dintptr = valptr + globalMaxNumNonzeros;
  intptr = (int *) dintptr;
    
  for (i=0; i<NumImportIDs; i++) {
    ToBlockRow = GRID64(ImportLIDs[i]);// FIXME long long
    assert((intptr[0])==ToBlockRow); // Sanity check
    RowDim = RowMap().ElementSize(ImportLIDs[i]);
    assert((intptr[1])==RowDim); // Sanity check
    NumBlockEntries = intptr[2];
    BlockIndices = intptr + 3; 
    ColDims = BlockIndices + globalMaxNumBlockEntries;
    if (CombineMode==Add) {
      if (StaticGraph() || IndicesAreLocal()) {
  // Replace any current values
  EPETRA_CHK_ERR(BeginSumIntoGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
      else {
  // Insert values
  EPETRA_CHK_ERR(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
    }
    else if (CombineMode==Insert) {
      if (StaticGraph() || IndicesAreLocal()) {
  // Replace any current values
  EPETRA_CHK_ERR(BeginReplaceGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
      else {
  // Insert values
  EPETRA_CHK_ERR(BeginInsertGlobalValues(ToBlockRow, NumBlockEntries, BlockIndices));
      }
    }
    // Now extract each block entry into send buffer
    values = valptr;
    for (j=0; j<NumBlockEntries; j++) {
      int LDA = RowDim;
      int ColDim = ColDims[j];
      SubmitBlockEntry(values, LDA, RowDim, ColDim);
      values += (LDA*ColDim);
    }
    EndSubmitEntries(); // Done with this block row
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + globalMaxNumNonzeros;
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
int Epetra_VbrMatrix::BlockMap2PointMap(const Epetra_BlockMap & BlockMap,
          Epetra_Map * & PointMap) const
{
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
    int StartID = BlockMap.GID64(i)*MaxElementSize;
    int ElementSize = BlockMap.ElementSize(i);
    for (int j=0; j<ElementSize; j++) PtMyGlobalElements[curID++] = StartID+j;
  }
  assert(curID==PtNumMyElements); // Sanity test

  PointMap = new Epetra_Map(-1, PtNumMyElements, PtMyGlobalElements, BlockMap.IndexBase(), BlockMap.Comm());// FIXME long long

  if (PtNumMyElements>0) delete [] PtMyGlobalElements;

  if (!BlockMap.PointSameAs(*PointMap)) {EPETRA_CHK_ERR(-1);} // Maps not compatible
  return(0);
}

//=========================================================================
int Epetra_VbrMatrix::UpdateOperatorXY(const Epetra_MultiVector& X,
               const Epetra_MultiVector& Y) const
{
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
int Epetra_VbrMatrix::Apply(const Epetra_MultiVector& X,
          Epetra_MultiVector& Y) const
{
  if (!Epetra_VbrMatrix::UseTranspose()) {
    EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Epetra_VbrMatrix::DoMultiply(Epetra_VbrMatrix::UseTranspose(), *OperatorX_, *OperatorY_));
  } 
  else { // Swap role of OperatorX_ and OperatorY_ to remain compatible with domain and range spaces.
    EPETRA_CHK_ERR(UpdateOperatorXY(Y,X)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Epetra_VbrMatrix::DoMultiply(Epetra_VbrMatrix::UseTranspose(), *OperatorY_, *OperatorX_));
  }
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::ApplyInverse(const Epetra_MultiVector& X,
           Epetra_MultiVector& Y) const
{
  if (!Epetra_VbrMatrix::UseTranspose()) {
    EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(DoSolve(UpperTriangular(), Epetra_VbrMatrix::UseTranspose(), NoDiagonal(), *OperatorX_, *OperatorY_));
  }
  else { // Swap role of OperatorX_ and OperatorY_ to remain compatible with domain and range spaces.
    EPETRA_CHK_ERR(UpdateOperatorXY(Y,X)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(DoSolve(UpperTriangular(), Epetra_VbrMatrix::UseTranspose(), NoDiagonal(), *OperatorY_, *OperatorX_));
  }
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::Multiply(bool TransA, const Epetra_MultiVector& X,
          Epetra_MultiVector& Y) const
{
  if (!TransA) {
    EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Epetra_VbrMatrix::DoMultiply(TransA, *OperatorX_, *OperatorY_));
  } 
  else { // Swap role of OperatorX_ and OperatorY_ to remain compatible with domain and range spaces.
    EPETRA_CHK_ERR(UpdateOperatorXY(Y,X)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(Epetra_VbrMatrix::DoMultiply(TransA, *OperatorY_, *OperatorX_));
  }
  return(0);
}
//=========================================================================
int Epetra_VbrMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,
           Epetra_MultiVector& Y) const
{
  if (!Trans) {
    EPETRA_CHK_ERR(UpdateOperatorXY(X,Y)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(DoSolve(Upper, Trans, UnitDiagonal, *OperatorX_, *OperatorY_));
  }
  else { // Swap role of OperatorX_ and OperatorY_ to remain compatible with domain and range spaces.
    EPETRA_CHK_ERR(UpdateOperatorXY(Y,X)); // Update X and Y vector whose maps are compatible with the Vbr Matrix
    EPETRA_CHK_ERR(DoSolve(Upper, Trans, UnitDiagonal, *OperatorY_, *OperatorX_));
  }
  return(0);
}
//=========================================================================
void Epetra_VbrMatrix::Print(ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (MyPID==0) {
  os <<  "\nNumber of Global Block Rows  = "; os << NumGlobalBlockRows64(); os << endl;
  os <<    "Number of Global Block Cols  = "; os << NumGlobalBlockCols64(); os << endl;
  os <<    "Number of Global Block Diags = "; os << NumGlobalBlockDiagonals64(); os << endl;
  os <<    "Number of Global Blk Entries = "; os << NumGlobalBlockEntries64(); os << endl;
  os <<    "Global Max Num Block Entries = "; os << GlobalMaxNumBlockEntries(); os << endl;
  os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows64(); os << endl;
  os <<    "Number of Global Cols        = "; os << NumGlobalCols64(); os << endl;
  os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals64(); os << endl;
  os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros64(); os << endl;
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
      Epetra_SerialDenseMatrix ** Entries1;
      int RowDim1, NumBlockEntries1;
      int i, j;

      if (MyPID==0) {
  os.width(8);
  os <<  "   Processor ";
  os.width(10);
  os <<  "   Block Row Index ";
  os.width(10);
  os <<  "   Block Col Index \n";
  os.width(20);
  os <<  "   values     ";
  os << endl;
      }
      for (i=0; i<NumBlockRows1; i++) {
  int BlockRow1 = GRID64(i); // Get global row number// FIXME long long
  ExtractGlobalBlockRowPointers(BlockRow1, MaxNumBlockEntries1, RowDim1, 
              NumBlockEntries1, BlockIndices1,
              Entries1);
  
  for (j = 0; j < NumBlockEntries1 ; j++) {   
    os.width(8);
    os <<  MyPID ; os << "    ";  
    os.width(10);
    os <<  BlockRow1 ; os << "    ";  
    os.width(10);
    os <<  BlockIndices1[j]; os << "    " << endl;
    os.width(20);

          if (Entries1[j] == 0) {
            os << "Block Entry == NULL"<<endl;
            continue;
          }

    Epetra_SerialDenseMatrix entry_view(View, Entries1[j]->A(), Entries1[j]->LDA(),
           RowDim1, Entries1[j]->N());
          os << entry_view; os << "    ";
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

#endif // EPETRA_NO_32BIT_GLOBAL_INDICES
