//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector<int> & RowStencil,
        int rowIndex,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, vector< vector<int> >(1,RowStencil), vector<int>(1,rowIndex), GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_int_( vector< vector<int> >(1,RowStencil) ),
    RowIndices_int_( vector<int>(1,rowIndex) ),
    ROffset_(BlockUtility::CalculateOffset64(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseGraph.ColMap()))
{
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector<long long> & RowStencil,
        long long rowIndex,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, vector< vector<long long> >(1,RowStencil), vector<long long>(1,rowIndex), GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_LL_( vector< vector<long long> >(1,RowStencil) ),
    RowIndices_LL_( vector<long long>(1,rowIndex) ),
    ROffset_(BlockUtility::CalculateOffset64(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseGraph.ColMap()))
{
}
#endif

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_int_( RowStencil ),
    RowIndices_int_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset64(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseGraph.ColMap()))
{
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<long long> > & RowStencil,
        const vector<long long> & RowIndices,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_LL_( RowStencil ),
    RowIndices_LL_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset64(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseGraph.ColMap()))
{
}
#endif

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, LocalBlockGraph, GlobalComm )) ),
    BaseGraph_( BaseGraph ),
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    RowStencil_int_( ),
    RowIndices_int_( ),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    RowStencil_LL_( ),
    RowIndices_LL_( ),
#endif
    ROffset_(BlockUtility::CalculateOffset64(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseGraph.ColMap()))
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt() && LocalBlockGraph.RowMap().GlobalIndicesInt())
    BlockUtility::GenerateRowStencil(LocalBlockGraph, RowIndices_int_, RowStencil_int_);
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && LocalBlockGraph.RowMap().GlobalIndicesLongLong())
    BlockUtility::GenerateRowStencil(LocalBlockGraph, RowIndices_LL_, RowStencil_LL_);
  else
#endif
    throw "EpetraExt::BlockCrsMatrix::BlockCrsMatrix: Error, Global indices unknown.";
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseMatrix, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( Copy, BaseMatrix.RowMatrixRowMap(), 1 ), //Junk to satisfy constructor
    RowStencil_int_( RowStencil ),
    RowIndices_int_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset64(BaseMatrix.RowMatrixRowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseMatrix.RowMatrixColMap()))
{
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<long long> > & RowStencil,
        const vector<long long> & RowIndices,
        const Epetra_Comm & GlobalComm  )
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseMatrix, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( Copy, BaseMatrix.RowMatrixRowMap(), 1 ), //Junk to satisfy constructor
    RowStencil_LL_( RowStencil ),
    RowIndices_LL_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset64(BaseMatrix.RowMatrixRowMap())),
    COffset_(BlockUtility::CalculateOffset64(BaseMatrix.RowMatrixColMap()))
{
}
#endif

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix( const BlockCrsMatrix & Matrix )
  : Epetra_CrsMatrix( dynamic_cast<const Epetra_CrsMatrix &>( Matrix ) ),
    BaseGraph_( Matrix.BaseGraph_ ),
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    RowStencil_int_( Matrix.RowStencil_int_ ),
    RowIndices_int_( Matrix.RowIndices_int_ ),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    RowStencil_LL_( Matrix.RowStencil_LL_ ),
    RowIndices_LL_( Matrix.RowIndices_LL_ ),
#endif
    ROffset_( Matrix.ROffset_ ),
    COffset_( Matrix.COffset_ )
{
}

//==============================================================================
BlockCrsMatrix::~BlockCrsMatrix()
{
}

//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TLoadBlock(const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col)
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine copies entries of a BaseMatrix into big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // this->val[i+rowOffset][j+ColOffset] = BaseMatrix.val[i][j]

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices_local(MaxIndices);
  vector<int_type> Indices_global(MaxIndices);
  vector<double> vals(MaxIndices);
  int NumIndices;
  int ierr=0;

  for (int i=0; i<BaseMap.NumMyElements(); i++) {
    BaseMatrix.ExtractMyRowCopy( i, MaxIndices, NumIndices, &vals[0], &Indices_local[0] );

    // Convert to BlockMatrix Global numbering scheme
    for( int l = 0; l < NumIndices; ++l )
       Indices_global[l] = ColOffset +  (int_type) BaseColMap.GID64(Indices_local[l]);

    int_type BaseRow = (int_type) BaseMap.GID64(i);
    ierr = this->ReplaceGlobalValues(BaseRow + RowOffset, NumIndices, &vals[0], &Indices_global[0]);
    if (ierr != 0) std::cout << "WARNING BlockCrsMatrix::LoadBlock ReplaceGlobalValues err = " << ierr <<
            "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices_global[0] << std::endl;

  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::LoadBlock(const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt() && BaseMatrix.RowMatrixRowMap().GlobalIndicesInt())
        return TLoadBlock<int>(BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::LoadBlock: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::LoadBlock(const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && BaseMatrix.RowMatrixRowMap().GlobalIndicesLongLong())
        return TLoadBlock<long long>(BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::LoadBlock: Global indices not long long";
}
#endif

//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TSumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col)
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine copies entries of a BaseMatrix into big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // this->val[i+rowOffset][j+ColOffset] = BaseMatrix.val[i][j]

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices_local(MaxIndices);
  vector<int_type> Indices_global(MaxIndices);
  vector<double> vals(MaxIndices);
  int NumIndices;
  int ierr=0;

  for (int i=0; i<BaseMap.NumMyElements(); i++) {
    BaseMatrix.ExtractMyRowCopy( i, MaxIndices, NumIndices, &vals[0], &Indices_local[0] );

    // Convert to BlockMatrix Global numbering scheme
    for( int l = 0; l < NumIndices; ++l ) {
       Indices_global[l] = ColOffset +  (int_type) BaseColMap.GID64(Indices_local[l]);
       vals[l] *= alpha;
    }

    int_type BaseRow = (int_type) BaseMap.GID64(i);
    ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices, &vals[0], &Indices_global[0]);
    if (ierr != 0) {
      std::cout << "WARNING BlockCrsMatrix::SumIntoBlock SumIntoGlobalValues "
        "err = " << ierr << std::endl << "\t  Row " << BaseRow + RowOffset <<
        "Col start" << Indices_global[0] << std::endl;
    }
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt() && BaseMatrix.RowMatrixRowMap().GlobalIndicesInt())
        return TSumIntoBlock<int>(alpha, BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::SumIntoBlock: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && BaseMatrix.RowMatrixRowMap().GlobalIndicesLongLong())
        return TSumIntoBlock<long long>(alpha, BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::SumIntoBlock: Global indices not long long";
}
#endif

//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TSumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col)
{
  int_type RowOffset = Row * ROffset_;
  int_type ColOffset = Col * COffset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine copies entries of a BaseMatrix into big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // this->val[i+rowOffset][j+ColOffset] = BaseMatrix.val[i][j]

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices_local(MaxIndices);
  vector<int_type> Indices_global(MaxIndices);
  vector<double> vals(MaxIndices);
  int NumIndices;
  int ierr=0;

  for (int i=0; i<BaseMap.NumMyElements(); i++) {
    BaseMatrix.ExtractMyRowCopy( i, MaxIndices, NumIndices, &vals[0], &Indices_local[0] );

    // Convert to BlockMatrix Global numbering scheme
    for( int l = 0; l < NumIndices; ++l ) {
       Indices_global[l] = ColOffset +  (int_type) BaseColMap.GID64(Indices_local[l]);
       vals[l] *= alpha;
    }

    int_type BaseRow = (int_type) BaseMap.GID64(i);
    ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices, &vals[0], &Indices_global[0]);
    if (ierr != 0) {
      std::cout << "WARNING BlockCrsMatrix::SumIntoBlock SumIntoGlobalValues "
                << "err = " << ierr << std::endl
                << "\t  Row " << BaseRow + RowOffset
                << " Col start" << Indices_global[0] << std::endl;
    }
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt() && BaseMatrix.RowMatrixRowMap().GlobalIndicesInt())
        return TSumIntoGlobalBlock<int>(alpha, BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::SumIntoGlobalBlock: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && BaseMatrix.RowMatrixRowMap().GlobalIndicesLongLong())
        return TSumIntoGlobalBlock<long long>(alpha, BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::SumIntoGlobalBlock: Global indices not long long";
}
#endif


//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TBlockSumIntoGlobalValues(const int_type BaseRow, int NumIndices,
     double* vals, const int_type* Indices, const int_type Row, const int_type Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

  // Convert to BlockMatrix Global numbering scheme
  vector<int_type> OffsetIndices(NumIndices);
  for( int l = 0; l < NumIndices; ++l ) OffsetIndices[l] = Indices[l] + ColOffset;

  int ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices,
                                   vals, &OffsetIndices[0]);

  if (ierr != 0) {
    std::cout << "WARNING BlockCrsMatrix::BlockSumIntoGlobalValues err = "
              << ierr << std::endl << "\t  Row " << BaseRow + RowOffset
              << " Col start" << Indices[0] << std::endl;
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockSumIntoGlobalValues(const int BaseRow, int NumIndices,
     double* vals, const int* Indices, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt())
    return TBlockSumIntoGlobalValues<int>(BaseRow, NumIndices, vals, Indices, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockSumIntoGlobalValues: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockSumIntoGlobalValues(const long long BaseRow, int NumIndices,
     double* vals, const long long* Indices, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong())
    return TBlockSumIntoGlobalValues<long long>(BaseRow, NumIndices, vals, Indices, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockSumIntoGlobalValues: Global indices not long long";
}
#endif

//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TBlockReplaceGlobalValues(const int_type BaseRow, int NumIndices,
     double* vals, const int_type* Indices, const int_type Row, const int_type Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

  // Convert to BlockMatrix Global numbering scheme
  vector<int_type> OffsetIndices(NumIndices);
  for( int l = 0; l < NumIndices; ++l ) OffsetIndices[l] = Indices[l] + ColOffset;

  int ierr = this->ReplaceGlobalValues(BaseRow + RowOffset, NumIndices,
                                       vals, &OffsetIndices[0]);

  if (ierr != 0) {
    std::cout << "WARNING BlockCrsMatrix::BlockReplaceGlobalValues err = "
              << ierr << "\n\t  Row " << BaseRow + RowOffset << "Col start"
              << Indices[0] << std::endl;
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockReplaceGlobalValues(const int BaseRow, int NumIndices,
     double* vals, const int* Indices, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt())
    return TBlockReplaceGlobalValues<int>(BaseRow, NumIndices, vals, Indices, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockReplaceGlobalValues: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockReplaceGlobalValues(const long long BaseRow, int NumIndices,
     double* vals, const long long* Indices, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong())
    return TBlockReplaceGlobalValues<long long>(BaseRow, NumIndices, vals, Indices, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockReplaceGlobalValues: Global indices not long long";
}
#endif

//==============================================================================
template<typename int_type>
void BlockCrsMatrix::TBlockExtractGlobalRowView(const int_type BaseRow,
                                               int& NumEntries,
                                               double*& vals,
                                               const int_type Row,
                                               const int_type Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

  // Get the whole row
  int ierr = this->ExtractGlobalRowView(BaseRow + RowOffset, NumEntries,
                                        vals);

  // Adjust for just this block column
  vals += ColOffset;
  NumEntries -= ColOffset;

  if (ierr != 0) {
    std::cout << "WARNING BlockCrsMatrix::BlockExtractGlobalRowView err = "
              << ierr << "\n\t  Row " << BaseRow + RowOffset
              << " Col " << Col+ColOffset << std::endl;
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockExtractGlobalRowView(const int BaseRow, int& NumEntries,
     double*& vals, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt())
    return TBlockExtractGlobalRowView<int>(BaseRow, NumEntries, vals, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockExtractGlobalRowView: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::BlockExtractGlobalRowView(const long long BaseRow, int& NumEntries,
     double*& vals, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong())
    return TBlockExtractGlobalRowView<long long>(BaseRow, NumEntries, vals, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::BlockExtractGlobalRowView: Global indices not long long";
}
#endif

//==============================================================================

template<typename int_type>
void BlockCrsMatrix::TExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int_type Row, const int_type Col)
{
  std::vector<int_type>& RowIndices_ = TRowIndices<int_type>();
  std::vector< std::vector<int_type> >& RowStencil_ = TRowStencil<int_type>();
  int_type RowOffset = RowIndices_[(std::size_t)Row] * ROffset_;
  int_type ColOffset = (RowIndices_[(std::size_t)Row] + RowStencil_[(std::size_t)Row][(std::size_t)Col]) * COffset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  //const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine extracts entries of a BaseMatrix from a big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // BaseMatrix.val[i][j] = this->val[i+rowOffset][j+ColOffset]

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int_type> Indices(MaxIndices);
  vector<double> vals(MaxIndices);
  int NumIndices;
  int_type indx,icol;
  double* BlkValues;
  int *BlkIndices;
  int BlkNumIndices;
  int ierr=0;
  (void) ierr; // Forestall compiler warning for unused variable.

  for (int i=0; i<BaseMap.NumMyElements(); i++) {

    // Get pointers to values and indices of whole block matrix row
    int_type BaseRow = (int_type) BaseMap.GID64(i);
    int myBlkBaseRow = this->RowMatrixRowMap().LID(BaseRow + RowOffset);
    ierr = this->ExtractMyRowView(myBlkBaseRow, BlkNumIndices, BlkValues, BlkIndices);

    NumIndices = 0;
    // Grab columns with global indices in correct range for this block
    for( int l = 0; l < BlkNumIndices; ++l ) {
       icol = (int_type) this->RowMatrixColMap().GID64(BlkIndices[l]);
       indx = icol - ColOffset;
       if (indx >= 0 && indx < COffset_) {
         Indices[NumIndices] = indx;
         vals[NumIndices] = BlkValues[l];
         NumIndices++;
       }
    }

    //Load this row into base matrix
    BaseMatrix.ReplaceGlobalValues(BaseRow, NumIndices, &vals[0], &Indices[0]);
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockCrsMatrix::ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int Row, const int Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesInt() && BaseMatrix.RowMatrixRowMap().GlobalIndicesInt())
        return TExtractBlock<int>(BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::ExtractBlock: Global indices not int";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockCrsMatrix::ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const long long Row, const long long Col)
{
  if(Epetra_CrsMatrix::RowMatrixRowMap().GlobalIndicesLongLong() && BaseMatrix.RowMatrixRowMap().GlobalIndicesLongLong())
        return TExtractBlock<long long>(BaseMatrix, Row, Col);
  else
    throw "EpetraExt::BlockCrsMatrix::ExtractBlock: Global indices not long long";
}
#endif

} //namespace EpetraExt

