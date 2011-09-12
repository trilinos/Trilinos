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

#include <EpetraExt_BTF_LinearProblem.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LinearProblem.h>

#include <set>

using std::vector;
using std::map;
using std::set;

#define MATTRANS_F77 F77_FUNC(mattrans,MATTRANS)
#define GENBTF_F77   F77_FUNC(genbtf,GENBTF)
                                                                                                  
extern "C" {
extern void MATTRANS_F77( int*, int*, int*, int*, int*, int* );
extern void GENBTF_F77( int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int* );
}

namespace EpetraExt {

LinearProblem_BTF::
~LinearProblem_BTF()
{
  deleteNewObjs_();
}

void
LinearProblem_BTF::
deleteNewObjs_()
{
  if( NewProblem_ ) delete NewProblem_;

  if( NewMatrix_ ) delete NewMatrix_;

  if( NewLHS_ ) delete NewLHS_;
  if( NewRHS_ ) delete NewRHS_;

  if( NewMap_ ) delete NewMap_;

  for( int i = 0; i < Blocks_.size(); ++i )
    for( int j = 0; j < Blocks_[i].size(); ++j )
      delete Blocks_[i][j];
}

LinearProblem_BTF::NewTypeRef
LinearProblem_BTF::
operator()( OriginalTypeRef orig )
{
  changedLP_ = false;

  //check if there is a previous analysis and if it is valid
  if( &orig == origObj_ && NewProblem_ )
  {
    int * indices;
    double * values;
    int currCnt;
    int numRows = OrigMatrix_->NumMyRows();

    for( int i = 0; i < numRows; ++i )
      if( ZeroElements_[i].size() )
      {
        int loc = 0;
        OrigMatrix_->ExtractMyRowView( i, currCnt, values, indices );
        for( int j = 0; j < currCnt; ++j )
          if( ZeroElements_[i].count( indices[j] ) )
          {
            if( values[j] > threshold_ ) changedLP_ = true;
            ++loc;
          }
      }

//    changedLP_ = true;

    //if changed, dump all the old stuff and start over
    if( changedLP_ )
      deleteNewObjs_();
    else
      return *newObj_;
  }
    
  origObj_ = &orig;

  OrigMatrix_ = dynamic_cast<Epetra_CrsMatrix*>(orig.GetMatrix());
  OrigLHS_ = orig.GetLHS();
  OrigRHS_ = orig.GetRHS();

  if( OrigMatrix_->RowMap().DistributedGlobal() )
  { cout << "FAIL for Global!\n"; abort(); }
  if( OrigMatrix_->IndicesAreGlobal() )
  { cout << "FAIL for Global Indices!\n"; abort(); }

  int n = OrigMatrix_->NumMyRows();
  int nnz = OrigMatrix_->NumMyNonzeros();

//  cout << "Orig Matrix:\n";
//  cout << *OrigMatrix_ << endl;

  //create std CRS format
  //also create graph without zero elements
  vector<int> ia(n+1,0);
  int maxEntries = OrigMatrix_->MaxNumEntries();
  vector<int> ja_tmp(maxEntries);
  vector<double> jva_tmp(maxEntries);
  vector<int> ja(nnz);
  int cnt;

  const Epetra_BlockMap & OldRowMap = OrigMatrix_->RowMap();
  const Epetra_BlockMap & OldColMap = OrigMatrix_->ColMap();
  Epetra_CrsGraph strippedGraph( Copy, OldRowMap, OldColMap, 0 );
  ZeroElements_.resize(n);

  for( int i = 0; i < n; ++i )
  {
    ZeroElements_[i].clear();
    OrigMatrix_->ExtractMyRowCopy( i, maxEntries, cnt, &jva_tmp[0], &ja_tmp[0] );
    ia[i+1] = ia[i];
    for( int j = 0; j < cnt; ++j )
    {
      if( fabs(jva_tmp[j]) > threshold_ )
        ja[ ia[i+1]++ ] = ja_tmp[j];
      else
        ZeroElements_[i].insert( ja_tmp[j] );
    }

    int new_cnt = ia[i+1] - ia[i];
    strippedGraph.InsertMyIndices( i, new_cnt, &ja[ ia[i] ] );
  }
  nnz = ia[n];
  strippedGraph.FillComplete();

  if( verbose_ > 2 )
  {
    cout << "Stripped Graph\n";
    cout << strippedGraph;
  }

  vector<int> iat(n+1,0);
  vector<int> jat(nnz);
  for( int i = 0; i < n; ++i )
    for( int j = ia[i]; j < ia[i+1]; ++j )
      ++iat[ ja[j]+1 ];
  for( int i = 0; i < n; ++i )
    iat[i+1] += iat[i];
  for( int i = 0; i < n; ++i )
    for( int j = ia[i]; j < ia[i+1]; ++j )
      jat[ iat[ ja[j] ]++ ] = i;
  for( int i = 0; i < n; ++i )
    iat[n-i] = iat[n-i-1];
  iat[0] = 0;

  //convert to Fortran indexing
  for( int i = 0; i < n+1; ++i ) ++ia[i];
  for( int i = 0; i < nnz; ++i ) ++ja[i];
  for( int i = 0; i < n+1; ++i ) ++iat[i];
  for( int i = 0; i < nnz; ++i ) ++jat[i];

  if( verbose_ > 2 )
  {
    cout << "-----------------------------------------\n";
    cout << "CRS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << i+1 << ": " << ia[i+1] << ": ";
      for( int j = ia[i]-1; j<ia[i+1]-1; ++j )
        cout << " " << ja[j];
      cout << endl;
    }
    cout << "-----------------------------------------\n";
  }

  if( verbose_ > 2 )
  {
    cout << "-----------------------------------------\n";
    cout << "CCS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << i+1 << ": " << iat[i+1] << ": ";
      for( int j = iat[i]-1; j<iat[i+1]-1; ++j )
        cout << " " << jat[j];
      cout << endl;
    }
    cout << "-----------------------------------------\n";
  }

  vector<int> w(10*n);

  vector<int> rowperm(n);
  vector<int> colperm(n);

  //horizontal block
  int nhrows, nhcols, hrzcmp;
  //square block
  int nsrows, sqcmpn;
  //vertial block
  int nvrows, nvcols, vrtcmp;

  vector<int> rcmstr(n+1);
  vector<int> ccmstr(n+1);

  int msglvl = 0;
  int output = 6;

  GENBTF_F77( &n, &n, &iat[0], &jat[0], &ia[0], &ja[0], &w[0],
          &rowperm[0], &colperm[0], &nhrows, &nhcols,
          &hrzcmp, &nsrows, &sqcmpn, &nvrows, &nvcols, &vrtcmp,
          &rcmstr[0], &ccmstr[0], &msglvl, &output );

  //convert back to C indexing
  for( int i = 0; i < n; ++i )
  {
    --rowperm[i];
    --colperm[i];
  }
  for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
  {
    --rcmstr[i];
    --ccmstr[i];
  }

  if( verbose_ > 0 )
  {
    cout << "-----------------------------------------\n";
    cout << "BTF Output\n";
    cout << "-----------------------------------------\n";
//    cout << "RowPerm and ColPerm\n";
//    for( int i = 0; i<n; ++i )
//      cout << rowperm[i] << "\t" << colperm[i] << endl;
    if( hrzcmp )
    {
      cout << "Num Horizontal: Rows, Cols, Comps\n";
      cout << nhrows << "\t" << nhcols << "\t" << hrzcmp << endl;
    }
    cout << "Num Square: Rows, Comps\n";
    cout << nsrows << "\t" << sqcmpn << endl;
    if( vrtcmp )
    {
      cout << "Num Vertical: Rows, Cols, Comps\n";
      cout << nvrows << "\t" << nvcols << "\t" << vrtcmp << endl;
    }
//    cout << "Row, Col of upper left pt in blocks\n";
//    for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
//      cout << i << " " << rcmstr[i] << " " << ccmstr[i] << endl;
//    cout << "-----------------------------------------\n";
  }

  if( hrzcmp || vrtcmp )
  {
    cout << "FAILED! hrz cmp's:" << hrzcmp << " vrtcmp: " << vrtcmp << endl;
    exit(0);
  }

  rcmstr[sqcmpn] = n;

  //convert rowperm to OLD->NEW
  //reverse ordering of permutation to get upper triangular
  vector<int> rowperm_t( n );
  vector<int> colperm_t( n );
  for( int i = 0; i < n; ++i )
  {
    rowperm_t[i] = rowperm[i];
    colperm_t[i] = colperm[i];
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  OldGlobalElements_.resize(n);
  OldRowMap.MyGlobalElements( &OldGlobalElements_[0] );

  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newDomainElements[ i ] = OldGlobalElements_[ rowperm_t[i] ];
    newRangeElements[ i ] = OldGlobalElements_[ colperm_t[i] ];
  }

  //Setup New Block Map Info
  Blocks_.resize( sqcmpn );
  BlockDim_.resize( sqcmpn );
  for( int i = 0; i < sqcmpn; ++i )
  {
    BlockDim_[i] = rcmstr[i+1]-rcmstr[i];
    for( int j = rcmstr[i]; j < rcmstr[i+1]; ++j )
    {
      BlockRowMap_[ newDomainElements[j] ] = i;
      SubBlockRowMap_[ newDomainElements[j] ] = j-rcmstr[i];
      BlockColMap_[ newRangeElements[j] ] = i;
      SubBlockColMap_[ newRangeElements[j] ] = j-rcmstr[i];
    }
  }

  if( verbose_ > 2 )
  {
/*
    cout << "Block Mapping!\n";
    cout << "--------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << "Row: " << newDomainElements[i] << " " << BlockRowMap_[newDomainElements[i]] << " " <<
              SubBlockRowMap_[newDomainElements[i]] << "\t" << "Col: " << newRangeElements[i] << " " <<
              BlockColMap_[newRangeElements[i]] << " " << SubBlockColMap_[newRangeElements[i]] << endl;
    }
    for( int i = 0; i < sqcmpn; ++i )
      cout << "BlockDim: " << i << " " << BlockDim_[i] << endl;
    cout << "--------------------------\n";
*/
    int MinSize = 1000000000;
    int MaxSize = 0;
    for( int i = 0; i < sqcmpn; ++i )
    {
      if( MinSize > BlockDim_[i] ) MinSize = BlockDim_[i];
      if( MaxSize < BlockDim_[i] ) MaxSize = BlockDim_[i];
    }
    cout << "Min Block Size: " << MinSize << " " << "Max Block Size: " << MaxSize << endl;
  }

  vector<int> myBlockElements( sqcmpn );
  for( int i = 0; i < sqcmpn; ++i ) myBlockElements[i] = i;
  NewMap_ = new Epetra_BlockMap( sqcmpn, sqcmpn, &myBlockElements[0], &BlockDim_[0], 0, OldRowMap.Comm() );

  if( verbose_ > 2 )
  {
    cout << "New Block Map!\n";
    cout << *NewMap_;
  }

  //setup new graph
  vector< set<int> > crsBlocks( sqcmpn );
  BlockCnt_.resize( sqcmpn );
  int maxLength = strippedGraph.MaxNumIndices();
  vector<int> sIndices( maxLength );
  int currLength;
  for( int i = 0; i < n; ++i )
  {
    strippedGraph.ExtractGlobalRowCopy( OldGlobalElements_[i], maxLength, currLength, &sIndices[0] );
    for( int j = 0; j < currLength; ++j )
      crsBlocks[ BlockRowMap_[ OldGlobalElements_[i] ] ].insert( BlockColMap_[ sIndices[j] ] );
  }

  for( int i = 0; i < sqcmpn; ++i )
  {
    BlockCnt_[i] = crsBlocks[i].size();
    Blocks_[i].resize( BlockCnt_[i] );
  }

  NewBlockRows_.resize( sqcmpn );
  for( int i = 0; i < sqcmpn; ++i )
  {
    NewBlockRows_[i] = vector<int>( crsBlocks[i].begin(), crsBlocks[i].end() );
    for( int j = 0; j < BlockCnt_[i]; ++j )
    {
      Blocks_[i][j] = new Epetra_SerialDenseMatrix();
      Blocks_[i][j]->Shape( BlockDim_[i], BlockDim_[ NewBlockRows_[i][j] ] );
    }
  }

  //put data in Blocks_ and new LHS and RHS
  NewLHS_ = new Epetra_MultiVector( *NewMap_, 1 );
  NewRHS_ = new Epetra_MultiVector( *NewMap_, 1 );

  maxLength = OrigMatrix_->MaxNumEntries();
  vector<int> indices( maxLength );
  vector<double> values( maxLength );
  for( int i = 0; i < n; ++i )
  {
    int BlockRow = BlockRowMap_[ OldGlobalElements_[i] ];
    int SubBlockRow = SubBlockRowMap_[ OldGlobalElements_[i] ];
    OrigMatrix_->ExtractGlobalRowCopy( OldGlobalElements_[i], maxLength, currLength, &values[0], &indices[0] );
    for( int j = 0; j < currLength; ++j )
    {
      int BlockCol = BlockColMap_[ indices[j] ];
      int SubBlockCol = SubBlockColMap_[ indices[j] ];
      for( int k = 0; k < BlockCnt_[BlockRow]; ++k )
        if( BlockCol == NewBlockRows_[BlockRow][k] )
        {
          if( values[j] > threshold_ )
          {
//          cout << BlockRow << " " << SubBlockRow << " " << BlockCol << " " << SubBlockCol << endl;
//	    cout << k << endl;
	    (*(Blocks_[BlockRow][k]))(SubBlockRow,SubBlockCol) = values[j];
            break;
          }
          else
            ZeroElements_[i].erase( OrigMatrix_->RowMap().LID( indices[j] ) );
        }
    }

//    NewLHS_->ReplaceGlobalValue( BlockCol, SubBlockCol, 0, (*OrigLHS_)[0][i] );
    NewRHS_->ReplaceGlobalValue( BlockRow, SubBlockRow, 0, (*OrigRHS_)[0][i] );
  }

  if( verbose_ > 2 )
  {
    cout << "Zero Elements: \n";
    cout << "--------------\n";
    int cnt = 0;
    for( int i = 0; i < n; ++i )
    {
      set<int>::iterator iterSI = ZeroElements_[i].begin();
      set<int>::iterator endSI = ZeroElements_[i].end();
      for( ; iterSI != endSI; ++iterSI )
      {
        cout << " " << *iterSI;
        ++cnt;
      }
      cout << endl;
    }
    cout << "ZE Cnt: " << cnt << endl;
    cout << "--------------\n";
  }

  //setup new matrix
  NewMatrix_ = new Epetra_VbrMatrix( View, *NewMap_, &BlockCnt_[0] );
  for( int i = 0; i < sqcmpn; ++i )
  {
    NewMatrix_->BeginInsertGlobalValues( i, BlockCnt_[i], &(NewBlockRows_[i])[0] );
    for( int j = 0; j < BlockCnt_[i]; ++j )
      NewMatrix_->SubmitBlockEntry( *(Blocks_[i][j]) );
    NewMatrix_->EndSubmitEntries();
  }
  NewMatrix_->FillComplete();

  if( verbose_ > 2 )
  {
    cout << "New Block Matrix!\n";
    cout << *NewMatrix_;
    cout << "New Block LHS!\n";
    cout << *NewLHS_;
    cout << "New Block RHS!\n";
    cout << *NewRHS_;
  }

  //create new LP
  NewProblem_ = new Epetra_LinearProblem( NewMatrix_, NewLHS_, NewRHS_ );
  newObj_ = NewProblem_;

  if( verbose_ ) cout << "-----------------------------------------\n";

  return *newObj_;
}

bool
LinearProblem_BTF::
fwd()
{
  //zero out matrix
  int NumBlockRows = BlockDim_.size();
  for( int i = 0; i < NumBlockRows; ++i )
  {
    int NumBlocks = BlockCnt_[i];
    for( int j = 0; j < NumBlocks; ++j )
    {
      int Size = BlockDim_[i] * BlockDim_[ NewBlockRows_[i][j] ];
      double * p = Blocks_[i][j]->A();
      for( int k = 0; k < Size; ++k ) *p++ = 0.0;
    }
  }

  int maxLength = OrigMatrix_->MaxNumEntries();
  int n = OldGlobalElements_.size();
  int currLength;
  vector<int> indices( maxLength );
  vector<double> values( maxLength );
  for( int i = 0; i < n; ++i )
  {
    int BlockRow = BlockRowMap_[ OldGlobalElements_[i] ];
    int SubBlockRow = SubBlockRowMap_[ OldGlobalElements_[i] ];
    OrigMatrix_->ExtractGlobalRowCopy( OldGlobalElements_[i], maxLength, currLength, &values[0], &indices[0] );
    for( int j = 0; j < currLength; ++j )
    {
      if( fabs(values[j]) > threshold_ )
      {
        int BlockCol = BlockColMap_[ indices[j] ];
        int SubBlockCol = SubBlockColMap_[ indices[j] ];
	for( int k = 0; k < BlockCnt_[BlockRow]; ++k )
          if( BlockCol == NewBlockRows_[BlockRow][k] )
          {
	    (*(Blocks_[BlockRow][k]))(SubBlockRow,SubBlockCol) = values[j];
            break;
          }
      }
    }

    NewRHS_->ReplaceGlobalValue( BlockRow, SubBlockRow, 0, (*OrigRHS_)[0][i] );
  }

/*
  //fill matrix
  int sqcmpn = BlockDim_.size();
  for( int i = 0; i < sqcmpn; ++i )
  {
    NewMatrix_->BeginReplaceGlobalValues( i, NewBlockRows_[i].size(), &(NewBlockRows_[i])[0] );
    for( int j = 0; j < NewBlockRows_[i].size(); ++j )
      NewMatrix_->SubmitBlockEntry( Blocks_[i][j]->A(), Blocks_[i][j]->LDA(), Blocks_[i][j]->M(), Blocks_[i][j]->N() );
    NewMatrix_->EndSubmitEntries();
  }
*/

  return true;
}

bool
LinearProblem_BTF::
rvs()
{
  //copy data from NewLHS_ to OldLHS_;
  int rowCnt = OrigLHS_->MyLength();
  for( int i = 0; i < rowCnt; ++i )
  {
    int BlockCol = BlockColMap_[ OldGlobalElements_[i] ];
    int SubBlockCol = SubBlockColMap_[ OldGlobalElements_[i] ];
    (*OrigLHS_)[0][i] = (*NewLHS_)[0][ NewMap_->FirstPointInElement(BlockCol) + SubBlockCol ];
  }

  return true;
}

} //namespace EpetraExt
