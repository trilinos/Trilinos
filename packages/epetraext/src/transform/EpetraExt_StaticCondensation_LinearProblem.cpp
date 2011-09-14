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

#include <EpetraExt_StaticCondensation_LinearProblem.h>

#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <vector>
#include <map>
#include <set>

namespace EpetraExt {

LinearProblem_StaticCondensation::
~LinearProblem_StaticCondensation()
{
  if( Exporter_ ) delete Exporter_;

  if( NewProblem_ ) delete NewProblem_;
  if( NewRHS_ ) delete NewRHS_;
  if( NewLHS_ ) delete NewLHS_;
  if( NewMatrix_ ) delete NewMatrix_;
  if( NewGraph_ ) delete NewGraph_;
  if( NewRowMap_ ) delete NewRowMap_;
  if( NewColMap_ ) delete NewColMap_;

  if( ULHS_ ) delete ULHS_;
  if( RLHS_ ) delete RLHS_;
  if( LLHS_ ) delete LLHS_;
  if( URHS_ ) delete URHS_;
  if( RRHS_ ) delete RRHS_;
  if( LRHS_ ) delete LRHS_;

  if( UUMatrix_ ) delete UUMatrix_;
  if( URMatrix_ ) delete URMatrix_;
  if( ULMatrix_ ) delete ULMatrix_;
  if( RRMatrix_ ) delete RRMatrix_;
  if( RLMatrix_ ) delete RLMatrix_;
  if( LLMatrix_ ) delete LLMatrix_;

  if( UUGraph_ ) delete UUGraph_;
  if( URGraph_ ) delete URGraph_;
  if( ULGraph_ ) delete ULGraph_;
  if( RRGraph_ ) delete RRGraph_;
  if( RLGraph_ ) delete RLGraph_;
  if( LLGraph_ ) delete LLGraph_;

  if( UExporter_ ) delete UExporter_;
  if( RExporter_ ) delete RExporter_;
  if( LExporter_ ) delete LExporter_;

  if( UMap_ ) delete UMap_;
  if( RMap_ ) delete RMap_;
  if( LMap_ ) delete LMap_;
}

LinearProblem_StaticCondensation::NewTypeRef
LinearProblem_StaticCondensation::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int ierr = 0;

  OldMatrix_ = dynamic_cast<Epetra_CrsMatrix*>( orig.GetMatrix() );
  OldGraph_  = &OldMatrix_->Graph();
  OldRHS_    = orig.GetRHS();
  OldLHS_    = orig.GetLHS();
  OldRowMap_ = &OldMatrix_->RowMap();

  const Epetra_Comm & CommObj = OldRowMap_->Comm();

  if( !OldMatrix_ ) ierr = -2;
  if( !OldRHS_ )    ierr = -3;
  if( !OldLHS_ )    ierr = -4;

  if( OldRowMap_->DistributedGlobal() ) ierr = -5;
  if( degree_ != 1 ) ierr = -6;

  int NRows = OldGraph_->NumMyRows();
  int IndexBase = OldRowMap_->IndexBase();

  vector<int> ColNZCnt( NRows );
  vector<int> CS_RowIndices( NRows );
  map<int,int> RS_Map;
  map<int,int> CS_Map;

  int NumIndices;
  int * Indices;
  for( int i = 0; i < NRows; ++i )
  {
    ierr = OldGraph_->ExtractMyRowView( i, NumIndices, Indices );

    for( int j = 0; j < NumIndices; ++j )
    {
      ++ColNZCnt[ Indices[j] ];
      CS_RowIndices[ Indices[j] ] = i;
    }

    if( NumIndices == 1 ) RS_Map[i] = Indices[0];
  }

  if( verbose_ )
  {
    cout << "-------------------------\n";
    cout << "Row Singletons\n";
    for( map<int,int>::iterator itM = RS_Map.begin(); itM != RS_Map.end(); ++itM )
      cout << (*itM).first << "\t" << (*itM).second << endl;
    cout << "Col Counts\n";
    for( int i = 0; i < NRows; ++i )
      cout << i << "\t" << ColNZCnt[i] << "\t" << CS_RowIndices[i] << endl;
    cout << "-------------------------\n";
  }

  set<int> RS_Set;
  set<int> CS_Set;

  for( int i = 0; i < NRows; ++i )
    if( ColNZCnt[i] == 1 )
    {
      int RowIndex = CS_RowIndices[i];
      if( RS_Map.count(i) && RS_Map[i] == RowIndex )
      {
        CS_Set.insert(i);
        RS_Set.insert( RowIndex );
      }
    }

  if( verbose_ )
  {
    cout << "-------------------------\n";
    cout << "Singletons: " << CS_Set.size() << endl;
    set<int>::iterator itRS = RS_Set.begin();
    set<int>::iterator itCS = CS_Set.begin();
    set<int>::iterator endRS = RS_Set.end();
    set<int>::iterator endCS = CS_Set.end();
    for( ; itRS != endRS; ++itRS, ++itCS )
      cout << *itRS << "\t" << *itCS << endl;
    cout << "-------------------------\n";
  }

  //build new maps
  int NSingletons = CS_Set.size();
  int NReducedRows = NRows - 2* NSingletons; 
  vector<int> ReducedIndices( NReducedRows );
  vector<int> CSIndices( NSingletons );
  vector<int> RSIndices( NSingletons );
  int Reduced_Loc = 0;
  int RS_Loc = 0;
  int CS_Loc = 0;
  for( int i = 0; i < NRows; ++i )
  {
    int GlobalIndex = OldRowMap_->GID(i);
    if     ( RS_Set.count(i) ) RSIndices[RS_Loc++] = GlobalIndex;
    else if( CS_Set.count(i) ) CSIndices[CS_Loc++] = GlobalIndex;
    else                       ReducedIndices[Reduced_Loc++] = GlobalIndex;
  }

  vector<int> RowPermutedIndices( NRows );
  copy( RSIndices.begin(), RSIndices.end(), RowPermutedIndices.begin() );
  copy( ReducedIndices.begin(), ReducedIndices.end(), RowPermutedIndices.begin() + NSingletons );
  copy( CSIndices.begin(), CSIndices.end(), RowPermutedIndices.begin() + NReducedRows + NSingletons );

  vector<int> ColPermutedIndices( NRows );
  copy( CSIndices.begin(), CSIndices.end(), ColPermutedIndices.begin() );
  copy( ReducedIndices.begin(), ReducedIndices.end(), ColPermutedIndices.begin() + NSingletons );
  copy( RSIndices.begin(), RSIndices.end(), ColPermutedIndices.begin() + NReducedRows + NSingletons );

  UMap_ = new Epetra_Map( NSingletons, NSingletons, &RSIndices[0], IndexBase, CommObj );
  RMap_ = new Epetra_Map( NReducedRows, NReducedRows, &ReducedIndices[0], IndexBase, CommObj );
  LMap_ = new Epetra_Map( NSingletons, NSingletons, &CSIndices[0], IndexBase, CommObj );

  NewRowMap_ = new Epetra_Map( NRows, NRows, &RowPermutedIndices[0], IndexBase, CommObj );
  NewColMap_ = new Epetra_Map( NRows, NRows, &ColPermutedIndices[0], IndexBase, CommObj );

  //Construct Permuted System
  Exporter_ = new Epetra_Export( *OldRowMap_, *NewRowMap_ );

  NewRHS_ = new Epetra_MultiVector( *NewRowMap_, OldRHS_->NumVectors() );
  NewRHS_->Export( *OldRHS_, *Exporter_, Insert );

  NewLHS_ = new Epetra_MultiVector( *NewRowMap_, OldLHS_->NumVectors() );
  NewLHS_->Export( *OldLHS_, *Exporter_, Insert );

  NewGraph_ = new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 );
  NewGraph_->Export( *OldGraph_, *Exporter_, Insert );
  NewGraph_->FillComplete();
  cout << "Permuted Graph:\n" << *NewGraph_;

  NewMatrix_ = new Epetra_CrsMatrix( Copy, *NewGraph_ );
  NewMatrix_->Export( *OldMatrix_, *Exporter_, Insert );
  NewMatrix_->FillComplete();
  cout << "Permuted Matrix:\n" << *NewMatrix_;

  UExporter_ = new Epetra_Export( *OldRowMap_, *UMap_ );
  RExporter_ = new Epetra_Export( *OldRowMap_, *RMap_ );
  LExporter_ = new Epetra_Export( *OldRowMap_, *LMap_ );
cout << "UExporter:\n" << *UExporter_;
cout << "RExporter:\n" << *RExporter_;
cout << "LExporter:\n" << *LExporter_;

  ULHS_ = new Epetra_MultiVector( *LMap_, OldLHS_->NumVectors() );
  ULHS_->Export( *OldLHS_, *LExporter_, Insert );
  cout << "ULHS:\n" << *ULHS_;

  RLHS_ = new Epetra_MultiVector( *RMap_, OldLHS_->NumVectors() );
  RLHS_->Export( *OldLHS_, *RExporter_, Insert );
  cout << "RLHS:\n" << *RLHS_;

  LLHS_ = new Epetra_MultiVector( *UMap_, OldLHS_->NumVectors() );
  LLHS_->Export( *OldLHS_, *UExporter_, Insert );
  cout << "LLHS:\n" << *LLHS_;

  URHS_ = new Epetra_MultiVector( *UMap_, OldRHS_->NumVectors() );
  URHS_->Export( *OldRHS_, *UExporter_, Insert );
  cout << "URHS:\n" << *URHS_;

  RRHS_ = new Epetra_MultiVector( *RMap_, OldRHS_->NumVectors() );
  RRHS_->Export( *OldRHS_, *RExporter_, Insert );
  cout << "RRHS:\n" << *RRHS_;

  LRHS_ = new Epetra_MultiVector( *LMap_, OldRHS_->NumVectors() );
  LRHS_->Export( *OldRHS_, *LExporter_, Insert );
  cout << "LRHS:\n" << *LRHS_;

  UUGraph_ = new Epetra_CrsGraph( Copy, *UMap_, *LMap_, 0 );
  UUGraph_->Export( *OldGraph_, *UExporter_, Insert );
  UUGraph_->FillComplete( LMap_, UMap_ );
  cout << "UUGraph:\n" << *UUGraph_;

  UUMatrix_ = new Epetra_CrsMatrix( Copy, *UUGraph_ );
  UUMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  UUMatrix_->FillComplete();
  cout << "UUMatrix:\n" << *UUMatrix_;

  URGraph_ = new Epetra_CrsGraph( Copy, *UMap_, *RMap_, 0 );
  URGraph_->Export( *OldGraph_, *UExporter_, Insert );
  URGraph_->FillComplete( RMap_, UMap_ );
  cout << "URGraph:\n" << *URGraph_;

  URMatrix_ = new Epetra_CrsMatrix( Copy, *URGraph_ );
  URMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  URMatrix_->FillComplete();
  cout << "URMatrix:\n" << *URMatrix_;

  ULGraph_ = new Epetra_CrsGraph( Copy, *UMap_, *UMap_, 0 );
  ULGraph_->Export( *OldGraph_, *UExporter_, Insert );
  ULGraph_->FillComplete();
  cout << "ULGraph:\n" << *ULGraph_;

  ULMatrix_ = new Epetra_CrsMatrix( Copy, *ULGraph_ );
  ULMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  ULMatrix_->FillComplete();
  cout << "ULMatrix:\n" << *ULMatrix_;

  RRGraph_ = new Epetra_CrsGraph( Copy, *RMap_, *RMap_, 0 );
  RRGraph_->Export( *OldGraph_, *RExporter_, Insert );
  RRGraph_->FillComplete();
  cout << "RRGraph:\n" << *RRGraph_;

  RRMatrix_ = new Epetra_CrsMatrix( Copy, *RRGraph_ );
  RRMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  RRMatrix_->FillComplete();
  cout << "RRMatrix:\n" << *RRMatrix_;

  RLGraph_ = new Epetra_CrsGraph( Copy, *RMap_, *UMap_, 0 );
  RLGraph_->Export( *OldGraph_, *RExporter_, Insert );
  RLGraph_->FillComplete( UMap_, RMap_ );
  cout << "RLGraph:\n" << *RLGraph_;

  RLMatrix_ = new Epetra_CrsMatrix( Copy, *RLGraph_ );
  RLMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  RLMatrix_->FillComplete();
  cout << "RLMatrix:\n" << *RLMatrix_;

  LLGraph_ = new Epetra_CrsGraph( Copy, *LMap_, *UMap_, 0 );
  LLGraph_->Export( *OldGraph_, *LExporter_, Insert );
  LLGraph_->FillComplete( UMap_, LMap_ );
  cout << "LLGraph:\n" << *LLGraph_;

  LLMatrix_ = new Epetra_CrsMatrix( Copy, *LLGraph_ );
  LLMatrix_->Export( *OldMatrix_, *LExporter_, Insert );
  LLMatrix_->FillComplete();
  cout << "LLMatrix:\n" << *LLMatrix_;

  if( verbose_ )
  {
    cout << "Full System Characteristics:" << endl
         << "----------------------------" << endl
         << " Dimension                   = " << NRows << endl
         << " NNZs                        = " << OldMatrix_->NumGlobalNonzeros() << endl
         << " Max Num Row Entries         = " << OldMatrix_->GlobalMaxNumEntries() << endl << endl
         << "Reduced System Characteristics:" << endl
         << " Dimension                   = " << NReducedRows << endl
         << " NNZs                        = " << RRMatrix_->NumGlobalNonzeros() << endl
         << " Max Num Row Entries         = " << RRMatrix_->GlobalMaxNumEntries() << endl
         << "Singleton Info:" << endl
         << " Num Singleton                 = " << NSingletons << endl
         << "Ratios:" << endl
         << " % Reduction in Dimension    = " << 100.0*(NRows-NReducedRows)/NRows << endl
         << " % Reduction in NNZs         = " << (OldMatrix_->NumGlobalNonzeros()-RRMatrix_->NumGlobalNonzeros())/100.0 << endl
         << "-------------------------------" << endl;
  }

  NewProblem_ = new Epetra_LinearProblem( RRMatrix_, RLHS_, RRHS_ );

  newObj_ = NewProblem_;

  cout << "done with SC\n";

  return *NewProblem_;
}

bool
LinearProblem_StaticCondensation::
fwd()
{
  if( verbose_ ) cout << "LP_SC::fwd()\n";
  if( verbose_ ) cout << "LP_SC::fwd() : Exporting to New LHS\n";
  ULHS_->Export( *OldLHS_, *LExporter_, Insert );
  RLHS_->Export( *OldLHS_, *RExporter_, Insert );
  LLHS_->Export( *OldLHS_, *UExporter_, Insert );

  if( verbose_ ) cout << "LP_SC::fwd() : Exporting to New RHS\n";
  URHS_->Export( *OldRHS_, *UExporter_, Insert );
  RRHS_->Export( *OldRHS_, *RExporter_, Insert );
  LRHS_->Export( *OldRHS_, *LExporter_, Insert );

  UUMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  URMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  ULMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  RRMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  RLMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  LLMatrix_->Export( *OldMatrix_, *LExporter_, Insert );

  Epetra_Vector LLDiag( *LMap_ );
  LLMatrix_->ExtractDiagonalCopy( LLDiag );
  Epetra_Vector LLRecipDiag( *LMap_ );
  LLRecipDiag.Reciprocal( LLDiag );

  if( verbose_ ) cout << "LP_SC::fwd() : Forming LLHS\n";
  LLDiag.Multiply( 1.0, LLRecipDiag, *LRHS_, 0.0 );
  int LSize = LMap_->NumMyElements();
  for( int i = 0; i < LSize; ++i ) (*LLHS_)[0][i] = LLDiag[i];

  if( verbose_ ) cout << "LP_SC::fwd() : Updating RRHS\n";
  Epetra_Vector RUpdate( *RMap_ );
  RLMatrix_->Multiply( false, *LLHS_, RUpdate );
  RRHS_->Update( -1.0, RUpdate, 1.0 );

  if( verbose_ ) cout << "LP_SC::fwd() : Updating URHS\n";
  Epetra_Vector UUpdate( *UMap_ );
  ULMatrix_->Multiply( false, *LLHS_, UUpdate );
  URHS_->Update( -1.0, UUpdate, 1.0 );
  
  return true;
}

bool
LinearProblem_StaticCondensation::
rvs()
{
  if( verbose_ ) cout << "LP_SC::rvs()\n";
  if( verbose_ ) cout << "LP_SC::rvs() : Updating URHS\n";
  Epetra_Vector UUpdate( *UMap_ );
  URMatrix_->Multiply( false, *RLHS_, UUpdate );
  URHS_->Update( -1.0, UUpdate, 1.0 );

  Epetra_Vector UUDiag( *UMap_ );
  UUMatrix_->ExtractDiagonalCopy( UUDiag );
  Epetra_Vector UURecipDiag( *UMap_ );
  UURecipDiag.Reciprocal( UUDiag );

  if( verbose_ ) cout << "LP_SC::rvs() : Forming ULHS\n";
  UUDiag.Multiply( 1.0, UURecipDiag, *URHS_, 0.0 );
  int USize = UMap_->NumMyElements();
  for( int i = 0; i < USize; ++i ) (*ULHS_)[0][i] = UUDiag[i];

  if( verbose_ ) cout << "LP_SC::rvs() : Importing to Old LHS\n";
  OldLHS_->Import( *ULHS_, *LExporter_, Insert );
  OldLHS_->Import( *RLHS_, *RExporter_, Insert );
  OldLHS_->Import( *LLHS_, *UExporter_, Insert );

  return true;
}

} // namespace EpetraExt

