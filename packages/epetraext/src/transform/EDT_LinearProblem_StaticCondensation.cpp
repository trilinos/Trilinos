
#include <EDT_LinearProblem_StaticCondensation.h>

#include <vector>
#include <map>
#include <set>

#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

EpetraExt::LinearProblem_StaticCondensation::~LinearProblem_StaticCondensation()
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

EpetraExt::LinearProblem_StaticCondensation::NewTypePtr EpetraExt::LinearProblem_StaticCondensation::operator()
	( EpetraExt::LinearProblem_StaticCondensation::OriginalTypeRef original )
{
  int ierr = 0;

  OldMatrix_ = dynamic_cast<Epetra_CrsMatrix*>( original.GetMatrix() );
  OldGraph_  = &OldMatrix_->Graph();
  OldRHS_    = original.GetRHS();
  OldLHS_    = original.GetLHS();
  OldRowMap_ = &OldMatrix_->RowMap();

  const Epetra_Comm & CommObj = OldRowMap_->Comm();

  if( !OldMatrix_ ) ierr = -2;
  if( !OldRHS_ )    ierr = -3;
  if( !OldLHS_ )    ierr = -4;

  if( OldRowMap_->DistributedGlobal() ) ierr = -5;
  if( degree_ != 1 ) ierr = -6;

  int NRows = OldGraph_->NumMyRows();
  int IndexBase = OldRowMap_->IndexBase();

  std::vector<int> ColNZCnt( NRows );
  std::vector<int> CS_RowIndices( NRows );
  std::map<int,int> RS_Map;
  std::map<int,int> CS_Map;

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
    for( std::map<int,int>::iterator itM = RS_Map.begin(); itM != RS_Map.end(); ++itM )
      cout << itM->first << "\t" << itM->second << endl;
    cout << "Col Counts\n";
    for( int i = 0; i < NRows; ++i )
      cout << i << "\t" << ColNZCnt[i] << "\t" << CS_RowIndices[i] << endl;
    cout << "-------------------------\n";
  }

  std::set<int> RS_Set;
  std::set<int> CS_Set;

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
    std::set<int>::iterator itRS = RS_Set.begin();
    std::set<int>::iterator itCS = CS_Set.begin();
    std::set<int>::iterator endRS = RS_Set.end();
    std::set<int>::iterator endCS = CS_Set.end();
    for( ; itRS != endRS; ++itRS, ++itCS )
      cout << *itRS << "\t" << *itCS << endl;
    cout << "-------------------------\n";
  }

  //build new maps
  int NSingletons = CS_Set.size();
  int NReducedRows = NRows - 2* NSingletons; 
  std::vector<int> ReducedIndices( NReducedRows );
  std::vector<int> CSIndices( NSingletons );
  std::vector<int> RSIndices( NSingletons );
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

  std::vector<int> RowPermutedIndices( NRows );
  copy( RSIndices.begin(), RSIndices.end(), RowPermutedIndices.begin() );
  copy( ReducedIndices.begin(), ReducedIndices.end(), RowPermutedIndices.begin() + NSingletons );
  copy( CSIndices.begin(), CSIndices.end(), RowPermutedIndices.begin() + NReducedRows + NSingletons );

  std::vector<int> ColPermutedIndices( NRows );
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
  NewGraph_->TransformToLocal();
  cout << "Permuted Graph:\n" << *NewGraph_;

  NewMatrix_ = new Epetra_CrsMatrix( Copy, *NewGraph_ );
  NewMatrix_->Export( *OldMatrix_, *Exporter_, Insert );
  NewMatrix_->TransformToLocal();
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
  UUGraph_->TransformToLocal();
  cout << "UUGraph:\n" << *UUGraph_;

  UUMatrix_ = new Epetra_CrsMatrix( Copy, *UUGraph_ );
  UUMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  UUMatrix_->TransformToLocal();
  cout << "UUMatrix:\n" << *UUMatrix_;

  URGraph_ = new Epetra_CrsGraph( Copy, *UMap_, *RMap_, 0 );
  URGraph_->Export( *OldGraph_, *UExporter_, Insert );
  URGraph_->TransformToLocal();
  cout << "URGraph:\n" << *URGraph_;

  URMatrix_ = new Epetra_CrsMatrix( Copy, *URGraph_ );
  URMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  URMatrix_->TransformToLocal();
  cout << "URMatrix:\n" << *URMatrix_;

  ULGraph_ = new Epetra_CrsGraph( Copy, *UMap_, *UMap_, 0 );
  ULGraph_->Export( *OldGraph_, *UExporter_, Insert );
  ULGraph_->TransformToLocal();
  cout << "ULGraph:\n" << *ULGraph_;

  ULMatrix_ = new Epetra_CrsMatrix( Copy, *ULGraph_ );
  ULMatrix_->Export( *OldMatrix_, *UExporter_, Insert );
  ULMatrix_->TransformToLocal();
  cout << "ULMatrix:\n" << *ULMatrix_;

  RRGraph_ = new Epetra_CrsGraph( Copy, *RMap_, *RMap_, 0 );
  RRGraph_->Export( *OldGraph_, *RExporter_, Insert );
  RRGraph_->TransformToLocal();
  cout << "RRGraph:\n" << *RRGraph_;

  RRMatrix_ = new Epetra_CrsMatrix( Copy, *RRGraph_ );
  RRMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  RRMatrix_->TransformToLocal();
  cout << "RRMatrix:\n" << *RRMatrix_;

  RLGraph_ = new Epetra_CrsGraph( Copy, *RMap_, *UMap_, 0 );
  RLGraph_->Export( *OldGraph_, *RExporter_, Insert );
  RLGraph_->TransformToLocal();
  cout << "RLGraph:\n" << *RLGraph_;

  RLMatrix_ = new Epetra_CrsMatrix( Copy, *RLGraph_ );
  RLMatrix_->Export( *OldMatrix_, *RExporter_, Insert );
  RLMatrix_->TransformToLocal();
  cout << "RLMatrix:\n" << *RLMatrix_;

  LLGraph_ = new Epetra_CrsGraph( Copy, *LMap_, *UMap_, 0 );
  LLGraph_->Export( *OldGraph_, *LExporter_, Insert );
  LLGraph_->TransformToLocal();
  cout << "LLGraph:\n" << *LLGraph_;

  LLMatrix_ = new Epetra_CrsMatrix( Copy, *LLGraph_ );
  LLMatrix_->Export( *OldMatrix_, *LExporter_, Insert );
  LLMatrix_->TransformToLocal();
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

  return NewProblem_;
}

bool EpetraExt::LinearProblem_StaticCondensation::fwd()
{
  ULHS_->Export( *OldLHS_, *LExporter_, Insert );
  RLHS_->Export( *OldLHS_, *RExporter_, Insert );
  LLHS_->Export( *OldLHS_, *UExporter_, Insert );

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

  LLDiag.Multiply( 'N', 'N', 1.0, LLRecipDiag, *LRHS_, 0.0 );
  int LSize = LMap_->NumMyElements();
  for( int i = 0; i < LSize; ++i ) (*LLHS_)[0][i] = LLDiag[i];

  Epetra_Vector RUpdate( *RMap_ );
  RLMatrix_->Multiply( false, *LLHS_, RUpdate );
  RRHS_->Update( -1.0, RUpdate, 1.0 );

  Epetra_Vector UUpdate( *UMap_ );
  ULMatrix_->Multiply( false, *LLHS_, UUpdate );
  URHS_->Update( -1.0, UUpdate, 1.0 );
  
  return true;
}

bool EpetraExt::LinearProblem_StaticCondensation::rvs()
{
  Epetra_Vector UUpdate( *UMap_ );
  URMatrix_->Multiply( false, *RLHS_, UUpdate );
  URHS_->Update( -1.0, UUpdate, 1.0 );

  Epetra_Vector UUDiag( *UMap_ );
  UUMatrix_->ExtractDiagonalCopy( UUDiag );
  Epetra_Vector UURecipDiag( *UMap_ );
  UURecipDiag.Reciprocal( UUDiag );

  UUDiag.Multiply( 'N', 'N', 1.0, UURecipDiag, *URHS_, 0.0 );
  int USize = UMap_->NumMyElements();
  for( int i = 0; i < USize; ++i ) (*ULHS_)[0][i] = UUDiag[i];

  OldLHS_->Import( *ULHS_, *LExporter_, Insert );
  OldLHS_->Import( *RLHS_, *RExporter_, Insert );
  OldLHS_->Import( *LLHS_, *UExporter_, Insert );

  return true;
}

