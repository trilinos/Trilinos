
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

namespace Epetra_Transform {

LinearProblem_StaticCondensation::~LinearProblem_StaticCondensation()
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

std::auto_ptr<Epetra_LinearProblem> LinearProblem_StaticCondensation::operator()
	( const Epetra_LinearProblem & original )
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

  NewProblem_ = new Epetra_LinearProblem( NewMatrix_, NewLHS_, NewRHS_ );

  UExporter_ = new Epetra_Export( *OldRowMap_, *UMap_ );
  RExporter_ = new Epetra_Export( *OldRowMap_, *RMap_ );
  LExporter_ = new Epetra_Export( *OldRowMap_, *LMap_ );
cout << "UExporter:\n" << *UExporter_;
cout << "RExporter:\n" << *RExporter_;
cout << "LExporter:\n" << *LExporter_;

  ULHS_ = new Epetra_MultiVector( *LMap_, OldLHS_->NumVectors() );
  ULHS_->Export( *OldLHS_, *UExporter_, Insert );
  cout << "ULHS:\n" << *ULHS_;

  RLHS_ = new Epetra_MultiVector( *RMap_, OldLHS_->NumVectors() );
  RLHS_->Export( *OldLHS_, *RExporter_, Insert );
  cout << "RLHS:\n" << *RLHS_;

  LLHS_ = new Epetra_MultiVector( *UMap_, OldLHS_->NumVectors() );
  LLHS_->Export( *OldLHS_, *LExporter_, Insert );
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
//         << " NNZs                        = " << ReducedMatrix->NumGlobalNonzeros() << endl
//         << " Max Num Row Entries         = " << ReducedMatrix->GlobalMaxNumEntries() << endl
         << "Singleton Info:" << endl
         << " Num Singleton                 = " << NSingletons << endl
         << "Ratios:" << endl
         << " % Reduction in Dimension    = " << 100.0*(NRows-NReducedRows)/NRows << endl
//         << " % Reduction in NNZs         = " << (Matrix.NumGlobalNonzeros()-ReducedMatrix->NumGlobalNonzeros())/100.0 << endl
         << "-------------------------------" << endl;
  }

  return std::auto_ptr<Epetra_LinearProblem>(NewProblem_);

#if 0
  Epetra_CrsMatrix * oMatrix = dynamic_cast<Epetra_CrsMatrix*>( original.GetMatrix() );
  if     ( !oMatrix )           EPETRA_CHK_ERR(-2);
  else if( !original.GetRHS() ) EPETRA_CHK_ERR(-3);
  else if( !original.GetLHS() ) EPETRA_CHK_ERR(-4);

  Epetra_Map & RowMap    = original.RowMap();
  Epetra_Map & DomainMap = original.DomainMap();
  Epetra_Map & ColumnMap = original.ColumnMap();

  Epetra_CrsMatrix & Matrix = *oMatrix;
  Epetra_CrsGraph & Graph = Matrix.Graph();

  Epetra_IntVector LocalColProfiles( ColumnMap );
  Epetra_IntVector ColProfiles( DomainMap );
  Epetra_IntVector LocalColSingletonRows( ColumnMap );
  Epetra_IntVector ColSingletonRows( DomainMap );

  int NRows = Graph.NumMyRows();
  int NCols = Graph..NumMyCols();
  int IndexBase = RowMap.IndexBase();

  vector<int> LocalRowIDs( NCols, -1 );
  vector<int> EliminatedRowIDs( NRows );
  vector<int> EliminatedColIDs( NCols );
  int NEliminatedRows = 0;
  int NEliminatedCols = 0;

  int NSingletons = 0;

  //Phase 1: Collect local singleton row and column info

  int NumIndices;
  int * Indices;
  int ColIndex;
  for( int i = 0; i < NRows; ++i )
  {
    EPETRA_CHK_ERR( Graph.ExtractMyRowView( i, NumIndices, Indices ) );

    for( int j = 0; j < NumIndices; ++j )
    {
      ColIndex = Indices[j];
      ++LocalColProfiles[ ColIndex ]; //Inc col cnt

      //Record row id for current col, to be used for elimination of singleton col
      LocalRowIDs[ ColIndex ] = i;
    }

    //Check for singletion rows and add 1 to cols to be eliminated
    if( NumIndices == 1 )
    {
      ColIndex = Indices[0];
      ++LocalColSingletonRows[ ColIndex ];
      EliminatedRowIDs[ NEliminatedRows++ ] = RowMap.GID(i);
      EliminatedColIDs[ NEliminatedCols++ ] = ColumnMap.GID( ColIndex );

      ++NSingletons;
    }
  }

  //Phase 2: Generate global singleton info
  // -Gather column entrie counts to determine singletons for global system
  // -Distribute column and row singleton info

  // Make copy fo LocalColProfiles for use in detecting cols that disappear locally
  Epetra_IntVector LocalColProfilesCopy( LocalColProfiles );

  if( ColumnMap.DistributedGlobal() )
  {
    EPETRA_CHK_ERR( ColProfiles.Export( LocalColProfiles, Matrix.Importer(), Add ) );
    EPETRA_CHK_ERR( LocalColProfiles.Import( ColProfiles, Matrix.Importer(), Insert ) );
    EPETRA_CHK_ERR( ColSingletonRows.Export( LocalColSingletonRows, Matrix.Importer(), Add ) );
    EPETRA_CHK_ERR( LocalColSingletonRows.Import( ColSingletonRows, Matrix.Importer(), Insert ) );
  }

  vector<bool> RowsWithSingletonCols( NRows, false );

  int NColSingletons = 0;

  //Count singleton cols not already counted with singleton rows
  for( int j = 0; j < NCols; ++j )
  {
    int i = LocalRowIDs[j];
    int ColID = ColumnMap.GID(j);

    if( LocalColProfiles[j] == 1 )
    {
      if( i == -1 ) //if no row entry locally, just eliminate column
        EliminatedColIDs[ NEliminatedCols++ ] = ColID;
      else //Check if column already eliminated by row check above
      {
        if( Matrix.NumMyEntries(i) != 1 )
        {
          if( RowsWithColSingletons[i] ) EPETRA_CHK_ERR(-2); //Row had 2 column singletons, GAK!
          RowsWithColSingletons[i] = true;
          EliminatedRowIDs[ NEliminatedRows++ ] = RowMap.GID(i);
          EliminatedColIDs[ NEliminatedCols++ ] = ColID;
          ++NColSingletons;

          //Keep track of columns associated with a deleted row
          //If eventually all column entries are deleted, the column should also be deleted
          EPETRA_CHK_ERR( Graph.ExtractMyRowView( i, NumIndices, Indices ) );
          for( int jj = 0; jj < NumIndices; ++jj ) --LocalColProfilesCopy[ Indices[jj] ];
        }
      }
    }
    else if( LocalColSingletonRows[j] == 1 ) //Check if other proc eliminated this column
    {
      if( i != -1 ) //If no local row entry, just eliminate column
      {
        if( Matrix.NumMyEntries(i) != 1 && !RowMap.MyLID(j) )
          EliminatedColIDs[ NEliminatedCols++ ] = ColID;
      }
    }
    else if( LocalColSingletonRows[j] > 1 )
      EPETRA_CHK_ERR( -3 ); //Column has 2 row singletons, GAK!
  }

  //Phase 3: Setup arrays for construction
  vector<int> ColSingletonRowLIDs( NColSingletons );
  vector<int> ColSingletonColLIDs( NColSingletons );
  vector<double> ColSingletonPivots( NColSingletons );

  //Register singleton columns not already counted as singleton rows
  //Check if columns disappear due to all assoc. rows being eliminated
  int NColSingletonsTmp = 0;

    if( LocalColProfiles[j] == 1 && LocalRowIDs[j] != -1 )
    {
      if( Matrix.NumMyEntries( LocalRowIDs[j] ) != 1 )
      {
        ColSingletonRowLIDs[ NColSingletonTmp ] = LocalRowIDs[j];
        ColSingletonColLIDs[ NColSingletonTmp++ ] = Matrix.NumMyEntries( LocalRowIDs[j] );
      }
    }
    else if( !LocalColProfilesCopy[j] && LocalColSingletonRows[j] != 1 && LocalRowIDs[j] != -1 )
    // if column has no row entry on this proc, just eliminate
    { 
      if( Matrix.NumMyEntries( LocalRowIDs[j] ) != 1 )
        EliminatedColIDs[ NEliminatedCols++ ] = ColumnMap.GID(j);
    }

  assert( NumColSingletonsTmp == NumColSingletons ); //Sanity Check

  //Sort Row LIDs
  int * CSCLTmp = &ColSingletonColLIDs[0];
  Epetra_Util().Sort( true, NColSingletons_, &ColSingletonRowLIDs[0], 0, 0, 1, &CSCLTmp );

  Epetra_Map RowEliminateMap( -1, NEliminateRows, &EliminatedRowIDs[0], IndexBase, ColumnMap.Comm() );
  Epetra_Map ColEliminateMap( -1, NEliminateCols, &EliminatedColIDs[0], IndexBase, ColumnMap.Comm() );

  int NReducedRows = NRows - NEliminatedRows;
  int NReducedCols = NCols - NEliminatedCols;

  vector<int> ReducedRows( NReducedRows );
  vector<int> ReducedCols( NReducedCols );

  int NReducedRowsTmp = 0;
  int NReducedColsTmp = 0;

  //Register reduced rows by determining if each row of matrix is in RowEliminateMap
  for( int i = 0; i < NRows; ++i )
    if( !RowEliminateMap.MyGID( RowMap.GID(i) ) )
      ReducedRows[NReducedRowsTmp++] = RowMap.GID(i);

  //Register reduced cols by determining if each col of matrix is in ColEliminateMap
  for( int i = 0; i < NCols; ++i )
    if( !ColEliminateMap.MyGID( ColumnMap.GID(i) ) )
      ReducedCols[NReducedColsTmp++] = ColumnMap.GID(i);

  assert( NReducedRowsTmp == NReducedRows );
  assert( NReducedColsTmp == NReducedCols );

  //Cannot handle nonsymmetric elimination
  for( int i = 0; i < NEliminatedRows; ++i )
    if( !ColEliminateMap.MyGID( RowEliminateMap.GID(i) ) ) EPETRA_CHK_ERR(-4);

  //Construct Reduced Matrix Maps
  Epetra_Map * ReducedRowMap
    = Epetra_Map( -1, NReducedRows, &ReducedRows[0], IndexBase, ColumnMap.Comm() );
  Epetra_Map * ReducedDomainMap
    = Epetra_Map( -1, NReducedRows, &ReducedRows[0], IndexBase, ColumnMap.Comm() );

  Epetra_Import LHSImporter( *ReducedDomainMap, DomainMap );
  Epetra_Import RHSImporter( *ReducedRowMap, RowMap );

  //Construction of Reduced Problem
  if( ReducedMatrix ) delete ReducedMatrix;
  Epetra_CrsMatrix * ReducedMatrix = new Epetra_CrsMatrix( Copy, *ReducedRowMap, 0 );

  //Temporary for X values due to explicit elimination
  Epetra_MultiVector ExportX( ColumnMap, original.GetLHS()->NumVectors() );

  //Access to RHS and LHS
  Epetra_MultiVector * rhs = original.GetRHS();
  Epetra_MultiVector * lhs = original.GetLHS();

  int MaxNumEntries = Matrix.MaxNumEntries();
  vector<int> ReducedIndices( MaxNumEntries );
  vector<double> ReducedValues( MaxNumEntries );

  int ColSingletonCnt = 0;
  double Values;
  for( int i = 0; i < NRows; ++i )
  {
    int cGRID = RowMap.GID(i);
    EPETRA_CHK_ERR( Matrix.ExtractMyRowView( i, NumEntries, Values, Indices ) );
    if( ReducedRowMap->MyGID( cGRID ) )
    {
      int ReducedNumEntries = 0;
      for( int j = 0; j < NumEntries; ++j )
      {
        int ColIndex = ColumnMap.GID( Indices[j] );
        if( !ColEliminateMap.MyGID( ColIndex ) )
        {
          ReducedIndices[ ReducedNumEntries ] = ColIndex;
          ReducedValues[ ReducedNumEntries++ ] = Values[j];
        }
      }
      EPETRA_CHK_ERR( ReducedMatrix->InsertGlobalValues( cGRID, ReducedNumEntries, &ReducedValues[0], &ReducedIndices[0] ) );
    }
    else if( NumEntries == 1 ) //singleton row, explicit elimination
    {
      double pivot = Values[0];
      if( !pivot ) EPETRA_CHK_ERR(-1);
      int indX = Indices[0];
      int NVecs = ExportX.NumVectors();
      for( int j = 0; j < NVecs; ++j ) ExportX[j][indX] = (*rhs)[j][i]/pivot;
    }
    else //singleton column, setup for post-solve step
    {
      assert( i == ColSingletonROWLIDs[ColSingletonCnt] ); //Sanity Check
      int targetCol = ColSingletonColLIDs[ColSingletonCnt];
      for( int j = 0; j < NumEntries; ++j )
        if( Indices[j] == targetCol )
        {
          double pivot = Values[j];
          if( !pivot ) EPETRA_CHK_ERR(-2);
          ColSingletonPivots_[ColSingletonCnt++] = pivot;
          break;
        }
    }
  }

  assert( ColSingletonCnt == NColSingletons );  //Sanity Check

  EPETRA_CHK_ERR( ReducedMatrix->TransformToLocal( ReducedDomainMap, ReducedRowMap );

  //Contruct Reduced LHS (Import from Original LHS)
  if( ReducedLHS ) delete ReducedLHS;
  ReducedLHS = new Epetra_MultiVector( *ReducedDomainMap, NumVectors );
  EPETRA_CHK_ERR( ReducedLHS->Import( *lhs, *LHSImporter, Insert ) );
  lhs->PutScalar( 0.0 ); //Clear LHS for insert from reduced solve
  
  //Construct Reduced RHS (Include Explicit Eliminaton of Rows)
  Epetra_MultiVector tmpX( DomainMap, NumVectors );
  Epetra_MultiVector tmpB( RowMap, NumVectors );

  EPETRA_CHK_ERR( tmpX.Export( ExportX, Importer, Add ) );
  EPETRA_CHK_ERR( lhs->Export( ExportX, Importer, Add ) );

  EPETRA_CHK_ERR( Matrix->Multiply( false, tmpX, tmpB ) );

  EPETRA_CHK_ERR( tmpB.Update( 1.0, *rhs, -1.0 ) );

  if( ReducedRHS ) delete ReducedRHS;
  ReducedRHS = new Epetra_MultiVector( *ReducedRowMap, NumVectors );
  EPETRA_CHK_ERR( ReducedRHS->Import( tmpB, *RHSImporter, Insert ) );

  if( ReducedProblem ) delete ReducedProblem;
  ReducedProblem = new Epetra_LinearProblem( ReducedMatrix, ReducedLHS, ReducedRHS );
  
  if( vebose_ )
  {
    cout << "Full System Characteristics:" << endl
         << "----------------------------" << endl
         << " Dimension                   = " << Matrix.NumGlobalRows() << endl
         << " NNZs                        = " << Matrix.NumGlobalNonzeros() << endl
         << " Max Num Row Entries         = " << Matrix.GlobaMaxNumEntries() << endl << endl
         << "Reduced System Characteristics:" << endl
         << " Dimension                   = " << ReducedMatrix->NumGlobalRows() << endl
         << " NNZs                        = " << ReducedMatrix->NumGlobalNonzeros() << endl
         << " Max Num Row Entries         = " << ReducedMatrix->GlobalMaxNumEntries() << endl
         << "Singleton Info:" << endl
         << " Num Singleton Rows          = " << NRowSingletons << endl
         << " Num Singleton Cols          = " << NColSingletons << endl
         << "Ratios:" << endl
         << " % Reduction in Dimension    = " << (Matrix.NumGlobalRows()-ReducedMatrix->NumGlobalRows())/100.0 << endl;
         << " % Reduction in NNZs         = " << (Matrix.NumGlobalNonzeros()-ReducedMatrix->NumGlobalNonzeros())/100.0 << endl;
         << "-------------------------------" << endl;
  }
#endif

}

bool LinearProblem_StaticCondensation::fwd()
{
  NewMatrix_->Export( *OldMatrix_, *Exporter_, Insert );
  NewRHS_->Export( *OldRHS_, *Exporter_, Insert );
  NewLHS_->Export( *OldLHS_, *Exporter_, Insert );

  return true;
}

bool LinearProblem_StaticCondensation::rvs()
{
  Epetra_MultiVector * tmpLHS = const_cast<Epetra_MultiVector*>(OldLHS_);
  tmpLHS->Import( *NewLHS_, *Exporter_, Insert );

  return true;
}

}
