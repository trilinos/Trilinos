//  #define DEBUG
//  #define USE_STL_SORT
#define USE_LOCAL

#ifndef USE_LOCAL
#ifndef USE_STL_SORT
At present, either USE_LOCAL or USE_STL_SORT is required
#endif
#endif

  /* Copyright (2003) Sandia Corportation. Under the terms of Contract 
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

#include "Amesos_Dscpack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#ifdef USE_STL_SORT
#include <algorithm>
#endif
#ifdef DEBUG
#include "Comm_assert_equal.h"
#endif

  //=============================================================================
  Amesos_Dscpack::Amesos_Dscpack(const Epetra_LinearProblem &prob, const AMESOS::Parameter::List &ParameterList ) {


  Problem_ = &prob ; 
  A_and_LU_built = false ; 
  FirstCallToSolve_ = true ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Dscpack::~Amesos_Dscpack(void) {

  if ( MyDscRank>=0 && A_and_LU_built ) { 
    //    cout << " At the top of the Amesos_Dscpack destructor" << endl ; 
    DSC_FreeAll( MyDSCObject ) ; 
    DSC_Close0( MyDSCObject ) ; 
    DSC_End( MyDSCObject ) ; 
    //    cout << " At the bottom of the Amesos_Dscpack destructor" << endl ; 
  }

  //    cout << " At the REAL bottom of the Amesos_Dscpack destructor" << endl ; 
}


int Amesos_Dscpack::PerformSymbolicFactorization() {
  bool factor = true; 

  vector <int> Replicates;
  vector <int> Ap;
  vector <int> Ai;

  bool CheckExtraction = false;    //  Set to true to force extraction for unit test

  Epetra_RowMatrix *RowMatrixA = 
    dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

#include "Epetra_CrsMatrix.h"

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
  Epetra_CrsMatrix *ExtractCrsMatrixA = 0;
#endif
  Epetra_CrsMatrix *Phase2Mat = 0 ;
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;
#if 0
  //
  //  The following lines allow us time to attach the debugger
  //
  int hatever;
  if ( iam == 0 )  cin >> hatever ; 
#endif
  Comm.Barrier();

  //
  //  Step 1)  Convert the matrix to an Epetra_CrsMatrix
  //
  //  If RowMatrixA is not a CrsMatrix, i.e. the dynamic cast fails, 
  //  extract a CrsMatrix from the RowMatrix.
  //
  if ( CastCrsMatrixA != 0 && ! CheckExtraction ) { 
    Phase2Mat = CastCrsMatrixA ; 
  } else {
    assert( false ) ;
  }


  const Epetra_Map &Phase2Matmap = Phase2Mat->RowMap() ; 

  //
  //  Step 2)  Coalesce the matrix onto process 0
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int IsLocal = ( Phase2Matmap.NumMyElements() == 
		  Phase2Matmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 

  Epetra_CrsMatrix *Phase3Mat = 0 ;

  int NumGlobalElements_ = Phase2Matmap.NumGlobalElements() ;
  //  Create a serial map in case we end up needing it 
  //  If it is created inside the else block below it would have to
  //  be with a call to new().
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;
  Epetra_Map SerialMap( NumGlobalElements_, NumMyElements_, 0, Comm );
  Epetra_CrsMatrix SerialCrsMatrixA(Copy, SerialMap, 0);


  if ( IsLocal==1 ) {
    Phase3Mat = Phase2Mat ;
  } else {

    Epetra_Export export_to_serial( Phase2Matmap, SerialMap);

    SerialCrsMatrixA.Export( *Phase2Mat, export_to_serial, Add ); 
    
    SerialCrsMatrixA.TransformToLocal() ; 
    Phase3Mat = &SerialCrsMatrixA ;

  }
  Comm.Barrier() ; 


  int numrows = Phase3Mat->NumGlobalRows();
  int numentries = Phase3Mat->NumGlobalNonzeros();

  assert( numrows == Phase3Mat->NumGlobalCols() );
  

  //
  //  Step 4)  Create a replicated map and matrix (someday we won't need this)
  //
  int * AllIDs = new int[numrows];
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  const Epetra_Map &Phase3Matmap = Phase3Mat->RowMap() ; 
  Epetra_Map ReplicatedMap( -1, numrows, AllIDs, 0, Comm);
  Epetra_Import importer( ReplicatedMap, Phase3Matmap );

  int nArows = Phase3Mat->NumGlobalRows() ; 
  int nAcols = Phase3Mat->NumGlobalCols() ; 

  Epetra_Import ImportToReplicated( ReplicatedMap, Phase2Matmap);

  Epetra_CrsMatrix Phase5Mat(Copy, ReplicatedMap, 0);
  EPETRA_CHK_ERR( Phase5Mat.Import( *Phase3Mat, importer, Insert) );
  EPETRA_CHK_ERR( Phase5Mat.TransformToLocal() ) ; 

  if ( factor ) { 
    //
    //  Step 6) Convert the matrix to Ap, Ai
    //
    Replicates.resize( numrows );
    for( int i = 0 ; i < numrows; i++ ) Replicates[i] = 1; 
    Ap.resize( numrows+1 );
    Ai.resize( EPETRA_MAX( numrows, numentries) ) ; 

    int NumEntriesPerRow ;
    double *RowValues = 0 ;
    int *ColIndices = 0 ;
    int Ai_index = 0 ; 
    int MyRow = -13 ;
    for ( MyRow = 0; MyRow <numrows; MyRow++ ) {
      int status = Phase5Mat.ExtractMyRowView( MyRow, NumEntriesPerRow, RowValues, ColIndices ) ;
      assert( status == 0 ) ; 
      Ap[MyRow] = Ai_index ; 
      for ( int j = 0; j < NumEntriesPerRow; j++ ) { 
	Ai[Ai_index] = ColIndices[j] ; 
	Ai_index++;
      }
    }
    assert( numrows == MyRow );
    assert( Ai_index == numentries ) ; 
    Ap[ numrows ] = Ai_index ; 
  }

  //
  //  Step 7)  Call Dscpack
  //  
  //  There is nothing special about -13.  It could be interpreted as meaning
  //  unitialized.
  //
  MyDSCObject = DSC_Begin() ; 
  int OrderCode = 2;
  vector<double> MyANonZ;
  int numprocs  = Comm.NumProc() ;            
  
  if ( factor ) { 
    DscNumProcs = -13 ; 
    MyDscRank = -13 ; 
    NumGlobalCols = -13 ; 
    NumLocalStructs = -13 ; 
    NumLocalCols = -13 ; 
    NumLocalNonz = 0 ; 
    GlobalStructNewColNum = 0 ; 
    GlobalStructNewNum = 0 ;  
    GlobalStructOwner = 0 ; 
    LocalStructOldNum = 0 ; 
    
    DscNumProcs = 1 ; 
    NumGlobalCols = 0 ; 
    assert( numprocs == Comm.NumProc() ) ; 

    int maxprocs = EPETRA_MIN( numprocs, 
			       DSC_Analyze( numrows, &Ap[0], &Ai[0], &Replicates[0] ) ) ; 
    while ( DscNumProcs * 2 <= maxprocs ) DscNumProcs *= 2 ;
    
    DSC_Open0( MyDSCObject, DscNumProcs, &MyDscRank, MPIC ) ; 

    NumLocalCols = 0 ; // This is for those processes not in the Dsc grid
    if ( MyDscRank >= 0 ) { 
      assert( iam == MyDscRank ) ; 
      EPETRA_CHK_ERR( DSC_Order ( MyDSCObject, OrderCode, numrows, &Ap[0], &Ai[0], 
				  &Replicates[0], &NumGlobalCols, &NumLocalStructs, 
				  &NumLocalCols, &NumLocalNonz, 
				  &GlobalStructNewColNum, &GlobalStructNewNum, 
				  &GlobalStructOwner, &LocalStructOldNum ) ) ; 
      assert( NumGlobalCols == numrows ) ; 
      assert( NumLocalCols == NumLocalStructs ) ; 
    }

    if ( MyDscRank >= 0 ) { 
      int TotalMemory, MaxSingleBlock; 

      const int Limit = 5000000 ;  //  Memory Limit set to 5 Terabytes 
      EPETRA_CHK_ERR( DSC_SFactor ( MyDSCObject, &TotalMemory, 
				    &MaxSingleBlock, Limit, DSC_LBLAS3, DSC_DBLAS2 ) ) ; 

    }        //     if ( MyDscRank >= 0 ) 

    A_and_LU_built = true; 
  } else {  // if ( factor)
    assert( numprocs == Comm.NumProc() ) ; 
  }  //End else if ( factor ) 

}

int Amesos_Dscpack::PerformNumericFactorization() {


  bool factor = true; 
  //
  //  I am going to put these here until I determine that I need them in 
  //  DscpackOO.h 
  //

  vector <int> Replicates;

  bool CheckExtraction = false;    //  Set to true to force extraction for unit test

  Epetra_RowMatrix *RowMatrixA = 
    dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

#include "Epetra_CrsMatrix.h"

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
  Epetra_CrsMatrix *ExtractCrsMatrixA = 0;
#endif
  Epetra_CrsMatrix *Phase2Mat = 0 ;
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;
#if 0
  //
  //  The following lines allow us time to attach the debugger
  //
  int hatever;
  if ( iam == 0 )  cin >> hatever ; 
#endif
  Comm.Barrier();

  //
  //  Step 1)  Convert the matrix to an Epetra_CrsMatrix
  //
  //  If RowMatrixA is not a CrsMatrix, i.e. the dynamic cast fails, 
  //  extract a CrsMatrix from the RowMatrix.
  //
  if ( CastCrsMatrixA != 0 && ! CheckExtraction ) { 
    Phase2Mat = CastCrsMatrixA ; 
  } else {
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
    ExtractCrsMatrixA = new Epetra_CrsMatrix( *RowMatrixA ) ; 

    Phase2Mat = ExtractCrsMatrixA ; 
#ifdef DEBUG
    //  DEBUG only works with DSCPACK1.0Ken and after you copy Comm_assert_equal.h over from the test directory
    if ( CheckExtraction ) 
      assert( CrsMatricesAreIdentical( CastCrsMatrixA, ExtractCrsMatrixA ) ) ; 
#endif
#else
    assert( false ) ;
#endif
  }


  const Epetra_Map &Phase2Matmap = Phase2Mat->RowMap() ; 

  //
  //  Step 2)  Coalesce the matrix onto process 0
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int IsLocal = ( Phase2Matmap.NumMyElements() == 
		  Phase2Matmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 
#ifdef DEBUG
  assert( Comm_assert_equal( &Comm, IsLocal ) );
#endif

  Epetra_CrsMatrix *Phase3Mat = 0 ;

  int NumGlobalElements_ = Phase2Matmap.NumGlobalElements() ;
  //  Create a serial map in case we end up needing it 
  //  If it is created inside the else block below it would have to
  //  be with a call to new().
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;
  Epetra_Map SerialMap( NumGlobalElements_, NumMyElements_, 0, Comm );
  Epetra_CrsMatrix SerialCrsMatrixA(Copy, SerialMap, 0);


  if ( IsLocal==1 ) {
    Phase3Mat = Phase2Mat ;
  } else {

    Epetra_Export export_to_serial( Phase2Matmap, SerialMap);

    SerialCrsMatrixA.Export( *Phase2Mat, export_to_serial, Add ); 
    
    SerialCrsMatrixA.TransformToLocal() ; 
    Phase3Mat = &SerialCrsMatrixA ;

  }
  Comm.Barrier() ; 


  int numrows = Phase3Mat->NumGlobalRows();
  int numentries = Phase3Mat->NumGlobalNonzeros();

  assert( numrows == Phase3Mat->NumGlobalCols() );
  

  //
  //  Step 4)  Create a replicated map and matrix (someday we won't need this)
  //
  int * AllIDs = new int[numrows];
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  const Epetra_Map &Phase3Matmap = Phase3Mat->RowMap() ; 
  Epetra_Map ReplicatedMap( -1, numrows, AllIDs, 0, Comm);
  Epetra_Import importer( ReplicatedMap, Phase3Matmap );
  int nArows = Phase3Mat->NumGlobalRows() ; 
  int nAcols = Phase3Mat->NumGlobalCols() ; 

  Epetra_Import ImportToReplicated( ReplicatedMap, Phase2Matmap);


  Epetra_CrsMatrix Phase5Mat(Copy, ReplicatedMap, 0);
  EPETRA_CHK_ERR( Phase5Mat.Import( *Phase3Mat, importer, Insert) );
  EPETRA_CHK_ERR( Phase5Mat.TransformToLocal() ) ; 

  //
  //  Step 7)  Call Dscpack
  //  
  int OrderCode = 2;
  vector<double> MyANonZ;
  int numprocs  = Comm.NumProc() ;            
  
  if ( factor ) { 

    assert( numprocs == Comm.NumProc() ) ; 

    Epetra_Map DscMap( numrows, NumLocalCols, LocalStructOldNum, 0, Comm ) ;

    //
    //  Import from the CrsMatrix (KEN GXX - can we go straight from the RowMatrix?)
    //
    Epetra_Import ImportToDsc( DscMap, Phase2Matmap );

    Epetra_CrsMatrix DscMat(Copy, DscMap, 0);
    EPETRA_CHK_ERR( DscMat.Import( *Phase2Mat, ImportToDsc, Insert) );
    EPETRA_CHK_ERR( DscMat.TransformToLocal() ) ; 
#if 0
    cout << endl << " Here is Phase2Mat : " << endl ; 
    Phase2Mat->Print( cout ) ; 
    cout << endl << " Here is DSCMAT : " << endl ; 
    DscMat.Print( cout ) ; 
#endif

    assert( MyDscRank >= 0 || NumLocalNonz == 0 ) ;
    assert( MyDscRank >= 0 || NumLocalCols == 0 ) ;
    assert( MyDscRank >= 0 || NumGlobalCols == 0  ) ; 
    MyANonZ.resize( NumLocalNonz ) ; 
    int NonZIndex = 0 ; 
    int num_my_row_entries  = -13 ; 

    int max_num_entries = DscMat.MaxNumEntries() ; 
    vector<int> col_indices( max_num_entries ) ; 
    vector<double> mat_values( max_num_entries ) ; 
    assert( NumLocalCols == DscMap.NumMyElements() ) ;
    vector<int> my_global_elements( NumLocalCols ) ; 
    EPETRA_CHK_ERR( DscMap.MyGlobalElements( &my_global_elements[0] ) ) ;

    vector<int> GlobalStructOldColNum( NumGlobalCols ) ; 
      
    typedef pair<int, double> Data;

    vector<Data> sort_array(NumGlobalCols);  // This is a gross 
    vector<int>  sort_indices(NumGlobalCols);  // This is a gross 
    // over-estimate of the max size of this array.  Ken work Fix this GXX


    for ( int i = 0; i < NumLocalCols ; i++ ) { 
      assert( my_global_elements[i] == LocalStructOldNum[i] ) ; 
      int num_entries_this_row; 
      //  USE_LOCAL and not USE_LOCAL both work
      //  #define USE_LOCAL
#ifdef USE_LOCAL
      EPETRA_CHK_ERR( DscMat.ExtractMyRowCopy( i, max_num_entries, num_entries_this_row, 
					       &mat_values[0], &col_indices[0] ) ) ; 
#else
      EPETRA_CHK_ERR( DscMat.ExtractGlobalRowCopy( DscMat.GRID(i), max_num_entries, num_entries_this_row, 
						   &mat_values[0], &col_indices[0] ) ) ; 
#endif
      //
      //  I would like to confirm num_entries against something provided by 
      //  Dscpack, but I can't tell what.
      //	assert ( num_entries_this_row ==   
      int OldRowNumber =  LocalStructOldNum[i] ;
      assert( GlobalStructOwner[ OldRowNumber ] != -1 ) ; 

      int NewRowNumber = GlobalStructNewColNum[ my_global_elements[ i ] ] ; 
      assert( numprocs > 1 || NewRowNumber == i ) ; 

      //
      //  Now we have to sort the column elements 
      //
      for ( int j = 0; j < num_entries_this_row; j++ ) { 
#ifdef USE_LOCAL
	sort_array[j].first = GlobalStructNewColNum[ DscMat.GCID( col_indices[j])] ; 
	sort_indices[j] =  GlobalStructNewColNum[ DscMat.GCID( col_indices[j])] ; 
#else
	sort_array[j].first = GlobalStructNewColNum[ col_indices[j] ]; 
#endif
	sort_array[j].second = mat_values[j] ; 
      }
#ifdef USE_STL_SORT
      sort(&sort_array[0], &sort_array[num_entries_this_row]);
#else
      double **DoubleCompanions = new double*[2] ;
      *DoubleCompanions = &mat_values[0] ; 
      Epetra_Util sorter;
      sorter.Sort( true, num_entries_this_row, &sort_indices[0],
		   1, DoubleCompanions, 0, 0 ) ;
      delete[] DoubleCompanions; 
#endif

      for ( int j = 0; j < num_entries_this_row; j++ ) { 
#ifdef USE_STL_SORT
	int NewColNumber = sort_array[j].first ; 
	if ( NewRowNumber <= NewColNumber ) MyANonZ[ NonZIndex++ ] = sort_array[j].second ; 
#else
	int NewColNumber = sort_indices[j] ; 
	if ( NewRowNumber <= NewColNumber ) MyANonZ[ NonZIndex++ ] = mat_values[j] ; 
#endif
#ifndef USE_LOCAL
	assert( NonZIndex <= NumLocalNonz );
#endif
      }
    }

    if ( MyDscRank >= 0 ) { 
      int TotalMemory, MaxSingleBlock; 


      const int SchemeCode = 1; 
#ifndef USE_LOCAL
      assert( NonZIndex == NumLocalNonz );
#endif

      EPETRA_CHK_ERR( DSC_NFactor ( MyDSCObject, SchemeCode, &MyANonZ[0], 
				    DSC_LLT,  DSC_LBLAS3, DSC_DBLAS2 ) ) ;

    }        //     if ( MyDscRank >= 0 ) 

  } else {  // if ( factor)
    assert( numprocs == Comm.NumProc() ) ; 
  }  //End else if ( factor ) 
}




bool Amesos_Dscpack::MatrixShapeOK() const { 
  bool OK =  GetProblem()->IsOperatorSymmetric() ;

  //
  //  The following test is redundant.  I have left it here in case the 
  //  IsOperatorSymmetric test turns out not to be reliable.
  //
  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Dscpack::SymbolicFactorization() {
  assert( false ) ;   // Not implemented yet ;
}

int Amesos_Dscpack::NumericFactorization() {
  assert( false ) ;   // Not implemented yet ;
}

//
//  Solve() uses several intermediate matrices to convert the input matrix
//  to one that we can pass to the Sparse Direct Solver
//
//  Epetra_RowMatrix *RowMatrixA - The input matrix
//
int Amesos_Dscpack::Solve() { 
  //  int Solve() { 

  PerformSymbolicFactorization();

  PerformNumericFactorization();

  Epetra_RowMatrix *RowMatrixA = 
    dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
  Epetra_CrsMatrix *ExtractCrsMatrixA = 0;
#endif
  Epetra_CrsMatrix *Phase2Mat = 0 ;
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;
  Comm.Barrier();

  //
  //  Step 1)  Convert the matrix to an Epetra_CrsMatrix
  //
  //  If RowMatrixA is not a CrsMatrix, i.e. the dynamic cast fails, 
  //  extract a CrsMatrix from the RowMatrix.
  //
  if ( CastCrsMatrixA != 0  ) { 
    Phase2Mat = CastCrsMatrixA ; 
  } else {
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
    ExtractCrsMatrixA = new Epetra_CrsMatrix( *RowMatrixA ) ; 

    Phase2Mat = ExtractCrsMatrixA ; 
#else
    assert( false ) ;
#endif
  }


  const Epetra_Map &Phase2Matmap = Phase2Mat->RowMap() ; 

  //
  //  Step 2)  Coalesce the matrix onto process 0
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int IsLocal = ( Phase2Matmap.NumMyElements() == 
		  Phase2Matmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 
#ifdef DEBUG
  assert( Comm_assert_equal( &Comm, IsLocal ) );
#endif

  Epetra_CrsMatrix *Phase3Mat = 0 ;

  int NumGlobalElements_ = Phase2Matmap.NumGlobalElements() ;
  //  Create a serial map in case we end up needing it 
  //  If it is created inside the else block below it would have to
  //  be with a call to new().
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;
  Epetra_Map SerialMap( NumGlobalElements_, NumMyElements_, 0, Comm );
  Epetra_CrsMatrix SerialCrsMatrixA(Copy, SerialMap, 0);


  if ( IsLocal==1 ) {
    Phase3Mat = Phase2Mat ;
  } else {

    Epetra_Export export_to_serial( Phase2Matmap, SerialMap);

    SerialCrsMatrixA.Export( *Phase2Mat, export_to_serial, Add ); 
    
    SerialCrsMatrixA.TransformToLocal() ; 
    Phase3Mat = &SerialCrsMatrixA ;

  }
  Comm.Barrier() ; 


  int numrows = Phase3Mat->NumGlobalRows();
  int numentries = Phase3Mat->NumGlobalNonzeros();

  assert( numrows == Phase3Mat->NumGlobalCols() );
  

  //
  //  Step 4)  Create a replicated map and matrix (someday we won't need this)
  //
  int * AllIDs = new int[numrows];
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  const Epetra_Map &Phase3Matmap = Phase3Mat->RowMap() ; 
  Epetra_Map ReplicatedMap( -1, numrows, AllIDs, 0, Comm);
  Epetra_Import importer( ReplicatedMap, Phase3Matmap );
  //
  //  Step 5)  Convert vector b to a replicated vector
  //
  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  int nArows = Phase3Mat->NumGlobalRows() ; 
  int nAcols = Phase3Mat->NumGlobalCols() ; 

#ifdef ONE_VECTOR_ONLY
  assert( vecX->NumVectors() == 1 ) ; 
  assert( vecB->NumVectors() == 1 ) ; 

  Epetra_Vector *vecXvector = dynamic_cast<Epetra_Vector*>(vecX) ; 
  Epetra_Vector *vecBvector = dynamic_cast<Epetra_Vector*>(vecB) ; 

  assert( vecXvector != 0 ) ; 
  assert( vecBvector != 0 ) ; 

  Epetra_Vector vecXreplicated( ReplicatedMap ) ; 
  Epetra_Vector vecBreplicated( ReplicatedMap ) ; 
#else
  Epetra_MultiVector *vecXvector = (vecX) ; 
  Epetra_MultiVector *vecBvector = (vecB) ; 

  Epetra_MultiVector vecXreplicated( ReplicatedMap, nrhs ) ; 
  Epetra_MultiVector vecBreplicated( ReplicatedMap, nrhs ) ; 
#endif


  Epetra_Import ImportToReplicated( ReplicatedMap, Phase2Matmap);

  vecXreplicated.Import( *vecXvector, ImportToReplicated, Insert ) ;
  vecBreplicated.Import( *vecBvector, ImportToReplicated, Insert ) ;

  assert( nArows == vecXreplicated.MyLength() ) ; 
  assert( nAcols == vecBreplicated.MyLength() ) ;

  double *bValues ;
  double *xValues ;
  int bLda, xLda ; 

  assert( vecBreplicated.ExtractView( &bValues, &bLda ) == 0 )  ; 
  assert( vecXreplicated.ExtractView( &xValues, &xLda ) == 0 ) ; 

  Epetra_CrsMatrix Phase5Mat(Copy, ReplicatedMap, 0);
  EPETRA_CHK_ERR( Phase5Mat.Import( *Phase3Mat, importer, Insert) );
  EPETRA_CHK_ERR( Phase5Mat.TransformToLocal() ) ; 

  //
  //  Step 7)  Call Dscpack
  //  
  int OrderCode = 2;
  vector<double> MyANonZ;
  int numprocs  = Comm.NumProc() ;            
  
  int ldb = numrows ; 

  //  amesos_test.exe DSCPACK Diagonal.mtx 0 0 -1 0 1e-15 1e-15 still dumps even if we disable the following code block
  if ( MyDscRank >= 0 ) {
    EPETRA_CHK_ERR( DSC_InputRhsGlobalVec ( MyDSCObject, bValues, numrows ) ) ;
    EPETRA_CHK_ERR( DSC_Solve ( MyDSCObject ) ) ; 
    vector<double> Dsc_outputs( numrows ) ; 
    vector<int> Dsc_indices( numrows) ;
    EPETRA_CHK_ERR( DSC_GetGlobalSolution ( MyDSCObject, &Dsc_indices[0], &Dsc_outputs[0] ) ) ; 
    
    for ( int i =0 ; i<numrows; i++ ) 
      xValues[Dsc_indices[i]] = Dsc_outputs[i] ; 
  }

  if ( iam == 0 ) assert( MyDscRank >= 0 ) ; // Make sure that process 0 has valid data for xValues

  Comm.Broadcast( &xValues[0], numrows, 0 ) ; 

  //
  //  Step 8)  Convert vector x back to a distributed vector
  //
  //  This is an ugly hack - it should be cleaned up someday
  //
  for (int i = 0 ; i < numrows; i++ ) { 
    int lid[1000] ; 
    lid[0] = Phase2Matmap.LID( i ) ; 
    if ( lid[0] >= 0 ) { 
      for (int j = 0 ; j < nrhs; j++ ) { 
	vecXvector->ReplaceMyValue(  lid[0], j, xValues[i + ldb * j]  ) ; 
      }
    }
  }
  return(0) ; 
}
