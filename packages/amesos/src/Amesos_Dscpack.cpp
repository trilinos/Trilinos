// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//  As of July 1st, USE_STL_SORT and USE_LOCAL both work together or separately
//  (But you have to set at least one)
//  #define USE_STL_SORT
#define USE_LOCAL

#ifndef USE_LOCAL
#ifndef USE_STL_SORT
At present, either USE_LOCAL or USE_STL_SORT is required
#endif
#endif

#include "Amesos_Dscpack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#ifdef USE_STL_SORT
#include <algorithm>
#endif

  //=============================================================================
  Amesos_Dscpack::Amesos_Dscpack(const Epetra_LinearProblem &prob, const AMESOS::Parameter::List &ParameterList ) : DscGraph_(0), SymbolicFactorizationOK_(false), NumericFactorizationOK_(false)  {


  Problem_ = &prob ; 
  A_and_LU_built = false ; 
  FirstCallToSolve_ = true ; 
  ParameterList_ = &ParameterList ; 
  MyDSCObject = DSC_Begin() ; 
}

//=============================================================================
Amesos_Dscpack::~Amesos_Dscpack(void) {

  if ( MyDscRank>=0 && A_and_LU_built ) { 
    DSC_FreeAll( MyDSCObject ) ; 
    DSC_Close0( MyDSCObject ) ; 
    DSC_End( MyDSCObject ) ; 
  }

  delete DscGraph_;  // This might not exist, is it dangerous to delete it?
}

int Amesos_Dscpack::ReadParameterList() {
  if (ParameterList_->isParameterSublist("Dscpack") ) {
    AMESOS::Parameter::List DscpackParams = ParameterList_->sublist("Dscpack") ;
  }  
  return 0;
}


int Amesos_Dscpack::PerformSymbolicFactorization() {


  vector <int> Replicates;
  vector <int> Ap;
  vector <int> Ai;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  const Epetra_Comm &Comm = RowMatrixA->Comm();
  int iam = Comm.MyPID() ;

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int numrows = CastCrsMatrixA->NumGlobalRows();
  int numentries = CastCrsMatrixA->NumGlobalNonzeros();
  assert( numrows == CastCrsMatrixA->NumGlobalCols() );
  

  //
  //  Create a replicated map and graph 
  //
  vector<int> AllIDs( numrows ) ; 
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  Epetra_Map ReplicatedMap( -1, numrows, &AllIDs[0], 0, Comm);
  Epetra_Import importer( ReplicatedMap, OriginalMap );

  Epetra_Import ImportToReplicated( ReplicatedMap, OriginalMap);
  Epetra_CrsGraph ReplicatedGraph( Copy, ReplicatedMap, 0 ); 
  EPETRA_CHK_ERR( ReplicatedGraph.Import( (CastCrsMatrixA->Graph()), importer, Insert) );
  EPETRA_CHK_ERR( ReplicatedGraph.TransformToLocal() ) ; 


    //
    //  Convert the matrix to Ap, Ai
    //
    Replicates.resize( numrows );
    for( int i = 0 ; i < numrows; i++ ) Replicates[i] = 1; 
    Ap.resize( numrows+1 );
    Ai.resize( EPETRA_MAX( numrows, numentries) ) ; 

    int NumEntriesPerRow ;
    int *ColIndices = 0 ;
    int Ai_index = 0 ; 
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      EPETRA_CHK_ERR( ReplicatedGraph.ExtractMyRowView( MyRow, NumEntriesPerRow, ColIndices ) );
      Ap[MyRow] = Ai_index ; 
      for ( int j = 0; j < NumEntriesPerRow; j++ ) { 
	Ai[Ai_index] = ColIndices[j] ; 
	Ai_index++;
      }
    }
    assert( Ai_index == numentries ) ; 
    Ap[ numrows ] = Ai_index ; 


  //
  //  Call Dscpack Symbolic Factorization
  //  
  int OrderCode = 2;
  vector<double> MyANonZ;
  int numprocs  = Comm.NumProc() ;            
  
    NumLocalNonz = 0 ; 
    GlobalStructNewColNum = 0 ; 
    GlobalStructNewNum = 0 ;  
    GlobalStructOwner = 0 ; 
    LocalStructOldNum = 0 ; 
    
    NumGlobalCols = 0 ; 

    /*
      Dscpack uses a number of processes that is a power of 2
     */
    DscNumProcs = 1 ; 
    while ( DscNumProcs * 2 <=EPETRA_MIN( numprocs, 
					  DSC_Analyze( numrows, &Ap[0], &Ai[0], 
						       &Replicates[0] ) ) )
      DscNumProcs *= 2 ;
    
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

    } 


  //    A_and_LU_built = true; 

  SymbolicFactorizationOK_ = true ; 
  return 0;
}

int Amesos_Dscpack::PerformNumericFactorization() {

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  const Epetra_Comm &Comm = RowMatrixA->Comm();
  int iam = Comm.MyPID() ;
  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int numrows = CastCrsMatrixA->NumGlobalRows();
  int numentries = CastCrsMatrixA->NumGlobalNonzeros();

  assert( numrows == CastCrsMatrixA->NumGlobalCols() );
  
  //
  //  Call Dscpack to perform Numeric Factorization
  //  
  vector<double> MyANonZ;
  Epetra_Map DscMap( numrows, NumLocalCols, LocalStructOldNum, 0, Comm ) ;
  
    //
    //  Import from the CrsMatrix
    //
    Epetra_Import ImportToDsc( DscMap, OriginalMap );

    Epetra_CrsMatrix DscMat(Copy, DscMap, 0);
    EPETRA_CHK_ERR( DscMat.Import( *CastCrsMatrixA, ImportToDsc, Insert) );
    EPETRA_CHK_ERR( DscMat.TransformToLocal() ) ; 

    //    cout << "Amesos_Dscpack.cpp:: DscMat = " << DscMat << endl ; 


    //    assert( DscGraph_ == 0 ) ; 
    if ( DscGraph_ ) delete DscGraph_ ; 
    DscGraph_ = new Epetra_CrsGraph ( DscMat.Graph() ); 

    assert( MyDscRank >= 0 || NumLocalNonz == 0 ) ;
    assert( MyDscRank >= 0 || NumLocalCols == 0 ) ;
    assert( MyDscRank >= 0 || NumGlobalCols == 0  ) ; 
    MyANonZ.resize( NumLocalNonz ) ; 
    int NonZIndex = 0 ; 


    int max_num_entries = DscMat.MaxNumEntries() ; 
    vector<int> col_indices( max_num_entries ) ; 
    vector<double> mat_values( max_num_entries ) ; 
    assert( NumLocalCols == DscMap.NumMyElements() ) ;
    vector<int> my_global_elements( NumLocalCols ) ; 
    EPETRA_CHK_ERR( DscMap.MyGlobalElements( &my_global_elements[0] ) ) ;

    vector<int> GlobalStructOldColNum( NumGlobalCols ) ; 
      
    typedef pair<int, double> Data; 
    vector<Data> sort_array(max_num_entries); 
    vector<int>  sort_indices(max_num_entries);

    for ( int i = 0; i < NumLocalCols ; i++ ) { 
      assert( my_global_elements[i] == LocalStructOldNum[i] ) ; 
      int num_entries_this_row; 
      //  USE_LOCAL and not USE_LOCAL both work
#ifdef USE_LOCAL
      EPETRA_CHK_ERR( DscMat.ExtractMyRowCopy( i, max_num_entries, num_entries_this_row, 
					       &mat_values[0], &col_indices[0] ) ) ; 
#else
      EPETRA_CHK_ERR( DscMat.ExtractGlobalRowCopy( DscMat.GRID(i), max_num_entries, num_entries_this_row, 
						   &mat_values[0], &col_indices[0] ) ) ; 
#endif
      int OldRowNumber =  LocalStructOldNum[i] ;
      assert( GlobalStructOwner[ OldRowNumber ] != -1 ) ; 

      int NewRowNumber = GlobalStructNewColNum[ my_global_elements[ i ] ] ; 

      //
      //  Sort the column elements 
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
      const int SchemeCode = 1; 
#ifndef USE_LOCAL
      assert( NonZIndex == NumLocalNonz );
#endif

      EPETRA_CHK_ERR( DSC_NFactor ( MyDSCObject, SchemeCode, &MyANonZ[0], 
				    DSC_LLT,  DSC_LBLAS3, DSC_DBLAS2 ) ) ;

    }        //     if ( MyDscRank >= 0 ) 


  NumericFactorizationOK_ = true ; 
  return 0;
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

 PerformSymbolicFactorization();

 NumericFactorizationOK_ = false; 
 return 0;
}

int Amesos_Dscpack::NumericFactorization() {

  if ( ! SymbolicFactorizationOK_ ) PerformSymbolicFactorization();

  PerformNumericFactorization();

  return 0;
}

//
//  Solve() uses several intermediate matrices to convert the input matrix
//  to one that we can pass to the Sparse Direct Solver
//
//  Epetra_RowMatrix *RowMatrixA - The input matrix
//
int Amesos_Dscpack::Solve() { 
  //  int Solve() { 

  if ( ! SymbolicFactorizationOK_ ) PerformSymbolicFactorization();

  if ( ! NumericFactorizationOK_ ) PerformNumericFactorization();

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;
  Comm.Barrier();

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  //
  //  Coalesce the matrix graph onto process 0
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPIC = comm1.Comm() ;

  int numrows = CastCrsMatrixA->NumGlobalRows();
  int numentries = CastCrsMatrixA->NumGlobalNonzeros();

  assert( numrows == CastCrsMatrixA->NumGlobalCols() );
  

  //
  //  Convert vector b to a vector in the form that DSCPACK needs it
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

  Epetra_MultiVector *vecBvector = (vecB) ; // Ken gxx- do we need vecBvector?
  Epetra_Map DscMap( numrows, NumLocalCols, LocalStructOldNum, 0, Comm ) ;


  double *dscmapXvalues ;
  int dscmapXlda ;
  Epetra_MultiVector dscmapX( DscMap, nrhs ) ; 
  assert( dscmapX.ExtractView( &dscmapXvalues, &dscmapXlda ) == 0 ) ; 
  assert( dscmapXlda == NumLocalCols ) ; 

  double *dscmapBvalues ;
  int dscmapBlda ;
  Epetra_MultiVector dscmapB( DscMap, nrhs ) ; 
  assert( dscmapB.ExtractView( &dscmapBvalues, &dscmapBlda ) == 0 ) ; 
  assert( dscmapBlda == NumLocalCols ) ; 

  Epetra_Import ImportOriginalToDsc( DscMap, OriginalMap );
  dscmapB.Import( *vecBvector, ImportOriginalToDsc, Insert ) ;


  vector<double> ValuesInNewOrder( NumLocalCols ) ; 

  if ( MyDscRank >= 0 ) {
    for ( int j =0 ; j < nrhs; j++ ) { 
      for ( int i = 0; i < NumLocalCols; i++ ) { 
	ValuesInNewOrder[i] = dscmapBvalues[ DscGraph_->LCID( LocalStructOldNum[i] ) +j*dscmapBlda ] ;
      }
      EPETRA_CHK_ERR( DSC_InputRhsLocalVec ( MyDSCObject, &ValuesInNewOrder[0], NumLocalCols ) ) ;
      EPETRA_CHK_ERR( DSC_Solve ( MyDSCObject ) ) ; 
      EPETRA_CHK_ERR( DSC_GetLocalSolution ( MyDSCObject, &ValuesInNewOrder[0], NumLocalCols ) ) ; 
      for ( int i = 0; i < NumLocalCols; i++ ) { 
	dscmapXvalues[ DscGraph_->LCID( LocalStructOldNum[i] ) +j*dscmapXlda ] = ValuesInNewOrder[i];
      }
    }
    
  }

  Epetra_Import ImportDscToOriginal( OriginalMap, DscMap );
  vecX->Import( dscmapX, ImportDscToOriginal, Insert ) ;

  return(0) ; 
}
