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
#include "Epetra_Time.h"
#ifdef USE_STL_SORT
#include <algorithm>
#endif

//=============================================================================
Amesos_Dscpack::Amesos_Dscpack(const Epetra_LinearProblem &prob ) : 
  SymbolicFactorizationOK_(false), 
  NumericFactorizationOK_(false),
  DscGraph_(0), 
  UseTranspose_(false), // Dscpack is only for symmetric systems
  DscNumProcs(-1), // will be set later
  PrintTiming_(false),
  PrintStatus_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  verbose_(1),
  debug_(0), // turn it on only when debugging is needed
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0),
  Time_(0),  
  ImportToSerial_(0),
  DscMap_(0),
  MaxProcs_(-1)
{  
  Problem_ = &prob ; 
  A_and_LU_built = false ; 
  FirstCallToSolve_ = true ; 
  MyDSCObject = DSC_Begin() ; 
}

//=============================================================================
Amesos_Dscpack::~Amesos_Dscpack(void) {

  if ( MyDscRank>=0 && A_and_LU_built ) { 
    DSC_FreeAll( MyDSCObject ) ; 
    DSC_Close0( MyDSCObject ) ; 
    DSC_End( MyDSCObject ) ; 
  }

  if( Time_ ) { delete Time_; Time_ = 0; }
  if( ImportToSerial_ ) { delete ImportToSerial_; ImportToSerial_ = 0; }
  if( DscMap_ ) { delete DscMap_; DscMap_ = 0; }
    
  if ( DscGraph_) delete DscGraph_;  // This might not exist, is it dangerous to delete it?

  // MS // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

}

int Amesos_Dscpack::SetParameters( Teuchos::ParameterList &ParameterList ) 
{

  if( debug_ == 1 ) cout << "Entering `SetParameters()' ..." << endl;

  // ========================================= //
  // retrive KLU's parameters from list.       //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

  // compute norms of some vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",false);

  // compute the true residual Ax-b after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",false);

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",1);

  // possible debug statements
  // 0 - no debug
  // 1 - debug
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get("DebugLevel",0);

  // define on how many processes should be used
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcs_ = ParameterList.get("MaxProcs",-1);
  
  // MS // NO DSCPACK-specify parameters at this point, uncomment
  // MS // as necessary
  /*
  if (ParameterList.isSublist("Dscpack") ) {
    Teuchos::ParameterList DscpackParams = ParameterList.sublist("Dscpack") ;
  }
  */
  
  return 0;
}


int Amesos_Dscpack::PerformSymbolicFactorization()
{
  
  if( debug_ == 1 ) cout << "Entering `PerformSymbolicFactorization()' ..." << endl;

  Time_->ResetStartTime();
  
  vector <int> Replicates;
  vector <int> Ap;
  vector <int> Ai;

  Epetra_RowMatrix *RowMatrixA = Problem_->GetMatrix();
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ;
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm());
  MPIC = comm1.Comm() ;

  int numrows = CastCrsMatrixA->NumGlobalRows();
  int numentries = CastCrsMatrixA->NumGlobalNonzeros();
  assert( numrows == CastCrsMatrixA->NumGlobalCols() );

  //
  //  Create a replicated map and graph 
  //
  vector<int> AllIDs( numrows ) ; 
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  Epetra_Map ReplicatedMap( -1, numrows, &AllIDs[0], 0, Comm());
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
  
  ConTime_ += Time_->ElapsedTime();

  Time_->ResetStartTime();
  
  //
  //  Call Dscpack Symbolic Factorization
  //  
  int OrderCode = 2;
  vector<double> MyANonZ;
  
  NumLocalNonz = 0 ; 
  GlobalStructNewColNum = 0 ; 
  GlobalStructNewNum = 0 ;  
  GlobalStructOwner = 0 ; 
  LocalStructOldNum = 0 ; 
  
  NumGlobalCols = 0 ; 
  
  // MS // Have to define the maximum number of processes to be used
  // MS // This is only a suggestion as Dscpack uses a number of processes that is a power of 2  

  int NumGlobalNonzeros = GetProblem()->GetMatrix()->NumGlobalNonzeros();
  int NumRows = GetProblem()->GetMatrix()->NumGlobalRows(); 

  // optimal value for MaxProcs == -1
  
  int OptNumProcs1 = 1+EPETRA_MAX( NumRows/10000, NumGlobalNonzeros/1000000 );
  OptNumProcs1 = EPETRA_MIN(Comm().NumProc(),OptNumProcs1 );

  // optimal value for MaxProcs == -2

  int OptNumProcs2 = (int)sqrt(1.0*Comm().NumProc());
  if( OptNumProcs2 < 1 ) OptNumProcs2 = 1;

  // fix the value of MaxProcs

  switch( MaxProcs_ ) {
  case -1:
    MaxProcs_ = OptNumProcs1;
    break;
  case -2:
    MaxProcs_ = OptNumProcs2;
    break;
  case -3:
    MaxProcs_ = Comm().NumProc();
    break;
  }

  // MS // here I continue with the old code...
  
  DscNumProcs = 1 ; 
  while ( DscNumProcs * 2 <=EPETRA_MIN( MaxProcs_, 
					DSC_Analyze( numrows, &Ap[0], &Ai[0], 
						     &Replicates[0] ) ) )
    DscNumProcs *= 2 ;
  
  DSC_Open0( MyDSCObject, DscNumProcs, &MyDscRank, MPIC ) ; 
  
  NumLocalCols = 0 ; // This is for those processes not in the Dsc grid
  if ( MyDscRank >= 0 ) { 
    assert( Comm().MyPID() == MyDscRank ) ; 
    EPETRA_CHK_ERR( DSC_Order ( MyDSCObject, OrderCode, numrows, &Ap[0], &Ai[0], 
				&Replicates[0], &NumGlobalCols, &NumLocalStructs, 
				&NumLocalCols, &NumLocalNonz, 
				&GlobalStructNewColNum, &GlobalStructNewNum, 
				&GlobalStructOwner, &LocalStructOldNum ) ) ; 
    assert( NumGlobalCols == numrows ) ; 
    assert( NumLocalCols == NumLocalStructs ) ; 
  }
  
  if ( MyDscRank >= 0 ) { 
    int MaxSingleBlock; 
    
    const int Limit = 5000000 ;  //  Memory Limit set to 5 Terabytes 
    EPETRA_CHK_ERR( DSC_SFactor ( MyDSCObject, &TotalMemory_, 
				  &MaxSingleBlock, Limit, DSC_LBLAS3, DSC_DBLAS2 ) ) ; 
    
  }
  
  //    A_and_LU_built = true; 
  
  SymbolicFactorizationOK_ = true ; 

  SymTime_ += Time_->ElapsedTime();

  return 0;

}

int Amesos_Dscpack::PerformNumericFactorization()
{

  if( debug_ == 1 ) cout << "Entering `PerformNumericFactorization()' ..." << endl;
  
  Time_->ResetStartTime();

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm());
  MPIC = comm1.Comm() ;

  int numrows = CastCrsMatrixA->NumGlobalRows();
  assert( numrows == CastCrsMatrixA->NumGlobalCols() );
  
  //
  //  Call Dscpack to perform Numeric Factorization
  //  
  vector<double> MyANonZ;
  DscMap_ = new Epetra_Map(numrows, NumLocalCols, LocalStructOldNum, 0, Comm());
  
  //
  //  Import from the CrsMatrix
  //
  if( ImportToSerial_ == 0 ) {
    ImportToSerial_ = new Epetra_Import( *DscMap_, OriginalMap );
    assert( ImportToSerial_ != 0 );
  }
  
  Epetra_CrsMatrix DscMat(Copy, *DscMap_, 0);
  EPETRA_CHK_ERR( DscMat.Import( *CastCrsMatrixA, *ImportToSerial_, Insert) );
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
  assert( NumLocalCols == DscMap_->NumMyElements() ) ;
  vector<int> my_global_elements( NumLocalCols ) ; 
  EPETRA_CHK_ERR( DscMap_->MyGlobalElements( &my_global_elements[0] ) ) ;
  
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

  NumTime_ += Time_->ElapsedTime();
  
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


int Amesos_Dscpack::SymbolicFactorization()
{
  if( debug_ == 1 ) cout << "Entering `SymbolicFactorization()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumSymbolicFact_++;
  
  PerformSymbolicFactorization();
  
  NumericFactorizationOK_ = false; 
  return 0;
}

int Amesos_Dscpack::NumericFactorization()
{

  if( debug_ == 1 ) cout << "Entering `NumericFactorization()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumNumericFact_++;
  
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
int Amesos_Dscpack::Solve()
{

  if( debug_ == 1 ) cout << "Entering `Solve()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumSolve_++;

  if ( ! SymbolicFactorizationOK_ ) PerformSymbolicFactorization();

  if ( ! NumericFactorizationOK_ ) PerformNumericFactorization();

  Time_->ResetStartTime();
  
  // MS // it was GetOperator, now is GetMatrix (only matrices can
  // MS // be factorized
  Epetra_RowMatrix *RowMatrixA = Problem_->GetMatrix();
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  //  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 

  // FIXME: what the hell I am doing here ??
  //  assert(CastCrsMatrixA!=0);

  Comm().Barrier();

  //...//  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm());
  MPIC = comm1.Comm() ;

  // MS // some checks on matrix size
  int numrows = RowMatrixA->NumGlobalRows();
  assert( numrows == RowMatrixA->NumGlobalCols() );

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
  //  Epetra_Map DscMap( numrows, NumLocalCols, LocalStructOldNum, 0, Comm() ) ;


  double *dscmapXvalues ;
  int dscmapXlda ;
  Epetra_MultiVector dscmapX( *DscMap_, nrhs ) ; 
  assert( dscmapX.ExtractView( &dscmapXvalues, &dscmapXlda ) == 0 ) ; 
  assert( dscmapXlda == NumLocalCols ) ; 

  double *dscmapBvalues ;
  int dscmapBlda ;
  Epetra_MultiVector dscmapB( *DscMap_, nrhs ) ; 
  assert( dscmapB.ExtractView( &dscmapBvalues, &dscmapBlda ) == 0 ) ; 
  assert( dscmapBlda == NumLocalCols ) ; 

  // MS // erase this, use Import allocated only once
  //....//  Epetra_Import ImportOriginalToDsc( DscMap, OriginalMap );
  //....//  dscmapB.Import( *vecBvector, ImportOriginalToDsc, Insert ) ;
  dscmapB.Import( *vecBvector, *ImportToSerial_, Insert ) ;

  VecTime_ += Time_->ElapsedTime();
  Time_->ResetStartTime();
  
  // MS // now solve the problem
  
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

  SolTime_ += Time_->ElapsedTime();

  // MS // use always the same Import/Export, avoid allocations
  // MS // add timing
  //....//  Epetra_Import ImportDscToOriginal( OriginalMap, DscMap );
  
  Time_->ResetStartTime();  
  vecX->Export( dscmapX, *ImportToSerial_, Insert ) ;
  VecTime_ += Time_->ElapsedTime();

  // MS // compute vector norms if required
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Dscpack : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }

  // MS // compute true residual if required
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for( int i=0 ; i<nrhs ; ++i ) {
      (Problem_->GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));

      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Dscpack : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }
  
  return(0) ; 
}

// ================================================ ====== ==== ==== == =

void Amesos_Dscpack::PrintStatus()
{

  if( Comm().MyPID() != 0  ) return;

  int numrows =  GetProblem()->GetMatrix()->NumGlobalRows();
  int numentries = GetProblem()->GetMatrix()->NumGlobalNonzeros();
  
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Dscpack : Matrix has " << numrows << " rows"
       << " and " << numentries << " nonzeros" << endl;
  cout << "Amesos_Dscpack : Nonzero elements per row = "
       << 1.0*numentries/numrows << endl;
  cout << "Amesos_Dscpack : Percentage of nonzero elements = "
       << 100.0*numentries/(pow(numentries,2.0)) << endl;
  cout << "Amesos_Dscpack : Available process(es) = " << Comm().NumProc() << endl;
  cout << "Amesos_Dscpack : Process(es) used = " << DscNumProcs
       << ", idle = " << Comm().NumProc() - DscNumProcs << endl;
  cout << "Amesos_Dscpack : Estimated total memory for factorization =  " 
       << TotalMemory_ << " Mbytes" << endl; 
  cout << "----------------------------------------------------------------------------" << endl;

  DSC_DoStats( MyDSCObject );
  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Dscpack::PrintTiming()
{
  if( Comm().MyPID() ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Dscpack : Time to convert matrix to DSCPACK format = "
       << ConTime_ << " (s)" << endl;
  cout << "Amesos_Dscpack : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Dscpack : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Dscpack : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Dscpack : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Dscpack : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Dscpack : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Dscpack : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
}
