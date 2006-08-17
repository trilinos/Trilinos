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

#include "Amesos_Dscpack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Teuchos_RefCountPtr.hpp"
#include <algorithm>

#define DBL_R_NUM
extern "C" {
#include "dscmain.h"
}

class Amesos_Dscpack_Pimpl {
public:
  DSC_Solver	MyDSCObject_;
} ;
//=============================================================================
Amesos_Dscpack::Amesos_Dscpack(const Epetra_LinearProblem &prob ) : 
  PrivateDscpackData_( rcp( new Amesos_Dscpack_Pimpl() ) ),
  DscNumProcs(-1), // will be set later
  MaxProcs_(-1)
{  
  Problem_ = &prob ; 
  A_and_LU_built = false ; 
  FirstCallToSolve_ = true ; 
  PrivateDscpackData_->MyDSCObject_ = DSC_Begin() ; 

  MyDscRank = -1 ; 
}

//=============================================================================
Amesos_Dscpack::~Amesos_Dscpack(void) 
{
  // MS // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

  if ( MyDscRank>=0 && A_and_LU_built ) { 
    DSC_FreeAll( PrivateDscpackData_->MyDSCObject_ ) ; 
    DSC_Close0( PrivateDscpackData_->MyDSCObject_ ) ; 
  }
  DSC_End( PrivateDscpackData_->MyDSCObject_ ) ; 
}

//=============================================================================
int Amesos_Dscpack::SetParameters( Teuchos::ParameterList &ParameterList) 
{
  // ========================================= //
  // retrive DSCPACK's parameters from list.   //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );
  
  // MS // NO DSCPACK-specify parameters at this point, uncomment
  // MS // as necessary
  /*
  if (ParameterList.isSublist("Dscpack") ) {
    Teuchos::ParameterList DscpackParams = ParameterList.sublist("Dscpack") ;
  }
  */
  
  return 0;
}

//=============================================================================
int Amesos_Dscpack::PerformSymbolicFactorization()
{
  ResetTime(0);
  ResetTime(1);

  MyPID_    = Comm().MyPID();
  NumProcs_ = Comm().NumProc();
  
  Epetra_RowMatrix *RowMatrixA = Problem_->GetMatrix();
  if (RowMatrixA == 0)
    AMESOS_CHK_ERR(-1);

  const Epetra_Map& OriginalMap = RowMatrixA->RowMatrixRowMap() ;
  const Epetra_MpiComm& comm1   = dynamic_cast<const Epetra_MpiComm &> (Comm());
  int numrows                   = RowMatrixA->NumGlobalRows();
  int numentries                = RowMatrixA->NumGlobalNonzeros();

  Teuchos::RefCountPtr<Epetra_CrsGraph> Graph;

  Epetra_CrsMatrix* CastCrsMatrixA = 
    dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA); 

  if (CastCrsMatrixA)
  {
    Graph = Teuchos::rcp(const_cast<Epetra_CrsGraph*>(&(CastCrsMatrixA->Graph())), false);
  }
  else
  {
    int MaxNumEntries = RowMatrixA->MaxNumEntries();
    Graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, OriginalMap, MaxNumEntries));

    vector<int>    Indices(MaxNumEntries);
    vector<double> Values(MaxNumEntries);

    for (int i = 0 ; i < RowMatrixA->NumMyRows() ; ++i)
    {
      int NumEntries;
      RowMatrixA->ExtractMyRowCopy(i, MaxNumEntries, NumEntries,
                                   &Values[0], &Indices[0]);

      for (int j = 0 ; j < NumEntries ; ++j)
        Indices[j] = RowMatrixA->RowMatrixColMap().GID(Indices[j]);

      int GlobalRow = RowMatrixA->RowMatrixRowMap().GID(i);
      Graph->InsertGlobalIndices(GlobalRow, NumEntries, &Indices[0]);
    }

    Graph->FillComplete();
  }

  //
  //  Create a replicated map and graph 
  //
  vector<int> AllIDs( numrows ) ; 
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  Epetra_Map      ReplicatedMap( -1, numrows, &AllIDs[0], 0, Comm());
  Epetra_Import   ReplicatedImporter(ReplicatedMap, OriginalMap);
  Epetra_CrsGraph ReplicatedGraph( Copy, ReplicatedMap, 0 ); 

  AMESOS_CHK_ERR(ReplicatedGraph.Import(*Graph, ReplicatedImporter, Insert));
  AMESOS_CHK_ERR(ReplicatedGraph.FillComplete());

  //
  //  Convert the matrix to Ap, Ai
  //
  vector <int> Replicates(numrows);
  vector <int> Ap(numrows + 1);
  vector <int> Ai(EPETRA_MAX(numrows, numentries));

  for( int i = 0 ; i < numrows; i++ ) Replicates[i] = 1; 
  
  int NumEntriesPerRow ;
  int *ColIndices = 0 ;
  int Ai_index = 0 ; 
  for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
    AMESOS_CHK_ERR( ReplicatedGraph.ExtractMyRowView( MyRow, NumEntriesPerRow, ColIndices ) );
    Ap[MyRow] = Ai_index ; 
    for ( int j = 0; j < NumEntriesPerRow; j++ ) { 
      Ai[Ai_index] = ColIndices[j] ; 
      Ai_index++;
    }
  }
  assert( Ai_index == numentries ) ; 
  Ap[ numrows ] = Ai_index ; 
  
  AddTime("matrix conversion", 0);

  ResetTime(0);

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
  OptNumProcs1 = EPETRA_MIN(NumProcs_,OptNumProcs1 );

  // optimal value for MaxProcs == -2

  int OptNumProcs2 = (int)sqrt(1.0 * NumProcs_);
  if( OptNumProcs2 < 1 ) OptNumProcs2 = 1;

  // fix the value of MaxProcs

  switch (MaxProcs_) 
  {
  case -1:
    MaxProcs_ = EPETRA_MIN(OptNumProcs1, NumProcs_);
    break;
  case -2:
    MaxProcs_ = EPETRA_MIN(OptNumProcs2, NumProcs_);
    break;
  case -3:
    MaxProcs_ = NumProcs_;
    break;
  }

#if 0
  if (MyDscRank>=0 && A_and_LU_built) { 
    DSC_ReFactorInitialize(PrivateDscpackData_->MyDSCObject);
  }
#endif
  //  if ( ! A_and_LU_built ) { 
  //    DSC_End( PrivateDscpackData_->MyDSCObject ) ; 
  //    PrivateDscpackData_->MyDSCObject = DSC_Begin() ;
  //  } 

  // MS // here I continue with the old code...
  
  AddTime("overhead", 1);

  DscNumProcs = 1 ; 
  int DscMax = DSC_Analyze( numrows, &Ap[0], &Ai[0], &Replicates[0] );

  while ( DscNumProcs * 2 <=EPETRA_MIN( MaxProcs_, DscMax ) )  DscNumProcs *= 2 ;
  
  MyDscRank = -1; 
  DSC_Open0( PrivateDscpackData_->MyDSCObject_, DscNumProcs, &MyDscRank, comm1.Comm()) ; 
  
  NumLocalCols = 0 ; // This is for those processes not in the Dsc grid
  if ( MyDscRank >= 0 ) { 
    assert( MyPID_ == MyDscRank ) ; 
    AMESOS_CHK_ERR( DSC_Order ( PrivateDscpackData_->MyDSCObject_, OrderCode, numrows, &Ap[0], &Ai[0], 
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
    AMESOS_CHK_ERR( DSC_SFactor ( PrivateDscpackData_->MyDSCObject_, &TotalMemory_, 
				  &MaxSingleBlock, Limit, DSC_LBLAS3, DSC_DBLAS2 ) ) ; 
    
  }
  
  //  A_and_LU_built = true;   // If you uncomment this, TestOptions fails
  
  AddTime("symbolic", 0);

  return(0);
}

//=============================================================================
int Amesos_Dscpack::PerformNumericFactorization()
{
  ResetTime(0);
  ResetTime(1);

  Epetra_RowMatrix* RowMatrixA = Problem_->GetMatrix();
  if (RowMatrixA == 0)
    AMESOS_CHK_ERR(-1);

  const Epetra_Map& OriginalMap = RowMatrixA->RowMatrixRowMap() ; 

  int numrows = RowMatrixA->NumGlobalRows();
  assert( numrows == RowMatrixA->NumGlobalCols() );
  
  //
  //  Call Dscpack to perform Numeric Factorization
  //  
  vector<double> MyANonZ;
#if 0
    if ( IsNumericFactorizationOK_ ) { 
      DSC_ReFactorInitialize(PrivateDscpackData_->MyDSCObject);
    }
#endif

  DscRowMap_ = Teuchos::rcp(new Epetra_Map(numrows, NumLocalCols, 
                                           LocalStructOldNum, 0, Comm()));

  if (DscRowMap_.get() == 0)
    AMESOS_CHK_ERR(-1);
  
  Importer_ = rcp(new Epetra_Import(DscRowMap(), OriginalMap));
    
  //
  //  Import from the CrsMatrix
  //
  Epetra_CrsMatrix DscMat(Copy, DscRowMap(), 0);
  AMESOS_CHK_ERR(DscMat.Import(*RowMatrixA, Importer(), Insert));
  AMESOS_CHK_ERR(DscMat.FillComplete()); 

  DscColMap_ = Teuchos::rcp(new Epetra_Map(DscMat.RowMatrixColMap()));

  assert( MyDscRank >= 0 || NumLocalNonz == 0 ) ;
  assert( MyDscRank >= 0 || NumLocalCols == 0 ) ;
  assert( MyDscRank >= 0 || NumGlobalCols == 0  ) ; 
  MyANonZ.resize( NumLocalNonz ) ; 
  int NonZIndex = 0 ;

  int max_num_entries = DscMat.MaxNumEntries() ; 
  vector<int> col_indices( max_num_entries ) ; 
  vector<double> mat_values( max_num_entries ) ; 
  assert( NumLocalCols == DscRowMap().NumMyElements() ) ;
  vector<int> my_global_elements( NumLocalCols ) ; 
  AMESOS_CHK_ERR(DscRowMap().MyGlobalElements( &my_global_elements[0] ) ) ;
  
  vector<int> GlobalStructOldColNum( NumGlobalCols ) ; 
  
  typedef pair<int, double> Data; 
  vector<Data> sort_array(max_num_entries); 
  vector<int>  sort_indices(max_num_entries);
  
  for ( int i = 0; i < NumLocalCols ; i++ ) { 
    assert( my_global_elements[i] == LocalStructOldNum[i] ) ; 
    int num_entries_this_row; 
#ifdef USE_LOCAL
    AMESOS_CHK_ERR( DscMat.ExtractMyRowCopy( i, max_num_entries, num_entries_this_row, 
					     &mat_values[0], &col_indices[0] ) ) ; 
#else
    AMESOS_CHK_ERR( DscMat.ExtractGlobalRowCopy( DscMat.GRID(i), max_num_entries, num_entries_this_row, 
						 &mat_values[0], &col_indices[0] ) ) ; 
#endif
    int OldRowNumber =  LocalStructOldNum[i] ;
    if (GlobalStructOwner[ OldRowNumber ] == -1)
      AMESOS_CHK_ERR(-1);
    
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
    sort(&sort_array[0], &sort_array[num_entries_this_row]);
    
    for ( int j = 0; j < num_entries_this_row; j++ ) { 
      int NewColNumber = sort_array[j].first ; 
      if ( NewRowNumber <= NewColNumber ) MyANonZ[ NonZIndex++ ] = sort_array[j].second ; 

#ifndef USE_LOCAL
      assert( NonZIndex <= NumLocalNonz );  // This assert can fail on non-symmetric matrices
#endif
    }
  }
  
  AddTime("overhead", 1);

  if ( MyDscRank >= 0 ) { 
    const int SchemeCode = 1; 
#ifndef USE_LOCAL
    assert( NonZIndex == NumLocalNonz );
#endif
    
    AMESOS_CHK_ERR( DSC_NFactor ( PrivateDscpackData_->MyDSCObject_, SchemeCode, &MyANonZ[0], 
				  DSC_LLT,  DSC_LBLAS3, DSC_DBLAS2 ) ) ;
    
  }        //     if ( MyDscRank >= 0 ) 
  
  IsNumericFactorizationOK_ = true ; 

  AddTime("numeric", 0);
  
  return(0);
}

bool Amesos_Dscpack::MatrixShapeOK() const 
{
  bool OK =  GetProblem()->IsOperatorSymmetric() ;

  //
  //  The following test is redundant.  I have left it here in case the 
  //  IsOperatorSymmetric test turns out not to be reliable.
  //
  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}

//=============================================================================
int Amesos_Dscpack::SymbolicFactorization()
{
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
  InitTime(Comm(), 2);

  AMESOS_CHK_ERR(PerformSymbolicFactorization());

  IsSymbolicFactorizationOK_ = true; 
  NumSymbolicFact_++;

  return(0);
}

//=============================================================================
int Amesos_Dscpack::NumericFactorization()
{
  IsNumericFactorizationOK_ = false;

  if (!IsSymbolicFactorizationOK_) 
    AMESOS_CHK_ERR(SymbolicFactorization());

  AMESOS_CHK_ERR(PerformNumericFactorization());

  IsNumericFactorizationOK_ = true;
  NumNumericFact_++;
  
  return(0);
}

//=============================================================================
int Amesos_Dscpack::Solve()
{
  if (IsNumericFactorizationOK_ == false) 
    AMESOS_CHK_ERR(NumericFactorization());

  ResetTime(0);
  ResetTime(1);
  
  Epetra_RowMatrix *RowMatrixA = Problem_->GetMatrix();
  if (RowMatrixA == 0)
    AMESOS_CHK_ERR(-1);

  // MS // some checks on matrix size
  if (RowMatrixA->NumGlobalRows() != RowMatrixA->NumGlobalCols())
    AMESOS_CHK_ERR(-1);

  //  Convert vector b to a vector in the form that DSCPACK needs it
  //
  Epetra_MultiVector* vecX = Problem_->GetLHS(); 
  Epetra_MultiVector* vecB = Problem_->GetRHS(); 

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1); // something wrong with input

  int NumVectors = vecX->NumVectors(); 
  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-2);

  double *dscmapXvalues ;
  int dscmapXlda ;
  Epetra_MultiVector dscmapX(DscRowMap(),NumVectors) ; 
  int ierr;
  AMESOS_CHK_ERR(dscmapX.ExtractView(&dscmapXvalues,&dscmapXlda));
  assert (dscmapXlda == NumLocalCols); 

  double *dscmapBvalues ;
  int dscmapBlda ;
  Epetra_MultiVector dscmapB(DscRowMap(), NumVectors ) ; 
  ierr = dscmapB.ExtractView( &dscmapBvalues, &dscmapBlda );
  AMESOS_CHK_ERR(ierr);
  assert( dscmapBlda == NumLocalCols ) ; 

  AMESOS_CHK_ERR(dscmapB.Import(*vecB, Importer(), Insert));

  AddTime("vector redistribution", 0);
  ResetTime(0);
  
  // MS // now solve the problem
  
  vector<double> ValuesInNewOrder( NumLocalCols ) ;

  AddTime("overhead", 1);

  if ( MyDscRank >= 0 ) {
    for ( int j =0 ; j < NumVectors; j++ ) { 
      for ( int i = 0; i < NumLocalCols; i++ ) { 
	ValuesInNewOrder[i] = dscmapBvalues[DscColMap().LID( LocalStructOldNum[i] ) +j*dscmapBlda ] ;
      }
      AMESOS_CHK_ERR( DSC_InputRhsLocalVec ( PrivateDscpackData_->MyDSCObject_, &ValuesInNewOrder[0], NumLocalCols ) ) ;
      AMESOS_CHK_ERR( DSC_Solve ( PrivateDscpackData_->MyDSCObject_ ) ) ; 
      AMESOS_CHK_ERR( DSC_GetLocalSolution ( PrivateDscpackData_->MyDSCObject_, &ValuesInNewOrder[0], NumLocalCols ) ) ; 
      for ( int i = 0; i < NumLocalCols; i++ ) { 
	dscmapXvalues[DscColMap().LID( LocalStructOldNum[i] ) +j*dscmapXlda ] = ValuesInNewOrder[i];
      }
    }
  }

  AddTime("solve", 0);
  ResetTime(0);
  ResetTime(1);

  vecX->Export( dscmapX, Importer(), Insert ) ;

  AddTime("vector redistribution", 0);

  if (ComputeTrueResidual_)
    ComputeTrueResidual(*(GetProblem()->GetMatrix()), *vecX, *vecB, 
                        false, "Amesos_Dscpack");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Dscpack");
  AddTime("overhead", 1);
  
  NumSolve_++;

  return(0) ; 
}

// ======================================================================
void Amesos_Dscpack::PrintStatus() const
{
  if (Problem_->GetOperator() != 0 && MyPID_ != 0)
  {
    string p = "Amesos_Dscpack : ";
    PrintLine();

    int n = GetProblem()->GetMatrix()->NumGlobalRows();
    int nnz = GetProblem()->GetMatrix()->NumGlobalNonzeros();

    cout << p << "Matrix has " << n << " rows"
         << " and " << nnz << " nonzeros" << endl;
    cout << p << "Nonzero elements per row = "
         << 1.0 *  nnz / n << endl;
    cout << p << "Percentage of nonzero elements = "
         << 100.0 * nnz /(pow(n,2.0)) << endl;
    cout << p << "Available process(es) = " << NumProcs_ << endl;
    cout << p << "Process(es) used = " << DscNumProcs
         << ", idle = " << NumProcs_ - DscNumProcs << endl;
    cout << p << "Estimated total memory for factorization =  " 
         << TotalMemory_ << " Mbytes" << endl; 
  }

  if ( MyDscRank >= 0 ) 
    DSC_DoStats( PrivateDscpackData_->MyDSCObject_ );

  if (!MyPID_)
    PrintLine();

  return;
}

// ====================================================================== 
void Amesos_Dscpack::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || MyPID_ != 0)
    return;

  double ConTime = GetTime("conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");
  double OveTime = GetTime("overhead");

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Dscpack : ";
  PrintLine();

  cout << p << "Time to convert matrix to DSCPACK format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = "
       << SymTime * NumSymbolicFact_ << " (s), avg = " << SymTime << " (s)" << endl;
  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime * NumNumericFact_ << " (s), avg = " << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime * NumSolve_ << " (s), avg = " << SolTime << " (s)" << endl;

  double tt = SymTime * NumSymbolicFact_ + NumTime * NumNumericFact_ + SolTime * NumSolve_;
  if (tt != 0)
  {
    cout << p << "Total time spent in Amesos = " << tt << " (s) " << endl;
    cout << p << "Total time spent in the Amesos interface = " << OveTime << " (s)" << endl;
    cout << p << "(the above time does not include DSCPACK time)" << endl;
    cout << p << "Amesos interface time / total time = " << OveTime / tt << endl;
  }

  PrintLine();

  return;
}
