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


// TO DO: use Stat structure to print statistics???
// allow users to specify usermap ???

#include "Amesos_Superludist.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
// #include "CrsMatrixTranspose.h"
#include "Epetra_Time.h"

int Superludist_NumProcRows( int NumProcs ) {
#ifdef TFLOP
  //  Else, parameter 6 of DTRSV CTNLU is incorrect 
  return 1;
#else
  int i;
  int NumProcRows ;
  for ( i = 1; i*i <= NumProcs; i++ ) 
    ;
  bool done = false ;
  for ( NumProcRows = i-1 ; done == false ; ) {
    int NumCols = NumProcs / NumProcRows ; 
    if ( NumProcRows * NumCols == NumProcs ) 
      done = true; 
    else 
      NumProcRows-- ; 
  }
  return NumProcRows;
#endif
}

  //=============================================================================
Amesos_Superludist::Amesos_Superludist(const Epetra_LinearProblem &prob ) :
    GridCreated_(0), 
    FactorizationDone_(0), 
    NumRows_(0), 
    NumGlobalNonzeros_(0), 
    UniformMap_(0), 
    UniformMatrix_(0) ,
    ExportToDist_(0),
    ImportToDistributed_(0),
    ImportBackToOriginal_(0),
    RowMatrixA_(0), 
    SuperluMat_(0),
    vecBdistributed_(0),
    vecXdistributed_(0),
    nprow_(0),
    npcol_(0),
    PrintNonzeros_(false),
    ColPerm_("MMD_AT_PLUS_A"),
    RowPerm_("NATURAL"),
    perm_c_(0),
    perm_r_(0),
    IterRefine_("NO"),
    ReplaceTinyPivot_(false),
    FactorizationOK_(false), 
    UseTranspose_(false),
    PrintStatus_(false),
    PrintTiming_(false),    
    verbose_(1),
    debug_(0),
    NumTime_(0.0),
    SolveTime_(0.0),
    NumSymbolicFact_(0),
    NumNumericFact_(0),
    NumSolve_(0)
{

  Problem_ = &prob ; 
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ; 
}

//=============================================================================
Amesos_Superludist::~Amesos_Superludist(void) {

  // print out the value of used parameters.
  // This is defaulted to no.

  if( PrintTiming_ ) PrintTiming();
  if( PrintStatus_ ) PrintStatus();
  
  if ( FactorizationDone_ ) {
    SUPERLU_FREE( SuperluA_.Store );
    ScalePermstructFree(&ScalePermstruct_);
    Destroy_LU(NumRows_, &grid_, &LUstruct_);
    LUstructFree(&LUstruct_);
    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
  }
  if ( GridCreated_ ) {
    superlu_gridexit(&grid_);
  }
  if (vecBdistributed_ ) delete vecBdistributed_ ; 
  if (vecXdistributed_ ) delete vecXdistributed_ ; 
  if (UniformMap_ ) delete UniformMap_ ; 
  if (UniformMatrix_ ) delete UniformMatrix_ ; 
  if (ExportToDist_) delete ExportToDist_ ;
  if (ImportToDistributed_) delete ImportToDistributed_;
  if (ImportBackToOriginal_) delete ImportBackToOriginal_;

}

//
//  SetParameters
//
//    Preconditions:
//       None
//
//    Postconditions:
//       MaxProcesses_, Redistribute_, AddZeroToDiag_, FactOption_ and 
//       ReuseSymbolic_ set according to ParameterList_ values.
//
int Amesos_Superludist::SetParameters( Teuchos::ParameterList &ParameterList ) {

  if( debug_ == 1 ) cout << "Entering `SetParameters()' ..." << endl;
  
  if( (int) &ParameterList == 0 ) return 0;

  if (ParameterList.isSublist("Superludist") ) {
    Teuchos::ParameterList SuperludistParams = ParameterList.sublist("Superludist") ;
  }  
  //
  //  We have to set these to their defaults here because user codes 
  //  are not guaranteed to have a "Superludist" parameter list.
  //
  Redistribute_ = true;
  AddZeroToDiag_ = false;
  FactOption_ = SamePattern_SameRowPerm ;
  ReuseSymbolic_ = false ; 

  MaxProcesses_ = - 1; 
  nprow_ = 0;
  npcol_ = 0;
  Equil_ = true;
  ColPerm_ = "MMD_AT_PLUS_A";
  perm_c_ = 0;
  RowPerm_ = "LargeDiag";
  perm_r_ = 0;
  IterRefine_ = "NO";
  ReplaceTinyPivot_ = true;
  
  PrintNonzeros_ = false;
  
  // Ken, I modified so that parameters are not added to the list if not present.
  // Is it fine for you ????
  
  // parameters for all packages

  if( ParameterList.isParameter("Redistribute") )
    Redistribute_ = ParameterList.get("Redistribute",Redistribute_);  

  if( ParameterList.isParameter("AddZeroToDiag") )
    AddZeroToDiag_ = ParameterList.get("AddZeroToDiag",AddZeroToDiag_); 
  // Ken, should we have a parameter like AddToDiag too ???
  /*
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);
  */

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // debug level:
  // 0 - no debug output
  // 1 - some output
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get("DebugLevel", 0);

  MaxProcesses_ = ParameterList.get("MaxProcs",MaxProcesses_);
  if( MaxProcesses_ > Comm().NumProc() ) MaxProcesses_ = Comm().NumProc();
  
  // parameters for Superludist only
  
  if (ParameterList.isSublist("Superludist") ) {
    Teuchos::ParameterList SuperludistParams = ParameterList.sublist("Superludist") ;
    if( SuperludistParams.isParameter("ReuseSymbolic") )
      ReuseSymbolic_ = SuperludistParams.get("ReuseSymbolic",ReuseSymbolic_);
    string FactOption;
    
    if( SuperludistParams.isParameter("Fact") )
      FactOption = SuperludistParams.get("Fact", "SamePattern_SameRowPerm");
    

    //
    //  Neither FactOption = DOFACT nor FactOption = FACTORED make any sense here.
    //    FactOption is only relevant if ReuseSympolic_ is set true.
    //        FactOption = DOFACT does not reuse the symbolic factorization, 
    //      the correct way not to reuse symbolic factorization is to 
    //      set reusesymbolic_ to false 
    //        FactOption == FACTORED causes the call to NumericFactorization 
    //      to have no affect.  
    //

    if( FactOption == "SamePattern_SameRowPerm" ) FactOption_ = SamePattern_SameRowPerm;
    else if( FactOption == "DOFACT" ) FactOption_ = DOFACT;
    else if( FactOption == "SamePattern" ) FactOption_ = SamePattern;
    else if( FactOption == "FACTORED" ) FactOption_ = FACTORED;

    /* I moved this on the non-superlu level 
    MaxProcesses_ = SuperludistParams.get("MaxProcesses",MaxProcesses_);
    if( MaxProcesses_ > Comm().NumProc() ) MaxProcesses_ = Comm().NumProc();
    */

    /* Ken, should we support these too ?? 
    if( SuperludistParams.isParameter("nprow") )
      nprow_ = SuperludistParams.get("nprow",MaxProcesses_);
    if( SuperludistParams.isParameter("npcol") ) 
      npcol_ = SuperludistParams.get("npcol",MaxProcesses_);
    */
    
    if(  SuperludistParams.isParameter("Equil") )
      Equil_ = SuperludistParams.get("Equil",true);

    if( SuperludistParams.isParameter("ColPerm") )
      ColPerm_ = SuperludistParams.get("ColPerm","MMD_AT_PLUS_A");

    if( ColPerm_ == "MY_PERMC" ) {
      if( SuperludistParams.isParameter("perm_c") )
	perm_c_ = SuperludistParams.get("perm_c",perm_c_);
    }

    if( SuperludistParams.isParameter("RowPerm") )
      RowPerm_ = SuperludistParams.get("RowPerm","LargeDiag");
    if( RowPerm_ == "MY_PERMR" ) {
      if( SuperludistParams.isParameter("perm_r") )
	perm_r_ = SuperludistParams.get("perm_r",perm_r_);
    }

    if( SuperludistParams.isParameter("IterRefine") )
      IterRefine_ = SuperludistParams.get("IterRefine","NO");

    if( SuperludistParams.isParameter("ReplaceTinyPivot") )
      ReplaceTinyPivot_ = SuperludistParams.get("ReplaceTinyPivot",true);

    if( SuperludistParams.isParameter("PrintNonzeros") )
      PrintNonzeros_ = SuperludistParams.get("PrintNonzeros",false);

  }
  
  return 0;
}

//
//  RedistributeA
//
//    Preconditions:
//
//    Postconditions:
//       nprow_, npcol_

//       UniformMap_
//       UniformMatrix_
//       ExportToDist_
//
int Amesos_Superludist::RedistributeA( ) {

  if( debug_ == 1 ) cout << "Entering `RedistributeA()' ..." << endl;
 
  const Epetra_Map &OriginalMap = RowMatrixA_->RowMatrixRowMap() ; 

  //  Establish the grid (nprow vs. npcol) 
  //  Note:  The simple heuristic below is untested and 
  //  will lead to a poor grid shape if ProcessNumHeuristic is large
  //  and nearly prime.
  //
  NumGlobalNonzeros_ = RowMatrixA_->NumGlobalNonzeros() ; 
  assert( NumRows_ == RowMatrixA_->NumGlobalRows() ) ; 

  const Epetra_Comm &Comm_ = RowMatrixA_->Comm();
  int NumberOfProcesses = Comm_.NumProc() ; 
  if ( MaxProcesses_ > 0 ) {
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses, MaxProcesses_ ) ; 
  }
  else {
    // Ken, what about this:
    // -1 ==> classical Amesos_Superludist approach
    // -2 ==> square root of number of processes (as we are doing in ML)
    if( MaxProcesses_ == -1 ) {
      int ProcessNumHeuristic = 1+EPETRA_MAX( NumRows_/10000, NumGlobalNonzeros_/1000000 );
      NumberOfProcesses = EPETRA_MIN( NumberOfProcesses,  ProcessNumHeuristic );
      nprow_ = Superludist_NumProcRows( NumberOfProcesses ) ; 
      npcol_ = NumberOfProcesses / nprow_ ;
      assert ( nprow_ * npcol_ == NumberOfProcesses ) ;
    }
    else if( MaxProcesses_ == -2 ) {
      nprow_ = (int) pow(1.0 * NumberOfProcesses, 0.26 );
      if( nprow_ < 1 ) nprow_ = 1;
      npcol_ = nprow_;
      NumberOfProcesses = npcol_ * nprow_;
    }
  }

  //
  //  Compute a cannonical uniform distribution: MyFirstElement - The
  //  first element which this processor would have
  //  NumExpectedElemetns - The number of elements which this
  //  processor would have
  //
  //
  //  This is the same distribution that Epetra_Map( NumRows_, 0, Comm ) 
  //  would give by default, when NumberOfProcesses = Comm_.NumProc() ; 
  //
  int m_per_p = NumRows_ / NumberOfProcesses ;
  int remainder = NumRows_ - ( m_per_p * NumberOfProcesses );
  int MyFirstElement = iam_ * m_per_p + EPETRA_MIN( iam_, remainder );
  int MyFirstNonElement = (iam_+1) * m_per_p + EPETRA_MIN( iam_+1, remainder );
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 

  if ( iam_ >= NumberOfProcesses ) {
    NumExpectedElements = 0 ; 
  }
  assert( NumRows_ ==  RowMatrixA_->NumGlobalRows() ) ; 
  if ( UniformMap_ ) delete( UniformMap_ ) ; 
  UniformMap_ = new Epetra_Map( NumRows_, NumExpectedElements, 0, Comm_ );
  if ( UniformMatrix_ ) delete( UniformMatrix_ ) ; 
  UniformMatrix_ = new Epetra_CrsMatrix(Copy, *UniformMap_, 0);
  if ( ExportToDist_ ) delete ExportToDist_;
  ExportToDist_ = new Epetra_Export( OriginalMap, *UniformMap_);
  if (ImportToDistributed_) delete ImportToDistributed_;
  ImportToDistributed_ = new Epetra_Import( *UniformMap_, OriginalMap);

  if (ImportBackToOriginal_) delete ImportBackToOriginal_;
  ImportBackToOriginal_ = new Epetra_Import( OriginalMap,*UniformMap_);
  UniformMatrix_->Export( *RowMatrixA_, *ExportToDist_, Add ); 

  if (AddZeroToDiag_ ) { 
    //
    //  Add 0.0 to each diagonal entry to avoid empty diagonal entries;
    //
    double zero = 0.0;
    UniformMatrix_->SetTracebackMode(0);
    for ( int i = 0 ; i < UniformMap_->NumGlobalElements(); i++ ) 
      if ( UniformMatrix_->LRID(i) >= 0 ) 
	UniformMatrix_->InsertGlobalValues( i, 1, &zero, &i ) ;
    UniformMatrix_->SetTracebackMode(1);
  }

  UniformMatrix_->FillComplete() ; 

  //  I can't imagine that the following line is necessary:
  //  Commented out:  Feb 28th, 2004
  //  UniformMatrix_->Export( *RowMatrixA_, *ExportToDist_, Insert ); 
  
  return 0;
}


//
//  Factor
//
//    Preconditions:
//       Problem_
//         ->GetOperator() must be a RowMatrix else return -1
//         If Redistribute_ is not set, ->GetOperator()->RowMatrixRowMap() 
//           must be a LinearMap() else return -2
//         ->GetOperator() must represent a square matrix, else return -3
//         ->GetOperator() must be non-singular (and factorizable by Superludist)
//           else positive return code.  Values less than n indicate the 
//           diagonal element which was zero.  Values greater than n indicate
//           a memory allocation failure.  
//
//    Postconditions:
//       The matrix specified by Problem_->Operator() will have been redistributed,
//         converted to the form needed by Superludist and factored.
//       nprow_, npcol_
//       UniformMap_
//       UniformMatrix_
//       ExportToDist_
//       SuperLUmat_
//       RowValuesV_
//       ColValuesV_
//       Ap_, Ai_, Aval_
//       SuperluA_
//       SuperLU internal data structures reflecting the LU factorization:
//         ScalePermstruct_
//         LUstructInit_
//
//   
int Amesos_Superludist::Factor( ) {

  if( debug_ == 1 ) cout << "Entering `Factor()' ..." << endl;
  
  //
  //  For now, if you change the shape of a matrix, you need to 
  //  create a new Amesos instance.
  //  
  //
  if ( NumRows_ != 0 &&   NumRows_ != RowMatrixA_->NumGlobalRows() ) 
    EPETRA_CHK_ERR(-5);
  NumRows_ = RowMatrixA_->NumGlobalRows() ; 

  //
  //  Set the matrix and grid shapes and populate the matrix SuperluMat
  //
  if ( Redistribute_ ) {

    RedistributeA() ; 
    SuperluMat_ = UniformMatrix_ ;

  } else {
    // compute unless the used specified them in parameter list

    //  Revision 1.7 Oct 29, 2003 has a detailed check for cannonical distribution
    if( ! ( RowMatrixA_->RowMatrixRowMap().LinearMap())  ) EPETRA_CHK_ERR(-2);

    const Epetra_Comm &Comm_ = RowMatrixA_->Comm();
    int numprocs = Comm_.NumProc() ; 
    nprow_ = Superludist_NumProcRows( numprocs ) ; 
    npcol_ = numprocs / nprow_ ;
    assert ( nprow_ * npcol_ == numprocs ) ; 

    SuperluMat_ = RowMatrixA_ ;
  }

  //
  //  Extract Ai_, Ap_ and Aval_ from SuperluMat_
  //
  int MyActualFirstElement = SuperluMat_->RowMatrixRowMap().MinMyGID() ; 
  int NumMyElements = SuperluMat_->NumMyRows() ; 
  int nnz_loc = SuperluMat_->NumMyNonzeros() ;
  Ap_.resize( NumMyElements+1 );
  Ai_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  Aval_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  
  int NzThisRow ;
  int Ai_index = 0 ; 
  int MyRow;
  int num_my_cols = SuperluMat_->NumMyCols() ; 
  double *RowValues;
  int *ColIndices;
  int MaxNumEntries_ = SuperluMat_->MaxNumEntries();

  Global_Columns_.resize( num_my_cols ) ; 
  for ( int i = 0 ; i < num_my_cols ; i ++ ) { 
    Global_Columns_[i] = SuperluMat_->RowMatrixColMap().GID( i ) ; 
  }
  Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SuperluMat_);
  if ( SuperluCrs != 0 ) {
    ColIndicesV_.resize(MaxNumEntries_);
    RowValuesV_.resize(MaxNumEntries_);
  }
  for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
    if ( SuperluCrs != 0 ) {
      EPETRA_CHK_ERR( SuperluCrs->
		      ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					ColIndices ) != 0 ) ;
    }
    else {
      EPETRA_CHK_ERR( SuperluMat_->
		      ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					NzThisRow, &RowValuesV_[0], 
					&ColIndicesV_[0] ) != 0 );
      RowValues =  &RowValuesV_[0];
      ColIndices = &ColIndicesV_[0];
    }
    Ap_[MyRow] = Ai_index ; 
    for ( int j = 0; j < NzThisRow; j++ ) { 
      Ai_[Ai_index] = Global_Columns_[ColIndices[j]] ; 
      Aval_[Ai_index] = RowValues[j] ; 
      Ai_index++;
    }
  }
  assert( NumMyElements == MyRow );
  Ap_[ NumMyElements ] = Ai_index ; 

  //
  //  Setup Superlu's grid 
  //
  const Epetra_Comm &Comm_ = RowMatrixA_->Comm();
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm_);
  MPI_Comm MPIC = comm1.Comm() ;

  //
  //  Bug:  If nprow_ and npcol_ have changed since grid_ was created, we 
  //  could have a problem here.
  //
  if ( ! GridCreated_ ) {
    GridCreated_ = true; 

    superlu_gridinit( MPIC, nprow_, npcol_, &grid_);
  }

  if ( FactorizationDone_ ) {
    SUPERLU_FREE( SuperluA_.Store );
    ScalePermstructFree(&ScalePermstruct_);
    Destroy_LU(NumRows_, &grid_, &LUstruct_);
    LUstructFree(&LUstruct_);
    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
  }

  //
  //  Only those processes in the grid participate from here on
  //
  if ( iam_ < nprow_ * npcol_ ) {
    //
    //  Set up Superlu's data structures
    //
#ifdef OLD_SUPERLUDIST
    set_default_options(&options_);
#else
    set_default_options_dist(&options_);
#endif
    if( PrintNonzeros_ ) options_.PrintStat != YES;
    else                 options_.PrintStat = NO;

    int numcols = RowMatrixA_->NumGlobalCols() ; 
    if( NumRows_ != numcols ) EPETRA_CHK_ERR(-3) ; 



    dCreate_CompRowLoc_Matrix_dist( &SuperluA_, NumRows_, numcols, 
				    nnz_loc, NumMyElements, MyActualFirstElement,
				    &Aval_[0], &Ai_[0], &Ap_[0], 
				    SLU_NR_loc, SLU_D, SLU_GE );

    FactorizationDone_ = true;   // i.e. clean up Superlu data structures in the destructor

    ScalePermstructInit(NumRows_, numcols, &ScalePermstruct_);
    LUstructInit(NumRows_, numcols, &LUstruct_);

    // stick options from ParameterList to options_ structure
    // Here we follow the same order of the SuperLU_dist 2.0 manual (pag 55/56)

    
    assert( options_.Fact == DOFACT );  
    options_.Fact = DOFACT ;       

    if( Equil_ ) options_.Equil = (yes_no_t)YES;
    else         options_.Equil = NO;

    if( ColPerm_ == "NATURAL" ) options_.ColPerm = NATURAL;
    else if( ColPerm_ == "MMD_AT_PLUS_A" ) options_.ColPerm = MMD_AT_PLUS_A;
    else if( ColPerm_ == "MMD_ATA" ) options_.ColPerm = MMD_ATA;
    else if( ColPerm_ == "COLAMD" ) options_.ColPerm = COLAMD;
    else if( ColPerm_ == "MY_PERMC" ) {
      options_.ColPerm = MY_PERMC;
      ScalePermstruct_.perm_c = perm_c_;
    }

    if( RowPerm_ == "NATURAL" ) options_.RowPerm = (rowperm_t)NATURAL;
    if( RowPerm_ == "LargeDiag" ) options_.RowPerm = LargeDiag;
    else if( ColPerm_ == "MY_PERMR" ) {
      options_.RowPerm = MY_PERMR;
      ScalePermstruct_.perm_r = perm_r_;
    }

    if( ReplaceTinyPivot_ ) options_.ReplaceTinyPivot = (yes_no_t)YES;
    else                    options_.ReplaceTinyPivot = NO;

    if( IterRefine_ == "NO" ) options_.IterRefine = (IterRefine_t)NO;
    else if( IterRefine_ == "DOUBLE" ) options_.IterRefine = DOUBLE;
    else if( IterRefine_ == "EXTRA" ) options_.IterRefine = EXTRA;

    if( PrintNonzeros_ ) options_.PrintStat = (yes_no_t)YES;
    else                 options_.PrintStat = NO;
    
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */

    //
    //  Factor A using Superludsit (via a call to pdgssvx)
    //
    int info ;
    double berr ;    //  Should be untouched
    double xValues;  //  Should be untouched
    int nrhs = 0 ;   //  Prevents forward and back solves
    int ldx = NumRows_;     //  Should be untouched

    pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues, ldx, nrhs, &grid_,
	    &LUstruct_, &SOLVEstruct_, &berr, &stat, &info);

    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
    EPETRA_CHK_ERR( info ) ; 

    PStatFree(&stat);
  }

  return 0;
}

//
//   Refactor - Refactor the matrix 
//
//     Preconditions:
//       The non-zero pattern of the matrix must not have changed since the 
//         previous call to Factor().  Refactor ensures that each process owns 
//         the same number of columns that it did on the previous call to Factor()
//         and returns -4 if a discrepancy is found.  However, that check does not
//         guarantee that no change was made to the non-zero structure of the matrix.
//       No call to SetParameters should be made between the call to Factor()
//         and the call to Refactor().  If the user does not call SetParameters, 
//         as they need never do, they are safe on this.
//
//     Postconditions:
//       The matrix specified by Problem_->Operator() will have been redistributed,
//         converted to the form needed by Superludist and factored.
//       Ai_, Aval_ 
//       SuperluA_
//       SuperLU internal data structures reflecting the LU factorization
//         ScalePermstruct_
//         LUstructInit_
//       
//     Performance notes:
//       Refactor does not allocate or de-allocate memory.
//         
//     Return codes:
//       -4 if we detect a change to the non-zero structure of the matrix.
//
int Amesos_Superludist::ReFactor( ) {

  if( debug_ == 1 ) cout << "Entering `ReFactor()' ..." << endl;
  
    //
    //  Update Ai_ and Aval_ (while double checking Ap_)
    //
    if ( Redistribute_ ) { 
      if( UniformMatrix_->Export( *RowMatrixA_, *ExportToDist_, Insert ) ) 
	EPETRA_CHK_ERR(-4); 
    }
    
    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SuperluMat_);
    
    double *RowValues;
    int *ColIndices;
    int MaxNumEntries_ = SuperluMat_->MaxNumEntries();
    int NumMyElements = SuperluMat_->NumMyRows() ; 
    int NzThisRow ;
    int Ai_index = 0 ; 
    int MyRow;
    int num_my_cols = SuperluMat_->NumMyCols() ; 
    for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	EPETRA_CHK_ERR( SuperluMat_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_,
					  NzThisRow, &RowValuesV_[0], 
					  &ColIndicesV_[0] ) != 0 );
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }

      if ( Ap_[MyRow] != Ai_index ) EPETRA_CHK_ERR(-4);
      for ( int j = 0; j < NzThisRow; j++ ) { 
	//  pdgssvx alters Ai_, so we have to set it again.
	Ai_[Ai_index] = Global_Columns_[ColIndices[j]];
	Aval_[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    if( Ap_[ NumMyElements ] != Ai_index ) EPETRA_CHK_ERR(-4); 


    if ( iam_ < nprow_ * npcol_ ) {
	
#ifdef OLD_SUPERLUDIST
      set_default_options(&options_);
#else
      set_default_options_dist(&options_);
#endif
      
      options_.Fact = FactOption_;
      SuperLUStat_t stat;
      PStatInit(&stat);    /* Initialize the statistics variables. */
      int info ;
      double berr ;    //  Should be untouched
      double xValues;  //  Should be untouched
      int nrhs = 0 ;   //  Prevents forward and back solves
      int ldx = NumRows_;     //  Should be untouched
      pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues, ldx, nrhs, &grid_,
	      &LUstruct_, &SOLVEstruct_, &berr, &stat, &info);
      PStatFree(&stat);
      EPETRA_CHK_ERR( info ) ;
    } 
  
  return 0;
}

bool Amesos_Superludist::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Superludist::SymbolicFactorization() {

  if( debug_ == 1 ) cout << "Entering `SymbolicFactorization()' ..." << endl;
  
  FactorizationOK_ = false ; 

  return 0;
}

int Amesos_Superludist::NumericFactorization() {
  
  if( debug_ == 1 ) cout << "Entering `NumericFactorizationA()' ..." << endl;

  NumNumericFact_++;
  
  Epetra_Time Time(Comm());
  
  RowMatrixA_ = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if( RowMatrixA_ == 0 ) EPETRA_CHK_ERR(-1); 

  const Epetra_Comm &Comm_ = RowMatrixA_->Comm();
  iam_ = Comm_.MyPID() ;

  if ( FactorizationOK_ && ReuseSymbolic_ ) {
    ReFactor() ; 
  }  else { 
    Factor() ; 
    FactorizationOK_ = true;   

  }

  NumTime_ += Time.ElapsedTime();
  
  return 0;
}


int Amesos_Superludist::Solve() { 

  if( debug_ == 1 ) cout << "Entering `Solve()' ..." << endl;

  NumSolve_++;
  
  Epetra_Time Time(Comm());
  
  //  NRformat_loc *Astore;
  //  Astore = (NRformat_loc *) SuperluA_.Store;

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  //
  //  Compute the number of right hands sides (and check that X and B have the same shape) 
  //
  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  //
  //  Set up vecXptr and vecBptr 
  //
  double *bValues ;
  double *xValues ;
  int ldb, ldx ; 


  Epetra_MultiVector* vecXptr; 
  Epetra_MultiVector* vecBptr; 

  if ( Redistribute_ ) {
    //
    //  I couldn't figure out how to change the number of vectors in a multivector,
    //  so if we need more vectors in the distributed versions of X and B, I 
    //  delete them and rebuild them.  Ugly, but it isn't all that likely that codes
    //  will change the number of right hand sides repeatedly.
    //
    if ( vecXdistributed_ != 0 ) {
      assert( vecBdistributed_ != 0 ) ;
      if ( vecXdistributed_->NumVectors() != nrhs ) {
	delete vecXdistributed_ ; 
	delete vecBdistributed_ ; 
	vecXdistributed_ = 0 ; 
	vecBdistributed_ = 0 ; 
      } else {
	assert(  vecBdistributed_->NumVectors() == nrhs ) ;
      }
    }
    if ( vecXdistributed_ == 0 ) {
      vecXdistributed_ = new Epetra_MultiVector( *UniformMap_, nrhs ) ; 
      vecBdistributed_ = new Epetra_MultiVector( *UniformMap_, nrhs ) ; 
    }

    vecXdistributed_->Import( *vecX, *ImportToDistributed_, Insert ) ;
    vecBdistributed_->Import( *vecB, *ImportToDistributed_, Insert ) ;

    vecXptr = vecXdistributed_ ; 
    vecBptr = vecBdistributed_ ; 
  } else {
    vecXptr = vecX ; 
    vecBptr = vecB ; 
  }


  int NumMyElements = vecBptr->MyLength(); 
  EPETRA_CHK_ERR( vecBptr->ExtractView( &bValues, &ldb ) )  ; 
  EPETRA_CHK_ERR( vecXptr->ExtractView( &xValues, &ldx ) ) ; 

  //
  //  pdgssvx returns x in b, so we copy b into x.  
  //
  for ( int j = 0 ; j < nrhs; j++ )
    for ( int i = 0 ; i < NumMyElements; i++ ) xValues[i+j*ldx] = bValues[i+j*ldb]; 
  
  /* Bail out if I do not belong in the grid. */
  if ( iam_ < nprow_ * npcol_ ) {
    int info ;
    vector<double>berr(nrhs);
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */
    
    assert( GridCreated_ ) ; 
    assert( FactorizationDone_ ) ; 
    options_.Fact = FACTORED ;       

    bool BlockSolve = false ; 
    if ( BlockSolve ) { 
      pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues[0], ldx, nrhs, &grid_,
	      &LUstruct_, &SOLVEstruct_, &berr[0], &stat, &info);
      EPETRA_CHK_ERR( info ) ;
    } else {
      for ( int j =0 ; j < nrhs; j++ ) { 
	pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues[j*ldx], ldx, 1, &grid_,
		&LUstruct_, &SOLVEstruct_, &berr[0], &stat, &info);
	if ( options_.SolveInitialized ) {
	  dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
	}
	EPETRA_CHK_ERR( info ) ;
      }
    } 
    
    PStatFree(&stat);
  }
  
  if ( Redistribute_ ) { 
    vecX->Import( *vecXptr, *ImportBackToOriginal_, Insert ) ;
  }

  SolveTime_ += Time.ElapsedTime();
  return; 
}

void Amesos_Superludist::PrintStatus() 
{
  if( Comm().MyPID() ) return 0;
  
  cout << "----------------------------------------------------------------------------" << endl;

  cout << "Amesos_Superludist : Redistribute = " << Redistribute_ << endl;
  cout << "Amesos_Superludist : # available processes = " << Comm().NumProc() << endl;
  cout << "Amesos_Superludist : # processes used in computation = " << nprow_*npcol_
       << " ( = " << nprow_ << "x" << npcol_ << ")" << endl;
  cout << "Amesos_Superludist : Equil = " << Equil_ << endl;
  cout << "Amesos_Superludist : ColPerm = " << ColPerm_ << endl;
  cout << "Amesos_Superludist : RowPerm = " << RowPerm_ << endl;
  cout << "Amesos_Superludist : IterRefine = " << IterRefine_ << endl;
  cout << "Amesos_Superludist : ReplaceTinyPivot = " << ReplaceTinyPivot_ << endl;
  cout << "Amesos_Superludist : AddZeroToDiag = " << AddZeroToDiag_ << endl;
  cout << "Amesos_Superludist : Redistribute = " << Redistribute_ << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  return;
}


// ================================================ ====== ==== ==== == =

void Amesos_Superludist::PrintTiming()
{
  if( Comm().MyPID() ) return;
  
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Superludist : Number of symbolic factorizaions = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Superludist : Number of numeric factorizaions = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Superludist : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Superludist : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Superludist : Time for solve = "
       << SolveTime_ << " (s), avg = " << SolveTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
   
  return;
}
