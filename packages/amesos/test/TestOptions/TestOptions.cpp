//
//  To run this under valgrind, try:
//  valgrind --suppressions=../Test_Basic/Suppressions --gen-suppressions=yes --leak-check=yes --show-reachable=yes ./TestOptions.exe -v
//
//  To run this with valgrind under mpirun, 
//  mpirun -np 2 valgrind --log-file=TestOpt.logfile --suppressions=../Test_Basic/Suppressions --gen-suppressions=yes --leak-check=yes --show-reachable=yes ./TestOptions.exe -v
//
//  test/scripts/daily/serial/TestMemoryLeaks[.exe] performs an automated test for memory leaks
//  using valgrind and this code.  To run TestMemoryLeaks, cd to test/TestOptions in the
//  build directory and type ../scripts/daily/serial/TestMemoryLeaks.exe.  The output is stored
//  in amesos/logLinux.txt.
//
//

//  TestOptions tests all options for each Amesos Class on a limited number 
//  of matrices.  
//

//  TestOneMatrix - Test one matrix 
//    - Distributed vs not distributed -
//    - Transpose vs not transposed -
//      TestAllClasses  - Test one matrix and one setting of distributed and transpose
//        - Calls TestOtherClasses (one for each Amesos class) and TestSuperludist 
//        TestOtherClasses
//        TestSuperludist
//        TestScalapack
//
//
//  Todo:
//    Write TestKlu, TestSuperlu, TestScalapack, TestUmfpack, TestDscpack
//    Enable tests of various parameter options
//    Make it test all four codes (DSCPACK, UMFPACK, SuperLU_DIST, KLU )
//    Valgrind it
//    Enable tests of transpose and distributed matrices  - DONE 
//    Enable FACTOR_B - DONE 
//

#include "Trilinos_Util.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#include "CrsMatrixTranspose.h"
#include "TestAllClasses.h"
#include <string>
#include "Teuchos_RCP.hpp"
#include "NewMatNewMap.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#if 0

#ifdef HAVE_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#else
#ifdef HAVE_VALGRIND_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#endif 
#endif 

#endif 

std::vector<std::string> AmesosClasses;

int NumAmesosClasses;

int CreateCrsMatrix( char *in_filename, const Epetra_Comm &Comm, 
		     Epetra_Map *& readMap,
		     const bool transpose, const bool distribute, 
		     bool& symmetric, Epetra_CrsMatrix *& Matrix ) {

  Epetra_CrsMatrix * readA = 0; 
  Epetra_Vector * readx = 0; 
  Epetra_Vector * readb = 0;
  Epetra_Vector * readxexact = 0;

  //
  //  This hack allows TestOptions to be run from either the test/TestOptions/ directory or from 
  //  the test/ directory (as it is in nightly testing and in make "run-tests")
  //
  FILE *in_file = fopen( in_filename, "r");

  char *filename;
  if (in_file == NULL ) 
    filename = &in_filename[1] ; //  Strip off ithe "." from
				 //  "../" and try again 
  else {
    filename = in_filename ;
    fclose( in_file );
  }

  symmetric = false ; 
  std::string FileName = filename ;

  int FN_Size = FileName.size() ; 
  std::string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  std::string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );

  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( filename, false, Comm, readMap, readA, readx, 
						      readb, readxexact) );
    symmetric = false; 
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( filename, true, Comm, readMap, readA, readx, 
							readb, readxexact) );
      symmetric = true; 
    } else {
      if (  LastFourBytes == ".mtx" ) { 
	EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( filename, Comm, readMap, 
							       readA, readx, readb, readxexact) );   
	in_file = fopen( filename, "r");
	assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 
	const int BUFSIZE = 800 ; 
	char buffer[BUFSIZE] ; 
	fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this std::string 
	std::string headerline1 = buffer;
#ifdef TFLOP
	if ( headerline1.find("symmetric") < BUFSIZE ) symmetric = true;
#else
	if ( headerline1.find("symmetric") != std::string::npos) symmetric = true; 

#endif
	fclose(in_file);

      } else {
	// Call routine to read in HB problem
	Trilinos_Util_ReadHb2Epetra( filename, Comm, readMap, readA, readx, 
						     readb, readxexact) ;
	if (  LastFourBytes == ".rsa" ) symmetric = true ; 
      }
    }
  }


  if ( readb )  delete readb;
  if ( readx ) delete readx;
  if ( readxexact ) delete readxexact;

  Epetra_CrsMatrix *serialA ; 
  Epetra_CrsMatrix *transposeA;

  int ierr = 0;

  if ( transpose ) {
    transposeA = new Epetra_CrsMatrix( Copy, *readMap, 0 );
    ierr = CrsMatrixTranspose( readA, transposeA );
    assert( ierr == 0 ); 
    serialA = transposeA ; 
    delete readA;
    readA = 0 ; 
  } else {
    serialA = readA ; 
  }

  assert( (void *) &serialA->Graph() ) ;
  assert( (void *) &serialA->RowMap() ) ;
  assert( serialA->RowMap().SameAs(*readMap) ) ; 

  if ( distribute ) { 
    // Create uniform distributed map
    Epetra_Map* mapPtr = 0;
    
    if(readMap->GlobalIndicesInt())
      mapPtr = new Epetra_Map((int) readMap->NumGlobalElements(), 0, Comm);
    else if(readMap->GlobalIndicesLongLong())
      mapPtr = new Epetra_Map(readMap->NumGlobalElements(), 0, Comm);
    else
      assert(false);
    
    Epetra_Map& DistMap = *mapPtr;
    
    // Create Exporter to distribute read-in matrix and vectors
    Epetra_Export exporter( *readMap, DistMap );
    
    Epetra_CrsMatrix *Amat = new Epetra_CrsMatrix( Copy, DistMap, 0 );
    Amat->Export(*serialA, exporter, Add);
    ierr = Amat->FillComplete();
    assert(ierr == 0);    
    
    Matrix = Amat; 
    //
    //  Make sure that deleting Amat->RowMap() will delete map 
    //
    //  Bug:  We can't manage to delete map his way anyway,
    //        and this fails on tranposes, so for now I just accept
    //        the memory loss.
    //    assert( &(Amat->RowMap()) == map ) ; 
    delete readMap; 
    readMap = 0 ; 
    delete serialA; 
    delete mapPtr;
  } else { 

    Matrix = serialA; 
  }


  return 0;
}

int TestErrors( const std::vector<bool> AmesosClassesInstalled, 
		   char *filename, 
#ifdef EPETRA_MPI
		   Epetra_MpiComm& Comm,
#else
		   Epetra_SerialComm& Comm,
#endif
		   const bool verbose, 
		   int &NumTests  ) {

  int NumErrors =0 ;
  double error = -13; 
  double residual = -13;


  for ( int iterTrans =0 ; iterTrans < 2; iterTrans++ ) {
    const bool transpose = iterTrans == 1 ; 
    
    const bool distribute = 1;

    const int iterRowindex = 0;
    const int iterColindex = 0 ;
    const int iterRangemap = 0 ;
    const int iterDomainmap = 0 ;
    const int iterDiagonalOpts = 0 ; 
    bool printit = true ; 
    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
    const int EpetraMatrixType = 0 ;
    bool symmetric = true;
    const int iterDist = 0 ; 
    
    Epetra_CrsMatrix *Amat = 0 ;
    Epetra_Map *readMap = 0 ;
    CreateCrsMatrix( filename, Comm, readMap, transpose, distribute, symmetric, Amat ) ;
    
    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << 
				" Creating matrix from " <<
				" filename = " << filename <<
				" symmetric = " << symmetric <<
				" distribute = " << distribute <<
				" iterRowindex = " << iterRowindex <<
				" iterColindex = " << iterColindex <<
				" iterRangemap = " << iterRangemap <<
				" iterDomainmap = " << iterDomainmap <<
				" EpetraMatrixType = " << EpetraMatrixType <<
				" iterDiagonalOpts = " << iterDiagonalOpts <<
				" transpose = "  << transpose 
				   << " iterDist = " << iterDist << std::endl ; 
    
    if ( iterDiagonalOpts )  Comm.SetTracebackMode(1);  // Turn off positive Epetra warnings (i.e. iniefficient code, such as memory re-allocation)
    
    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
    
    RCP<Epetra_CrsMatrix> Bmat = NewMatNewMap( *Amat, 
						       iterDiagonalOpts, 
						       iterRowindex,
						       iterColindex,
						       iterRangemap,
						       iterDomainmap
						       ) ; 
    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
    Comm.SetTracebackMode(2);
    
    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << 
				" filename = " << filename <<
				" symmetric = " << symmetric <<
				" distribute = " << distribute <<
				" iterRowindex = " << iterRowindex <<
				" iterColindex = " << iterColindex <<
				" iterRangemap = " << iterRangemap <<
				" iterDomainmap = " << iterDomainmap <<
				" EpetraMatrixType = " << EpetraMatrixType <<
				" iterDiagonalOpts = " << iterDiagonalOpts <<
				" transpose = "  << transpose << " iterDist = " << iterDist << std::endl ; 
    
    //
    //  This causes a failure in Amesos_Superludist:
    Epetra_CrsMatrix* Cmat = &*Bmat;
    //  Epetra_CrsMatrix* Cmat = Amat ;
    
    
    const int Level = 1; 
    const double MaxError = 1e-3;
    
    int NumTheseTests = 0 ; 
    if ( verbose ) {
      std::cout << " About to test  " << filename 
	   << __FILE__ << "::"  << __LINE__
	   << " EpetraMatrixType = " <<  EpetraMatrixType 
	   << (transpose?" transpose":"" ) 
	   << (distribute?" distribute":"" )
	   << std::endl ; 
    }
    int Error = TestAllClasses( AmesosClasses, EpetraMatrixType, 
				AmesosClassesInstalled, 
				Cmat, 
				transpose ,
				verbose, 
				symmetric, 
				Level,
				MaxError, 
				iterDiagonalOpts, 
				iterRowindex,
				iterColindex,
				iterRangemap,
				iterDomainmap,
				distribute,
				filename,
				error, 
				residual, 
				NumTheseTests ) ;
    NumTests += NumTheseTests ;
    NumErrors += Error ;
    //      BUG:  Memory leak 
    //      delete &(Amat->RowMap()) ; 
    if ( Amat ) delete Amat ; 
    if ( readMap ) delete readMap ; 
  }
  if ( verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 

  return NumErrors;
} 

int TestOneMatrix( const std::vector<bool> AmesosClassesInstalled, 
		   char *filename, 
#ifdef EPETRA_MPI
		   Epetra_MpiComm& Comm,
#else
		   Epetra_SerialComm& Comm,
#endif
		   //		   Epetra_Comm &Comm, 
		   const bool verbose, 
		   const bool PerformDiagonalTest, 
		   double Rcond,
		   int &NumTests  ) {

  int NumErrors =0 ;
  double error = -13; 
  double residual = -13;
  //  double errors[NumAmesosClasses];
  //  double residuals[NumAmesosClasses];
  //  for (int i = 0 ; i < NumAmesosClasses; i ++ ) errors[i] = residuals[i] = 0.0 ; 

  //#ifdef HAVE_AMESOS_UMFPACK
#if 0
  Epetra_CrsMatrix *Amat ;

  //
  //  Compute the reciprocal condition number using Amesos_UMFPACK via the Amesos interface
  //
  Epetra_Map *readMap = 0 ; 
  CreateCrsMatrix( filename, Comm, readMap, false, false, symmetric, Amat ) ;
  Teuchos::ParameterList ParamList ;
  Epetra_LinearProblem Problem;
  Amesos Afactory;

  Amesos_BaseSolver* Abase ; 
  Abase = Afactory.Create( "Amesos_Umfpack", Problem ) ; 
  if ( Abase == 0 ) {
    std::cerr << " AMESOS_UMFPACK is required for this test " << std::endl ;
    exit(13);
  }  ;

  //
  //  Factor A to compute Rcond = reciprocal condition number estimate
  //
  Problem.SetOperator( Amat );
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
  Amesos_Umfpack* UmfpackOperator = dynamic_cast<Amesos_Umfpack *> (Abase) ; 
  //  double Rcond = UmfpackOperator->GetRcond();

  int ind[1];
  double val[1];
  ind[0] = 0;
  val[0] = 1 ; 
  double AnormInf =  Amat->NormInf() ;
  if (verbose) std::cout << " norm(Amat) = " << AnormInf << std::endl; 
  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  AnormInf =  Amat->NormInf() ;
  if (verbose) std::cout << " norm(Amat) = " << AnormInf << std::endl; 


  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
  double Rcond1 = UmfpackOperator->GetRcond();

  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  AnormInf =  Amat->NormInf() ;
  if (verbose) std::cout << " norm(Amat) = " << AnormInf << std::endl; 
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
   double Rcond2 = UmfpackOperator->GetRcond();

  if (verbose) std::cout << " Rcond1 = " << Rcond1 << std::endl; 
  if (verbose) std::cout << " Rcond2 = " << Rcond2 << std::endl; 

  if ( readMap ) delete readMap ;
#else
  double Rcond1 = Rcond ;
  double Rcond2 = Rcond ;
#endif

  //
  //  Rowindex and Colindex control the maps and indices used to create the matrix
  //
  //  These tests are all disabled in TestAllClasses.cpp
  //
  const int RowindexMax = 3;   // bug should be three ( 1 based, 3 based, non contiguous )
  const int ColindexMax = 2;   // bug should be two:  ( row map, 4 based )

  //
  //  Rangemap and Domainmap control the Range and Domain maps used in the call to FillComplete
  //  If both are "no change", FillComplete is called with no parameters (i.e. without specifying maps)
  //  Otherwise, domain and range maps are specified in the call to FillComplete
  //
  //  These tests are all disabled in TestAllClasses.cpp
  //
  int RangemapMax = 4; // bug should be four:  ( no change, serial, bizarre dist, replicated )
  int DomainmapMax = 4; // bug should be four:  ( no change, serial, bizarre dist, replicated )  IRRELEVANT see ThisDomainMax  

  //
  //  DiagonalOpts controls whether diagonal elements are left alone,
  //  or removed from both the matrix and the underlying map
  //
  int DiagonalOptsMax = 2;   // should be two:  ( no change, elements missing from map )
  //
  //
  //
  int EpetraMatrixTypeMax = 3; // 0 = Epetra_CrsMatrix; 1 = Epetra_RowMatriw; 2 = StorageOptimized Epetra_CrsMatrix
  //
  //  No point in trying to run distributed memory tests on a serial run
  //
  int iterDistMax = 2;
  if ( Comm.NumProc() == 1 ) {
    iterDistMax = 1 ; 
    RangemapMax = 1 ; 
    DomainmapMax = 1 ; 
  }

  

  if (! PerformDiagonalTest ) DiagonalOptsMax = 1 ; 

  for ( int iterTrans =0 ; iterTrans < 2; iterTrans++ ) {
    bool transpose = iterTrans == 1 ; 
    
    for ( int iterDist =0; iterDist < iterDistMax; iterDist++ ) {  
      bool distribute = ( iterDist == 1 ); 
	
#if 1
      for ( int iterRowindex = 0 ; iterRowindex < RowindexMax; iterRowindex++ ) {
	for ( int iterColindex = 0 ; iterColindex < ColindexMax; iterColindex++ ) {
	  //
	  //  The current version of NewMatNewMap.cpp supports only trivial 
	  //  replicated maps, hence we do not allow any fancy indexing
	  //
	  int ThisRangemapMax = RangemapMax ;
              // Bug #1920 Amesos classes can't handle replicated domain or ranges  if ( iterRowindex > 0 || iterColindex > 0 ) 
	  ThisRangemapMax = EPETRA_MIN( 3, ThisRangemapMax );
	  int ThisDomainmapMax =  EPETRA_MIN( 3, ThisRangemapMax );  // Bug #1920 
	  for ( int iterRangemap = 0 ; iterRangemap < ThisRangemapMax; iterRangemap++ ) {
	    for ( int iterDomainmap = 0 ; iterDomainmap < ThisDomainmapMax; iterDomainmap++ ) {
	      for ( int iterDiagonalOpts = 0 ; iterDiagonalOpts < DiagonalOptsMax; iterDiagonalOpts++ ) { 
#else
		int iterRowindex = 0; { 
		  int iterColindex = 0 ; { 
		    int iterRangemap = 0 ; { 
		      int iterDomainmap = 0 ; {
			for ( int iterDiagonalOpts = 1 ; iterDiagonalOpts < DiagonalOptsMax; iterDiagonalOpts++ ) { 
#endif
			  const bool printit = false; 
		//  diagonal opts testing only works on distributed matrices whose row and column indices match 
		//  On a serial matrix, eliminate a column from the map makes the matrix singular
		//  If the row and column indices don't match, eliminating a column from the map is, typically, irrelevant

		if ( ( iterColindex == 0 && distribute ) || iterDiagonalOpts == 0 ) { 
		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
		  for ( int EpetraMatrixType = 0 ; EpetraMatrixType < EpetraMatrixTypeMax;  EpetraMatrixType++ ) {

		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
		    //  These tests presently take over 7 hours on some platforms.  But, I don't want to eliminate any category of tests
		    //  The following test will cull 90% of the tests and still cover every type of test and most combinations
		    Epetra_Util EU;
		    int RunTest[1] ;
		    RunTest[0] =  (EU.RandomDouble() > 0.8)?1:0 ; 
		    if ( iterRowindex == 0 &&
			 iterColindex == 0 &&
			 iterRangemap == 0 &&
			 iterDomainmap == 0 ) RunTest[0] = 1; 
		    Comm.Broadcast( RunTest, 1, 0 ) ; 
		    if ( RunTest[0] ) { 
		    //
		    //  We test only one level for different indexing or different Range and Domain maps
		    //  to avoid hassles of moving data from the domain space to the range space and back
		    //
		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
		    int MaxLevel = 3 ; 
		    if ( iterRowindex > 0 ) MaxLevel = 1 ; 
		    if ( iterColindex > 0 ) MaxLevel = 1 ; 
		    if ( iterRangemap > 0 ) MaxLevel = 1 ; 
		    if ( iterDomainmap > 0 ) MaxLevel = 1 ; 

		    bool symmetric = true;
		    
		    Epetra_CrsMatrix *Amat = 0 ;
		    Epetra_Map *readMap = 0 ;
		    CreateCrsMatrix( filename, Comm, readMap, transpose, distribute, symmetric, Amat ) ;
		    //		  assert( symmetric == false ) ; 
		    
		    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << 
				     " Creating matrix from " <<
				     " filename = " << filename <<
				     " symmetric = " << symmetric <<
				     " distribute = " << distribute <<
				     " iterRowindex = " << iterRowindex <<
				     " iterColindex = " << iterColindex <<
				     " iterRangemap = " << iterRangemap <<
				     " iterDomainmap = " << iterDomainmap <<
				     " EpetraMatrixType = " << EpetraMatrixType <<
				     " iterDiagonalOpts = " << iterDiagonalOpts <<
				     " transpose = "  << transpose << " iterDist = " << iterDist << std::endl ; 


		    if ( iterDiagonalOpts )  Comm.SetTracebackMode(1);  // Turn off positive Epetra warnings (i.e. iniefficient code, such as memory re-allocation)

		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 

		    RCP<Epetra_CrsMatrix> Bmat = NewMatNewMap( *Amat, 
								       iterDiagonalOpts, 
								       iterRowindex,
								       iterColindex,
								       iterRangemap,
								       iterDomainmap
								       ) ; 
		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
		    Comm.SetTracebackMode(2);

		    if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << 
				     " filename = " << filename <<
				     " symmetric = " << symmetric <<
				     " distribute = " << distribute <<
				     " iterRowindex = " << iterRowindex <<
				     " iterColindex = " << iterColindex <<
				     " iterRangemap = " << iterRangemap <<
				     " iterDomainmap = " << iterDomainmap <<
				     " EpetraMatrixType = " << EpetraMatrixType <<
				     " iterDiagonalOpts = " << iterDiagonalOpts <<
				     " transpose = "  << transpose << " iterDist = " << iterDist << std::endl ; 

		    //
		    //  This causes a failure in Amesos_Superludist:
		    Epetra_CrsMatrix* Cmat = &*Bmat;
		    //  Epetra_CrsMatrix* Cmat = Amat ;
		 

		    int Level ; 
		    double MaxError ;
		    if ( Rcond*Rcond1*Rcond2 > 1e-16 ) { 
		      Level = EPETRA_MIN( 3, MaxLevel );
		      MaxError = Rcond*Rcond1*Rcond2;
		    } else if  ( Rcond*Rcond1 > 1e-16 ) {
		      Level = EPETRA_MIN( 2, MaxLevel );
		      MaxError = Rcond*Rcond1;
		    } else {
		      Level = EPETRA_MIN( 1, MaxLevel );
		      MaxError = Rcond;
		    }

		    int NumTheseTests = 0 ; 
		    if ( verbose ) {
		      std::cout << " About to test  " << filename 
			   << __FILE__ << "::"  << __LINE__
			   << " EpetraMatrixType = " <<  EpetraMatrixType 
			   << (transpose?" transpose":"" ) 
			   << (distribute?" distribute":"" )
			   << std::endl ; 
		    }
		    if ( iterDiagonalOpts == 0 ) 
		      Comm.SetTracebackMode(2);
		    else
		      Comm.SetTracebackMode(1);  // In PerformOneSolveAndTest, MyMatWithDiag->ReplaceDiagonalValues may return 1 indicating that structurally non-zero elements were left untouched.

		    int Error = TestAllClasses( AmesosClasses, EpetraMatrixType, 
						AmesosClassesInstalled, 
						Cmat, 
						transpose ,
						verbose, 
						symmetric, 
						Level,
						MaxError, 
						iterDiagonalOpts, 
						iterRowindex,
						iterColindex,
						iterRangemap,
						iterDomainmap,
						distribute,
						filename,
						error, 
						residual, 
						NumTheseTests ) ;
		    NumTests += NumTheseTests ;
		    NumErrors += Error ;
		    if ( Comm.MyPID() == 0  && ( ( verbose && NumTheseTests ) || Error ) ) {
		      std::cout << " Tested  " << filename 
			   << __FILE__ << "::"  << __LINE__
			   << " EpetraMatrixType = " <<  EpetraMatrixType 
			   << (transpose?" transpose":"" ) 
			   << (distribute?" distribute":"" ) << " error = " 
			   << error 
			   << " residual = " 
			   << residual 
			   << std::endl ; 
		    }
		    //      BUG:  Memory leak 
		    //      delete &(Amat->RowMap()) ; 
		    if ( Amat ) delete Amat ; 
		    if ( readMap ) delete readMap ; 
#if 0
		    double relresidual = 
		      errors[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( errors[ (int) AMESOS_SUPERLUDIST], error ) ; 
		    residuals[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( residuals[ (int) AMESOS_SUPERLUDIST], residual ) ; 
		    NumErrors += ( residual > maxresidual ) ; 
#endif
		    }
		  }
		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
		}
		  if ( printit && verbose ) std::cout << __FILE__ << "::" << __LINE__ << std::endl ; 
	      }
	    }
	  }
	}
      }
    }
  }

  return NumErrors;
} 

#if 0
#define TEST_P(variable) { { \
                      if ( true ) { std::cerr << "AMESOS_PRINT " << # variable << "= " << variable << std::endl; };  }\
                   }


#define TEST_PRINT(variable) { { \
                      if ( true ) { std::cerr << "AMESOS_PRINT " << # variable  << "= " << variable <<  ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; };  }\
                   }

#endif

//
//  Usage:  TestOptions [-s] [-v] [-q]
//
//  -s = short
//  -v = verbose
//  -q = quiet 
//

int NextMain( int argc, char *argv[] ) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif


  bool verbose = false; 
  bool small = false ; 
  bool quiet = false ; 
  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 'v') ) 
    verbose = true ; 
  if ( argc >= 3 && (argv[2][0] == '-') &&  (argv[2][1] == 'v') ) 
    verbose = true ; 
  if ( argc >= 4 && (argv[3][0] == '-') &&  (argv[3][1] == 'v') ) 
    verbose = true ; 

  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 's') ) 
    small = true ; 
  if ( argc >= 3 && (argv[2][0] == '-') &&  (argv[2][1] == 's') ) 
    small = true ; 
  if ( argc >= 4 && (argv[3][0] == '-') &&  (argv[3][1] == 's') ) 
    small = true ; 

  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 'q') ) 
    quiet = true ; 
  if ( argc >= 3 && (argv[2][0] == '-') &&  (argv[2][1] == 'q') ) 
    quiet = true ; 
  if ( argc >= 4 && (argv[3][0] == '-') &&  (argv[3][1] == 'q') ) 
    quiet = true ; 


  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 'h') ) {
    std::cerr << "Usage TestOptions [-s] [-v] [-q] " << std::endl ; 
    std::cerr << "-v:  verbose  " << std::endl ; 
    std::cerr << "-s:  small  " << std::endl ; 
    std::cerr << "-q:  quiet  " << std::endl ; 
    exit(-1);
  }



#ifdef HAVE_AMESOS_KLU
  AmesosClasses.push_back( "Amesos_Klu" );
#endif

#if 1

#ifdef HAVE_AMESOS_PARAKLETE
  AmesosClasses.push_back( "Amesos_Paraklete" );
#endif





#ifdef HAVE_AMESOS_PARDISO
  //  bug #1915  
  //  bug #1998 - Enabling Amesos_Pardiso causes Amesos_Klu to fail - strange but true
  //  AmesosClasses.push_back( "Amesos_Pardiso" );
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  AmesosClasses.push_back( "Amesos_Superludist" );
#endif




#ifdef HAVE_AMESOS_LAPACK
  AmesosClasses.push_back( "Amesos_Lapack" );
#endif

#ifdef HAVE_AMESOS_SUPERLU
  AmesosClasses.push_back( "Amesos_Superlu" );
#endif

#ifdef HAVE_AMESOS_TAUCS
  AmesosClasses.push_back( "Amesos_Taucs" );
#endif

#ifdef HAVE_AMESOS_UMFPACK
  AmesosClasses.push_back( "Amesos_Umfpack" );
#endif



#ifdef HAVE_AMESOS_DSCPACK
  if ( ! quiet ) AmesosClasses.push_back( "Amesos_Dscpack" );         //  bug #1205 
#endif

#ifdef HAVE_AMESOS_MUMPS
  AmesosClasses.push_back( "Amesos_Mumps" );
#endif

#ifdef HAVE_AMESOS_SCALAPACK
  AmesosClasses.push_back( "Amesos_Scalapack" ) ;
#endif



#endif

  NumAmesosClasses = AmesosClasses.size();
  std::vector<bool> AmesosClassesInstalled( NumAmesosClasses );

  assert( NumAmesosClasses > 0 ) ; 


  if ( Comm.MyPID() != 0 ) verbose = false ; 
#if 0
  //
  //  Wait for a character to allow time to attach the debugger
  //
  if ( Comm.MyPID() == 0 ) {
    char what = 'a'; 
    while ( what == 'a' )    // I don't know why we need this while loop  at least on bewoulf
      std::cin >> what ; 
  }
  Comm.Barrier();

  std::cout << __FILE__ << "::" << __LINE__ << " Comm.MyPID() = "  << Comm.MyPID() << std::endl ; 
#endif



  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  Comm.SetTracebackMode(2);

#ifndef HAVE_AMESOS_EPETRAEXT
    if ( ! quiet && Comm.MyPID() == 0 ) 
      std::cout << "Amesos has been built without epetraext, capabilites requiring epetraext, such as reindexing and Amesos_Paraklete non-transpose solves, will not be tested" << std::endl ; 
#endif

  for (int i=0; i < NumAmesosClasses; i++ ) {


    Abase = Afactory.Create( &AmesosClasses[i][0], Problem ) ; 
    if ( Abase == 0 ) {
      if ( !quiet && Comm.MyPID() == 0  ) std::cout << AmesosClasses[i] << " not built in this configuration"  << std::endl ;
      AmesosClassesInstalled[i] = false;
    } else {
      if (  !quiet && Comm.MyPID() == 0  ) std::cout << " Testing " << AmesosClasses[i] << std::endl ;
      AmesosClassesInstalled[i] = true;
      Teuchos::ParameterList ParamList ;
      ParamList.set( "NoDestroy", true );    // Prevents Amesos_Mumps from deleting data
      Abase->SetParameters( ParamList );     // which causes Amesos_Mumps to crash  on this trivial instantiation 
    }
    delete Abase ; 
    }

  int result = 0 ; 
  int numtests = 0 ;

  //  ImpcolB.rua fails - the failure could be in the test code, in particular in NewMatNewMap.cpp
  //  result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/ImpcolB.rua", Comm, verbose, false, 1e-9 , numtests ) ;

  //
  //  Khead.triS remains non-singular even after a diagaonal element is removed from the map 
  //
  // Khead.triS fails on DSCPACK 
  result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/SuperLU.triU", Comm, verbose, false, 1e-6 , numtests ) ;

  //
  //  small is set by TestValgrind - keep testing to a minimum because execution time is so slow
  //  quiet is set by TestQuietAmesos - dscpack is not quiet at the moment, hence we can't test symmetric matrices
  //  in TestQuietAmesos
  //
  //  bug #1205 
  //
  if ( ! small ) {
  result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/MissingADiagonal.mtx", Comm, verbose, false, 1e-2 , numtests ) ;
    result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/Khead.triS", Comm, verbose, true, 1e-6 , numtests ) ;
    result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/bcsstk04.mtx", Comm, verbose, false, 1e-4 , numtests ) ;
    //
    //  The file reader for .rua files is not quiet 
    //
    if (! quiet) { 
      result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/Diagonal.mtx", Comm, verbose, false, 1e-1 , numtests ) ;
      if ( Comm.NumProc() == 1) {
	result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/662_bus_out.rsa", Comm, verbose, false, 1e-5 , numtests ) ;
	result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/SuperLU.rua", Comm, verbose, false, 1e-2 , numtests ) ;
	result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/ImpcolB.rua", Comm, verbose, false, 1e-6 , numtests ) ;
	//	result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/ImpcolB.rua", Comm, verbose, false, 1e-6 , numtests ) ;
      }
    }
  }

  // bug #2184 - Amesos_Klu fails to detect Structurally singular matrices 
  // This test, as modified in PerformOneSolveAndTest.cpp, tests the present
  // beahviour - i.e. that SymbolicFactorization() fails to detect the
  // structurally singular matrix, but that NumericFactorization() 
  // catches the singularity instead.  
  result+=TestErrors( AmesosClassesInstalled, (char *) "../Test_Basic/StructurallySingular.mtx",
		      Comm, verbose, numtests ) ;
  result+=TestErrors( AmesosClassesInstalled, (char *) "../Test_Basic/NumericallySingular.mtx",
		      Comm, verbose, numtests ) ;
 
  if ( ! quiet && Comm.MyPID() == 0 ) std::cout << result << " Tests failed "  << numtests << " Tests performed " << std::endl ; 

  if ( result == 0 && numtests > 0 ) {
    if (! quiet && Comm.MyPID() == 0)
      std::cout << std::endl << "TEST PASSED" << std::endl << std::endl;
  }
  else {
    if (Comm.MyPID() == 0)
      std::cout << std::endl << "TEST FAILED" << std::endl << std::endl;
    AMESOS_CHK_ERR( 1 ) ; 
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return result ; 
}

//
//  I put this in hoping that this would eliminate a bogus memory leak report 
//  from valgrind. 
//
int main( int argc, char *argv[] ) {
  int retval = NextMain( argc, argv ) ; 
  return retval ;
}
