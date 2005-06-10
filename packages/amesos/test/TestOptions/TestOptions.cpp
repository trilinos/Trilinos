//
//  To run this under valgrind, try:
//  valgrind --suppressions=../Test_Basic/Suppressions --gen-suppressions=yes --leak-check=yes --show-reachable=yes ./TestOptions.exe -v
//
//  To run this with valgrind under mpirun, 
//  mpirun -np 2 valgrind --logfile=TestOpt.logfile --suppressions=../Test_Basic/Suppressions --gen-suppressions=yes --leak-check=yes --show-reachable=yes ./TestOptions.exe -v
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
//  TestOptions - Calls TestOneMatrix for each of several matrices
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
#include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Amesos_Umfpack.h"
#include "CrsMatrixTranspose.h"
#include "TestAllClasses.h"
#include <string>
#include "Teuchos_RefCountPtr.hpp"
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

vector<string> AmesosClasses;

int NumAmesosClasses;

int CreateCrsMatrix( char *filename, const Epetra_Comm &Comm, 
		     Epetra_Map *& readMap,
		     const bool transpose, const bool distribute, 
		     bool& symmetric, Epetra_CrsMatrix *& Matrix ) {

  Epetra_CrsMatrix * readA = 0; 
  Epetra_Vector * readx = 0; 
  Epetra_Vector * readb = 0;
  Epetra_Vector * readxexact = 0;

  symmetric = false ; 
  string FileName = filename ;
  int FN_Size = FileName.size() ; 
  string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );

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
	FILE* in_file = fopen( filename, "r");
	assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 
	const int BUFSIZE = 800 ; 
	char buffer[BUFSIZE] ; 
	fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this string 
	string headerline1 = buffer;
#ifdef TFLOP
	if ( headerline1.find("symmetric") < BUFSIZE ) symmetric = true;
#else
	if ( headerline1.find("symmetric") != string::npos) symmetric = true; 

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

  if ( transpose ) {
    transposeA = new Epetra_CrsMatrix( Copy, *readMap, 0 );
    assert( CrsMatrixTranspose( readA, transposeA ) == 0 ); 
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
    Epetra_Map DistMap(readMap->NumGlobalElements(), 0, Comm);

    // Create Exporter to distribute read-in matrix and vectors
    Epetra_Export exporter( *readMap, DistMap );
    
    Epetra_CrsMatrix *Amat = new Epetra_CrsMatrix( Copy, DistMap, 0 );
    Amat->Export(*serialA, exporter, Add);
    assert(Amat->FillComplete()==0);    
    
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
  } else { 

    Matrix = serialA; 
  }


  return 0;
}

int TestOneMatrix( const vector<bool> AmesosClassesInstalled, 
		   char *filename, 
		   Epetra_Comm &Comm, 
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
    cerr << " AMESOS_UMFPACK is required for this test " << endl ;
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
  if (verbose) cout << " norm(Amat) = " << AnormInf << endl; 
  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  AnormInf =  Amat->NormInf() ;
  if (verbose) cout << " norm(Amat) = " << AnormInf << endl; 


  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
  double Rcond1 = UmfpackOperator->GetRcond();

  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  AnormInf =  Amat->NormInf() ;
  if (verbose) cout << " norm(Amat) = " << AnormInf << endl; 
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
   double Rcond2 = UmfpackOperator->GetRcond();

  if (verbose) cout << " Rcond1 = " << Rcond1 << endl; 
  if (verbose) cout << " Rcond2 = " << Rcond2 << endl; 

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
  const int RowindexMax = 3;   // should be three ( 1 based, 3 based, non contiguous )
  const int ColindexMax = 2;   // should be two:  ( row map, 4 based )

  //
  //  Rangemap and Domainmap control the Range and Domain maps used in the call to FillComplete
  //  If both are "no change", FillComplete is called with no parameters (i.e. without specifying maps)
  //  Otherwise, domain and range maps are specified in the call to FillComplete
  //
  //  These tests are all disabled in TestAllClasses.cpp
  //
  int RangemapMax = 3; // should be three:  ( no change, serial, bizarre dist )
  int DomainmapMax = 1; // should be three:  ( no change, serial, bizarre dist )

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

  for ( int iterTrans =0 ; iterTrans < 1; iterTrans++ ) {
    bool transpose = iterTrans == 1 ; 
    
    for ( int iterDist =0 ; iterDist < iterDistMax; iterDist++ ) {
      bool distribute = ( iterDist == 1 ); 
	
      //
      for ( int iterRowindex = 0 ; iterRowindex < RowindexMax; iterRowindex++ ) {
	for ( int iterColindex = 0 ; iterColindex < ColindexMax; iterColindex++ ) {
	  for ( int iterRangemap = 0 ; iterRangemap < RangemapMax; iterRangemap++ ) {
	    for ( int iterDomainmap = 0 ; iterDomainmap < DomainmapMax; iterDomainmap++ ) {
	      for ( int iterDiagonalOpts = 0 ; iterDiagonalOpts < DiagonalOptsMax; iterDiagonalOpts++ ) {
		//  diagonal opts testing only works on distributed matrices whose row and column indices match 
		//  On a serial matrix, eliminate a column from the map makes the matrix singular
		//  If the row and column indices don't match, eliminating a column from the map is, typically, irrelevant
		if ( ( iterColindex == 0 && distribute ) || iterDiagonalOpts == 0 ) { 
		  for ( int EpetraMatrixType = 0 ; EpetraMatrixType < EpetraMatrixTypeMax;  EpetraMatrixType++ ) {
		    if ( verbose ) cout << __FILE__ << "::" << __LINE__ << 
				     " distribute = " << distribute <<
				     " iterRowindex = " << iterRowindex <<
				     " iterColindex = " << iterColindex <<
				     " iterRangemap = " << iterRangemap <<
				     " iterDomainmap = " << iterDomainmap <<
				     " EpetraMatrixType = " << EpetraMatrixType <<
				     " iterDiagonalOpts = " << iterDiagonalOpts <<
				     " transpose = "  << transpose << " iterDist = " << iterDist << endl ; 

		    bool symmetric = true;
		    
		    Epetra_CrsMatrix *Amat = 0 ;
		    Epetra_Map *readMap = 0 ;
		    CreateCrsMatrix( filename, Comm, readMap, transpose, distribute, symmetric, Amat ) ;
		    //		  assert( symmetric == false ) ; 
		    
		    RefCountPtr<Epetra_CrsMatrix> Bmat = NewMatNewMap( *Amat, 
								       iterDiagonalOpts, 
								       iterRowindex,
								       iterColindex,
								       iterRangemap,
								       iterDomainmap
								       ) ; 
		    
		    //
		    //  This causes a failure in Amesos_Superludist:
		    Epetra_CrsMatrix* Cmat = &*Bmat;
		    //  Epetra_CrsMatrix* Cmat = Amat ;
		 
   
		    if ( Rcond*Rcond1*Rcond2 > 1e-16 ) 
		      { 
			NumErrors += TestAllClasses( AmesosClasses, EpetraMatrixType, 
						     AmesosClassesInstalled, 
						     Cmat, 
						     transpose, 
						     verbose, 
						     symmetric, 
						     3, 
						     Rcond*Rcond1*Rcond2, 
						     iterDiagonalOpts, 
						     iterRowindex,
						     iterColindex,
						     iterRangemap,
						     iterDomainmap,
						     distribute,
						     filename,
						     error, 
						     residual, 
						     NumTests ) ;
		      }
		    else if ( Rcond*Rcond1 > 1e-16 ) 
		      {
			NumErrors += TestAllClasses( AmesosClasses, EpetraMatrixType, 
						     AmesosClassesInstalled, 
						     Cmat, 
						     transpose, 
						     verbose, 
						     symmetric, 
						     2, 
						     Rcond*Rcond1, 
						     iterDiagonalOpts, 
						     iterRowindex,
						     iterColindex,
						     iterRangemap,
						     iterDomainmap,
						     distribute,
						     filename,
						     error, 
						     residual, 
						     NumTests ) ;
		      }
		    else
		      {
			NumErrors += TestAllClasses( AmesosClasses, EpetraMatrixType, 
						     AmesosClassesInstalled, 
						     Cmat, 
						     transpose, 
						     verbose, 
						     symmetric, 
						     1, 
						     Rcond, 
						     iterDiagonalOpts, 
						     iterRowindex,
						     iterColindex,
						     iterRangemap,
						     iterDomainmap,
						     distribute,
						     filename,
						     error, 
						     residual, 
						     NumTests ) ;
		      }
		    if ( verbose ) {
		      cout << " Testing  " << filename 
			   << (transpose?" transpose":"" ) 
			   << (distribute?" distribute":"" ) << " error = " 
			   << error 
			   << " residual = " 
			   << residual 
			   << endl ; 
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
                      if ( true ) { cerr << "AMESOS_PRINT " << # variable << "= " << variable << endl; };  }\
                   }


#define TEST_PRINT(variable) { { \
                      if ( true ) { cerr << "AMESOS_PRINT " << # variable  << "= " << variable <<  ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; };  }\
                   }

#endif

//
//  Usage:  TestOptions [-s] [-v]
//

int NextMain( int argc, char *argv[] ) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif


#ifdef HAVE_AMESOS_SCALAPACK
  AmesosClasses.push_back( "Amesos_Scalapack" ) ;
#endif

#if 1
#ifdef HAVE_AMESOS_KLU
  AmesosClasses.push_back( "Amesos_Klu" );
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  AmesosClasses.push_back( "Amesos_Superludist" );
#endif

#ifdef HAVE_AMESOS_PARDISO
  AmesosClasses.push_back( "Amesos_Pardiso" );
#endif

#ifdef HAVE_AMESOS_TAUCS
  AmesosClasses.push_back( "Amesos_Taucs" );
#endif

#ifdef HAVE_AMESOS_LAPACK
  AmesosClasses.push_back( "Amesos_Lapack" );
#endif

#ifdef HAVE_AMESOS_UMFPACK
  AmesosClasses.push_back( "Amesos_Umfpack" );
#endif
#ifdef HAVE_AMESOS_MUMPS
  AmesosClasses.push_back( "Amesos_Mumps" );
#endif

#ifdef HAVE_AMESOS_SUPERLU
  AmesosClasses.push_back( "Amesos_Superlu" );
#endif
#ifdef HAVE_AMESOS_DSCPACK
  //  This fails on my Fedora Core linux box
  //  AmesosClasses.push_back( "Amesos_Dscpack" );         //  bug #1205 
#endif
#endif

  NumAmesosClasses = AmesosClasses.size();

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
    cerr << "Usage TestOptions [-s] [-v] [-q] " << endl ; 
    cerr << "-v:  verbose  " << endl ; 
    cerr << "-s:  small  " << endl ; 
    cerr << "-q:  quiet  " << endl ; 
    exit(-1);
  }

  vector<bool> AmesosClassesInstalled( NumAmesosClasses );



  if ( Comm.MyPID() != 0 ) verbose = false ; 
  if ( Comm.MyPID() == 0 ) {
    char what; 
    //    cin >> what ; 
  }



  Teuchos::ParameterList ParamList ;
    ParamList.set( "DebugLevel", 1 );
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  Comm.SetTracebackMode(2);

  //  for (int i=0; i < 0; i++ ) {
  for (int i=0; i < NumAmesosClasses; i++ ) {

    Abase = Afactory.Create( &AmesosClasses[i][0], Problem ) ; 
    if ( Abase == 0 ) {
      if ( verbose ) cout << AmesosClasses[i] << " not built in this configuration"  << endl ;
      AmesosClassesInstalled[i] = false;
    } else {
      if ( verbose ) cout << " Testing " << AmesosClasses[i] << endl ;
      AmesosClassesInstalled[i] = true;
      Teuchos::ParameterList ParamList ;
      ParamList.set( "NoDestroy", true );    // Only affects Amesos_Mumps
      Abase->SetParameters( ParamList );
    }
    delete Abase ; 
    }

  int result = 0 ; 
  int numtests = 0 ;

  bool symmetric ; 

  result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/Diagonal.mtx", Comm, verbose, false, 1e-6 , numtests ) ;
  symmetric = false ; 
  result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/MissingADiagonal.mtx", Comm, verbose, false, 1e-6 , numtests ) ;
#if 0
  //
  //  TriDiagonal.mtx remains non-singular even after a diagaonal element is removed from the map 
  //
  result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/TriDiagonal.mtx", Comm, verbose, true, 1e-6 , numtests ) ;
  //  symmetric = false; // Bizarre bug causes the following assert to fail even though symmetric, when printed, is false.  When the print statement is in, this line is not necessary and the assert does NOT fail.  
  //  cout << __FILE__ << "::" << __LINE__ << " symmetric = " << symmetric << endl ; 
  
  
#if 0
      // Khead.triS fails on DSCPACK 
      result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/Khead.triS", Comm, verbose, true, 1e-6 , numtests ) ;
#endif

      //
      //  small is set by TestValgrind - keep testing to a minimum because execution time is so slow
      //  quiet is set by TestQuietAmesos - dscpack is not quiet at the moment, hence we can't test symmetric matrices
      //  in TestQuietAmesos
      //
      //  bug #1205 
      //
      if ( ! small && ! quiet ) {
	result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/bcsstk04.mtx", Comm, verbose, false, 1e-6 , numtests ) ;
	result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/662_bus_out.rsa", Comm, verbose, false, 1e-6 , numtests ) ;
	result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/SuperLU.rua", Comm, verbose, false, 1e-6 , numtests ) ;
	result += TestOneMatrix( AmesosClassesInstalled, "../Test_Basic/ImpcolB.rua", Comm, verbose, false, 1e-6 , numtests ) ;
      }

      //      result += TestOneMatrix( AmesosClassesInstalled, (char *) "../Test_Basic/SuperLU.triU", Comm, verbose, 1e-6 , numtests ) ;


#endif
  if ( verbose) cout << result << " Tests failed " ; 

  if (verbose ) cout << numtests << " Tests performed " << endl ; 


  if ( result == 0 && numtests > 0 ) {
    if (verbose && Comm.MyPID() == 0)
      cout << endl << "TEST PASSED" << endl << endl;
  }
  else {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST FAILED" << endl << endl;
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
