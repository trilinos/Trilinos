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
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#ifdef HAVE_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#else
#ifdef HAVE_VALGRIND_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#endif 
#endif 

string AS = "Amesos_Scalapack" ;
string AU = "Amesos_Umfpack" ;

vector<string> FactorySet( 7, ("Amesos_Scalapack",
				     "Amesos_Umfpack",
				     "Amesos_Mumps",
				     "Amesos_Superludist",
				     "Amesos_Klu",
				     "Amesos_Superlu",
				     "Amesos_Dscpack")   );
#if 0
  char* FactorySet[] = { "Amesos_Scalapack",
			 "Amesos_Umfpack",
			 "Amesos_Mumps",
			 "Amesos_Superludist",
			 "Amesos_Klu",
			 "Amesos_Superlu",
			 "Amesos_Dscpack" };
#endif


const int NumAmesosClasses = FactorySet.size();

int CreateCrsMatrix( char *filename, Epetra_Comm &Comm, 
		     bool transpose, bool distribute, 
		     bool& symmetric, Epetra_CrsMatrix *& Matrix ) {

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;

  assert( NumAmesosClasses == 7 ) ;
   
  symmetric = false ; 
  string FileName = filename ;
  int FN_Size = FileName.size() ; 
  string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
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
      }
    }
  }

  delete readb;
  delete readxexact;

  Epetra_CrsMatrix *serialA ; 

  if ( transpose ) {
    Epetra_CrsMatrix *transposeA = new Epetra_CrsMatrix( Copy, *readMap, 0 );
    assert( CrsMatrixTranspose( readA, transposeA ) == 0 ); 
    serialA = transposeA ; 
    delete readA;
  } else {
    serialA = readA ; 
  }

  assert( serialA->RowMap().SameAs(*readMap) ) ; 

  if ( distribute ) { 
    // Create uniform distributed map
    Epetra_Map *map = new Epetra_Map(readMap->NumGlobalElements(), 0, Comm);

    // Create Exporter to distribute read-in matrix and vectors
    Epetra_Export exporter( *readMap, *map);
    
    Epetra_CrsMatrix *Amat = new Epetra_CrsMatrix( Copy, *map, 0 );
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
    delete serialA; 
  } else { 

    Matrix = serialA; 
  }


  return 0;
}

int TestOneMatrix( 
		  char *filename, 
		  Epetra_Comm &Comm, 
		  bool verbose, 
		  bool symmetric, 
		  double Rcond,
		  int &NumTests  ) {

  if ( verbose ) cout << endl << endl << " Matrix = " << filename << endl ;

  int NumErrors =0 ;
  double error; 
  double residual;
  //  double errors[NumAmesosClasses];
  //  double residuals[NumAmesosClasses];
  //  for (int i = 0 ; i < NumAmesosClasses; i ++ ) errors[i] = residuals[i] = 0.0 ; 

#if COMPUTE_RCOND 
  Epetra_CrsMatrix *Amat ;

  //
  //  Compute the reciprocal condition number using Amesos_UMFPACK via the Amesos interface
  //
  bool symmetric; 
  CreateCrsMatrix( filename, Comm, false, false, symmetric, Amat ) ;
  Teuchos::ParameterList ParamList ;
  Epetra_LinearProblem Problem;
  Amesos Afactory;

  Amesos_BaseSolver* Abase ; 
  Abase = Afactory.Create( AMESOS_UMFPACK, Problem, ParamList ) ; 
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
  if (verbose) cout << " norm(Amat) = " << Amat->NormInf() << endl; 
  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  if (verbose) cout << " norm(Amat) = " << Amat->NormInf() << endl; 

  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
   double Rcond1 = UmfpackOperator->GetRcond();

  if ( Amat->MyGRID( 0 ) )
    Amat->SumIntoMyValues( 0, 1, val, ind ) ; 
  if (verbose) cout << " norm(Amat) = " << Amat->NormInf() << endl; 
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 
   double Rcond2 = UmfpackOperator->GetRcond();

  if (verbose) cout << " Rcond = " << Rcond << endl; 
  if (verbose) cout << " Rcond1 = " << Rcond1 << endl; 
  if (verbose) cout << " Rcond2 = " << Rcond2 << endl; 


#else
  double Rcond1 = Rcond ;
  double Rcond2 = Rcond ;
#endif
  for ( int iterTrans =0 ; iterTrans < 2; iterTrans++ ) {
    //  for ( int iterTrans =0 ; iterTrans < 1; iterTrans++ ) {
    bool transpose = iterTrans == 1 ; 
    
    for ( int iterDist =0 ; iterDist < 2; iterDist++ ) {
      //    for ( int iterDist =0 ; iterDist < 1; iterDist++ ) {
      bool distribute = ( iterDist == 1 ); 

      Epetra_CrsMatrix *Amat ;
      CreateCrsMatrix( filename, Comm, transpose, distribute, symmetric, Amat ) ;


      if ( Rcond*Rcond1*Rcond2 > 1e-16 ) 
	{ 
	  NumErrors += TestAllClasses( Amat, 
					transpose, 
					verbose, 
					symmetric, 
					3, 
					Rcond*Rcond1*Rcond2, 
					error, 
					residual, 
					NumTests ) ;
	}
      else if ( Rcond*Rcond1 > 1e-16 ) 
	{
	  NumErrors += TestAllClasses( Amat, 
					transpose, 
					verbose, 
					symmetric, 
					2, 
					Rcond*Rcond1, 
					error, 
					residual, 
					NumTests ) ;
	}
      else
	{
	  NumErrors += TestAllClasses( Amat, 
					transpose, 
					verbose, 
					symmetric, 
					1, 
					Rcond, 
					error, 
					residual, 
					NumTests ) ;
	}
      if ( verbose ) {
	cout << " Amesos_Superludist " << filename 
	<< (transpose?" transpose":"" ) 
	<< (distribute?" distribute":"" ) << " error = " 
	<< error 
	<< " residual = " 
	<< residual 
	<< endl ; 
      }
      //      BUG:  Memory leak 
      //      delete &(Amat->RowMap()) ; 
      delete Amat ; 
#if 0
      double relresidual = 
      errors[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( errors[ (int) AMESOS_SUPERLUDIST], error ) ; 
      residuals[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( residuals[ (int) AMESOS_SUPERLUDIST], residual ) ; 
      NumErrors += ( residual > maxresidual ) ; 
#endif
    }
  }

  return NumErrors;
} 


int main( int argc, char *argv[] ) {

  bool verbose = false; 
  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 'v') ) 
    verbose = true ; 

  bool Short = false; 
  if ( argc >= 2 && (argv[1][0] == '-') &&  (argv[1][1] == 's') ) 
    Short = true ; 
  if ( argc >= 3 && (argv[2][0] == '-') &&  (argv[2][1] == 's') ) 
    Short = true ; 

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if ( Comm.MyPID() != 0 ) verbose = false ; 

  Teuchos::ParameterList ParamList ;
    ParamList.set( "DebugLevel", 1 );
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos Afactory;

  //  for (int i=0; i < NumAmesosClasses; i++ ) {
  for (int i=0; i < 0; i++ ) {
    Abase = Afactory.Create( &FactorySet[i][0], Problem ) ; 
    if ( Abase == 0 ) {
      if ( verbose ) cout << FactorySet[i] << " not built in this configuration"  << endl ;
    } else {
      if ( verbose ) cout << " Testing " << FactorySet[i] << endl ;
    }
    delete Abase ; 
    }




  int result = 0 ; 
  int numtests = 0 ;

  bool symmetric ; 

  //  result += TestOneMatrix("Tri.triS", Comm, verbose, symmetric, 1e-1 , numtests ) ;
  //  result += TestOneMatrix("Tri2.triS", Comm, verbose, symmetric, 1e-5 , numtests ) ;
  //  result += TestOneMatrix("../bcsstk01.mtx", Comm, verbose, symmetric, 1e-6 , numtests ) ;
  result += TestOneMatrix("../ImpcolB.rua", Comm, verbose, symmetric, 1e-6 , numtests ) ;

  //
  //  This is really slow when run on valgrind, so we don't want to run 
  //  the following larger matrices when we are using valgrind.
  //
  //  This test is not foolproof - it is possible to have valgrind and not valgrind.h.
  //
#ifdef HAVE_VALGRIND 
  if ( ! RUNNING_ON_VALGRIND ) {
#endif

    if ( ! Short) { 
      result += TestOneMatrix("../bcsstk13.mtx", Comm, verbose, symmetric, 1e-6 , numtests ) ;
      //  result += TestOneMatrix("../bcsstk02.mtx", Comm, verbose, symmetric, 1e-6 , numtests ) ;
      result += TestOneMatrix("../bcsstk04.mtx", Comm, verbose, symmetric, 1e-6 , numtests ) ;
      result += TestOneMatrix("../bcsstk08.mtx", Comm, verbose, symmetric, 1e-6 , numtests ) ;

    }
#ifdef HAVE_VALGRIND 
  }
#endif

  if ( verbose) cout << result << " Tests failed " ; 

  if (verbose ) cout << numtests << " Tests performed " << endl ; 


#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

#ifdef HAVE_VALGRIND
  //
  //  If this is being run under valgrind, query valgrind to see if valgrind
  //  detected ay errors.
  //
  //  This does not catch memory leaks.  grep "loss" in the valgrind log files
  //  to look for memory leaks.  
  //
  if ( RUNNING_ON_VALGRIND ) { 
    if (verbose) cout <<  VALGRIND_COUNT_ERRORS << " valgrind errors " << endl; 
    result += VALGRIND_COUNT_ERRORS;
  }
#endif
  return result ; 
}
