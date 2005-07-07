#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "TestAllClasses.h"
#include "TestOtherClasses.h"
#include "TestSuperludist.h"
#include "TestScalapack.h"
#include "TestKlu.h"
 
int TestAllClasses( const vector<string> AmesosClasses,
		    int EpetraMatrixType,
		    const vector<bool> AmesosClassesInstalled,
		    Epetra_CrsMatrix *& Amat, 
		    const bool transpose, 
		    const bool verbose, 
		    const bool symmetric, 
		    const int Levels,
		    const double Rcond,
		    int Diagonal,
		    int ReindexRowMap,
		    int ReindexColMap,
		    int RangeMapType,
		    int DomainMapType,
		    bool distribute,
		    char *filename,
		    double &maxrelerror, 
		    double &maxrelresidual,
		    int &NumTests ) {

  bool RowMapEqualsColMap = ( ReindexColMap == 0 ) ; 

  string StringFilename = filename ; 
  bool MissingADiagonal = ( StringFilename == "../Test_Basic/MissingADiagonal.mtx" );
  if ( MissingADiagonal ) RowMapEqualsColMap = false ; // Bug #1405 - this turns off the AddToDiag test 

  const Epetra_Map& row_map = Amat->RowMap() ; 

  const int NumAmesosClasses = AmesosClasses.size();
  int errors = 0 ;

  for (int i=0; i < NumAmesosClasses; i++ ) {
    if ( AmesosClassesInstalled[i] ) { 
      if ( AmesosClasses[i] == "Amesos_Scalapack") { 
	bool RunScalapackTest = true;
	if ( ReindexRowMap || ReindexColMap ) RunScalapackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunScalapackTest = false ;   //  Bug #1403
	if ( RunScalapackTest && verbose) cout << " Testing SCALAPACK " << endl ; 
	if ( RunScalapackTest )	errors += TestScalapack( Amat, 
							 EpetraMatrixType,
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 maxrelerror, 
							 maxrelresidual, 
							 NumTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Umfpack" ) {
	bool RunUmfpackTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) && row_map.DistributedGlobal() ) 
	  RunUmfpackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunUmfpackTest = false ;   //  Bug #1403

	if ( RunUmfpackTest && verbose) cout << " Testing UMFPACK " << endl ; 
	
	if ( RunUmfpackTest ) errors += TestOtherClasses("Amesos_Umfpack",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Taucs" ) {
	bool RunTaucsTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) ) 
	  RunTaucsTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunTaucsTest = false ;   //  Bug #1403
	if ( MissingADiagonal ) RunTaucsTest = false ; // Bug #1449

	if ( RunTaucsTest && verbose) cout << " Testing TAUCS " << endl ; 
	
	if ( RunTaucsTest ) errors += TestOtherClasses("Amesos_Taucs",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Pardiso" ) {
	bool RunPardisoTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) && row_map.DistributedGlobal() ) 
	  RunPardisoTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunPardisoTest = false ;   //  Bug #1403

	if ( RunPardisoTest && verbose) cout << " Testing PARDISO " << endl ; 
	
	if ( RunPardisoTest ) errors += TestOtherClasses("Amesos_Pardiso",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Mumps" ) {
	bool RunMumpsTest = true;
	if ( ( ReindexRowMap || ReindexColMap ) ) 
	  RunMumpsTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunMumpsTest = false ;   //  Bug #1403
	if (  RunMumpsTest && verbose) cout << " Testing MUMPS " << endl ; 
	if ( MissingADiagonal ) RunMumpsTest = false ; // Bug #1435

	if ( RunMumpsTest ) errors += TestOtherClasses("Amesos_Mumps",
							 EpetraMatrixType,
						       Amat, 
						       transpose, 
						       verbose, 
						       Levels, 
						       Rcond, 
						       RowMapEqualsColMap,
						       maxrelerror, 
						       maxrelresidual, 
						       NumTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Klu" ) {
	bool RunKluTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && distribute )  //  Bug #969
	  RunKluTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunKluTest = false ;   //  Bug #1403
	if ( RunKluTest && verbose) cout << " Testing KLU " << endl ; 

	if ( RunKluTest ) errors += TestKlu( Amat, 
							 EpetraMatrixType,
					     transpose, 
					     verbose, 
					     Levels, 
					     Rcond, 
					     RowMapEqualsColMap, 
					     maxrelerror, 
					     maxrelresidual, 
					     NumTests ) ;
  
      } else if ( AmesosClasses[i] == "Amesos_Superlu" ) {
	bool RunSuperluTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && Amat->Comm().NumProc() > 1  )  //  Bug #969
	  RunSuperluTest = false ;   //  Bug #969
	if ( MissingADiagonal ) RunSuperluTest = false ; // Bug #1404
	RowMapEqualsColMap = false ; // Bug #1405 - this turns off the AddToDiag test 
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunSuperluTest = false ;   //  Bug #1403

	if ( RunSuperluTest ) {
	  if ( verbose) cout << " Testing SUPERLU " << endl ; 
	  errors += TestOtherClasses("Amesos_Superlu",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTests ) ;
	}
      } else if ( AmesosClasses[i] == "Amesos_Pastix" ) {
	bool RunPastixTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && Amat->Comm().NumProc() > 1  )  //  Bug #969
	  RunPastixTest = false ;   //  Bug #969
	if ( MissingADiagonal ) RunPastixTest = false ; // Bug #1404
	RowMapEqualsColMap = false ; // Bug #1405 - this turns off the AddToDiag test 
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunPastixTest = false ;   //  Bug #1403

	if ( RunPastixTest ) {
	  if ( verbose) cout << " Testing Pastix " << endl ; 
	  errors += TestOtherClasses("Amesos_Pastix",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTests ) ;
	}

      } else if ( AmesosClasses[i] == "Amesos_Paraklete" ) {
	bool RunParakleteTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && Amat->Comm().NumProc() > 1  )  //  Bug #969
	  RunParakleteTest = false ;   //  Bug #969
	if ( MissingADiagonal ) RunParakleteTest = false ; // Bug #1404
	RowMapEqualsColMap = false ; // Bug #1405 - this turns off the AddToDiag test 
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunParakleteTest = false ;   //  Bug #1403

	if ( RunParakleteTest ) {
	  if ( verbose) cout << " Testing Paraklete " << endl ; 
	  errors += TestOtherClasses("Amesos_Paraklete",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTests ) ;
	}

      } else if ( AmesosClasses[i] == "Amesos_Dscpack" ) {
	//
	//  A quick sanity check - make sure symmetric is the same on all processes
	//
	const int sym_int = symmetric?0:1 ; 
	int sym_int_out = sym_int; 
	Amat->Comm().Broadcast( &sym_int_out, 1, 0 ) ; 
	assert( sym_int == sym_int_out ) ; 

	bool RunDscpackTest = true;
	if ( symmetric ) RunDscpackTest = false ; 
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && distribute )  //  Bug #969
	  RunDscpackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunDscpackTest = false ;   //  Bug #1403
	if ( RunDscpackTest ) { 
	  if ( verbose) cout << " Testing DSCPACK " << endl ; 
    
	  errors += TestOtherClasses("Amesos_Dscpack",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTests ) ;
	} else {
	  if ( verbose ) cout << " DSCPACK not tested on unsymmetric matrices " 
			      << endl ; 
	}
    
      } else if ( AmesosClasses[i] == "Amesos_Superludist" ) {
	bool RunSuperludistTest = true;
	if ( transpose ) { 
	  if ( verbose ) cout << "Superludist does not support transpose " << endl ; 
	  RunSuperludistTest = false ;    // Bug #822
	}
	if ( ReindexRowMap || ReindexColMap ) RunSuperludistTest = false ;    //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunSuperludistTest = false ;   //  Bug #1403
	if ( RunSuperludistTest ) { 
	  if ( verbose) cout << " Testing Superludist " << endl ; 
  
	  errors += TestSuperludist(Amat, 
							 EpetraMatrixType,
				    transpose, 
				    verbose, 
				    Levels, 
				    Rcond, 
				    maxrelerror, 
				    maxrelresidual, 
				    NumTests ) ;
	}
      }
    }
  }
	  
  if ( verbose) cout << " TestAllClasses errors = " << errors << endl ; 

  return errors;
}

