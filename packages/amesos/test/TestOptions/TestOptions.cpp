//
//  TestOptions tests all options for each Amesos Class on a limited number 
//  of matrices.  
//

#include "Trilinos_Util.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Amesos_Factory.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "CrsMatrixTranspose.h"
#include <string>
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int CreateCrsMatrix( char *filename, Epetra_Comm &Comm, 
		     bool transpose, bool distribute, Epetra_CrsMatrix *& Matrix ) {

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
  double idiotic;
   
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
    } else {
      if (  LastFourBytes == ".mtx" ) { 
	EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( filename, Comm, readMap, 
							       readA, readx, readb, readxexact) );
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

  assert( serialA->RowMap() == readMap ) ; 

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
    assert( Amat->RowMap() == map ) ; 
    delete readMap; 
    delete serialA; 
  } else { 

    Matrix = serialA; 
  }



  
}

int TestOneMatrix( char *filename, Epetra_Comm &Comm, bool verbose, double MaxRelResidual ) {

  bool distribute;
  int errcode ;
  double errors[NumAmesosClasses];
  double residuals[NumAmesosClasses];
  for (int i = 0 ; i < NumAmesosClasses; i ++ ) errors[i] = residuals[i] = 0.0 ; 

  assert( (int) AMESOS_KLU == 0 ) ;
  assert( (int) AMESOS_DSCPACK == 3 ) ;

  Epetra_CrsMatrix *Amat ;
  for ( int iterTrans =0 ; iterTrans < 2; iterTrans++ ) {
    bool transpose = iterTrans == 1 ; 
    
    for ( int iterDist =0 ; iterDist < 2; iterDist++ ) {
      bool distribute = Dist == 1 ; 

      CreateCrsMatrix( filename, Comm, transpose, distribute, Amat ) ;

      TestSuperludist( Amat, transpose, error, residual ) ; 


      if ( verbose ) {
	cout << " Amesos_Superludist " << filename << (transpose?" transpose":"" ) 
	     << (distribute?" distribute":"" ) << " error = " << error << " residual = " 
	     << residual << endl ; 
      }
      double relresidual = 
      errors[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( errors[ (int) AMESOS_SUPERLUDIST], error ) ; 
      residuals[(int) AMESOS_SUPERLUDIST] = EPETRA_MAX( residuals[ (int) AMESOS_SUPERLUDIST], residual ) ; 
      NumFailures += ( residual > maxresidual ) ; 

    }
  }

    
  delete Amat->RowMap() ; 
  delete Amat ; 
  
} 


int main( int argc, char *argv[] ) {

  bool verbose = false; 
  if ( argc >= 2 && argv[1][0] == '-' &&  argv[1][0]t == 'v' ) verbose = true ; 

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if ( Comm.MyPID() != 0 ) verbose = false ; 

  AmesosClassType FactorySet[] = { AMESOS_KLU,   AMESOS_UMFPACK,
				  AMESOS_SUPERLUDIST,
				  AMESOS_DSCPACK }; 
  char *AmesosClassNames[] =  { "AMESOS_KLU",   "AMESOS_UMFPACK",
				  "AMESOS_SUPERLUDIST",
				  "AMESOS_DSCPACK" }; 

  AMESOS::Parameter::List ParamList ;
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos_Factory Afactory;

  assert(  sizeof(FactorySet)/sizeof(FactorySet[0]) == 
	   sizeof(AmesosClassNames)/sizeof(AmesosClassNames[0]) );

  for (int i=0; i < sizeof(FactorySet)/sizeof(FactorySet[0]); i++ ) {
    Abase = Afactory.Create( FactorySet[i], Problem, ParamList ) ; 
    if ( Abase == 0 ) {
      cout << AmesosClassNames[i] << " not built in this configuration"  << endl ;
    } else {
      cout << " Testing " << AmesosClassNames[i] << endl ;
    }
  }




  int result ; 
  result += TestOneMatrix("../bcsstk01.mtx", Comm ) ;
}
