// @HEADER
// ***********************************************************************
// 
//                Amesos: An Interface to Direct Solvers
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

#include "Amesos_ConfigDefs.h"
// This example needs Galeri to generate the linear system.
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

using namespace Teuchos;
using namespace Galeri;

int MyCreateCrsMatrix( char *in_filename, const Epetra_Comm &Comm, 
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
	FILE* in_file = fopen( filename, "r");
	assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 
	const int BUFSIZE = 800 ; 
	char buffer[BUFSIZE] ; 
	fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this string 
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






// ===================== //
// M A I N   D R I V E R //
// ===================== //
//
// This example compares all the available Amesos solvers
// for the solution of the same linear system. 
//
// The example can be run in serial and in parallel.
//
// Author: Marzio Sala, SNL 9214
// Last modified: Oct-05.

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);
  double TotalResidual = 0.0;

  // Create the Map, defined as a grid, of size nx x ny x nz,
  // subdivided into mx x my x mz cubes, each assigned to a 
  // different processor.

#if 0
  ParameterList GaleriList;
  GaleriList.set("nx", 4);
  GaleriList.set("ny", 4);
  GaleriList.set("nz", 4 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1);
  GaleriList.set("mz", Comm.NumProc());
  Epetra_Map* Map = CreateMap("Cartesian3D", Comm, GaleriList);
  
  // Create a matrix, in this case corresponding to a 3D Laplacian
  // discretized using a classical 7-point stencil. Please refer to 
  // the Galeri documentation for an overview of available matrices.
  // 
  // NOTE: matrix must be symmetric if DSCPACK is used.

  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace3D", Map, GaleriList);
#else
  bool transpose = false ; 
  bool distribute = false ; 
  bool symmetric ; 
  Epetra_CrsMatrix *Matrix = 0 ;
  Epetra_Map *Map = 0 ;
  CreateCrsMatrix( "ibm.triU", Comm, Map, transpose, distribute, &symmetric, Matrix ) ;




#endif

  // build vectors, in this case with 1 vector
  Epetra_MultiVector LHS(*Map, 1); 
  Epetra_MultiVector RHS(*Map, 1); 
    
  // create a linear problem object
  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  // use this list to set up parameters, now it is required
  // to use all the available processes (if supported by the
  // underlying solver). Uncomment the following two lines
  // to let Amesos print out some timing and status information.
  ParameterList List;
  List.set("PrintTiming",true);
  List.set("PrintStatus",true);
  List.set("MaxProcs",Comm.NumProc());

  std::vector<std::string> SolverType;
  SolverType.push_back("Amesos_Paraklete");

  Epetra_Time Time(Comm);
  
  // this is the Amesos factory object that will create 
  // a specific Amesos solver.
  Amesos Factory;

  // Cycle over all solvers.
  // Only installed solvers will be tested.
  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) 
  {
    // Check whether the solver is available or not
    if (Factory.Query(SolverType[i])) 
    {
      // 1.- set exact solution (constant vector)
      LHS.PutScalar(1.0);
 
      // 2.- create corresponding rhs
      Matrix->Multiply(false, LHS, RHS);
 
      // 3.- randomize solution vector
      LHS.Random();
 
      // 4.- create the amesos solver object
      Amesos_BaseSolver* Solver = Factory.Create(SolverType[i], Problem);
      assert (Solver != 0);

      Solver->SetParameters(List);

      // 5.- factorize and solve
      
      Time.ResetStartTime();
      AMESOS_CHK_ERR(Solver->SymbolicFactorization());

      // 7.- delete the object
      delete Solver;

    }
  }

  delete Matrix;
  delete Map;

  if (TotalResidual > 1e-9) 
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} // end of main()
