
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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

// Trilinos Tutorial
// -----------------
// Compute the lowest eigenvalue and the corresponding eigenvector
//
// Example of usagee:
// [msala:anasazi]> mpirun -np 4 ./ex1.exe -problem_type=laplace_2d -problem_size=100 -map_type=box
// This will create a matrix corresponding to a 2D Laplacian. The grid is decomposed into squares.

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[])
{
    
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize the command line parser
  CommandLineParser CLP(argc,argv);
  
  // initialize an Gallery object
  CrsMatrixGallery Gallery("", Comm);

  // add default values
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "100" ); 

  // initialize the gallery as specified in the command line
  Gallery.Set(CLP);

  // retrive pointers to matrix
  Epetra_CrsMatrix * RowA = Gallery.GetMatrix();

  // set Anasazi default values
  int block = 1;
  int length = 30;
  int nev = 1;
  double tol = 1.0e-14;
  // SM is the smallest eigenvalue, LM the largest one (in magnitude).
  // Please check the Anasazi documentation.
  string which="SM";
  int restarts = 100;
  //int step = 1;
  int step = restarts*length*block;

  // create an Anasazi vector, based on Epetra vector, and fill it with
  // defautl values
  const Epetra_Map * Map = &(RowA->DomainMap());
  Anasazi::PetraVec<double> ivec(*Map, block);
  ivec.MvRandom();

  // Call the ctor that calls the petra ctor for a matrix
  Anasazi::PetraMat<double> Amat(*RowA);	
  Anasazi::Eigenproblem<double> MyProblem(&Amat, &ivec);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, tol, nev, length, block, 
					       which, step, restarts);
	
  // Inform the solver that the problem is not symmetric (change if required)
  MyBlockArnoldi.setSymmetric(false);
  MyBlockArnoldi.setDebugLevel(0);

  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi.solve();
  
  // Obtain results directly
  double * resids = MyBlockArnoldi.getResiduals();
  double * evalr = MyBlockArnoldi.getEvals(); 
  double * evali = MyBlockArnoldi.getiEvals();

  // Retrieve eigenvectors
  Anasazi::PetraVec<double> evecr(*Map, nev);
  MyBlockArnoldi.getEvecs( evecr );
  Anasazi::PetraVec<double> eveci(*Map, nev);
  MyBlockArnoldi.getiEvecs( eveci );

  // Output results to screen
  if( Comm.MyPID() == 0 ) {
    for( int i=0 ; i<nev ; ++i ) 
      cout << "eval[" << i << "] = " << evalr[i] << " + i " << evali[i] << endl;
    MyBlockArnoldi.currentStatus();
  }
  
  // free memory
  if( resids ) delete [] resids;
  if( evalr ) delete [] evalr;
  if( evali ) delete [] evali;
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;

}
