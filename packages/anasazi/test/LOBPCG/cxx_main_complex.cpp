// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
//
// This test is for LOBPCG solving a standard (Ax=xl) complex Hermitian
// eigenvalue problem, using the LOBPCGSolMgr solver manager.
//
// The matrix used is from MatrixMarket:
// Name: MHD1280B: Alfven Spectra in Magnetohydrodynamics
// Source: Source: A. Booten, M.N. Kooper, H.A. van der Vorst, S. Poedts and J.P. Goedbloed University of Utrecht, the Netherlands
// Discipline: Plasma physics
// URL: http://math.nist.gov/MatrixMarket/data/NEP/mhd/mhd1280b.html
// Size: 1280 x 1280
// NNZ: 22778 entries

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"

#include "AnasaziMVOPTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "iohb.h"
#endif

// templated multivector and sparse matrix classes
#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
  int info = 0;
  bool boolret;
  int MyPID = 0;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &MyPID );
#endif

  bool testFailed;
  bool verbose = 0;
  std::string filename("mhd1280b.cua");
  std::string which("LM");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ScalarType;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ScalarType;
#else
  typedef double ScalarType;
  // no complex. quit with failure.
  if (verbose && MyPID == 0) {
    cout << "Not compiled with complex support." << endl;
    cout << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
#endif
  typedef ScalarTraits<ScalarType>                   SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Anasazi::MultiVec<ScalarType>               MV;
  typedef Anasazi::Operator<ScalarType>               OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
  ScalarType ONE  = SCT::one();
  ScalarType ZERO = SCT::zero();



  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }


#ifndef HAVE_ANASAZI_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  if (verbose && MyPID == 0) {
    cout << "End Result: TEST FAILED" << endl;	
  }
  return -1;
#endif

  Anasazi::ReturnType returnCode;

  // Get the data from the HB file
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  ScalarType *cvals;
  nnz = -1;
  info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
                              &colptr,&rowind,&dvals);
  if (info == 0 || nnz < 0) {
    if (verbose && MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  // Convert interleaved doubles to complex values
  cvals = new ScalarType[nnz];
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ScalarType(dvals[ii*2],dvals[ii*2+1]);
  }
  // Build the problem matrix
  RefCountPtr< MyBetterOperator<ScalarType> > A 
    = rcp( new MyBetterOperator<ScalarType>(dim,colptr,nnz,rowind,cvals) );

  // Create initial vectors
  int blockSize = 5;
  RefCountPtr<MyMultiVec<ScalarType> > ivec = rcp( new MyMultiVec<ScalarType>(dim,blockSize) );
  ivec->MvRandom();

  // Create eigenproblem
  int nev = 4;
  RefCountPtr<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > MyProblem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(A,ivec) );
  //
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian(true);
  //
  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;	
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }


  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }


  // Eigensolver parameters
  int maxIters = 450;
  MagnitudeType tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "Which", which );
  MyPL.set( "Full Ortho", true );
  MyPL.set( "Verbosity", verbosity );
  //
  // Create the solver manager
  Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    testFailed = true; 
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<MagnitudeType> evals = sol.Evals;
  RefCountPtr<MV> evecs = sol.Evecs;
  int nevecs = MVT::GetNumberVecs(*evecs);

  // Compute the direct residual
  std::vector<MagnitudeType> normV( nevecs );
  SerialDenseMatrix<int,ScalarType> T(nevecs,nevecs);
  for (int i=0; i<nevecs; i++) {
    T(i,i) = evals[i];
  }
  RefCountPtr<MV > Avecs = MVT::Clone( *evecs, nevecs );
  OPT::Apply( *A, *evecs, *Avecs );
  MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
  MVT::MvNorm( *Avecs, &normV );

  for (int i=0; i<nevecs; i++) {
    if ( SCT::magnitude(normV[i]/evals[i]) > 5.0e-5 ) {
      testFailed = true;
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}	
