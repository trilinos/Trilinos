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
// This test is for BlockKrylovSchur solving a standard (Ax=xl) complex Hermitian
// eigenvalue problem.
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

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

using namespace Teuchos;
using Tpetra::Platform;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;

int main(int argc, char *argv[]) 
{
  using std::cout;
  using std::endl;

  GlobalMPISession mpisess(&argc,&argv,&std::cout);

  int info = 0;
  bool boolret;
  int MyPID = 0;

  RCP<const Platform<int> > platform = Tpetra::DefaultPlatform<int>::getPlatform();
  RCP<Comm<int> > comm = platform->createComm();

  MyPID = rank(*comm);

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  bool shortrun = false;
  bool insitu = false;
  std::string which("LM");
  std::string filename("mhd1280b.cua");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;

  typedef std::complex<double>                ST;
  typedef ScalarTraits<ST>                   SCT;
  typedef SCT::magnitudeType                  MT;
  typedef MultiVector<int,ST>         MV;
  typedef Operator<int,ST>            OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  ST ONE  = SCT::one();

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }
  
  // Get the data from the HB file
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  if (MyPID == 0) {
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    colptr = NULL;
    rowind = NULL;
  }
  Teuchos::broadcast(*comm,0,&info);
  Teuchos::broadcast(*comm,0,&nnz);
  Teuchos::broadcast(*comm,0,&dim);
  if (info == 0 || nnz < 0) {
    if (verbose && MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  // create map
  Map<int> map(dim,0,*platform);
  RCP<CrsMatrix<int,ST> > K = rcp(new CrsMatrix<int,ST>(map));
  if (MyPID == 0) {
    // Convert interleaved doubles to complex values
    // HB format is compressed column. CrsMatrix is compressed row.
    const double *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        K->submitEntry(*rptr++ - 1,c,ST(dptr[0],dptr[1]));
        dptr += 2;
      }
    }
  }
  K->fillComplete();
  cout << *K << endl;

  // Create initial vectors
  int blockSize = 5;
  RCP<MultiVector<int,ST> > ivec = rcp( new MultiVector<int,ST>(map,blockSize) );
  ivec->random();

  // Create eigenproblem
  const int nev = 4;
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,ivec) );
  //
  // Inform the eigenproblem that the operator K is symmetric
  problem->setHermitian(false);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  boolret = problem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }

  // Eigensolver parameters
  int numBlocks;
  int maxRestarts;
  if (shortrun) {
    maxRestarts = 25;
    numBlocks = 5;
  }
  else {
    maxRestarts = 50;
    numBlocks = 10;
  }
  MT tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "In Situ Restarting", insitu );
  //
  // Create the solver manager
  Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> MySolverMgr(problem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && shortrun==false) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MT> normV( numev );
    SerialDenseMatrix<int,ST> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = sol.Evals[i].realpart;
    }
    RCP<MV> Kvecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *K, *evecs, *Kvecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
    MVT::MvNorm( *Kvecs, normV );
  
    os << "Direct residual norms computed in BlockKrylovSchurComplex_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);
      }
      os << std::setw(20) << sol.Evals[i].realpart << std::setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

  if (MyPID == 0) {
    // Clean up.
    free( dvals );
    free( colptr );
    free( rowind );
  }

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
