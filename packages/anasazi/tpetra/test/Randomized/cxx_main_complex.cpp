// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziRandomizedSolMgr.hpp"

#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Trilinos_Util_iohb.h>

template <class ST>
void print(std::vector<ST> const &a) {
  std::cout << "The vector elements are : ";

  for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
}

int main(int argc, char *argv[])
{
  using namespace Teuchos;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Tpetra::Map;
  using std::vector;
  using std::cout;
  using std::endl;

  typedef std::complex<double>                ST;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Tpetra::MultiVector<ST>             MV;
  typedef MV::global_ordinal_type             GO;
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  const ST ONE  = SCT::one();

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  int info = 0;
  const int MyPID = comm->getRank ();

  bool testFailed;
  bool verbose = true;
  bool debug = false;
  std::string ortho("ICGS");
  std::string filename("mhd1280b.mtx");
  int nev = 4;
  int nsteps = 50;
  int blockSize = 5;
  MT tol = 1.0e-1;
  int resFreq = 0;
  int orthoFreq = 1;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("ortho",&ortho,"Orthogonalization method (DGKS, ICGS, or SVQB)");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("nsteps",&nsteps,"Number of times to apply A. ('q' parameter)");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  cmdp.setOption("filename",&filename,"Filename for Matrix Market test matrix.");
  cmdp.setOption("orthoFreq",&orthoFreq,"How many iterations should pass before reorthogonalizing the basis (not including when computing residual norms).");
  cmdp.setOption("resFreq",&resFreq,"How many iterations should pass before computing the residuals (Orthogonalization of the basis is done here).");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if ( blockSize < nev ){
    if(MyPID == 0) std::cout << "Block size must be greater than or equal to num evals. Increasing block size to " << nev << std::endl;
    blockSize = nev;
  }


  RCP<CrsMatrix<ST>> A;
  RCP<const Tpetra::Map<> > map;
  std::string mtx = ".mtx", mm = ".mm", cua = ".cua";

  if(std::equal(mtx.rbegin(), mtx.rend(), filename.rbegin()) || std::equal(mm.rbegin(), mm.rend(), filename.rbegin())) {
    A = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(filename,comm);
    map = A->getDomainMap();
  } else if(std::equal(cua.rbegin(), cua.rend(), filename.rbegin())) {
    int dim,dim2,nnz;
    int rnnzmax;
    double *dvals;
    int *colptr,*rowind;
    nnz = -1;
    if (MyPID == 0) {
      info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
      // find maximum NNZ over all rows
      vector<int> rnnz(dim,0);
      for (int *ri=rowind; ri<rowind+nnz; ++ri) {
        ++rnnz[*ri-1];
      }
      rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
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
    Teuchos::broadcast(*comm,0,&rnnzmax);
    if (info == 0 || nnz < 0) {
      if (MyPID == 0) {
        cout << "Error reading '" << filename << "'" << endl
          << "End Result: TEST FAILED" << endl;
      }
      return -1;
    }
    // create map
    map = rcp (new Map<> (dim, 0, comm));
    A = rcp (new CrsMatrix<ST> (map, rnnzmax));
    if (MyPID == 0) {
      // Convert interleaved doubles to complex values
      // HB format is compressed column. CrsMatrix is compressed row.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          A->insertGlobalValues (static_cast<GO> (*rptr++ - 1), tuple<GO> (c), tuple (ST (dptr[0], dptr[1])));
          dptr += 2;
        }
      }
    }
    if (MyPID == 0) {
      // Clean up.
      free( dvals );
      free( colptr );
      free( rowind );
    }
    A->fillComplete();

  } else {
    if (MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  // Create initial vectors
  RCP<MV> randVecs = rcp(new MV(map,blockSize));
  randVecs->randomize();

  // Create eigenproblem
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    rcp (new Anasazi::BasicEigenproblem<ST,MV,OP> (A, randVecs));
  problem->setNEV (nev);
  // Inform the eigenproblem that you are done passing it information
  bool boolret = problem->setProblem ();
  if (! boolret) {
    if (MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
  if (verbose) {
    verbosity += Anasazi::IterationDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }

  // Eigensolver parameters
  //
  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", nsteps );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Orthogonalization", ortho ); 
  MyPL.set( "Orthogonalization Frequency", orthoFreq ); 
  MyPL.set( "Residual Frequency", orthoFreq ); 
  //
  // Create the solver manager
  Anasazi::Experimental::RandomizedSolMgr<ST,MV,OP> MySolverMgr(problem, MyPL);
  std::cout << "DEBUG: Created solver manager. Calling solve." << std::endl;
  //int numItsSlvr = MySolverMgr.getNumIters();
  //std::cout << "DEBUG: Num iters is : " << numItsSlvr << std::endl;

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    testFailed = true;
  }

  // Extract evects/evals from solution
  Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  std::cout << "Verifying Eigenvector residuals" << std::endl;
  // Verify Evec residuals by hand. 
  if (numev > 0) {
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MT> normV( numev );
    Teuchos::SerialDenseMatrix<int,ST> T (numev, numev);
    //std::cout << "About to try to access the evals. Bet it crashes here. Haha." << std::endl;
    for (int i = 0; i < numev; ++i) {
      T(i,i) = sol.Evals[i].realpart;
    }
    //std::cout << "If you see this, it didn't really crash there." << std::endl;
    RCP<MV> Avecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *A, *evecs, *Avecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
    MVT::MvNorm( *Avecs, normV );

    os << "Direct residual norms computed in Tpetra_RandomizedSolver_complex_test.exe" << endl
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
    if (MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

  if (testFailed) {
    if (MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}
