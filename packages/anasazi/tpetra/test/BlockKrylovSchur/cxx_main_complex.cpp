// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>
//I/O for Matrix Market files
#include <MatrixMarket_Tpetra.hpp>

using Tpetra::CrsMatrix;
using Tpetra::Map;
using std::vector;

int main(int argc, char *argv[])
{
#ifndef HAVE_TPETRA_COMPLEX_DOUBLE
#  error "Anasazi: This test requires Scalar = std::complex<double> to be enabled in Tpetra."
#else
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
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

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  bool success = false;

  const ST ONE = SCT::one ();


  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  const int MyPID = comm->getRank ();

  bool verbose = false;
  bool debug = false;
  bool insitu = false;
  bool herm = false;
  std::string which("LM");
  std::string filename;
  int nev = 4;
  int blockSize = 4;
  MT tol = 1.0e-6;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("herm","nonherm",&herm,"Solve Hermitian or non-Hermitian problem.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix (assumes non-Hermitian unless specified otherwise).");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;
  if (filename == "") filename = "mhd1280b.cua";

  if (MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  // Get the data from the HB file
  int dim,dim2,nnz;
  int rnnzmax;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;

  RCP<CrsMatrix<ST>> K;
  RCP<const Tpetra::Map<> > map;
  std::string mtx = ".mtx", mm = ".mm", cua = ".cua";

  if(std::equal(mtx.rbegin(), mtx.rend(), filename.rbegin()) || std::equal(mm.rbegin(), mm.rend(), filename.rbegin())) {
    K = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(filename,comm);
    map = K->getDomainMap();
  } else if(std::equal(cua.rbegin(), cua.rend(), filename.rbegin())) {

    int info = 0;
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
    K = rcp (new CrsMatrix<ST> (map, rnnzmax));
    if (MyPID == 0) {
      // Convert interleaved doubles to complex values
      // HB format is compressed column. CrsMatrix is compressed row.
      const double *dptr = dvals;
      const int *rptr = rowind;
      for (int c=0; c<dim; ++c) {
        for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
          K->insertGlobalValues (static_cast<GO> (*rptr++ - 1), tuple<GO> (c), tuple (ST (dptr[0], dptr[1])));
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
    K->fillComplete();
    // cout << *K << endl;
  } else {
    if (MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
        << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // Create initial vectors
  RCP<MV> ivec = rcp( new MV(map,blockSize) );
  ivec->randomize ();

  // Create eigenproblem
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,ivec) );
  //
  // Inform the eigenproblem that the operator K is symmetric
  problem->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  bool boolret = problem->setProblem();
  if (boolret != true) {
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
  int numBlocks = 8;
  int maxRestarts = 10;
  //
  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
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
  success = (returnCode == Anasazi::Converged);

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
    Teuchos::SerialDenseMatrix<int,ST> T (numev, numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = ST(sol.Evals[i].realpart,sol.Evals[i].imagpart);
    }
    RCP<MV> Kvecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *K, *evecs, *Kvecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
    MVT::MvNorm( *Kvecs, normV );

    os << "Direct residual norms computed in BlockKrylovSchurComplex_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(T(i,i)) != SCT::zero() ) {
        normV[i] = SCT::magnitude(normV[i]/T(i,i));
      }
      os << std::setw(20) << T(i,i) << std::setw(20) << normV[i] << endl;
      success = (normV[i] < tol);
    }
    if (MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

  if (MyPID==0) {
    if (success)
      cout << "End Result: TEST PASSED" << endl;
    else
      cout << "End Result: TEST FAILED" << endl;
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#endif // HAVE_TPETRA_COMPLEX_DOUBLE
}
