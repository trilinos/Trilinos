// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test is for the LOBPCG solver manager on a standard (Ax=xl) complex Hermitian
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
#include "AnasaziLOBPCGSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using std::vector;

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;

  typedef std::complex<double>                ST;
  typedef ScalarTraits<ST>                   SCT;
  typedef SCT::magnitudeType                  MT;
  typedef MultiVector<ST>                     MV;
  typedef MV::global_ordinal_type             GO;
  typedef Operator<ST>                        OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  const ST ONE  = SCT::one();

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  int info = 0;
  int MyPID = 0;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  MyPID = rank(*comm);

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  std::string filename("mhd1280b.cua");
  std::string which("LM");
  int nev = 4;
  int blockSize = 4;
  MT tol = 1.0e-6;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;

  if (MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  // Get the data from the HB file
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
  RCP<const Map<> > map = rcp (new Map<> (dim, 0, comm));
  RCP<CrsMatrix<ST> > K = rcp (new CrsMatrix<ST> (map, rnnzmax));
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
  int maxIters = 450;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "Full Ortho", true );
  //
  // Create the solver manager
  Anasazi::LOBPCGSolMgr<ST,MV,OP> MySolverMgr(problem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
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

    os << "Direct residual norms computed in Tpetra_LOBPCG_complex_test.exe" << endl
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
