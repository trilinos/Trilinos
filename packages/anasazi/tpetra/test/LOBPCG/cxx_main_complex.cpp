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

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

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
  typedef MultiVector<ST,int>                 MV;
  typedef Operator<ST,int>                    OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  const ST ONE  = SCT::one();

  GlobalMPISession mpisess(&argc,&argv,&std::cout);

  int info = 0;
  int MyPID = 0;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

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
  ST *cvals;
  int *offsets,*colinds;
  nnz = -1;
  if (MyPID == 0) {
    double *dvals;
    int *colptr, *rowind;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // Find number of non-zeros for each row
    vector<size_t> rnnz(dim,0);
    // Move through all row indices, adding the contribution to the appropriate row.
    // Remember, the file uses one-based indexing. We'll convert it to zero-based later.
    int curcol_0 = 0, currow_0;
    for (size_t curnnz_1=1; curnnz_1 <= nnz; ++curnnz_1) {
      // if colptr[c] <= curnnz_1 < colptr[c+1], then curnnz_1 belongs to column c
      // invariant: curcol_0 is the column for curnnz_1, i.e., curcol_0 is smallest number such that colptr[curcol_0+1] > curnnz_1
      while (colptr[curcol_0+1] <= curnnz_1) ++curcol_0;
      // entry curnnz_1 corresponds to (curcol_0, rowind[curnnz_1]) and (rowind[curnnz_1], curcol_0)
      // make sure not to count it twice
      ++rnnz[curcol_0];
    }
    const size_t totalnnz = std::accumulate( rnnz.begin(), rnnz.end(), 0 );
    // mark the maximum nnz per row, used to allocated data for the crsmatrix
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
    // allocate row structure and fill it
    cvals           = new ST[totalnnz];
    int *newoffs    = new int[dim+1];
    int *newinds    = new int[totalnnz];
    const double *curdval = dvals;
    // set up pointers
    newoffs[0] = 0;
    for (size_t row=1; row != dim+1; ++row) {
      newoffs[row] = newoffs[row-1] + rnnz[row-1];
    }
    // reorganize data from column oriented to row oriented, duplicating symmetric part as well
    // rowind are one-based; account for that, and convert them to zero-based.
    // use nnz as the number of entries add per row thus far
    std::fill( rnnz.begin(), rnnz.end(), 0 );
    curcol_0 = 0;
    for (size_t curnnz_1=1; curnnz_1 <= nnz; ++curnnz_1) {
      // if colptr[c] <= curnnz_1 < colptr[c+1], then curnnz_1 belongs to column c
      // invariant: curcol_0 is the column for curnnz_1, i.e., curcol_0 is smallest number such that colptr[curcol_0+1] > curnnz_1
      while (colptr[curcol_0+1] <= curnnz_1) ++curcol_0;
      // entry curnnz_1 corresponds to (rowind[curnnz_1], curcol_0)
      currow_0 = rowind[curnnz_1-1] - 1;
      // add (currow_0, curcol_0)
      size_t newnnz = newoffs[currow_0] + rnnz[currow_0];
      cvals[newnnz] = ST(curdval[0],curdval[1]);
      curdval += 2;
      newinds [newnnz] = curcol_0;
      ++rnnz[currow_0];
    }
    const size_t totalnnz2 = std::accumulate( rnnz.begin(), rnnz.end(), 0 );
    TEST_FOR_EXCEPT( totalnnz2 != totalnnz );
    // free the original data, point to new dada
    delete [] dvals;
    delete [] colptr; colptr = NULL;
    delete [] rowind; rowind = NULL;
    offsets = newoffs;
    colinds = newinds;
  }
  else {
    // address uninitialized data warnings
    cvals = NULL;
    offsets = NULL;
    colinds = NULL;
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
  RCP<const Map<int> > map = rcp( new Map<int>(dim,0,comm) );
  RCP<CrsMatrix<ST,int> > K = rcp( new CrsMatrix<ST,int>(map,rnnzmax) );
  if (MyPID == 0) {
    for (size_t row=0; row < dim; ++row) {
      const size_t nE = offsets[row+1] - offsets[row];
      if (nE > 0) {
        // add row to matrix
        K->insertGlobalValues( row, arrayView<const int>(colinds+offsets[row], nE), 
                                    arrayView<const  ST>(  cvals+offsets[row], nE));
      }
    }
  }
  if (MyPID == 0) {
    // Clean up.
    free( cvals );
    free( colinds );
    free( offsets );
  }
  K->fillComplete(Tpetra::DoOptimizeStorage);

  // Create initial vectors
  RCP<MV> ivec = rcp( new MV(map,blockSize) );
  MVT::MvRandom(*ivec);

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
