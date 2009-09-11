// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero. 
// The problem is solver for multiple scalar types, and timings are reported.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>

// #define PRINT_MATRIX

using namespace Teuchos;
using namespace Belos;
using Tpetra::global_size_t;
using Tpetra::DefaultPlatform;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Vector;
using Tpetra::Map;
using std::endl;
using std::cout;
using std::string;
using std::setw;
using std::vector;
using Teuchos::arrayView;


bool proc_verbose = false, reduce_tol, precond = true, dumpdata = false;
RCP<Map<int> > vmap;
ParameterList mptestpl; 
size_t rnnzmax;
int mptestdim, numrhs; 
double *dvals; 
const int *offsets, *colinds;
int mptestmypid = 0;
int mptestnumimages = 1;
int maxIters = 100;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<LinearProblem<Scalar,MultiVector<Scalar,int>,Operator<Scalar,int> > > buildProblem() {
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar,int>         OP;
  typedef MultiVector<Scalar,int>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;
  RCP<CrsMatrix<Scalar,int> > A = rcp(new CrsMatrix<Scalar,int>(vmap,rnnzmax));
  Array<Scalar> vals(rnnzmax);
  if (mptestmypid == 0) {
    for (size_t row=0; row < mptestdim; ++row) {
      const size_t nE = offsets[row+1] - offsets[row];
      if (nE > 0) {
        // convert row values from double to Scalar
        std::copy( dvals+offsets[row], dvals+offsets[row]+nE, vals.begin() );
        // add row to matrix
        A->insertGlobalValues( row, arrayView<const int>(colinds+offsets[row], nE), vals(0,nE) );
      }
    }
  }
  // distribute matrix data to other nodes
  A->fillComplete(Tpetra::DoOptimizeStorage);
  // Create initial MV and solution MV
  RCP<MV> B, X;
  X = rcp( new MV(vmap,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MV(vmap,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );
  // Construct a linear problem instance with zero initial MV
  RCP<LinearProblem<Scalar,MV,OP> > problem = rcp( new LinearProblem<Scalar,MV,OP>(A,X,B) );
  problem->setLabel(Teuchos::typeName(SCT::one()));
  TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters) {
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar,int>         OP;
  typedef MultiVector<Scalar,int>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;

  const Scalar ONE  = SCT::one();
  mptestpl.set<MT>( "Convergence Tolerance", ltol );         // Relative convergence tolerance requested
  mptestpl.set( "Explicit Residual Test", true );
  mptestpl.set( "Timer Label", typeName(ONE) );         // Set timer label to discern between the two solvers.

  if (mptestmypid==0) cout << "Testing Scalar == " << typeName(ONE) << endl;

  RCP<LinearProblem<Scalar,MV,OP> > problem;
  Time btimer("Build Timer"), ctimer("Construct Timer"), stimer("Solve Timer");
  ReturnType ret;
  if (mptestmypid==0) cout << "Building problem..." << endl;
  { 
    TimeMonitor localtimer(btimer);
    problem = buildProblem<Scalar>();
  }
  RCP<SolverManager<Scalar,MV,OP> > solver;
  if (mptestmypid==0) cout << "Constructing solver..." << endl; 
  {
    TimeMonitor localtimer(ctimer);
    solver = rcp(new BlockGmresSolMgr<Scalar,MV,OP>( problem, rcp(&mptestpl,false) ));
  }
  if (mptestmypid==0) cout << "Solving problem..." << endl;
  {
    TimeMonitor localtimer(stimer);
    try {
      ret = solver->solve();
    }
    catch (std::exception &e) {
      cout << "Caught exception: " << endl << e.what() << endl;
      ret = Unconverged;
    }
  }
  if (solver->isLOADetected()) {
    if (proc_verbose) cout << "Solver reports Loss of Accuracy. Computing residual anyway." << endl;
    // set Converged and check below.
    ret = Converged;
  }
  numIters = solver->getNumIters();
  times[0] = btimer.totalElapsedTime();
  times[1] = ctimer.totalElapsedTime();
  times[2] = stimer.totalElapsedTime();

  bool badRes = false;
  if (ret == Converged) {
    if (proc_verbose) cout << endl;
    if (mptestmypid==0) cout << "Computing residuals..." << endl;
    //
    // Compute actual residuals.
    //
    RCP<const OP> A = problem->getOperator();
    RCP<MV> X = problem->getLHS();
    RCP<const MV> B = problem->getRHS();
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(vmap, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -ONE, resid, ONE, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) cout << "Actual Residuals (normalized)" << endl;
    for ( int i=0; i<numrhs; i++) {
      MT actRes = actual_resids[i]/rhs_norm[i];
      if (proc_verbose) cout << "Problem "<<i<<" : \t"<< actRes <<endl;
      if (actRes > ltol) badRes = true;
    }
    if (proc_verbose) cout << endl;
  }
  if (ret != Converged || badRes) return false;
  return true;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  numrhs = 1;      // total number of right-hand sides to solve for
  int blocksize = 1;   // blocksize used by solver
  int numblocks = 5;   // blocksize used by solver
  string filename("bcsstk17.rsa");
  double tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by the solver.");
  cmdp.setOption("num-blocks",&numblocks,"Number of blocks in the Krylov basis.");
  cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
  cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");
  cmdp.setOption("num-iters",&maxIters,"Number of iterations.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  mptestnumimages = size(*comm);
  mptestmypid = rank(*comm);
  proc_verbose = ( verbose && (mptestmypid==0) );
  if (proc_verbose) cout << Belos_Version() << endl << endl;
  if (mptestmypid != 0) dumpdata = 0;

  //
  // Get the data (double) from the HB file and build the Map,Matrix
  //
  int nnz, info;
  nnz = -1;
  if (mptestmypid == 0) {
    int dim2;
    int *colptr, *rowind;
    info = readHB_newmat_double(filename.c_str(),&mptestdim,&dim2,&nnz,&colptr,&rowind,&dvals);
#ifdef PRINT_MATRIX
    std::cout << "\n\nColumn format: " << std::endl;
    for (int c=0; c < mptestdim; ++c) {
      for (int o=colptr[c]; o != colptr[c+1]; ++o) {
        std::cout << rowind[o-1]-1 << ", " << c << ", " << dvals[o-1] << std::endl;
      }
    }
#endif
    // Find number of non-zeros for each row
    vector<size_t> rnnz(mptestdim,0);
    // Move through all row indices, adding the contribution to the appropriate row
    // Skip the diagonals, we'll catch them below on the column pass.
    // Remember, the file uses one-based indexing. We'll convert it to zero-based later.
    int curcol_0 = 0, currow_0;
    for (size_t curnnz_1=1; curnnz_1 <= nnz; ++curnnz_1) {
      // if colptr[c] <= curnnz_1 < colptr[c+1], then curnnz_1 belongs to column c
      // invariant: curcol_0 is the column for curnnz_1, i.e., curcol_0 is smallest number such that colptr[curcol_0+1] > curnnz_1
      while (colptr[curcol_0+1] <= curnnz_1) ++curcol_0;
      // entry curnnz_1 corresponds to (curcol_0, rowind[curnnz_1]) and (rowind[curnnz_1], curcol_0)
      // make sure not to count it twice
      ++rnnz[curcol_0];
      currow_0 = rowind[curnnz_1-1] - 1;
      if (curcol_0 != currow_0) {
        ++rnnz[currow_0];
      }
    }
    const size_t totalnnz = std::accumulate( rnnz.begin(), rnnz.end(), 0 );
    // mark the maximum nnz per row, used to allocated data for the crsmatrix
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
    // allocate row structure and fill it
    double *newdvals  = new double[totalnnz];
    int *newoffs = new int[mptestdim+1];
    int *newinds    = new int[totalnnz];
    // set up pointers
    newoffs[0] = 0;
    for (size_t row=1; row != mptestdim+1; ++row) {
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
      // entry curnnz_1 corresponds to (curcol_0, rowind[curnnz_1]) and (rowind[curnnz_1], curcol_0)
      // it must be added twice if curcol_0 != rowind[curnnz_1]
      currow_0 = rowind[curnnz_1-1] - 1;
      // add (currow_0, curcol_0)
      const double curval = dvals[curnnz_1-1];
      size_t newnnz = newoffs[currow_0] + rnnz[currow_0];
      newdvals[newnnz] = curval;
      newinds [newnnz] = curcol_0;
      ++rnnz[currow_0];
      if (curcol_0 != currow_0) {
        newnnz = newoffs[curcol_0] + rnnz[curcol_0];
        newdvals[newnnz] = curval;
        newinds [newnnz] = currow_0;
        ++rnnz[curcol_0];
      }
    }
    const size_t totalnnz2 = std::accumulate( rnnz.begin(), rnnz.end(), 0 );
    TEST_FOR_EXCEPT( totalnnz2 != totalnnz );
    // free the original data, point to new dada
    delete [] dvals;
    dvals = newdvals;
    delete [] colptr; colptr = NULL;
    delete [] rowind; rowind = NULL;
    offsets = newoffs;
    colinds = newinds;
#ifdef PRINT_MATRIX
    std::cout << "\n\nRow format: " << std::endl;
    for (int r=0; r < mptestdim; ++r) {
      for (int o=offsets[r]; o != offsets[r+1]; ++o) {
        std::cout << r << ", " << colinds[o] << ", " << dvals[o] << std::endl;
      }
    }
#endif
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    offsets = NULL;
    colinds = NULL;
  }
  broadcast(*comm,0,&info);
  broadcast(*comm,0,&nnz);
  broadcast(*comm,0,&mptestdim);
  broadcast(*comm,0,&rnnzmax);
  if (info == 0 || nnz < 0) {
    if (mptestmypid == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // create map
  vmap = rcp(new Map<int>(static_cast<global_size_t>(mptestdim),static_cast<int>(0),comm));
  //
  mptestpl.set<int>( "Maximum Restarts", maxIters/blocksize + 1);
  mptestpl.set<int>( "Maximum Iterations", maxIters);
  mptestpl.set<int>( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  mptestpl.set<int>( "Num Blocks", numblocks);
  int verbLevel = Errors + Warnings;
  if (debug) {
    verbLevel |= Debug;
  }
  if (verbose) {
    verbLevel |= TimingDetails + FinalSummary + StatusTestDetails;
  }
  mptestpl.set( "Verbosity", verbLevel );
  if (verbose) {
    if (frequency > 0) {
      mptestpl.set( "Output Frequency", frequency );
    }
  }

  //
  // **********Print out information about problem*******************
  //
  if (mptestmypid==0) {
    cout << "Filename: " << filename << endl;
    cout << "Dimension of matrix: " << mptestdim << endl;
    cout << "Number of nonzeros: " << nnz << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Number of blocks used by solver: " << numblocks << endl;
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter);

  // done with the matrix data now; delete it
  if (mptestmypid == 0) {
    delete [] dvals; dvals = NULL;
    delete [] colinds; colinds = NULL;
    delete [] offsets; offsets = NULL;
  }

  if (mptestmypid==0) {
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << setw(12) << "float"   << "     " << setw(10) <<  ftime[0] << "     " << setw(9) <<  ftime[1] << "     " << setw(10) <<  ftime[2] << "     " << setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << setw(12) << "double"  << "     " << setw(10) <<  dtime[0] << "     " << setw(9) <<  dtime[1] << "     " << setw(10) <<  dtime[2] << "     " << setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
  }

  if (dumpdata) {
    std::ofstream fout("data.dat",std::ios_base::app | std::ios_base::out);
    fout << setw(80) << "filename" << setw(14) << "numimages"       << setw(14) << "blocksize" << setw(14) << "numblocks" << setw(14) << "dim" << setw(14) << "nnz" << setw(8) << "Scalar"   << setw(12) <<  "build"  << setw(12) <<  "init"   << setw(12) <<  "solve"  << setw(12) <<  "numiter" << endl;
    fout << setw(80) << filename   << setw(14) <<  mptestnumimages  << setw(14) <<  blocksize  << setw(14) << numblocks   << setw(14) << mptestdim   << setw(14) << nnz   << setw(8) << "float"    << setw(12) <<  ftime[0] << setw(12) <<  ftime[1] << setw(12) <<  ftime[2] << setw(12) <<  fiter     << endl;
    fout << setw(80) << filename   << setw(14) <<  mptestnumimages  << setw(14) <<  blocksize  << setw(14) << numblocks   << setw(14) << mptestdim   << setw(14) << nnz   << setw(8) << "double"   << setw(12) <<  dtime[0] << setw(12) <<  dtime[1] << setw(12) <<  dtime[2] << setw(12) <<  diter     << endl;
    fout.close();
  }

  bool allpass = fpass && dpass;

  if (!allpass) {
    if (mptestmypid==0) cout << "\nEnd Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (mptestmypid==0) cout << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}
