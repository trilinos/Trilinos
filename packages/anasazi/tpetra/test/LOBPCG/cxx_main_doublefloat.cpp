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
// This test is for the LOBPCG solver manager on a standard (Ax=xl) symmetric
// eigenvalue problem, with a diagonal preconditioner. The test can be run in either
// float or double.
//
// This driver reads the problem from a Harwell-Boeing (HB) file.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"

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
using namespace Anasazi;
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

// TODO: adjust this test to use CrsMatrix<double> regardless of the scalar type for MultiVector
//       this requires the ability to apply those operators to those MultiVectors, which Tpetra
//       currently cannot do.

bool proc_verbose = false, reduce_tol, precond = true, dumpdata = false;
RCP<Map<int> > vmap;
ParameterList mptestpl; 
size_t rnnzmax;
int mptestdim, mptestnev, blockSize;
double *dvals; 
const int *offsets, *colinds;
int mptestmypid = 0;
string mptestwhich("LR");

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<Eigenproblem<Scalar,MultiVector<Scalar,int>,Operator<Scalar,int> > > buildProblem() {
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
  // Create initial MV
  RCP<MV> X;
  X = rcp( new MV(vmap,blockSize) );
  MVT::MvRandom( *X );
  // Construct a linear problem instance with zero initial MV
  RCP<Eigenproblem<Scalar,MV,OP> > problem = rcp( new BasicEigenproblem<Scalar,MV,OP>(A,X) );
  problem->setHermitian(true);
  problem->setNEV(mptestnev);
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
  mptestpl.set<MT>( "Locking Tolerance", 0.1*ltol );         // Relative convergence tolerance requested

  if (mptestmypid==0) cout << "Testing Scalar == " << typeName(ONE) << endl;

  RCP<Eigenproblem<Scalar,MV,OP> > problem;
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
    solver = rcp(new LOBPCGSolMgr<Scalar,MV,OP>( problem, mptestpl));
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
    Eigensolution<Scalar,MV> solution = problem->getSolution();
    int numsol = solution.numVecs;
    RCP<const Operator<Scalar,int> > A = problem->getOperator();
    RCP<MultiVector<Scalar,int> > X = solution.Evecs;
    RCP<MultiVector<Scalar,int> > R = rcp(new MultiVector<Scalar,int>(*X));
    RCP<MultiVector<Scalar,int> > XL= rcp(new MultiVector<Scalar,int>(*X));
    A->apply(*X,*R);  // R = A*X
    Array<MT> L( numsol );
    for (int i=0; i<numsol; ++i) {
      L[i] = solution.Evals[i].realpart;
    }
    XL->scale(L);
    R->update(-1.0,*XL,1.0);
    Array<MT> resnorms( numsol );
    R->norm2(resnorms);
    for (int i=0; i<numsol; ++i) {
      if (L[i] != 0.0) {
        resnorms[i] /= SCT::magnitude(L[i]);
      }
    }
    if (proc_verbose) cout << setw(16) << "Lambda" << setw(16) << "Residuals" << endl;
    for ( int i=0; i<numsol; ++i) {
      if (proc_verbose) cout << setw(16) << L[i] << setw(16) << resnorms[i] << endl;
      if (resnorms[i] > ltol) badRes = true;
    }
    if (proc_verbose) cout << endl;
  }
  if (ret != Converged || badRes) return false;
  return true;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) 
{
  GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, debug = false;
  mptestnev = 1;             // number of eigenpairs to find
  blockSize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  string filename("bcsstk14.hb");
  double tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("which",&mptestwhich,"Which eigenvalues to compute.");
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("nev",&mptestnev,"Number of eigenpairs to find.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blockSize,"Block size to be used by the CG solver.");
  cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
  cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }

  mptestmypid = rank(*comm);
  proc_verbose = ( verbose && (mptestmypid==0) );
  if (proc_verbose) cout << Anasazi_Version() << endl << endl;

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
  // create the parameter list
  if (maxiters == -1) {
    maxiters = mptestdim/blockSize - 1; // maximum number of iterations to run
  }
  //
  mptestpl.set( "Which", mptestwhich );
  mptestpl.set( "Block Size", blockSize );              // Blocksize to be used by iterative solver
  mptestpl.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  int verbLevel = Errors + Warnings;
  if (debug) {
    verbLevel += Debug;
  }
  if (verbose) {
    verbLevel += TimingDetails + FinalSummary;
  }
  mptestpl.set( "Verbosity", verbLevel );

  //
  // **********Print out information about problem*******************
  //
  if (mptestmypid==0) {
    cout << "Filename: " << filename << endl;
    cout << "Dimension of matrix: " << mptestdim << endl;
    cout << "Number of nonzeros: " << nnz << endl;
    cout << "NEV: " << mptestnev << endl;
    cout << "Block size used by solver: " << blockSize << endl;
    cout << "Max number of iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter);
  /*
#ifdef HAVE_TEUCHOS_QD
  unsigned int old_cw; fpu_fix_start(&old_cw); // fix floating point rounding for intel processors; required for proper implementation
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double ddtime[3]; int dditer; bool ddpass = runTest<dd_real>(ltol,ddtime,dditer);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double qdtime[3]; int qditer; bool qdpass = runTest<qd_real>(ltol,qdtime,qditer);
  fpu_fix_end(&old_cw);                        // restore previous floating point rounding
#endif
  */

  // done with the matrix data now; delete it
  if (mptestmypid == 0) {
    delete [] dvals; dvals = NULL;
    delete [] colinds; colinds = NULL;
    delete [] offsets; offsets = NULL;
  }

  if (mptestmypid==0) {
    cout << endl << endl;
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << setw(12) << "float"   << "     " << setw(10) <<  ftime[0] << "     " << setw(9) <<  ftime[1] << "     " << setw(10) <<  ftime[2] << "     " << setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << setw(12) << "double"  << "     " << setw(10) <<  dtime[0] << "     " << setw(9) <<  dtime[1] << "     " << setw(10) <<  dtime[2] << "     " << setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
    /*
#ifdef HAVE_TEUCHOS_QD
    cout << setw(12) << "dd_real" << "     " << setw(10) << ddtime[0] << "     " << setw(9) << ddtime[1] << "     " << setw(10) << ddtime[2] << "     " << setw(9) << dditer << "     " << (ddpass ? "pass" : "fail") << endl;
    cout << setw(12) << "qd_real" << "     " << setw(10) << qdtime[0] << "     " << setw(9) << qdtime[1] << "     " << setw(10) << qdtime[2] << "     " << setw(9) << qditer << "     " << (qdpass ? "pass" : "fail") << endl;
#endif
    */
  }

  if (dumpdata) {
    std::ofstream fout("data.dat",std::ios_base::app | std::ios_base::out);
    fout << setw(14) << nnz << setw(8) << "float"   << setw(12) <<  ftime[0] << setw(12) <<  ftime[1] << setw(12) <<  ftime[2] << setw(12) <<  fiter << endl;
    fout << setw(14) << nnz << setw(8) << "double"  << setw(12) <<  dtime[0] << setw(12) <<  dtime[1] << setw(12) <<  dtime[2] << setw(12) <<  diter << endl;
    /*
#ifdef HAVE_TEUCHOS_QD
    fout << setw(14) << nnz << setw(8) << "dd_real" << setw(12) << ddtime[0] << setw(12) << ddtime[1] << setw(12) << ddtime[2] << setw(12) << dditer << endl;
    fout << setw(14) << nnz << setw(8) << "qd_real" << setw(12) << qdtime[0] << setw(12) << qdtime[1] << setw(12) << qdtime[2] << setw(12) << qditer << endl;
#endif
    */
    fout.close();
  }

  bool allpass = fpass && dpass;
  /*
#ifdef HAVE_TEUCHOS_QD
  allpass = allpass && ddpass && qdpass;
#endif
  */

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
