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
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#ifdef HAVE_TEUCHOS_QD
  #include <qd/dd_real.h>
  #include <qd/qd_real.h>
  #include <qd/fpu.h>
#endif

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
RCP<MultiVector<float,int> > actX;
ParameterList mptestpl; 
size_t rnnzmax;
int dim, numrhs; 
ArrayRCP<float> fvals;
ArrayRCP<const int> offsets, colinds;
int mptestmypid = 0;

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
    for (int row=0; row < dim; ++row) {
      const size_t nE = offsets[row+1] - offsets[row];
      if (nE > 0) {
        // convert row values from double to Scalar
        std::copy( fvals+offsets[row], fvals+offsets[row]+nE, vals.begin() );
        // add row to matrix
        A->insertGlobalValues( row, colinds(offsets[row], nE), vals(0,nE) );
      }
    }
  }
  // distribute matrix data to other nodes
  A->fillComplete(Tpetra::DoOptimizeStorage);
  // Create initial MV and solution MV
  RCP<MV> B, X;
  X = rcp( new MV(vmap,numrhs) );
  const size_t myLen = X->getLocalLength();
  // Set LHS to actX
  {
    ArrayRCP<ArrayRCP<     Scalar> >    Xvals =    X->get2dViewNonConst();
    ArrayRCP<ArrayRCP<const float> > actXvals = actX->get2dView();
    for (int j=0; j<numrhs; ++j) {
      std::copy(actXvals[j], actXvals[j] + myLen, Xvals[j]);
    }
  }
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
bool runTest(double ltol, double times[], int &numIters, RCP<MultiVector<Scalar,int> > &X) {
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar,int>         OP;
  typedef MultiVector<Scalar,int>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;

  const Scalar ONE  = SCT::one();
  mptestpl.set<MT>( "Convergence Tolerance", (MT) (ltol) );         // Relative convergence tolerance requested
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
    solver = rcp(new PseudoBlockCGSolMgr<Scalar,MV,OP>( problem, rcp(&mptestpl,false) ));
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
    RCP<const OP> A = problem->getOperator();
    X = problem->getLHS();
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
int main(int argc, char *argv[]) 
{
  GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  numrhs = 1;      // total number of right-hand sides to solve for
  int blocksize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  string filename("bcsstk14.hb");
  double tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
  cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
  cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  mptestmypid = rank(*comm);
  proc_verbose = ( verbose && (mptestmypid==0) );
  if (proc_verbose) cout << Belos_Version() << endl << endl;

  //
  // Get the data (double) from the HB file and build the Map,Matrix
  //
  int nnz, info;
  nnz = -1;
  if (mptestmypid == 0) {
    int dim2;
    double *dvals;
    int *colptr, *rowind;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // truncate data into float
    fvals = arcp<float>(nnz);
    std::copy( dvals, dvals+nnz, fvals.begin() );
    free(dvals); dvals = NULL;
    // Find number of non-zeros for each row
    vector<size_t> rnnz(dim,0);
    // Move through all row indices, adding the contribution to the appropriate row
    // Skip the diagonals, we'll catch them below on the column pass.
    // Remember, the file uses one-based indexing. We'll convert it to zero-based later.
    int curcol_0 = 0, currow_0;
    for (int curnnz_1=1; curnnz_1 <= nnz; ++curnnz_1) {
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
    ArrayRCP<float> newfvals = arcp<float>(totalnnz);
    ArrayRCP<int> newoffs = arcp<int>(dim+1);
    ArrayRCP<int> newinds = arcp<int>(totalnnz);
    // set up pointers
    newoffs[0] = 0;
    for (size_t row=1; row != static_cast<size_t>(dim+1); ++row) {
      newoffs[row] = newoffs[row-1] + rnnz[row-1];
    }
    // reorganize data from column oriented to row oriented, duplicating symmetric part as well
    // rowind are one-based; account for that, and convert them to zero-based.
    // use nnz as the number of entries add per row thus far
    std::fill( rnnz.begin(), rnnz.end(), 0 );
    curcol_0 = 0;
    for (int curnnz_1=1; curnnz_1 <= nnz; ++curnnz_1) {
      // if colptr[c] <= curnnz_1 < colptr[c+1], then curnnz_1 belongs to column c
      // invariant: curcol_0 is the column for curnnz_1, i.e., curcol_0 is smallest number such that colptr[curcol_0+1] > curnnz_1
      while (colptr[curcol_0+1] <= curnnz_1) ++curcol_0;
      // entry curnnz_1 corresponds to (curcol_0, rowind[curnnz_1]) and (rowind[curnnz_1], curcol_0)
      // it must be added twice if curcol_0 != rowind[curnnz_1]
      currow_0 = rowind[curnnz_1-1] - 1;
      // add (currow_0, curcol_0)
      const float curval = fvals[curnnz_1-1];
      size_t newnnz = newoffs[currow_0] + rnnz[currow_0];
      newfvals[newnnz] = curval;
      newinds [newnnz] = curcol_0;
      ++rnnz[currow_0];
      if (curcol_0 != currow_0) {
        newnnz = newoffs[curcol_0] + rnnz[curcol_0];
        newfvals[newnnz] = curval;
        newinds [newnnz] = currow_0;
        ++rnnz[curcol_0];
      }
    }
    const size_t totalnnz2 = std::accumulate( rnnz.begin(), rnnz.end(), 0 );
    TEST_FOR_EXCEPT( totalnnz2 != totalnnz );
    // free the original data, point to new dada
    fvals = newfvals;
    offsets = newoffs;
    colinds = newinds;
  }
  broadcast(*comm,0,&info);
  broadcast(*comm,0,&nnz);
  broadcast(*comm,0,&dim);
  broadcast(*comm,0,&rnnzmax);
  if (info == 0 || nnz < 0) {
    if (mptestmypid == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // create map
  vmap = rcp(new Map<int>(static_cast<global_size_t>(dim),static_cast<int>(0),comm));
  // create a consistent LHS
  actX = rcp(new MultiVector<float,int>(vmap,numrhs));
  MultiVecTraits<float,MultiVector<float,int> >::MvRandom(*actX);
  // create the parameter list
  if (maxiters == -1) {
    maxiters = dim/blocksize - 1; // maximum number of iterations to run
  }
  //
  mptestpl.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  mptestpl.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
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
    cout << "Dimension of matrix: " << dim << endl;
    cout << "Number of nonzeros: " << nnz << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Max number of CG iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  RCP<MultiVector<float,int> > Xf;
  RCP<MultiVector<double,int> > Xd;
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter,Xf);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter,Xd);
#ifdef HAVE_TEUCHOS_QD
  RCP<MultiVector<dd_real,int> > Xdd;
  RCP<MultiVector<qd_real,int> > Xqd;
  unsigned int old_cw; fpu_fix_start(&old_cw); // fix floating point rounding for intel processors; required for proper implementation
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double ddtime[3]; int dditer; bool ddpass = runTest<dd_real>(ltol,ddtime,dditer,Xdd);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double qdtime[3]; int qditer; bool qdpass = runTest<qd_real>(ltol,qdtime,qditer,Xqd);
  fpu_fix_end(&old_cw);                        // restore previous floating point rounding
#endif

  // test the relative error between Xf and Xd
  if (mptestmypid==0) {
    cout << "Relative error |xdi - xfi|/|xdi|" << endl;
    cout << "--------------------------------" << endl;
  }
  if (Xd != Teuchos::null && Xf != Teuchos::null) {
    for (int j=0; j<numrhs; ++j) {
      // compute |xdj-xfj|/|xdj|
      const size_t myLen = Xd->getLocalLength();
      RCP<const Vector<float,int> > xfj = Xf->getVector(j);
      RCP<Vector<double,int> > xdj = Xd->getVectorNonConst(j);
      double mag = xdj->norm2();
      ArrayRCP<double> dval = xdj->get1dViewNonConst();
      ArrayRCP<const float> fval = xfj->get1dView();
      for (size_t i=0; i < myLen; ++i) {
        dval[i] -= (double)fval[i];
      }
      fval = Teuchos::null;
      dval = Teuchos::null;
      double err = xdj->norm2() / mag;
      if (mptestmypid==0) {
        cout << err << endl;
      }
    }
  }
  Xf = Teuchos::null;
  Xd = Teuchos::null;

  // done with the matrix data now; delete it
  fvals = Teuchos::null;
  offsets = Teuchos::null;
  colinds = Teuchos::null;

  if (mptestmypid==0) {
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << setw(12) << "float"   << "     " << setw(10) <<  ftime[0] << "     " << setw(9) <<  ftime[1] << "     " << setw(10) <<  ftime[2] << "     " << setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << setw(12) << "double"  << "     " << setw(10) <<  dtime[0] << "     " << setw(9) <<  dtime[1] << "     " << setw(10) <<  dtime[2] << "     " << setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
#ifdef HAVE_TEUCHOS_QD
    cout << setw(12) << "dd_real" << "     " << setw(10) << ddtime[0] << "     " << setw(9) << ddtime[1] << "     " << setw(10) << ddtime[2] << "     " << setw(9) << dditer << "     " << (ddpass ? "pass" : "fail") << endl;
    cout << setw(12) << "qd_real" << "     " << setw(10) << qdtime[0] << "     " << setw(9) << qdtime[1] << "     " << setw(10) << qdtime[2] << "     " << setw(9) << qditer << "     " << (qdpass ? "pass" : "fail") << endl;
#endif
  }

  if (dumpdata) {
    std::ofstream fout("data.dat",std::ios_base::app | std::ios_base::out);
    fout << setw(14) << nnz << setw(8) << "float"   << setw(12) <<  ftime[0] << setw(12) <<  ftime[1] << setw(12) <<  ftime[2] << setw(12) <<  fiter << endl;
    fout << setw(14) << nnz << setw(8) << "double"  << setw(12) <<  dtime[0] << setw(12) <<  dtime[1] << setw(12) <<  dtime[2] << setw(12) <<  diter << endl;
#ifdef HAVE_TEUCHOS_QD
    fout << setw(14) << nnz << setw(8) << "dd_real" << setw(12) << ddtime[0] << setw(12) << ddtime[1] << setw(12) << ddtime[2] << setw(12) << dditer << endl;
    fout << setw(14) << nnz << setw(8) << "qd_real" << setw(12) << qdtime[0] << setw(12) << qdtime[1] << setw(12) << qdtime[2] << setw(12) << qditer << endl;
#endif
    fout.close();
  }

  bool allpass = fpass && dpass;
#ifdef HAVE_TEUCHOS_QD
  allpass = allpass && ddpass && qdpass;
#endif

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


