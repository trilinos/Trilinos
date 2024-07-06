// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero. 
// The problem is solver for multiple scalar types, and timings are reported.
//
// NOTE: No preconditioner is used in this case. 

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DiagPrecond.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

using namespace Teuchos;
using namespace Belos;

using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::DiagPrecond;
using Tpetra::MultiVector;
using Tpetra::Vector;
using Tpetra::Map;

using std::endl;
using std::cout;
using std::string;
using std::setw;
using std::vector;

using Teuchos::tuple;

bool proc_verbose = false, reduce_tol, precond = true, dumpdata = false;
RCP<Map<int> > vmap;
RCP<MultiVector<float,int> > actX;
ParameterList mptestpl; 
int rnnzmax;
int dim, numrhs; 
float *fvals;
int *colptr, *rowind;
int mptestmypid = 0;

///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<LinearProblem<Scalar,MultiVector<Scalar,int>,Operator<Scalar,int> > > buildProblem()
{
  using SCT = typename ScalarTraits<Scalar>;
  using MT = typename SCT::magnitudeType;
  using OP = typename Operator<Scalar,int>;
  using MV = typename MultiVector<Scalar,int>;
  using OPT = typename OperatorTraits<Scalar,MV,OP>;
  using MVT = typename MultiVecTraits<Scalar,MV>;

  RCP<CrsMatrix<Scalar,int> > A = rcp(new CrsMatrix<Scalar,int>(*vmap,rnnzmax));
  if (mptestmypid == 0) {
    // HB format is compressed column. CrsMatrix is compressed row.
    const MT *fptr = fvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        A->insertGlobalValues(*rptr-1,tuple(c),tuple<Scalar>(*fptr));
        A->insertGlobalValues(c,tuple(*rptr-1),tuple<Scalar>(*fptr));
        ++rptr;
        ++fptr;
      }
    }
  }
  // Distribute matrix data to other nodes
  A->fillComplete();
  // Create initial MV and solution MV
  RCP<MV> B, X;
  X = rcp( new MV(*vmap,numrhs) );
  // Set LHS to actX
  {
    typename MultiVector<Scalar,int>::double_pointer Xvals = X->extractView2D();
    typename MultiVector<Scalar,int>::const_double_pointer actXvals = actX->extractConstView2D();
    for (Teuchos_Ordinal j=0; j<numrhs; ++j) {
      std::copy(actXvals[j], actXvals[j]+actX->myLength(), Xvals[j]);
    }
  }
  B = rcp( new MV(*vmap,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );

  // Construct a linear problem instance with zero initial MV
  RCP<LinearProblem<Scalar,MV,OP> > problem = rcp( new LinearProblem<Scalar,MV,OP>(A,X,B) );
  problem->setLabel(Teuchos::typeName(SCT::one()));
  // Diagonal preconditioner
  if (precond) {
    Vector<Scalar,int> diags(A->getRowMap());
    A->getLocalDiagCopy(diags);
    for (Teuchos_Ordinal i=0; i<vmap->getNumMyEntries(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(diags[i] <= SCT::zero(), std::runtime_error,"Matrix is not positive-definite: " << diags[i]);
      diags[i] = SCT::one() / diags[i];
    }
    RCP<Operator<Scalar,int> > P = rcp(new DiagPrecond<Scalar,int>(diags));
    problem->setRightPrec(P);
  }
  TEUCHOS_TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}

///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters, RCP<MultiVector<Scalar,int> > &X) 
{

  using SCT = typename ScalarTraits<Scalar>;
  using MT = typename SCT::magnitudeType;
  using OP = typename Operator<Scalar,int>;
  using MV = typename MultiVector<Scalar,int>;
  using OPT = typename OperatorTraits<Scalar,MV,OP>;
  using MVT = typename MultiVecTraits<Scalar,MV>;

  const Scalar ONE  = SCT::one();
  mptestpl.set<MT>( "Convergence Tolerance", ltol );         // Relative convergence tolerance requested
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
    solver = rcp(new PseudoBlockCGSolMgr<Scalar,MV,OP>( problem, rcpFromRef(mptestpl) ));
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
    // Compute actual residuals.
    RCP<const OP> A = problem->getOperator();
    X = problem->getLHS();
    RCP<const MV> B = problem->getRHS();
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(*vmap, numrhs);
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
int main(int argc, char *argv[]) 
{
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // Get test parameters from command-line processor
  bool verbose = false, debug = false;
  bool printX = false;
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
  cmdp.setOption("print-x","no-print-x",&printX,"Print solution multivectors.");
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

  // Get the data (double) from the HB file and build the Map,Matrix
  int nnz, info;
  nnz = -1;
  if (mptestmypid == 0) {
    int dim2;
    double *dvals;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // Truncate data into float
    fvals = new float[nnz];
    std::copy( dvals, dvals+nnz, fvals );
    free(dvals);
    // Find maximum NNZ over all rows
    vector<int> rnnz(dim,0);
    for (int *ri=rowind; ri<rowind+nnz; ++ri) {
      ++rnnz[*ri-1];
    }
    for (int c=0; c<dim; ++c) {
      rnnz[c] += colptr[c+1]-colptr[c];
    }
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
  }
  else {
    // Address uninitialized data warnings
    fvals = NULL;
    colptr = NULL;
    rowind = NULL;
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
  // Create map
  vmap = rcp(new Map<int>(dim,0,comm));
  // Create a consistent LHS
  actX = rcp(new MultiVector<float,int>(*vmap,numrhs));
  MultiVecTraits<float,MultiVector<float,int> >::MvRandom(*actX);
  // Create the parameter list
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

  // Print out information about problem
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

  // Run tests for different scalar types
  RCP<MultiVector<float,int> > Xf;
  RCP<MultiVector<double,int> > Xd;
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter,Xf);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter,Xd);
// #ifdef HAVE_TEUCHOS_QD
//   unsigned int old_cw; fpu_fix_start(&old_cw); // fix floating point rounding for intel processors; required for proper implementation
//   ltol = (reduce_tol ? ltol*ltol : ltol);
//   double ddtime[3]; int dditer; bool ddpass = runTest<dd_real>(ltol,ddtime,dditer);
//   ltol = (reduce_tol ? ltol*ltol : ltol);
//   double qdtime[3]; int qditer; bool qdpass = runTest<qd_real>(ltol,qdtime,qditer);
//   fpu_fix_end(&old_cw);                        // restore previous floating point rounding
// #endif

  if (printX) {
    // doesn't compie; fix
    // RCP<FancyOStream> fos = Teuchos::getFancyOStream(rcp(&std::cout,false));
    // if (mptestmypid==0) cout << "X (float): " << endl;
    // Xf->describe(fos,Teuchos::VERB_EXTREME);
    // if (mptestmypid==0) cout << "X (double): " << endl;
    // Xd->describe(fos,Teuchos::VERB_EXTREME);
  }

  // Test the relative error between Xf and Xd
  if (mptestmypid==0) {
    cout << "Relative error |xdi - xfi|/|xdi|" << endl;
    cout << "--------------------------------" << endl;
  }
  for (int j=0; j<numrhs; ++j) {
    // Compute |xdj-xfj|/|xdj|
    RCP<Vector<float,int> > xfj = (*Xf)(j);
    RCP<Vector<double,int> > xdj = (*Xd)(j);
    double mag = xdj->norm2();
    ArrayView<double> dval;
    ArrayView<const float> fval;
    xdj->extractView1D(dval);
    xfj->extractConstView1D(fval);
    for (int i=0; i<xdj->myLength(); ++i) {
      dval[i] -= (double)fval[i];
    }
    double err = xdj->norm2() / mag;
    if (mptestmypid==0) {
      cout << err << endl;
    }
  }
  Xf = Teuchos::null;
  Xd = Teuchos::null;

  // Done with the matrix data now; delete it
  if (mptestmypid == 0) {
    free( fvals );
    free( colptr );
    free( rowind );
  }

  if (mptestmypid==0) {
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << setw(12) << "float"   << "     " << setw(10) <<  ftime[0] << "     " << setw(9) <<  ftime[1] << "     " << setw(10) <<  ftime[2] << "     " << setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << setw(12) << "double"  << "     " << setw(10) <<  dtime[0] << "     " << setw(9) <<  dtime[1] << "     " << setw(10) <<  dtime[2] << "     " << setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
// #ifdef HAVE_TEUCHOS_QD
//     cout << setw(12) << "dd_real" << "     " << setw(10) << ddtime[0] << "     " << setw(9) << ddtime[1] << "     " << setw(10) << ddtime[2] << "     " << setw(9) << dditer << "     " << (ddpass ? "pass" : "fail") << endl;
//     cout << setw(12) << "qd_real" << "     " << setw(10) << qdtime[0] << "     " << setw(9) << qdtime[1] << "     " << setw(10) << qdtime[2] << "     " << setw(9) << qditer << "     " << (qdpass ? "pass" : "fail") << endl;
// #endif
  }

  if (dumpdata) {
    std::ofstream fout("data.dat",std::ios_base::app | std::ios_base::out);
    fout << setw(14) << nnz << setw(8) << "float"   << setw(12) <<  ftime[0] << setw(12) <<  ftime[1] << setw(12) <<  ftime[2] << setw(12) <<  fiter << endl;
    fout << setw(14) << nnz << setw(8) << "double"  << setw(12) <<  dtime[0] << setw(12) <<  dtime[1] << setw(12) <<  dtime[2] << setw(12) <<  diter << endl;
// #ifdef HAVE_TEUCHOS_QD
//     fout << setw(14) << nnz << setw(8) << "dd_real" << setw(12) << ddtime[0] << setw(12) << ddtime[1] << setw(12) << ddtime[2] << setw(12) << dditer << endl;
//     fout << setw(14) << nnz << setw(8) << "qd_real" << setw(12) << qdtime[0] << setw(12) << qdtime[1] << setw(12) << qdtime[2] << setw(12) << qditer << endl;
// #endif
    fout.close();
  }

  bool allpass = fpass && dpass;
// #ifdef HAVE_TEUCHOS_QD
//   allpass = allpass && ddpass && qdpass;
// #endif

  if (!allpass) {
    if (mptestmypid==0) cout << "\nEnd Result: TEST FAILED" << endl;	
    return -1;
  }
  
  // Default return value
  if (mptestmypid==0) cout << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}
