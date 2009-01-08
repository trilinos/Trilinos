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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_ScalarTraits.hpp>

#ifdef HAVE_TEUCHOS_QD
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <qd/fpu.h>
#endif

using namespace Teuchos;
using namespace Belos;
using Tpetra::DefaultPlatform;
using Tpetra::Platform;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using std::endl;
using std::cout;
using std::string;

bool proc_verbose = false, reduce_tol;
RCP<Map<int> > vmap;
ParameterList pl; 
int dim, numrhs; 
double *dvals; 
int *colptr, *rowind;
int MyPID = 0;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<LinearProblem<Scalar,MultiVector<int,Scalar>,Operator<int,Scalar> > > buildProblem()
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<int,Scalar>         OP;
  typedef MultiVector<int,Scalar>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;
  RCP<CrsMatrix<int,Scalar> > A = rcp(new CrsMatrix<int,Scalar>(*vmap));
  if (MyPID == 0) {
    // HB format is compressed column. CrsMatrix is compressed row.
    const double *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        A->submitEntry(*rptr-1,c,*dptr);
        A->submitEntry(c,*rptr-1,*dptr);
        ++rptr;
        ++dptr;
      }
    }
  }
  // distribute matrix data to other nodes
  A->fillComplete();
  // Create initial MV and solution MV
  RCP<MV> B, X;
  X = rcp( new MV(*vmap,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MV(*vmap,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );
  // Construct an unpreconditioned linear problem instance with zero initial MV
  RCP<LinearProblem<Scalar,MV,OP> > problem = rcp( new LinearProblem<Scalar,MV,OP>(A,X,B) );
  TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters) 
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<int,Scalar>         OP;
  typedef MultiVector<int,Scalar>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;

  const Scalar ONE  = SCT::one();
  pl.set<MT>( "Convergence Tolerance", ltol );         // Relative convergence tolerance requested

  if (MyPID==0) cout << "Testing Scalar == " << typeName(ONE) << endl;

  RCP<LinearProblem<Scalar,MV,OP> > problem;
  Time btimer("Build Timer"), ctimer("Construct Timer"), stimer("Solve Timer");
  ReturnType ret;
  if (MyPID==0) cout << "Building problem..." << endl;
  { 
    TimeMonitor localtimer(btimer);
    problem = buildProblem<Scalar>();
  }
  RCP<SolverManager<Scalar,MV,OP> > solver;
  if (MyPID==0) cout << "Constructing solver..." << endl; 
  {
    TimeMonitor localtimer(ctimer);
    solver = rcp(new PseudoBlockCGSolMgr<Scalar,MV,OP>( problem, rcp(&pl,false) ));
  }
  if (MyPID==0) cout << "Solving problem..." << endl;
  {
    TimeMonitor localtimer(stimer);
    ret = solver->solve();
  }
  numIters = solver->getNumIters();
  times[0] = btimer.totalElapsedTime();
  times[1] = ctimer.totalElapsedTime();
  times[2] = stimer.totalElapsedTime();

  bool badRes = false;
  if (ret == Converged) {
    if (proc_verbose) cout << endl;
    if (MyPID==0) cout << "Computing residuals..." << endl;
    //
    // Compute actual residuals.
    //
    RCP<const OP> A = problem->getOperator();
    RCP<MV> X = problem->getLHS();
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
///////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) 
{
  GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Platform<int> > platform = DefaultPlatform<int>::getPlatform();
  RCP<Comm<int> > comm = platform->createComm();

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
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  MyPID = rank(*comm);
  proc_verbose = ( verbose && (MyPID==0) );
  if (proc_verbose) cout << Belos_Version() << endl << endl;

  //
  // Get the data (double) from the HB file and build the Map,Matrix
  //
  int nnz, info;
  nnz = -1;
  if (MyPID == 0) {
    int dim2;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    colptr = NULL;
    rowind = NULL;
  }
  broadcast(*comm,0,&info);
  broadcast(*comm,0,&nnz);
  broadcast(*comm,0,&dim);
  if (info == 0 || nnz < 0) {
    if (MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }

  // create map
  vmap = rcp(new Map<int>(dim,0,*platform));
  // create the parameter list
  if (maxiters == -1) {
    maxiters = dim/blocksize - 1; // maximum number of iterations to run
  }
  //
  pl.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  pl.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  int verbLevel = Errors + Warnings;
  if (debug) {
    verbLevel += Debug;
  }
  if (verbose) {
    verbLevel += TimingDetails + FinalSummary + StatusTestDetails;
  }
  pl.set( "Verbosity", verbLevel );
  if (verbose) {
    if (frequency > 0) {
      pl.set( "Output Frequency", frequency );
    }
  }

  //
  // **********Print out information about problem*******************
  //
  if (MyPID==0) {
    cout << "Dimension of matrix: " << dim << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Max number of CG iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }

  // run tests for different scalar types
  double ltol = tol;
  double ftime[3]; int fiter; bool fpass = runTest<float>(ltol,ftime,fiter);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double dtime[3]; int diter; bool dpass = runTest<double>(ltol,dtime,diter);
#ifdef HAVE_TEUCHOS_QD
  unsigned int old_cw; fpu_fix_start(&old_cw); // fix floating point rounding for intel processors; required for proper implementation
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double ddtime[3]; int dditer; bool ddpass = runTest<dd_real>(ltol,ddtime,dditer);
  ltol = (reduce_tol ? ltol*ltol : ltol);
  double qdtime[3]; int qditer; bool qdpass = runTest<qd_real>(ltol,qdtime,qditer);
  fpu_fix_end(&old_cw);                        // restore previous floating point rounding
#endif

  // done with the matrix data now; delete it
  if (MyPID == 0) {
    free( dvals );
    free( colptr );
    free( rowind );
  }

  if (MyPID==0) {
    cout << "Scalar field     Build time     Init time     Solve time     Num Iters     Test Passsed" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << std::setw(12) << "float"   << "     " << std::setw(10) <<  ftime[0] << "     " << std::setw(9) <<  ftime[1] << "     " << std::setw(10) <<  ftime[2] << "     " << std::setw(9) <<  fiter << "     " << ( fpass ? "pass" : "fail") << endl;
    cout << std::setw(12) << "double"  << "     " << std::setw(10) <<  dtime[0] << "     " << std::setw(9) <<  dtime[1] << "     " << std::setw(10) <<  dtime[2] << "     " << std::setw(9) <<  diter << "     " << ( dpass ? "pass" : "fail") << endl;
#ifdef HAVE_TEUCHOS_QD
    cout << std::setw(12) << "dd_real" << "     " << std::setw(10) << ddtime[0] << "     " << std::setw(9) << ddtime[1] << "     " << std::setw(10) << ddtime[2] << "     " << std::setw(9) << dditer << "     " << (ddpass ? "pass" : "fail") << endl;
    cout << std::setw(12) << "qd_real" << "     " << std::setw(10) << qdtime[0] << "     " << std::setw(9) << qdtime[1] << "     " << std::setw(10) << qdtime[2] << "     " << std::setw(9) << qditer << "     " << (qdpass ? "pass" : "fail") << endl;
#endif
  }

  bool allpass = fpass && dpass;
#ifdef HAVE_TEUCHOS_QD
  allpass = allpass && ddpass && qdpass;
#endif

  if (!allpass) {
    if (MyPID==0) cout << "\nEnd Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (MyPID==0) cout << "\nEnd Result: TEST PASSED" << endl;
  return 0;
}


