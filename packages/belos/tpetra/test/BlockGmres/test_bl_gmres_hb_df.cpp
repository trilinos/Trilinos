//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER
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
#include "BelosBlockGmresSolMgr.hpp"

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>

using namespace Teuchos;
using namespace Belos;
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
using Teuchos::tuple;

bool proc_verbose = false, reduce_tol, precond = true, dumpdata = false;
RCP<Map<int> > vmap;
ParameterList mptestpl; 
int rnnzmax;
int mptestdim, numrhs; 
double *dvals; 
int *colptr, *rowind;
int mptestmypid = 0;
int mptestnumimages = 1;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar> 
RCP<LinearProblem<Scalar,MultiVector<Scalar,int>,Operator<Scalar,int> > > buildProblem()
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType  MT;
  typedef Operator<Scalar,int>         OP;
  typedef MultiVector<Scalar,int>      MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;
  RCP<CrsMatrix<Scalar,int> > A = rcp(new CrsMatrix<Scalar,int>(vmap,rnnzmax));
  if (mptestmypid == 0) {
    // HB format is compressed column. CrsMatrix is compressed row.
    const double *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<mptestdim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        A->insertGlobalValues(*rptr-1,tuple(c),tuple<Scalar>(*dptr));
        if (c != *rptr -1) {
          A->insertGlobalValues(c,tuple(*rptr-1),tuple<Scalar>(*dptr));
        }
        ++rptr;
        ++dptr;
      }
    }
  }
  // distribute matrix data to other nodes
  A->fillComplete();
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
  // diagonal preconditioner
  // if (precond) {
  //   Vector<Scalar,int> diags(A->getRowMap());
  //   A->getLocalDiagCopy(diags);
  //   for (Teuchos_Ordinal i=0; i<vmap->getNumMyEntries(); ++i) {
  //     TEUCHOS_TEST_FOR_EXCEPTION(diags[i] <= SCT::zero(), std::runtime_error,"Matrix is not positive-definite: " << diags[i]);
  //     diags[i] = SCT::one() / diags[i];
  //   }
  //   RCP<Operator<Scalar,int> > P = rcp(new DiagPrecond<Scalar,int>(diags));
  //   problem->setRightPrec(P);
  // }
  TEUCHOS_TEST_FOR_EXCEPT(problem->setProblem() == false);
  return problem;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Scalar>
bool runTest(double ltol, double times[], int &numIters) 
{
  typedef ScalarTraits<Scalar>         SCT;
  typedef typename SCT::magnitudeType   MT;
  typedef Operator<Scalar,int> 		OP;
  typedef MultiVector<Scalar,int>       MV;
  typedef OperatorTraits<Scalar,MV,OP> OPT;
  typedef MultiVecTraits<Scalar,MV>    MVT;

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
int main(int argc, char *argv[]) 
{
  GlobalMPISession mpisess(&argc,&argv,&cout);

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

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
    info = readHB_newmat_double(filename.c_str(),&mptestdim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // find maximum NNZ over all rows
    vector<int> rnnz(mptestdim,0);
    for (int *ri=rowind; ri<rowind+nnz; ++ri) {
      ++rnnz[*ri-1];
    }
    for (int c=0; c<mptestdim; ++c) {
      rnnz[c] += colptr[c+1]-colptr[c];
    }
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    colptr = NULL;
    rowind = NULL;
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
  vmap = rcp(new Map<int>(mptestdim,0,comm));
  //
  mptestpl.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  mptestpl.set( "Num Blocks", numblocks);
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
    free( dvals );
    free( colptr );
    free( rowind );
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


