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
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

using namespace Teuchos;
using Tpetra::DefaultPlatform;
using Tpetra::Platform;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;
using std::endl;
using std::cout;

int main(int argc, char *argv[]) {

  typedef double                           ST;
  typedef ScalarTraits<ST>                SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<int,ST>         OP;
  typedef Tpetra::MultiVector<int,ST>      MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  GlobalMPISession mpisess(&argc,&argv,&cout);

  const ST one  = SCT::one();
  const ST zero = SCT::zero();	

  int info = 0;
  int MyPID = 0;

  RCP<const Platform<int> > platform = DefaultPlatform<int>::getPlatform();
  RCP<Comm<int> > comm = platform->createComm();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, proc_verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;      // total number of right-hand sides to solve for
  int blocksize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  std::string filename("bcsstk14.hb");
  MT tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
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

  if (proc_verbose) {
    std::cout << Belos::Belos_Version() << std::endl << std::endl;
  }

  //
  // Get the data from the HB file and build the Map,Matrix
  //
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  if (MyPID == 0) {
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
  Map<int> map(dim,0,*platform);
  RCP<CrsMatrix<int,ST> > A = rcp(new CrsMatrix<int,ST>(map));
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
  if (MyPID == 0) {
    // Clean up.
    free( dvals );
    free( colptr );
    free( rowind );
  }
  // distribute matrix data to other nodes
  A->fillComplete();

  // Create initial vectors
  RCP<MultiVector<int,ST> > B, X;
  X = rcp( new MultiVector<int,ST>(map,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MultiVector<int,ST>(map,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->globalLength();
  if (maxiters == -1) {
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  }
  //
  ParameterList belosList;
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  int verbLevel = Belos::Errors + Belos::Warnings;
  if (debug) {
    verbLevel += Belos::Debug;
  }
  if (verbose) {
    verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
  }
  belosList.set( "Verbosity", verbLevel );
  if (verbose) {
    if (frequency > 0) {
      belosList.set( "Output Frequency", frequency );
    }
  }
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  // *******************************************************************
  // *************Start the block CG iteration***********************
  // *******************************************************************
  //
  Belos::BlockCGSolMgr<ST,MV,OP> solver( rcp(&problem,false), rcp(&belosList,false) );

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of CG iterations: " << maxiters << std::endl; 
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver.solve();
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MultiVector<int,ST> resid(map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -one, resid, one, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  }
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
    if (proc_verbose) {
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
    }
    if (actRes > tol) badRes = true;
  }

  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
  //
} // end test_bl_cg_hb.cpp
