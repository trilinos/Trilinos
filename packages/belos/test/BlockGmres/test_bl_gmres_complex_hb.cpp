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
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_BELOS_TRIUTILS
#include "Trilinos_Util_iohb.h"
#endif

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) {
  //
#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef std::complex<double> ST;
#else
  std::cout << "Not compiled with std::complex support." << std::endl;
  std::cout << "End Result: TEST FAILED" << std::endl;
  return EXIT_FAILURE;
#endif

  typedef ScalarTraits<ST>                 SCT;
  typedef SCT::magnitudeType                MT;
  typedef Belos::MultiVec<ST>               MV;
  typedef Belos::Operator<ST>               OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
  ST one  = SCT::one();
  ST zero = SCT::zero();

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  int MyPID = session.getRank();
  //
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool success = false;
  bool verbose = false;
  try {
    int info = 0;
    bool norm_failure = false;
    bool proc_verbose = false;
    bool pseudo = false;   // use pseudo block GMRES to solve this linear system.
    int frequency = -1;  // how often residuals are printed by solver
    int blocksize = 1;
    int numrhs = 1;
    int maxrestarts = 15;
    int length = 50;
    std::string filename("mhd1280b.cua");
    MT tol = 1.0e-5;  // relative residual tolerance

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("pseudo","regular",&pseudo,"Use pseudo-block GMRES to solve the linear systems.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("num-restarts",&maxrestarts,"Maximum number of restarts allowed for the GMRES solver.");
    cmdp.setOption("blocksize",&blocksize,"Block size used by GMRES.");
    cmdp.setOption("subspace-length",&length,"Maximum dimension of block-subspace used by GMRES solver.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return EXIT_FAILURE;
    }

    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */
    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

#ifndef HAVE_BELOS_TRIUTILS
    std::cout << "This test requires Triutils. Please configure with --enable-triutils." << std::endl;
    if (MyPID==0) {
      std::cout << "End Result: TEST FAILED" << std::endl;
    }
    return EXIT_FAILURE;
#endif

    // Get the data from the HB file
    int dim,dim2,nnz;
    MT *dvals;
    int *colptr,*rowind;
    ST *cvals;
    nnz = -1;
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
        &colptr,&rowind,&dvals);
    if (info == 0 || nnz < 0) {
      if (MyPID==0) {
        std::cout << "Error reading '" << filename << "'" << std::endl;
        std::cout << "End Result: TEST FAILED" << std::endl;
      }
      return EXIT_FAILURE;
    }
    // Convert interleaved doubles to std::complex values
    cvals = new ST[nnz];
    for (int ii=0; ii<nnz; ii++) {
      cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
    }
    // Build the problem matrix
    RCP< MyBetterOperator<ST> > A
      = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,cvals) );
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    int maxits = dim/blocksize; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Num Blocks", length );                 // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxits );         // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // Construct the right-hand side and solution multivectors.
    // NOTE:  The right-hand side will be constructed such that the solution is
    // a vectors of one.
    //
    RCP<MyMultiVec<ST> > soln = rcp( new MyMultiVec<ST>(dim,numrhs) );
    RCP<MyMultiVec<ST> > rhs = rcp( new MyMultiVec<ST>(dim,numrhs) );
    MVT::MvRandom( *soln );
    OPT::Apply( *A, *soln, *rhs );
    MVT::MvInit( *soln, zero );
    //
    //  Construct an unpreconditioned linear problem instance.
    //
    RCP<Belos::LinearProblem<ST,MV,OP> > problem =
      rcp( new Belos::LinearProblem<ST,MV,OP>( A, soln, rhs ) );
    bool set = problem->setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    //
    // *******************************************************************
    // *************Start the block Gmres iteration***********************
    // *******************************************************************
    //
    Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > solver;
    if (pseudo)
      solver = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( problem, Teuchos::rcp(&belosList,false) ) );
    else
      solver = Teuchos::rcp( new Belos::BlockGmresSolMgr<ST,MV,OP>( problem, Teuchos::rcp(&belosList,false) ) );
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << dim << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of Gmres iterations: " << maxits << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver->solve();
    //
    // Compute actual residuals.
    //
    RCP<MyMultiVec<ST> > temp = rcp( new MyMultiVec<ST>(dim,numrhs) );
    OPT::Apply( *A, *soln, *temp );
    MVT::MvAddMv( one, *rhs, -one, *temp, *temp );
    std::vector<MT> norm_num(numrhs), norm_denom(numrhs);
    MVT::MvNorm( *temp, norm_num );
    MVT::MvNorm( *rhs, norm_denom );
    for (int i=0; i<numrhs; ++i) {
      if (proc_verbose)
        std::cout << "Relative residual "<<i<<" : " << norm_num[i] / norm_denom[i] << std::endl;
      if ( norm_num[i] / norm_denom[i] > tol ) {
        norm_failure = true;
      }
    }

    // Clean up.
    delete [] dvals;
    delete [] colptr;
    delete [] rowind;
    delete [] cvals;

    success = ret==Belos::Converged && !norm_failure;
    if (success) {
      if (proc_verbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_bl_gmres_complex_hb.cpp
