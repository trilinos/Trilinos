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

#include "simple_test_helpers.hpp"
#include "Belos_NonblockGmres.hpp"
#include "Belos_LinearProblem.hpp"
#include "Belos_NativeNormStatusTest.hpp"
#include "Belos_ResidualNormStatusTest.hpp"
#include "Belos_SummaryOutputterStatusTest.hpp"
#include "TSFCoreSolversGMRESSolver.hpp"
#include "TSFCoreSolversSummaryOutputter.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"

int main(int argc, char *argv[]) {

	typedef double Scalar;
	typedef Teuchos::ScalarTraits<Scalar>  ST;
	typedef ST::magnitudeType ScalarMag;
	using Teuchos::CommandLineProcessor;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  Teuchos::Time timer("");
  bool success = true;
  bool verbose = true;

  try {
		SimpleLinearProblemGenerator<Scalar> linearProblemGenerator;
		//
		// A) Read from the command line
		//
		bool            dumpAll      = false;
		int             maxNumIters  = 16;
		int             maxKylovDim  = 16;
		bool            useOldGmres  = false;
		int             outputPrec   = 4;
		CommandLineProcessor  clp(false); // Don't throw exceptions
		linearProblemGenerator.setupCLP(&verbose,&dumpAll,&clp);
		clp.setOption( "max-num-iters", &maxNumIters, "The maximum number of iterative solver iterations allowed." );
		clp.setOption( "max-krylov-dim", &maxKylovDim, "Maximum dimension of the Krylov subspace before a restart is performed" );
		clp.setOption( "use-old-gmres", "no-use-old-gmres", &useOldGmres, "Use the old GMRES solver or not" );
		clp.setOption( "output-prec", &outputPrec, "Number of digits to output" );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
		//
		if(verbose) {
			std::cout << "\n*** Simple test program for Belos::NonblockGmres solver\n";
			std::cout << std::setprecision(outputPrec) << std::scientific;
			std::cout << "\nEchoing options ...\n";
			clp.printHelpMessage("nbl_pgmres_tsfcore_simple.exe",std::cout);
		}
		//
		// B) Setup linear problem complete with operator, RHS, LHS
		// status test.
		//
		std::vector<ScalarMag> tols;
		Belos::LinearProblem<Scalar> lpup;
		linearProblemGenerator.setupLinearProblemEtc(
			rcp(&std::cout,false),verbose,dumpAll
			,&tols,&lpup
			);
		//
		// C) Solve the linear system using old TSFCore::Solvers::GMRESSolver
		//
		if(useOldGmres) {
			if(verbose)
				std::cout << "\nSolving linear system with TSFCore::Solvers::GMRESSolver ...\n";
			TSFCore::Solvers::GMRESSolver<Scalar>  oldGmresSolver;
			oldGmresSolver.set_out(rcp(new std::ofstream("oldGmresSolver.out")));
			oldGmresSolver.dump_all(dumpAll);
			TSFCore::Solvers::SummaryOutputter<Scalar> sumOut(rcp(&std::cout,false)," ");
			assign( &*lpup.getLhs(), ST::zero() );
			oldGmresSolver.solve(
				*lpup.getOperator().op(),lpup.getOperator().defaultTrans()
				,*lpup.getRhs(),&*lpup.getLhs(),ST::one(),maxNumIters
				,&sumOut,NULL,TSFCore::NOTRANS,NULL,TSFCore::NOTRANS
				);
			const Belos::LinearProblemState<Scalar> &lps = lpup;
			if(!checkResidual(
					 lps.getOperator(),lps.getRhs(),lps.getLhs()
					 ,verbose,tols,std::cout)
				) success = false;
		}
		//
		// D) Solve the linear system using the new Belos::NonblockGmres
		//
		if(verbose)
			std::cout << "\nSolving linear system with Belos::NonblockGmres ...\n";
		Belos::NonblockGmres<Scalar> newGmresSolver(maxKylovDim);
		newGmresSolver.set_out(rcp(new std::ofstream("newGmresSolver.out")));
		newGmresSolver.dump_all(dumpAll);
		assign( &*lpup.getLhs(), ST::zero() );
		newGmresSolver.setProblem(rcp(&lpup,false));
		Belos::IterateReturn newGmresSolverReturn = newGmresSolver.iterate(maxNumIters);
		if(verbose) {
			std::cout
				<< "\nBelos::NonblockGmres::iterate() returned:"
				<< "\n  return.iterateTermination  = " << toString(newGmresSolverReturn.iterateTermination)
				<< "\n  return.numCumulativeIter() = " << newGmresSolverReturn.numCumulativeIter
				<< std::endl;
		}
		// E) Compute actual residual norms relative to B: R = A*X - B
		const Belos::LinearProblemState<Scalar> &lps = lpup;
		if(!checkResidual(
				 lps.getOperator(),lps.getRhs(),lps.getLhs()
				 ,verbose,tols,std::cout)
			) success = false;
  }
  catch ( const std::exception &excpt ) {
    if (verbose)
      std::cerr << "*** Caught a standard exception : "<< excpt.what() << std::endl;
    success = false;
  }
  catch ( ... ) {
    if (verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }
  if ( success ) {
    std::cout<< "\n***************** The test PASSED !!!********************"<<std::endl;
    return 0;
  }
  std::cout<< "\n********************The test FAILED!!! ********************"<<std::endl;
  return 1; 

}

