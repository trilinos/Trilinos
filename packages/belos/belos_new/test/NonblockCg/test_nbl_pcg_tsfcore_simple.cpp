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

#include "../NonblockGmres/simple_test_helpers.hpp"
#include "Belos_NonblockCg.hpp"
#include "Belos_LinearProblem.hpp"
#include "Belos_NativeNormStatusTest.hpp"
#include "Belos_ResidualNormStatusTest.hpp"
#include "Belos_SummaryOutputterStatusTest.hpp"
#include "TSFCoreSolversCGSolver.hpp"
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
		bool            useOldCg     = false;
		int             outputPrec   = static_cast<int>(ST::t())/2;
		CommandLineProcessor  clp(false); // Don't throw exceptions
		linearProblemGenerator.setupCLP(&verbose,&dumpAll,&clp);
		clp.setOption( "max-num-iters", &maxNumIters, "The maximum number of iterative solver iterations allowed." );
		clp.setOption( "use-old-cg", "no-use-old-cg", &useOldCg, "Use the old CG solver or not" );
		clp.setOption( "output-prec", &outputPrec, "Number of digits to output" );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
		//
		if(verbose) {
			std::cout << "\n*** Simple test program for Belos::NonblockCg solver\n";
			std::cout << std::setprecision(outputPrec) << std::scientific;
			std::cout << "\nEchoing options ...\n";
			clp.printHelpMessage("nbl_pcg_tsfcore_simple.exe",std::cout);
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
		// C) Solve the linear system using old TSFCore::Solvers::CGSolver
		//
		if(useOldCg) {
			if(verbose)
				std::cout << "\nSolving linear system with TSFCore::Solvers::CGSolver ...\n";
			TSFCore::Solvers::CGSolver<Scalar>  oldCgSolver;
			oldCgSolver.set_out(rcp(new std::ofstream("oldCgSolver.out")));
			oldCgSolver.dump_all(dumpAll);
			TSFCore::Solvers::SummaryOutputter<Scalar> sumOut(rcp(&std::cout,false)," ");
			assign( &*lpup.getLhs(), ST::zero() );
			oldCgSolver.solve(
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
		// D) Solve the linear system using the new Belos::NonblockCg
		//
		if(verbose)
			std::cout << "\nSolving linear system with Belos::NonblockCg ...\n";
		Belos::NonblockCg<Scalar> newCgSolver;
		newCgSolver.set_out(rcp(new std::ofstream("newCgSolver.out")));
		newCgSolver.dump_all(dumpAll);
		assign( &*lpup.getLhs(), ST::zero() );
		newCgSolver.setProblem(rcp(&lpup,false));
		Belos::IterateReturn newCgSolverReturn = newCgSolver.iterate(maxNumIters);
		if(verbose) {
			std::cout
				<< "\nBelos::NonblockCg::iterate() returned:"
				<< "\n  return.iterateTermination  = " << toString(newCgSolverReturn.iterateTermination)
				<< "\n  return.numCumulativeIter() = " << newCgSolverReturn.numCumulativeIter
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

