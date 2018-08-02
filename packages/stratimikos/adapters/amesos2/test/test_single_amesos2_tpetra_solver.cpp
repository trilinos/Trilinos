/*
// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "test_single_amesos2_tpetra_solver.hpp"

#include "Thyra_Amesos2LinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_ParameterList.hpp"

bool Thyra::test_single_amesos2_tpetra_solver(
  const std::string                       matrixFile
  ,const int                              numRhs
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *amesos2LOWSFPL
  ,Teuchos::FancyOStream                  *out_arg
  ,const Teuchos::RCP<const Teuchos::Comm<int> >& comm
  )
{
  using Teuchos::rcp;
  using Teuchos::OSTab;
  bool result, success = true;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(out_arg,false);

  try {

    if(out.get()) {
      *out << "\n***"
           << "\n*** Testing Thyra::BelosLinearOpWithSolveFactory (and Thyra::BelosLinearOpWithSolve)"
           << "\n***\n"
           << "\nEchoing input options:"
           << "\n  matrixFile             = " << matrixFile
           << "\n  numRhs                 = " << numRhs
           << "\n  numRandomVectors       = " << numRandomVectors
           << "\n  maxFwdError            = " << maxFwdError
           << "\n  maxResid               = " << maxResid
           << "\n  showAllTests           = " << showAllTests
           << "\n  dumpAll                = " << dumpAll
           << std::endl;
    }

    if(out.get()) *out << "\nA) Reading in a tpetra matrix A from the file \'"<<matrixFile<<"\' ...\n";

    using Scalar = double;
    using LOWS = Thyra::Amesos2LinearOpWithSolve<Scalar>;
    using MAT = typename LOWS::MAT;
    using LOWSF = Thyra::Amesos2LinearOpWithSolveFactory<Scalar>;
    auto A_tpetra = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(matrixFile, comm);
    auto domain_thyra = Thyra::tpetraVectorSpace<Scalar>(A_tpetra->getDomainMap());
    auto range_thyra = Thyra::tpetraVectorSpace<Scalar>(A_tpetra->getRangeMap());
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> domain_thyra_const = domain_thyra;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> range_thyra_const = range_thyra;
    Teuchos::RCP<const Tpetra::Operator<Scalar>> A_tpetra_const = A_tpetra;

    auto A_thyra = Thyra::constTpetraLinearOp(range_thyra_const, domain_thyra_const, A_tpetra_const);

    if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A_thyra,Teuchos::VERB_EXTREME);

    if(out.get()) *out << "\nB) Creating an Amesos2LinearOpWithSolveFactory object opFactory ...\n";

    Teuchos::RCP<LinearOpWithSolveFactoryBase<Scalar> >
      lowsFactory;
    {
      Teuchos::RCP<LOWSF>
        amesos2LowsFactory = Teuchos::rcp(new LOWSF());
      lowsFactory = amesos2LowsFactory;
    }

    if(out.get()) {
      *out << "\nlowsFactory.getValidParameters():\n";
      lowsFactory->getValidParameters()->print(OSTab(out).o(),0,true,false);
      *out << "\namesos2LOWSFPL before setting parameters:\n";
      amesos2LOWSFPL->print(OSTab(out).o(),0,true);
    }

    lowsFactory->setParameterList(Teuchos::rcp(amesos2LOWSFPL,false));

    if(out.get()) {
      *out << "\namesos2LOWSFPL after setting parameters:\n";
      amesos2LOWSFPL->print(OSTab(out).o(),0,true);
    }

    if(out.get()) *out << "\nC) Creating a Amesos2LinearOpWithSolve object nsA from A ...\n";

    Teuchos::RCP<LinearOpWithSolveBase<Scalar> > nsA = lowsFactory->createOp();
    Thyra::initializeOp<Scalar>(*lowsFactory, A_thyra, nsA.ptr());

    if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

    LinearOpTester<Scalar> linearOpTester;
    linearOpTester.check_adjoint(false);
    linearOpTester.num_rhs(numRhs);
    linearOpTester.num_random_vectors(numRandomVectors);
    linearOpTester.set_all_error_tol(maxFwdError);
    linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
    linearOpTester.show_all_tests(showAllTests);
    linearOpTester.dump_all(dumpAll);
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpTester.check(*nsA,Teuchos::Ptr<Teuchos::FancyOStream>(out.get()));
    if(!result) success = false;

    if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";

    LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
    linearOpWithSolveTester.num_rhs(numRhs);
    linearOpWithSolveTester.turn_off_all_tests();
    linearOpWithSolveTester.check_forward_default(true);
    linearOpWithSolveTester.check_forward_residual(true);
    linearOpWithSolveTester.check_adjoint_default(false);
    linearOpWithSolveTester.check_adjoint_residual(false);
    linearOpWithSolveTester.set_all_solve_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_error_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_warning_tol(1e+1*maxResid);
    linearOpWithSolveTester.forward_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.forward_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.adjoint_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.adjoint_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.show_all_tests(showAllTests);
    linearOpWithSolveTester.dump_all(dumpAll);
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(out.get()) {
      *out << "\namesos2LOWSFPL after solving:\n";
      amesos2LOWSFPL->print(OSTab(out).o(),0,true);
    }

  }
  catch( const std::exception &excpt ) {
    if(out.get()) *out << std::flush;
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }

  return success;

}
