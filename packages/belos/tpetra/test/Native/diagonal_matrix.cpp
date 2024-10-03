// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "KokkosBlas1_mult.hpp"
#include "Teuchos_ParameterList.hpp"
#include <iostream>

namespace { // (anonymous)

struct CommandLineOptions {
  std::string solverName {"TPETRA CG"};
  // mfh 14 Aug 2018: Most of these CG solvers take 29 or 30 iterations
  // on this problem.  We add five iterations to allow for rounding
  // error.
  int maxAllowedNumIters {35};
  bool verbose {true};
};
CommandLineOptions commandLineOptions;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions (true);
  clp.setOption ("solver", &commandLineOptions.solverName,
		 "Name of the solver to test.  Belos::SolverFactory::create "
		 "must accept this string.  Protect with double quotes if it "
		 "has spaces: e.g., \"TPETRA CG PIPELINE\".");
  clp.setOption ("maxNumIters", &commandLineOptions.maxAllowedNumIters,
		 "Maximum number of iterations");
  clp.setOption ("verbose", "quiet", &commandLineOptions.verbose,
		 "Whether to print verbose output");
}

template<class TpetraOperatorType = Tpetra::Operator<>>
class TestDiagonalOperator : public TpetraOperatorType {
public:
  using multivector_type = Tpetra::MultiVector<
    typename TpetraOperatorType::scalar_type,
    typename TpetraOperatorType::local_ordinal_type,
    typename TpetraOperatorType::global_ordinal_type,
    typename TpetraOperatorType::node_type>;
  using map_type = Tpetra::Map<
    typename TpetraOperatorType::local_ordinal_type,
    typename TpetraOperatorType::global_ordinal_type,
    typename TpetraOperatorType::node_type>;

private:
  using device_type = typename map_type::device_type;
  using val_type = typename TpetraOperatorType::scalar_type;
  using STS = Teuchos::ScalarTraits<val_type>;

public:
  using mag_type = typename Teuchos::ScalarTraits<val_type>::magnitudeType;

  TestDiagonalOperator (const Teuchos::RCP<const map_type>& map,
                        const mag_type minSingularValue) :
    map_ (map)
  {
    using dev_view_type = Kokkos::View<mag_type*, device_type>;
    using host_view_type = typename dev_view_type::HostMirror;
    using SC = typename TpetraOperatorType::scalar_type;
    using LO = typename TpetraOperatorType::local_ordinal_type;

    const LO lclNumRows = map_.is_null () ?
      static_cast<LO> (0) :
      static_cast<LO> (map_->getLocalNumElements ());
    dev_view_type diag_d ("diag", lclNumRows);
    host_view_type diag_h = Kokkos::create_mirror_view (diag_d);

    if (lclNumRows != 0) {
      const SC ONE = STS::one ();
      const SC exponent = ONE / static_cast<mag_type> (lclNumRows - 1);
      const SC scalingFactor = ::pow (minSingularValue, exponent);
      diag_h(0) = ONE;
      for (LO lclRow = 1; lclRow < lclNumRows; ++lclRow) {
        diag_h(lclRow) = diag_h(lclRow-1) * scalingFactor;
      }
    }
    Kokkos::deep_copy (diag_d, diag_h);
    diag_ = diag_d;
  }

  Teuchos::RCP<const map_type> getDomainMap () const {
    return map_;
  }

  Teuchos::RCP<const map_type> getRangeMap () const {
    return map_;
  }

  void
  apply (const multivector_type& X,
         multivector_type& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         val_type alpha = Teuchos::ScalarTraits<val_type>::one (),
         val_type beta = Teuchos::ScalarTraits<val_type>::zero ()) const
  {
    using ISC = typename multivector_type::impl_scalar_type;
    TEUCHOS_TEST_FOR_EXCEPTION
      (mode == Teuchos::CONJ_TRANS && STS::isComplex,
       std::logic_error, "Conjugate transpose case not implemented.");

    const size_t numVecs = X.getNumVectors ();
    for (size_t j = 0; j < numVecs; ++j) {
      auto X_j = X.getVector (j);
      auto Y_j = Y.getVectorNonConst (j);

      auto X_j_lcl_2d = X_j->getLocalViewDevice (Tpetra::Access::ReadOnly);
      auto X_j_lcl = Kokkos::subview (X_j_lcl_2d, Kokkos::ALL (), 0);
      auto Y_j_lcl_2d = Y_j->getLocalViewDevice (Tpetra::Access::ReadWrite);
      auto Y_j_lcl = Kokkos::subview (Y_j_lcl_2d, Kokkos::ALL (), 0);

      KokkosBlas::mult (static_cast<ISC> (beta), Y_j_lcl,
                        static_cast<ISC> (alpha), X_j_lcl, diag_);
    }
  }

private:
  Teuchos::RCP<const map_type> map_;
  Kokkos::View<const mag_type*, device_type> diag_;
};

void
testSolver (Teuchos::FancyOStream& out,
	    bool& success,
	    const std::string& solverName,
	    const int maxAllowedNumIters,
	    const bool verbose)
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::endl;
  using map_type = Tpetra::Map<>;
  using MV = Tpetra::MultiVector<>;
  using OP = Tpetra::Operator<>;
  using SC = MV::scalar_type;
  using GO = map_type::global_ordinal_type;
  using mag_type = MV::mag_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using STM = Teuchos::ScalarTraits<mag_type>;
  // The Teuchos unit test framework likes to capture output to 'out',
  // and not print anything until the test is done.  This can hinder
  // debugging.  If the test crashes without useful output, try
  // setting this to 'true'.  That will change 'myOut' from an alias
  // to 'out', into a wrapper for std::cerr.
  constexpr bool debug = false;

  const SC ZERO = STS::zero ();
  const SC ONE = STS::one ();

  RCP<FancyOStream> myOutPtr =
    debug ? getFancyOStream (rcpFromRef (std::cerr)) : rcpFromRef (out);
  Teuchos::FancyOStream& myOut = *myOutPtr;

  myOut << "Test \"native\" Tpetra version of solver \"" << solverName << "\""
	<< endl;
  Teuchos::OSTab tab1 (out);

  myOut << "Create the linear system to solve" << endl;

  auto comm = Tpetra::getDefaultComm ();
  const GO gblNumRows = 10000;
  const GO indexBase = 0;
  RCP<const map_type> map (new map_type (gblNumRows, indexBase, comm));
  const mag_type minSingVal = 0.1;
  RCP<TestDiagonalOperator<OP> > A =
    rcp (new TestDiagonalOperator<> (map, minSingVal));

  MV X (A->getDomainMap (), 1);
  MV B (A->getRangeMap (), 1);
  B.randomize ();
  // Get ready for next solve by resetting initial guess to zero.
  X.putScalar (ZERO);

  myOut << "Create solver instance using Belos::SolverFactory" << endl;

  RCP<Belos::SolverManager<SC, MV, OP> > solver;
  try {
    Belos::SolverFactory<SC, MV, OP> factory;
    solver = factory.create (solverName, Teuchos::null);
  }
  catch (std::exception& e) {
    myOut << "*** FAILED: Belos::SolverFactory::create threw an exception: "
        << e.what () << endl;
    success = false;
    return;
  }

  TEST_ASSERT( solver.get () != nullptr );
  if (solver.get () == nullptr) {
    myOut << "Belos::SolverFactory returned a null solver." << endl;
    return;
  }

  myOut << "Set parameters" << endl;
  RCP<ParameterList> params = parameterList ("Belos");
  params->set ("Verbosity", verbose ? 1 : 0);
  try {
    solver->setParameters (params);
  }
  catch (std::exception& e) {
    myOut << "*** FAILED: setParameters threw an exception: "
        << e.what () << endl;
    success = false;
    return;
  }
  catch (...) {
    myOut << "*** FAILED: setParameters threw an exception "
      "not a subclass of std::exception." << endl;
    success = false;
    return;
  }

  myOut << "Set up the linear system to solve" << endl;
  auto lp = rcp (new Belos::LinearProblem<SC, MV, OP> (A, rcpFromRef (X),
                                                       rcpFromRef (B)));
  lp->setProblem ();

  myOut << "Solve the linear system" << endl;
  solver->setProblem (lp);
  const Belos::ReturnType belosResult = solver->solve ();

  myOut << "Belos solver wrapper result: "
	<< (belosResult == Belos::Converged ? "Converged" : "Unconverged")
	<< endl
	<< "Number of iterations: " << solver->getNumIters ()
	<< endl;

  TEST_ASSERT( solver->getNumIters () <= maxAllowedNumIters );
  TEST_ASSERT( belosResult == Belos::Converged );

  myOut << "Check the explicit residual norm(s)" << endl;

  // Get the tolerance that the solver actually used.
  const mag_type tol = [&] () {
      const char tolParamName[] = "Convergence Tolerance";
      auto pl = solver->getCurrentParameters ();
      if (! pl->isType<mag_type> (tolParamName)) {
	pl = solver->getValidParameters ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (! pl->isType<mag_type> (tolParamName), std::logic_error,
       "Solver lacks \"" << tolParamName << "\" parameter, in either "
	 "getCurrentParameters() or getValidParameters().");
      return pl->get<mag_type> (tolParamName);
    } ();

  MV X_copy (X, Teuchos::Copy);
  MV R (B.getMap (), B.getNumVectors ());
  A->apply (X_copy, R);
  R.update (ONE, B, -ONE);
  Teuchos::Array<mag_type> R_norms (R.getNumVectors ());
  R.norm2 (R_norms ());
  Teuchos::Array<mag_type> B_norms (B.getNumVectors ());
  B.norm2 (B_norms ());

  for (size_t j = 0; j < R.getNumVectors (); ++j) {
    const mag_type relResNorm = (B_norms[j] == STM::zero ()) ?
      R_norms[j] :
      R_norms[j] / B_norms[j];
    myOut << "Column " << (j+1) << " of " << R.getNumVectors ()
	  << ": Absolute residual norm: " << R_norms[j]
	  << ", Relative residual norm: " << relResNorm
	  << endl;
    TEST_ASSERT( relResNorm <= tol );
  }
  myOut << endl;
}

TEUCHOS_UNIT_TEST( TpetraNativeSolvers, Diagonal )
{
  testSolver (out, success, commandLineOptions.solverName,
	      commandLineOptions.maxAllowedNumIters,
	      commandLineOptions.verbose);
}

} // namespace (anonymous)

int main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
}
