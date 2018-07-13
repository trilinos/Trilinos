// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Amesos2_Factory.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
// Define typedefs and macros for testing over all Tpetra types.
// They work whether or not ETI is enabled.
#include "TpetraCore_ETIHelperMacros.h"

// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
namespace Amesos2 {
namespace Details {
  // FIXME (mfh 21 Aug 2015) NONE of the commented-out things work.
  //
  // extern void __attribute__((weak)) registerLinearSolverFactory ();
  // void __attribute__((weak)) registerLinearSolverFactory ();
  // #pragma weak registerLinearSolverLibrary

  extern void registerLinearSolverFactory ();

} // namespace Details
} // namespace Amesos2
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

namespace {
// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
  TEUCHOS_STATIC_SETUP()
  {
    // if (Amesos2::Details::registerLinearSolverFactory == NULL) {
    //   std::cout << "-- Amesos2::Details::registerLinearSolverFactory is NULL" << std::endl;
    // } else {
    //   Amesos2::Details::registerLinearSolverFactory ();
    // }

    Amesos2::Details::registerLinearSolverFactory ();
  }
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

  // Create a very simple square test matrix.  We use the identity
  // matrix here.  The point of this test is NOT to exercise the
  // solver; it's just to check that its LinearSolverFactory can
  // create working solvers.  Amesos2 has more rigorous tests for each
  // of its solvers.
  template<class SC, class LO, class GO, class NT>
  Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
  createTestMatrix (Teuchos::FancyOStream& out,
                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Tpetra::global_size_t gblNumRows)
  {
    using Teuchos::Comm;
    using Teuchos::tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
    typedef Tpetra::Map<LO,GO,NT> map_type;
    //typedef Teuchos::ScalarTraits<SC> STS; // unused

    Teuchos::OSTab tab0 (out);
    out << "Create test matrix with " << gblNumRows << " row(s)" << endl;
    Teuchos::OSTab tab1 (out);

    TEUCHOS_TEST_FOR_EXCEPTION
      ( gblNumRows == 0, std::invalid_argument, "gblNumRows = 0");

    const GO indexBase = 0;
    RCP<const map_type> rowMap (new map_type (gblNumRows, indexBase, comm));
    // The simple matrix with just row and column map the same and a diagonal
    // causes trouble for lot of reasons. One example being the parmetis call
    // from Superlu_dist which will fail with a null array, and a diagonal
    // matrix results in a null array as parmetis doesn't like self edges.
    // We will use a simple tridiagonal matrix instead as our test. -Siva
    RCP<const map_type> colMap = rowMap;
    const size_t maxNumEntPerRow = 3;
    const SC two = static_cast<SC> (2.0);
    const SC minusOne = static_cast<SC> (-1.0);
    RCP<MAT> A (new MAT (rowMap, maxNumEntPerRow, Tpetra::DynamicProfile));

    if (rowMap->getNodeNumElements () != 0) {
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement(lclRow);
        if (gblRow == rowMap->getMinAllGlobalIndex())
        {
          A->insertGlobalValues (gblRow,
               tuple<GO> (gblRow, gblRow + 1),
               tuple<SC> (two, minusOne));
        }
        else if ( gblRow == rowMap->getMaxAllGlobalIndex())
        {
          A->insertGlobalValues (gblRow,
               tuple<GO> (gblRow - 1, gblRow),
               tuple<SC> (minusOne, two));
        }
        else
        {
          A->insertGlobalValues (gblRow,
               tuple<GO> (gblRow - 1, gblRow, gblRow+1),
               tuple<SC> (minusOne, two, minusOne));
        }
      }
    }

    RCP<const map_type> domMap = rowMap;
    RCP<const map_type> ranMap = rowMap;
    A->fillComplete (domMap, ranMap);
    return A;
  }

  // Create a very simple square test linear system (matrix,
  // right-hand side(s), and exact solution(s).  We use the tridiagonal
  // matrix here.  The point of this test is NOT to exercise the
  // solver; it's just to check that its LinearSolverFactory can
  // create working solvers.  Amesos2 has more rigorous tests for each
  // of its solvers.
  template<class SC, class LO, class GO, class NT>
  void
  createTestProblem (Teuchos::FancyOStream& out,
                     Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> >& X,
                     Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
                     Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> >& B,
                     const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t gblNumRows,
                     const size_t numVecs)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;

    A = createTestMatrix<SC, LO, GO, NT> (out, comm, gblNumRows);
    X = rcp (new MV (A->getDomainMap (), numVecs));
    B = rcp (new MV (A->getRangeMap (), numVecs));

    X->putScalar (STS::one ());
    A->apply (*X, *B);
  }

  template<class SC, class LO, class GO, class NT>
  void
  testSolver (Teuchos::FancyOStream& out,
              bool& success,
              Tpetra::MultiVector<SC, LO, GO, NT>& X,
              const Teuchos::RCP<const Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
              const Tpetra::MultiVector<SC, LO, GO, NT>& B,
              const Tpetra::MultiVector<SC, LO, GO, NT>& X_exact,
              const std::string& solverName)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::Operator<SC,LO,GO,NT> OP;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;
    typedef typename MV::mag_type mag_type;
    typedef Teuchos::ScalarTraits<mag_type> STM;

    Teuchos::OSTab tab0 (out);
    out << "Test solver \"" << solverName << "\" from Amesos2 package" << endl;
    Teuchos::OSTab tab1 (out);

    // NDE: Beginning changes towards passing parameter list to shylu basker
    // for controlling various parameters per test, matrix, etc.

    Teuchos::ParameterList amesos2_paramlist;
    if ( solverName == "ShyLUBasker" || solverName == "shylubasker" ) {
      amesos2_paramlist.setName("Amesos2");
      Teuchos::ParameterList & shylubasker_paramlist = amesos2_paramlist.sublist("Basker");

      shylubasker_paramlist.set("num_threads", 1,
          "Number of threads");
      shylubasker_paramlist.set("pivot", false,
          "Should not pivot");
      shylubasker_paramlist.set("pivot_tol", .0001,
          "Tolerance before pivot, currently not used");
      shylubasker_paramlist.set("symmetric", false,
          "Should Symbolic assume symmetric nonzero pattern");
      shylubasker_paramlist.set("realloc" , false,
          "Should realloc space if not enough");
      shylubasker_paramlist.set("verbose", false,
          "Information about factoring");
      shylubasker_paramlist.set("verbose_matrix", false,
          "Give Permuted Matrices");
      shylubasker_paramlist.set("matching", true,
          "Use WC matching (Not Supported)");
      shylubasker_paramlist.set("matching_type", 0,
          "Type of WC matching (Not Supported)");
      shylubasker_paramlist.set("btf", true,
          "Use BTF ordering");
      shylubasker_paramlist.set("amd_btf", true,
          "Use AMD on BTF blocks (Not Supported)");
      shylubasker_paramlist.set("amd_dom", true,
          "Use CAMD on ND blocks (Not Supported)");
      shylubasker_paramlist.set("transpose", false,
          "Solve the transpose A");
    }

    RCP<Trilinos::Details::LinearSolver<MV, OP, mag_type> > solver;
    try {
      solver = Trilinos::Details::getLinearSolver<MV, OP, mag_type> ("Amesos2", solverName);
      solver->setParameters(Teuchos::rcpFromRef(amesos2_paramlist));

      //set parameter list here as well
    } catch (std::exception& e) {
      out << "*** FAILED: getLinearSolver threw an exception: " << e.what () << endl;
      success = false;
      return;
    }
    TEST_EQUALITY_CONST( solver.is_null (), false );
    if (solver.is_null ()) {
      out << "*** FAILED to create solver \"" << solverName
          << "\" from Amesos2 package" << endl;
      success = false;
      return;
    }

    out << "Set matrix" << endl;
    solver->setMatrix (A);

    out << "Compute symbolic and numeric factorization" << endl;
    solver->symbolic ();
    solver->numeric ();

    out << "Solve AX=B for X" << endl;
    X.putScalar (STS::zero ());
    solver->solve (X, B);

    out << "Check accuracy of result, using known solution" << endl;
    const size_t numVecs = B.getNumVectors ();
    Teuchos::Array<mag_type> X_norms (numVecs);
    Teuchos::Array<mag_type> X_exact_norms (numVecs);
    X.norm2 (X_norms ());
    X_exact.norm2 (X_exact_norms ());

    const mag_type tol = STM::squareroot (STM::eps ());
    TEST_COMPARE_FLOATING_ARRAYS( X_norms, X_exact_norms, tol );
  }

  //
  // The actual unit test.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SolverFactory, Solve, SC, LO, GO, NT )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::TypeNameTraits;
    using std::endl;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;

#if ! defined(TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION)
    out << "LinearSolverFactory run-time registration disabled; "
      "not running test" << endl;
    return;
#endif // NOT TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const Tpetra::global_size_t gblNumRows = comm->getSize () * 10;
    const size_t numVecs = 3;

    out << "Create test problem" << endl;
    RCP<MV> X_exact, B;
    RCP<MAT> A;
    createTestProblem (out, X_exact, A, B, comm, gblNumRows, numVecs);

    RCP<MV> X = rcp (new MV (X_exact->getMap (), numVecs));

    //X_exact->describe(out, Teuchos::VERB_EXTREME);
    //B->describe(out, Teuchos::VERB_EXTREME);


    const int numSolvers = 10;
    const char* solverNames[10] = {"shylubasker", "basker", "klu2", "superlu_dist",
                                  "superlu_mt", "superlu", "pardiso_mkl",
                                  "lapack", "mumps", "amesos2_cholmod"};
    // The number of solvers that Amesos2::create actually supports,
    // for the current MV and MAT types.  If it doesn't support _any_
    // of the solvers, we consider this test to have failed.
    int numSolversSupported = 0;
    for (int k = 0; k < numSolvers; ++k) {
      const std::string solverName (solverNames[k]);

      // Don't actually test the solver unless it is supported.
      if (Amesos2::query (solverName)) {
        // Unfortunately, Amesos2::query() doesn't tell us whether the
        // solver supports the given combination of template
        // parameters.  However, we can use Amesos2::create to tell us
        // this.  If it throws or returns null, the combination is
        // (presumably) not supported.
        bool skip = false;
        RCP<Amesos2::Solver<MAT, MV> > solver;
        try {
          solver = Amesos2::create<MAT, MV> (solverName, A, X, B);

          // NDE: Beginning changes towards passing parameter list to shylu basker
          // for controlling various parameters per test, matrix, etc.

          Teuchos::ParameterList amesos2_paramlist;
          if ( solverName == "shylubasker" ) {
            amesos2_paramlist.setName("Amesos2");
            Teuchos::ParameterList & shylubasker_paramlist = amesos2_paramlist.sublist("Basker");

            shylubasker_paramlist.set("num_threads", 1,
                "Number of threads");
            shylubasker_paramlist.set("pivot", false,
                "Should not pivot");
            shylubasker_paramlist.set("pivot_tol", .0001,
                "Tolerance before pivot, currently not used");
            shylubasker_paramlist.set("symmetric", false,
                "Should Symbolic assume symmetric nonzero pattern");
            shylubasker_paramlist.set("realloc" , false,
                "Should realloc space if not enough");
            shylubasker_paramlist.set("verbose", false,
                "Information about factoring");
            shylubasker_paramlist.set("verbose_matrix", false,
                "Give Permuted Matrices");
            shylubasker_paramlist.set("matching", true,
                "Use WC matching (Not Supported)");
            shylubasker_paramlist.set("matching_type", 0,
                "Type of WC matching (Not Supported)");
            shylubasker_paramlist.set("btf", true,
                "Use BTF ordering");
            shylubasker_paramlist.set("amd_btf", true,
                "Use AMD on BTF blocks (Not Supported)");
            shylubasker_paramlist.set("amd_dom", true,
                "Use CAMD on ND blocks (Not Supported)");
            shylubasker_paramlist.set("transpose", false,
                "Solve the transpose A");

            solver->setParameters(Teuchos::rcpFromRef(amesos2_paramlist));
          }
        }
        catch (...) {
          out << "Amesos2::create threw an exception for solverName = \"" <<
            solverName << "\", MAT = " << TypeNameTraits<MAT>::name () <<
            ", and MV = " << TypeNameTraits<MV>::name () << ".  As a result, "
            "we will skip attempting to create this solver using Trilinos::"
            "Details::getLinearSolver." << endl;
          skip = true;
        }
        if (solver.is_null ()) {
          out << "Amesos2::create returned null for solverName = \"" <<
            solverName << "\", MAT = " << TypeNameTraits<MAT>::name () <<
            ", and MV = " << TypeNameTraits<MV>::name () << ".  As a result, "
            "we will skip attempting to create this solver using Trilinos::"
            "Details::getLinearSolver." << endl;
          skip = true;
        }
        if (! skip) {
          ++numSolversSupported;
          bool testSuccess = true;
          testSolver<SC, LO, GO, NT> (out, testSuccess, *X, A, *B, *X_exact, solverName);
          success = success && testSuccess;
        }
      }
    }

    if (numSolversSupported == 0) {
      success = false;
      out << "*** Amesos2::create doesn't actually support any solvers for "
          << "SC = " << TypeNameTraits<SC>::name ()
          << " , LO = " << TypeNameTraits<LO>::name ()
          << ", GO = " << TypeNameTraits<GO>::name ()
          << ", and NT = " << TypeNameTraits<NT>::name () << "." << endl;
    }
    TEST_ASSERT ( success == true);
  }

  // Define typedefs that make the Tpetra macros work.
  TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that instantiates the unit test
#define LCLINST( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( SolverFactory, Solve, SC, LO, GO, NT )

  // Instantiate the unit test.
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

} // namespace (anonymous)
