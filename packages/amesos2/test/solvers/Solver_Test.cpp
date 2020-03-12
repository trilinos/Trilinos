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

/**
 * \file   Solver_Test.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Wed May 25 12:17:25 2011
 *
 * \brief  Tests Amesos2 solver interfaces using various matrix/vector
 *         objects, scalar/ordinal types, and input matrices.  Test
 *         parameters are specified by an input XML file.
 */

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp> // For reading matrix-market files

#include "Amesos2.hpp"          // includes everything from Amesos2

#include "KokkosBlas.hpp"

// #ifdef HAVE_TPETRA_INST_INT_INT
#if defined(HAVE_AMESOS2_EPETRA) && defined(HAVE_AMESOS2_EPETRAEXT)
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_CrsMatrixIn.h>
#endif  // HAVE_AMESOS2_EPETRAEXT
//#endif

using std::string;

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ptrInArg;
using Teuchos::outArg;
using Teuchos::ETransp;
using Teuchos::TRANS;
using Teuchos::NO_TRANS;
using Teuchos::CONJ_TRANS;
using Teuchos::ParameterList;
using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::Array;
using Teuchos::ArrayView;

/*
 * An example input xml file can be found in the default: "solvers_test.xml"
 */

// TODO: flush out the timers
RCP<Time> total_timer = TimeMonitor::getNewTimer("total time");
RCP<Teuchos::FancyOStream> fos; // for global output

// For global output during solution comparison
RCP<Teuchos::FancyOStream> compare_fos;

int verbosity = 0;
std::string filedir = "../matrices/";

/**
 * If \c true, then the solution routine will perform multiple solves
 * using a matrix.
 */
bool multiple_solves = true;

/**
 * If \c true, then the solution routine will at some point reperform
 * numeric factorization with a numerically different, but
 * symbolically same, matrix.
 */
bool refactor = false;

/*
 * Takes the given parameter list and performs the test solves that it describes.
 *
 * The top-level list is expected to describe the parameters for a single
 *
 * \return whether all tests for this matrix passed
 */
bool do_mat_test(const ParameterList& parameters);

bool test_mat_with_solver(const string& mm_file,
                          const string& solver_name,
                          const ParameterList& test_params,
                          ParameterList solve_params);

#ifdef HAVE_TPETRA_INST_INT_INT
#if defined(HAVE_AMESOS2_EPETRA) && defined(HAVE_AMESOS2_EPETRAEXT)
/*
 * Tests a matrix solve with the given solver on the matrix found in
 * the named file.  EPETRA_RUNS is a parameter list that describes the
 * runs that should be performed for the named solver with the Epetra
 * linear algebra objects.  As opposed to the tpetra tests, this
 * parameter list cannot let you specify ordinal or scalar types
 * (because this is not supported by epetra), but does let you specify
 * run-specific solver parameters.  If a non-list parameter entry is
 * found in EPETRA_RUNS, then a single test run with the default
 * parameters will be executed following all other epetra test runs.
 *
 * Instantiates the matrix as an Epetra_CrsMatrix, so types will be <int,int,double>
 *
 * This test is run for every non-complex matrix, with the idea that
 * all solvers will support the <int,int,double> types.
 *
 * Note: the input matrix must be square!
 *
 * \return Whether all tests for this matrix with epetra objects passed
 */
bool test_epetra(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& epetra_runs,
                 ParameterList solve_params);
#endif
#endif //HAVE_TPETRA_INST_INT_INT

/*
 * Tests a matrix solve with the matrix found in the named file using
 * the solver named in SOLVER_NAME.  TPETRA_RUNS is a parameter list
 * describing what test runs should be performed for tpetra objects.
 * Each run can describe the ordinal and scalar types that should be
 * tested, as well as any solver parameters specific for that run.  An
 * example can be found in \c solvers_test.xml.
 *
 * Instantiates the matrix as an Tpetra::CrsMatrix with the given template types.
 *
 * \return Whether all tests for this matrix with tpetra objects passed
*/
bool test_tpetra(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& tpetra_runs,
                 ParameterList solve_params);

/*
 * Will use Kokkos CrsMatrix so that any memory space can be applied as the
 * source memory space.
 *
 * \return Whether all tests for this matrix with kokkos objects passed
*/
bool test_kokkos(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& kokkos_runs,
                 ParameterList solve_params);

typedef Tpetra::Map<>::node_type DefaultNode;

int main(int argc, char*argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  TimeMonitor TotalTimer(*total_timer);

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int root = 0;

  string xml_file("solvers_test.xml"); // default xml file
  string src_memory_space_name("Undefined"); // default src memory space (no special testing)
  bool allprint = false;
  Teuchos::CommandLineProcessor cmdp;
  cmdp.setDocString("A test driver for Amesos2 solvers.  It reads parameters\n"
                    "from a given (or default) xml file which describes:\n"
                    " - What test matrices to use, for each\n"
                    " - Which solvers to use, for each\n"
                    " - Which linear alg objects (currently Tpetra and Epetra), and finally, for each\n"
                    " - a list of test runs, with optional parameters\n"
                    "An example can be found in `solvers_test.xml'");

  cmdp.setOption("xml-params", &xml_file, "XML Parameters file");
  cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  cmdp.setOption("filedir", &filedir, "Directory to search for matrix files");
  cmdp.setOption("verbosity", &verbosity, "Set verbosity level of output");
  cmdp.setOption("multiple-solves","single-solve", &multiple_solves, "Perform multiple solves with different RHS arguments");
  cmdp.setOption("refactor","no-refactor", &refactor, "Recompute L and U using a numerically different matrix at some point");
  try{
    cmdp.parse(argc,argv);
  } catch (const Teuchos::CommandLineProcessor::HelpPrinted& hp) {
    return EXIT_SUCCESS;        // help was printed, exit gracefully.
  }

  // set up output streams based on command-line parameters
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if( !allprint ) fos->setOutputToRootOnly( root );

  Teuchos::oblackholestream blackhole;
  if( verbosity > 3 ){
    compare_fos = fos;
  } else {
    compare_fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(blackhole));
  }

  //Read the contents of the xml file into a ParameterList.
  if( verbosity > 0 ){
    *fos << "Every proc reading parameters from xml_file: "
         << xml_file << std::endl;
  }
  ParameterList test_params =
    Teuchos::ParameterXMLFileReader(xml_file).getParameters();

  // Check the parameterlist for the presence of any of the other params
  if( test_params.isParameter("all-print") ){
    allprint = test_params.get<bool>("all-print");
  }
  if( test_params.isParameter("filedir") ){
    filedir = test_params.get<string>("filedir");
  }
  if( test_params.isParameter("verbosity") ){
    verbosity = test_params.get<int>("verbosity");
  }


  // Go through the input parameters and execute tests accordingly.
  bool success = true;
  ParameterList::ConstIterator mat_it;
  for( mat_it = test_params.begin(); mat_it != test_params.end(); ++mat_it ){
    if( test_params.entry(mat_it).isList() ){ // each matrix entry must be a list
      success &= do_mat_test(Teuchos::getValue<ParameterList>(test_params.entry(mat_it)));
    } else {
      *fos << "unexpected non-list entry in xml input, ignoring..." << std::endl;
    }
  }

  // The summary table is very verbose
  if( verbosity > 3 ){
    TimeMonitor::summarize();
  }

  // This output is used to indicate a passed test, the test framework
  // will parse for it.
  *fos << std::endl << "End Result: ";
  if( success ){
    *fos << "TEST PASSED";
  } else {
    *fos << "TEST FAILED";
  }
  *fos << std::endl;

  return 0;
}

bool do_mat_test(const ParameterList& parameters)
{
  bool success = true;

  const string mm_file = parameters.name();
  if( verbosity > 0 ){
    *fos << "Test matrix " << mm_file << " ... " << std::flush;
    if( verbosity > 1) *fos << std::endl;
  }

  // This parameter is ignored for now
  bool complex = false;
  if( parameters.isParameter("complex") ){
    if( ! parameters.isType<bool>("complex") ){
      *fos << "'complex' parameter type should be bool! ignoring..." << std::endl;
    } else {
      complex = parameters.get<bool>("complex");
    }
  }
  (void) complex; // forestall warning for set but unused variable

  ParameterList solve_params("Amesos2");
  if( parameters.isSublist("all_solver_params") ){
    solve_params.setParameters( parameters.sublist("all_solver_params") );
  }

  ParameterList::ConstIterator solver_it;
  for( solver_it = parameters.begin(); solver_it != parameters.end(); ++solver_it ){
    if( ! parameters.entry(solver_it).isUsed() ){
      string solver_name = parameters.name(solver_it);

      if ( ! parameters.entry(solver_it).isList() ){
        *fos << "no support for default test runs yet.  "
             << "Ignoring the " << solver_name << " solver..."
             << std::endl;
      } else {
        if( Amesos2::query(solver_name) ){
          // then we have support for this solver

          if( verbosity > 1) *fos << "  | with " << solver_name << " : " << std::endl;

          ParameterList test_params = Teuchos::getValue<ParameterList>(parameters.entry(solver_it));
          bool solver_success = test_mat_with_solver(mm_file, solver_name, test_params, solve_params);

          if( verbosity > 1 ){
            *fos << "  ";
            if( solver_success ) *fos << "+ ";
            else *fos << "- ";
            *fos << "Testing with " << solver_name << " ";
            if( solver_success ) *fos << "passed";
            else *fos << "failed!";
            *fos << std::endl;
          }

          success &= solver_success;
        } else {
          if( verbosity > 1 ) *fos << "  - " << solver_name << " not enabled, skipping..." << std::endl;
        }
      }
    }
  }

  if( verbosity > 0 ){
    if( verbosity > 1 ){
      *fos << mm_file << " ";
    }
    if( !success ) *fos << "failed!";
    else *fos << "passed";

    *fos << std::endl;
  }

  return( success );
}

bool
test_mat_with_solver (const string& mm_file,
                      const string& solver_name,
                      const ParameterList& test_params,
                      ParameterList solve_params)
{
  using std::cerr;
  using std::endl;
  bool success = true;

  if (test_params.isSublist("solver_params")) {
    solve_params.sublist (solver_name) = test_params.sublist ("solver_params");
  }

  ParameterList::ConstIterator object_it;
  for (object_it = test_params.begin(); object_it != test_params.end(); ++object_it) {
    if (! test_params.entry (object_it).isUsed ()) {
      const string object_name = test_params.name (object_it);

      // There may be a more flexible way (e.g. a `query' function) to do this check
      if (object_name == "epetra") {
        if (verbosity > 1) {
          *fos << "    Testing Epetra objects" << endl;
        }
#ifdef HAVE_TPETRA_INST_INT_INT
#if defined(HAVE_AMESOS2_EPETRA) && defined(HAVE_AMESOS2_EPETRAEXT)
        const ParameterList epetra_runs = Teuchos::getValue<ParameterList> (test_params.entry (object_it));
        const bool epetraSuccess = test_epetra (mm_file, solver_name, epetra_runs, solve_params);
        success &= epetraSuccess;
        if (verbosity > 1) {
          *fos << "      - Epetra test " << (epetraSuccess ? "succeeded" : "failed") << endl;
        }
#else
        if (verbosity > 1) {
          *fos << "    EpetraExt must be enabled for testing of Epetra objects.  Skipping this test."
               << endl;
        }
#endif // HAVE_AMESOS2_EPETRAEXT
#endif
      }
      else if (object_name == "kokkos") {
        if (verbosity > 1) {
          *fos << "    Testing Kokkos objects" << endl;
        }
        const ParameterList tpetra_runs = Teuchos::getValue<ParameterList> (test_params.entry (object_it));
        bool kokkosSuccess = test_kokkos(mm_file, solver_name, tpetra_runs, solve_params);

        if (verbosity > 1) {
          *fos << "      - Kokkos test " << (kokkosSuccess ? "succeeded" : "failed") << endl;
        }
        success &= kokkosSuccess;
      }
      else if (object_name == "tpetra") {
        if (verbosity > 1) {
          *fos << "    Testing Tpetra objects" << endl;
        }
        const ParameterList tpetra_runs = Teuchos::getValue<ParameterList> (test_params.entry (object_it));
        const bool tpetraSuccess = test_tpetra (mm_file, solver_name, tpetra_runs, solve_params);
        if (verbosity > 1) {
          *fos << "      - Tpetra test " << (tpetraSuccess ? "succeeded" : "failed") << endl;
        }
        success &= tpetraSuccess;
      }
      else if (object_name == "solver_params") {
        // FIXME (mfh 23 Jan 2014) Skip this sublist.  It's a quirk of
        // the test that the sublist of solver parameters is at the
        // same level as the lists for the different linear algebra
        // types.  We should fix this in the input parameter list
        // (test_params), but for now I'll just skip that sublist.
      }
      else {
        *fos << "    Linear algebra objects of type \"" << object_name << "\" are not supported." << endl
             << "    Here is the test_params parameter list:" << endl
             << test_params << endl
             << "    Here is the solve_params parameter list:" << endl
             << solve_params << endl;

      }
    }
  }

  return success;
}


// The generic version of this struct has an operator which always
// returns false.  Thus, by default, the test always fails (the
// solution is not correct).  Specializations of this class for
// specific (Multi)Vector types actually check the solution and may
// return true or false.
template <class Vector>
struct solution_checker {
  bool operator() (RCP<Vector> true_solution, RCP<Vector> solution) {
    return false;
  }
};

// Partial specialization of solution_checker for Tpetra::MultiVector.
template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct solution_checker<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> t_mv;
  bool operator()(RCP<t_mv> true_solution, RCP<t_mv> given_solution)
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_t;
    size_t num_vecs = true_solution->getNumVectors();

    Teuchos::Array<mag_t> ts_norms(num_vecs), gs_norms(num_vecs);
    true_solution->norm2(ts_norms());
    given_solution->norm2(gs_norms());

    return Teuchos::compareFloatingArrays(ts_norms, "true_solution",
                                          gs_norms, "given_solution",
                                          Teuchos::as<mag_t> (0.005), *compare_fos);
  }
};

template <typename Scalar,
          typename ExecutionSpace>
struct solution_checker<Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> > {
  typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> t_mv;
  bool operator()(RCP<t_mv> true_solution, RCP<t_mv> given_solution)
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_t;
    size_t num_vecs = true_solution->extent(1);
    Teuchos::Array<mag_t> ts_norms(num_vecs), gs_norms(num_vecs);
    for (size_t k = 0; k < num_vecs; ++k) {
      ts_norms[k] = KokkosBlas::nrm2_squared(
        Kokkos::subview(*true_solution, Kokkos::ALL, k));
      gs_norms[k] = KokkosBlas::nrm2_squared(
        Kokkos::subview(*given_solution, Kokkos::ALL, k));
    }

    return Teuchos::compareFloatingArrays(ts_norms, "true_solution",
                                          gs_norms, "given_solution",
                                          Teuchos::as<mag_t> (0.005), *compare_fos);
  }
};


#if defined(HAVE_AMESOS2_EPETRA) && defined(HAVE_AMESOS2_EPETRAEXT)
template <>
struct solution_checker<Epetra_MultiVector> {
  bool operator() (RCP<Epetra_MultiVector> true_solution, RCP<Epetra_MultiVector> given_solution)
  {
    int num_vecs = true_solution->NumVectors();

    Teuchos::Array<double> ts_norms(num_vecs), gs_norms(num_vecs);
    true_solution->Norm2(ts_norms.getRawPtr());
    given_solution->Norm2(gs_norms.getRawPtr());

    return Teuchos::compareFloatingArrays(ts_norms, "true_solution",
                                          gs_norms, "given_solution",
                                          0.005, *compare_fos);
  }
};
#endif


///////////////////////////
// Generic solve routine //
///////////////////////////

/* \param solver_name The name of the Amesos2 solver interface to use
 *            when solving.
 *
 * \param A1  an initialized matrix which should be the first subject
 *            of consideration for solution.
 *
 * \param A2  another pre-initialized matrix which will replace A1 at
 *            some point using the setA() method of an Amesos2 Solver,
 *            but should have the same symbolic structure as that of
 *            A1.
 *
 * \param Xhat A multivector which will be used to store potential
 *            solutions.
 *
 * \param x A non-empty array of solutions.  x[0] corresponds to the
 *           solution of A1 x = b[0], x[1] corresponds to the solution
 *           of A1 x = b[1], and so on.
 *
 * \param b   A non-empty array of right-hand-side multivectors.  The
 *            entire array will be stepped through twice.  During the
 *            first round, the first multivector in the array will be
 *            used to create an Amesos2 Solver object, and further
 *            rhs's will be substituted using the setB syntax.  During
 *            the second round, the Solver object will be created
 *            without a RHS, and the solve(X,B) method will be used.
 *
 * \param x2  True solution to A2 x2 = b2
 *
 * \param b2  RHS vector for use with A2
 *
 * \param num_vecs the number of vectors in each multivector of x and b
 *
 * \param solve_params parameters to give the the solver
 *            post-construction.
 *
 * The following solution scheme will be used:
 *
 *   symb -> num -> solve [ -> solve -> solve -> ... ] [ -> num -> solve ]
 *
 * We will also test using the shortcuts:
 *
 *   solve [ -> solve -> ... ] [ -> solve ]
 *
 * Where the last solve is done after replacing A2 with A1 as
 * indicated by the global `refactor' option.
 *
 * Although perhaps not necessary, but for completeness, we will test
 * incomplete use of an Amesos2 Solver object.  That is, a Solver
 * object is created, but is never used completely to solution of a
 * system.
 *
 *  - solver is created and then promptly destroyed
 *  - solver is created, computes numeric factorization, and is then
 *    destroyed after retrieving L and U statistics
 */
template <class Matrix, class Vector>
bool
do_solve_routine(const string& solver_name,
                 const RCP<Matrix> A1,
                 const RCP<Matrix> A2,
                 const RCP<Vector> Xhat,
                 const ArrayView<const RCP<Vector> > x,
                 const ArrayView<const RCP<Vector> > b,
                 const RCP<Vector> x2,
                 const RCP<Vector> b2,
                 const size_t num_vecs,
                 ParameterList solve_params)
{
  using std::endl;
  typedef typename ArrayView<const RCP<Vector> >::iterator rhs_it_t;
  bool success = true;          // prove me wrong!

  solution_checker<Vector> checker;
  RCP<Amesos2::Solver<Matrix,Vector> > solver;

  int phase = Amesos2::CLEAN;   // start with no computation

  while (phase != Amesos2::SOLVE) {
    // and X and B will never be needed, so we create the solver without them
    solver = Amesos2::create<Matrix,Vector> (solver_name, A1);
    //JDB: We should really use the parameters the user gives
    solver->setParameters( rcpFromRef(solve_params) );
    switch (phase) {
    case Amesos2::CLEAN:
      break;
    case Amesos2::PREORDERING:
      solver->preOrdering ();
      break;
    case Amesos2::SYMBFACT:
      solver->symbolicFactorization ();
      break;
    case Amesos2::NUMFACT:
      solver->numericFactorization ();
    }
    ++phase;
  }

  enum ESolveStyle {
    SOLVE_VERBOSE,
    SOLVE_XB,
    SOLVE_SHORT
  };

  // declare as 'int' to allow incrementing
  int style = SOLVE_VERBOSE;

  size_t count = 0;
  while (style <= SOLVE_SHORT) {
    rhs_it_t rhs_it = b.begin();
    rhs_it_t x_it = x.begin();
    count++;
    // Create our solver according to the current style
    switch (style) {
    case SOLVE_VERBOSE:
      solver = Amesos2::create<Matrix,Vector>(solver_name, A1, Xhat, *rhs_it);
      break;
    case SOLVE_XB:
      solver = Amesos2::create<Matrix,Vector>(solver_name, A1);
      break;
    case SOLVE_SHORT:
      solver = Amesos2::create<Matrix,Vector>(solver_name, A1);
      break;

    }
    if(num_vecs < count)
      {}


    solver->setParameters( rcpFromRef(solve_params) );

    // Do first solve according to our current style
    switch( style ){
    case SOLVE_VERBOSE:
      solver->preOrdering();
      solver->symbolicFactorization();
      solver->numericFactorization();
      solver->solve();
      break;
    case SOLVE_XB:
      solver->preOrdering();
      solver->symbolicFactorization();
      solver->numericFactorization();
      solver->solve(outArg(*Xhat), ptrInArg(**rhs_it));
      break;
    case SOLVE_SHORT:
      solver->solve(outArg(*Xhat), ptrInArg(**rhs_it));
      break;
    }

    success &= checker(*x_it, Xhat);
    if (! success) {
      return success; // bail out early if necessary
    }

    if( multiple_solves ){
      rhs_it_t rhs_end = b.end();
      for( ; rhs_it != rhs_end; ++rhs_it, ++x_it ){
        switch( style ){
        case SOLVE_VERBOSE:
          solver->setB(*rhs_it);
          solver->solve();
          break;
        case SOLVE_XB:
        case SOLVE_SHORT:
          solver->solve(outArg(*Xhat), ptrInArg(**rhs_it));
          break;
        }
        success &= checker(*x_it, Xhat);
        if( !success ) return success; // bail out early if necessary
      }
    }

    if( refactor ){


      // Keep the symbolic factorization from A1
      solver->setA(A2, Amesos2::SYMBFACT);

      switch( style ){
      case SOLVE_VERBOSE:
        solver->numericFactorization();
        solver->setB(b2);
        solver->solve();
        break;
      case SOLVE_XB:
        solver->numericFactorization();
        solver->solve(outArg(*Xhat), ptrInArg(*b2));
        break;
      case SOLVE_SHORT:
        solver->solve(outArg(*Xhat), ptrInArg(*b2));
        break;
      }

      success &= checker(x2, Xhat);
      if( !success ) return success; // bail out early if necessary
    }

    ++style;
  }

  return success;
}


#ifdef HAVE_TPETRA_INST_INT_INT
#if defined(HAVE_AMESOS2_EPETRA) && defined(HAVE_AMESOS2_EPETRAEXT)

//////////////////////////
//     Epetra Tests     //
//////////////////////////

bool do_epetra_test(const string& mm_file,
                    const string& solver_name,
                    ParameterList solve_params)
{
  using Teuchos::ScalarTraits;

  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;
  const size_t numVecs = 5;     // arbitrary number
  const size_t numRHS  = 5;     // also quite arbitrary

  bool transpose = solve_params.get<bool>("Transpose", false);

#ifdef HAVE_MPI
  const Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  const Epetra_SerialComm comm;
#endif

  if( verbosity > 2 ){
    *fos << std::endl << "      Reading matrix from " << mm_file << " ... " << std::flush;
  }
  std::string path = filedir + mm_file;
  MAT* A;
  int ret = EpetraExt::MatrixMarketFileToCrsMatrix(path.c_str(), comm, A, false, false);
  if( ret == -1 ){
    *fos << "error reading file from disk, aborting run." << std::endl;
    return( false );
  }
  if( verbosity > 2 ) *fos << "done" << std::endl;
  if( verbosity > 3 ) A->Print(std::cout);

  const Epetra_Map dmnmap = A->DomainMap();
  const Epetra_Map rngmap = A->RangeMap();

  RCP<MAT> A_rcp(A);
  RCP<MAT> A2;
  RCP<MV> x2, b2;
  RCP<MV> Xhat;
  if( transpose ){
    Xhat = rcp(new MV(dmnmap,numVecs));
    if( refactor ){
      x2 = rcp(new MV(dmnmap,numVecs));
      b2 = rcp(new MV(rngmap,numVecs));
    }
  } else {
    Xhat = rcp(new MV(rngmap,numVecs));
    if( refactor ){
      x2 = rcp(new MV(rngmap,numVecs));
      b2 = rcp(new MV(dmnmap,numVecs));
    }
  }
  Xhat->SetLabel("Xhat");

  Array<RCP<MV> > x(numRHS);
  Array<RCP<MV> > b(numRHS);
  for( size_t i = 0; i < numRHS; ++i ){
    if( transpose ){
      x[i] = rcp(new MV(dmnmap,numVecs));
      b[i] = rcp(new MV(rngmap,numVecs));
    } else {
      x[i] = rcp(new MV(rngmap,numVecs));
      b[i] = rcp(new MV(dmnmap,numVecs));
    }
    std::ostringstream xlabel, blabel;
    xlabel << "x[" << i << "]";
    blabel << "b[" << i << "]";
    x[i]->SetLabel(xlabel.str().c_str());
    b[i]->SetLabel(blabel.str().c_str());

    x[i]->Random();
    A->Multiply(transpose, *x[i], *b[i]);
  }

  if( refactor ){
    // There isn't a really nice way to get a deep copy of an entire
    // CrsMatrix, so we just read the file again.
    MAT* A2_ptr;
    ret = EpetraExt::MatrixMarketFileToCrsMatrix(path.c_str(), comm,
                                                 A2_ptr,
                                                 false, false);
    if( ret == -1 ){
      *fos << "error reading file from disk, aborting run." << std::endl;
      return( false );
    }
    A2.reset(A2_ptr);

    // perturb the values just a bit (element-wise square of first row)
    int l_fst_row_nnz = 0;
    A2->NumMyRowEntries(0, l_fst_row_nnz);
    Array<int> indices(l_fst_row_nnz);
    Array<double> values(l_fst_row_nnz);
    A2->ExtractMyRowCopy(0, l_fst_row_nnz, l_fst_row_nnz,
                         values.getRawPtr(), indices.getRawPtr());
    for( int i = 0; i < l_fst_row_nnz; ++i ){
      values[i] = values[i] * values[i];
    }
    A2->ReplaceMyValues(0, l_fst_row_nnz, values.getRawPtr(), indices.getRawPtr());

    A2->Multiply(transpose, *x2, *b2);
  } // else A2 is never read

  const bool result = do_solve_routine<MAT,MV>(solver_name, A_rcp, A2,
                                               Xhat, x(), b(), x2, b2,
                                               numRHS, solve_params);

  if (!result) {
    if( verbosity > 1 ){
      *fos << "failed!" << std::endl;
    }
    return( false );
  } else {
    if( verbosity > 1 ){
      *fos << "passed" << std::endl;
    }
    return( true );
  }
}

bool test_epetra(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& epetra_runs,
                 ParameterList solve_params)
{
  bool success = true;
  bool do_default = false;

  ParameterList::ConstIterator run_it;
  for( run_it = epetra_runs.begin(); run_it != epetra_runs.end(); ++run_it ){
    // an empty parameters list or any plain Parameter causes a default test for epetra
    if( epetra_runs.entry(run_it).isList() ){
      ParameterList run_list = Teuchos::getValue<ParameterList>(epetra_runs.entry(run_it));
      if( run_list.isSublist("solver_run_params") ){
        ParameterList solve_params_copy(solve_params);
        solve_params_copy.sublist(solver_name).setParameters( run_list.sublist("solver_run_params") );
        if( run_list.isSublist("amesos2_params") ){
          solve_params_copy.setParameters( run_list.sublist("amesos2_params") );
        }

        string run_name = epetra_runs.name(run_it);
        if( verbosity > 1 ){
          *fos << "    Doing epetra test run `" << run_name << "' ... " << std::flush;
        }
        success &= do_epetra_test(mm_file, solver_name, solve_params_copy);
      } else {
        do_default = true;
      }
    } else {
      do_default = true;
    }
  }

  // only do one default run
  if( do_default ){
    if( verbosity > 1 ){
      *fos << "    Doing epetra test default test run ... " << std::flush;
    }
    success &= do_epetra_test(mm_file, solver_name, solve_params);
  }

  return( success );
}
#endif  // HAVE_AMESOS2_EPETRAEXT
#endif // HAVE_TPETRA_INST_INT_INT
//////////////////////////
//     Tpetra Tests     //
//////////////////////////

template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal,
         typename Node>
bool do_tpetra_test_with_types(const string& mm_file,
                               const string& solver_name,
                               ParameterList solve_params)
{
  using Tpetra::Map;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Teuchos::Comm;
  using Teuchos::ScalarTraits;
  using std::endl;
  using std::flush;

  typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MAT;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  const size_t numVecs = 5;     // arbitrary number
  const size_t numRHS = 5;        // also arbitrary

  bool transpose = solve_params.get<bool>("Transpose", false);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  if (verbosity > 2) {
    *fos << endl << "      Reading matrix from " << mm_file << " ... " << flush;
  }
  std::string path = filedir + mm_file;
  RCP<MAT> A =
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile (path, comm);

  if (verbosity > 2) {
    *fos << "done" << endl;
    switch (verbosity) {
    case 6:
      A->describe (*fos, Teuchos::VERB_EXTREME); break;
    case 5:
      A->describe (*fos, Teuchos::VERB_HIGH); break;
    case 4:
      A->describe (*fos, Teuchos::VERB_LOW); break;
    }
  }

  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > dmnmap = A->getDomainMap();
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rngmap = A->getRangeMap();

  ETransp trans = transpose ? CONJ_TRANS : NO_TRANS;

  if (verbosity > 2) {
    *fos << endl << "      Creating right-hand side and solution vectors" << endl;
  }

  RCP<MAT> A2;
  RCP<MV> Xhat, x2, b2;
  if (transpose) {
    Xhat = rcp(new MV(dmnmap,numVecs));
    if (refactor) {
      x2 = rcp(new MV(dmnmap,numVecs));
      b2 = rcp(new MV(rngmap,numVecs));
    }
  } else {
    Xhat = rcp(new MV(rngmap,numVecs));
    if (refactor) {
      x2 = rcp(new MV(rngmap,numVecs));
      b2 = rcp(new MV(dmnmap,numVecs));
    }
  }
  Xhat->setObjectLabel("Xhat");

  Array<RCP<MV> > x(numRHS);
  Array<RCP<MV> > b(numRHS);
  for( size_t i = 0; i < numRHS; ++i ){
    if( transpose ){
      x[i] = rcp(new MV(dmnmap,numVecs));
      b[i] = rcp(new MV(rngmap,numVecs));
    } else {
      x[i] = rcp(new MV(rngmap,numVecs));
      b[i] = rcp(new MV(dmnmap,numVecs));
    }
    std::ostringstream xlabel, blabel;
    xlabel << "x[" << i << "]";
    blabel << "b[" << i << "]";
    x[i]->setObjectLabel(xlabel.str());
    b[i]->setObjectLabel(blabel.str());

    x[i]->randomize();
    A->apply(*x[i], *b[i], trans);
  }

  if (refactor) {
    if (verbosity > 2) {
      *fos << endl << "      Creating near-copy of matrix for refactor test" << endl;
    }

    // Make a deep copy of the entire CrsMatrix.
    // originalNode can be null; only needed for type deduction.
    A2 = rcp(new MAT(*A,Teuchos::Copy));

    // // There isn't a really nice way to get a deep copy of an entire
    // // CrsMatrix, so we just read the file again.
    // A2 = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path, comm);

    // perturb the values just a bit (element-wise square of first row)
    size_t l_fst_row_nnz = A2->getNumEntriesInLocalRow(0);
    Array<LocalOrdinal> indices(l_fst_row_nnz);
    Array<Scalar> values(l_fst_row_nnz);
    A2->getLocalRowCopy(0, indices, values, l_fst_row_nnz);
    for( size_t i = 0; i < l_fst_row_nnz; ++i ){
      values[i] = values[i] * values[i];
    }
    A2->resumeFill ();
    A2->replaceLocalValues (0, indices, values);
    A2->fillComplete (A->getDomainMap (), A->getRangeMap ());

    A2->apply (*x2, *b2, trans);
  } // else A2 is never read

  return do_solve_routine<MAT,MV>(solver_name, A, A2,
                                  Xhat, x(), b(), x2, b2,
                                  numRHS, solve_params);
}


bool test_tpetra(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& tpetra_runs,
                 ParameterList solve_params)
{
  bool success = true;

  typedef DefaultNode DN;

  ParameterList::ConstIterator run_it;
  for( run_it = tpetra_runs.begin(); run_it != tpetra_runs.end(); ++run_it ){
    if( tpetra_runs.entry(run_it).isList() ){
      ParameterList run_list = Teuchos::getValue<ParameterList>(tpetra_runs.entry(run_it));
      ParameterList solve_params_copy(solve_params);
      if( run_list.isSublist("solver_run_params") ){
        solve_params_copy.sublist(solver_name).setParameters( run_list.sublist("solver_run_params") );
      }
      if( run_list.isSublist("amesos2_params") ){
        solve_params_copy.setParameters( run_list.sublist("amesos2_params") );
      }

      string scalar, mag, lo, go, node;
      if( run_list.isParameter("Scalar") ){
        scalar = run_list.get<string>("Scalar");
        if( scalar == "complex" ){
#ifdef HAVE_TEUCHOS_COMPLEX
          // get the magnitude parameter
          if( run_list.isParameter("Magnitude") ){
            mag = run_list.get<string>("Magnitude");
          } else {
            *fos << "    Must provide a type for `Magnitude' when Scalar='complex', aborting run `"
                 << run_list.name() << "'..." << std::endl;
            continue;
          }
#else
          if( verbosity > 1 ){
            *fos << "    Complex support not enabled, skipping run `"
                 << run_list.name() << "'..." << std::endl;
          }
          continue;
#endif  // HAVE_TEUCHOS_COMPLEX
        }
      } else {
        *fos << "    Must provide a type for `Scalar', aborting run `"
             << run_list.name() << "'..." << std::endl;
        continue;
      }
      if( run_list.isParameter("LocalOrdinal") ){
        lo = run_list.get<string>("LocalOrdinal");
      } else {
        *fos << "    Must provide a type for `LocalOrdinal', aborting run `"
             << run_list.name() << "'..." << std::endl;
        continue;
      }
      if( run_list.isParameter("GlobalOrdinal") ){
        go = run_list.get<string>("GlobalOrdinal");
      } else {
        go = "default";
      }
      if( run_list.isParameter("Node") ){
        scalar = run_list.get<string>("Node");
      } else {
        node = "default";
      }

      string run_name = tpetra_runs.name(run_it);
      if( verbosity > 1 ){
        *fos << "    Doing tpetra test run `"
             << run_name << "' with"
             << " s=" << scalar;
        if( scalar == "complex" ){
          *fos << "(" << mag << ")";
        }
        *fos << " lo=" << lo
             << " go=" << go
             << " ... " << std::flush;
      }

      string timer_name = mm_file + "_" + scalar + "_" + lo + "_" + go;
      RCP<Time> timer = TimeMonitor::getNewTimer(timer_name);
      TimeMonitor LocalTimer(*timer);

      bool test_done = false;

#define AMESOS2_SOLVER_TPETRA_TEST(S,LO,GO,N)                           \
      test_done = true;                                                 \
      if (verbosity > 1) {                                              \
        *fos << std::endl                                               \
             << "Running test with types "                              \
             << "S=" << Teuchos::TypeNameTraits<S>::name ()             \
             << ", LO=" << Teuchos::TypeNameTraits<LO>::name ()         \
             << ", GO=" << Teuchos::TypeNameTraits<GO>::name ()         \
             << ", N=" << Teuchos::TypeNameTraits<N>::name ()           \
             << std::endl;                                              \
      }                                                                 \
      bool run_success = do_tpetra_test_with_types<S,LO,GO,N>(mm_file,solver_name, \
                                                              solve_params_copy); \
      if (verbosity > 1) {                                              \
        if (!run_success)                                               \
          *fos << "failed!" << std::endl;                               \
        else                                                            \
          *fos << "passed" << std::endl;                                \
      }                                                                 \
      success &= run_success


      using default_go_type = Tpetra::Map<>::global_ordinal_type;

      // I can't think of any better way to do the types as
      // specified in the parameter list at runtime but to branch
      // out on all possibilities, sorry...  Note: we're going to
      // ignore the `node' parameter for now
      if( scalar == "float" ){
#ifdef HAVE_TPETRA_INST_FLOAT
        if( lo == "int" ){
          if( go == "default" ){
            AMESOS2_SOLVER_TPETRA_TEST(float,int,default_go_type,DN);
          }
          else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
            AMESOS2_SOLVER_TPETRA_TEST(float,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
          }
          else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
            AMESOS2_SOLVER_TPETRA_TEST(float,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
          }
          else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
            AMESOS2_SOLVER_TPETRA_TEST(float,int,long long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG_LONG
          }
        }
        else if( lo == "long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
        else if( lo == "long long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
#else // NOT HAVE_TPETRA_INST_FLOAT
        *fos << "Scalar=float was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_FLOAT
      } // end scalar == "float"

      if( scalar == "double" ){
#ifdef HAVE_TPETRA_INST_DOUBLE
        if( lo == "int" ){
          if( go == "default" ){
            AMESOS2_SOLVER_TPETRA_TEST(double,int,default_go_type,DN);
          }
          else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
            AMESOS2_SOLVER_TPETRA_TEST(double,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
          }
          else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
            AMESOS2_SOLVER_TPETRA_TEST(double,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
          }
          else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INT_LONG_LONG
            AMESOS2_SOLVER_TPETRA_TEST(double,int,long long int,DN);
#else // NOT HAVE_TPETRA_INT_LONG_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INT_LONG_LONG
          }
        }
        else if( lo == "long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
        else if( lo == "long long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
#else // NOT HAVE_TPETRA_INST_DOUBLE
        *fos << "Scalar=double was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_DOUBLE
      } // end scalar == "double"

      if( scalar == "double double" ){
#ifdef HAVE_TPETRA_INST_DD_REAL
        if( lo == "int" ){
          if( go == "default" ){
            AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,default_go_type,DN);
          }
          else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
            AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
          }
          else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
            AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
          }
          else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INT_LONG_LONG
            AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,long long int,DN);
#else // NOT HAVE_TPETRA_INT_LONG_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INT_LONG_LONG
          }
        }
        else if( lo == "long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
        else if( lo == "long long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
#else // NOT HAVE_TPETRA_INST_DD_REAL
        *fos << "Scalar=dd_real was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_DD_REAL
      } // end scalar == "double double"

      if( scalar == "quad" || scalar == "quad double" ){
#ifdef HAVE_TPETRA_INST_QD_REAL
        if( lo == "int" ){
          if( go == "default" ){
            AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,default_go_type,DN);
          }
          else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
            AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
          }
          else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
            AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
          }
          else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INT_LONG_LONG
            AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,long long int,DN);
#else // NOT HAVE_TPETRA_INT_LONG_LONG
            *fos << "GO=" << go << " was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INT_LONG_LONG
          }
        }
        else if( lo == "long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
        else if( lo == "long long int" ){
          *fos << "Trilinos does not currently support LO=" << lo << std::endl;
        }
#else // NOT HAVE_TPETRA_INST_QD_REAL
        *fos << "Scalar=qd_real was not enabled at configure time" << std::endl;
#endif // HAVE_TPETRA_INST_QD_REAL
      } // end scalar == "quad double"

      if( scalar == "complex" ){
        if( mag == "float" ){
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
          typedef std::complex<float> cmplx;
          if( lo == "int" ){
            if( go == "default" ){
              AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,default_go_type,DN);
            }
            else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
              AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
            }
            else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
              AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
            }
            else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INT_LONG_LONG
              AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,long long int,DN);
#else // NOT HAVE_TPETRA_INT_LONG_LONG
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INT_LONG_LONG
            }
          }
          else if( lo == "long int" ){
            *fos << "Trilinos does not currently support LO=" << lo << std::endl;
          }
          else if( lo == "long long int" ){
            *fos << "Trilinos does not currently support LO=" << lo << std::endl;
          }
#else // NOT HAVE_TPETRA_INST_COMPLEX_FLOAT
          *fos << "Scalar=complex<float> was not enabled at configure time"
               << std::endl;
#endif // HAVE_TPETRA_INST_COMPLEX_FLOAT
        }

        if( mag == "double" ){
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
          typedef std::complex<double> cmplx_double;
          if( lo == "int" ){
            if( go == "default" ){
              AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,default_go_type,DN);
            }
            else if( go == "int" ){
#ifdef HAVE_TPETRA_INST_INT_INT
              AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,int,DN);
#else // NOT HAVE_TPETRA_INST_INT_INT
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INST_INT_INT
            }
            else if( go == "long int" ){
#ifdef HAVE_TPETRA_INST_INT_LONG
              AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,long int,DN);
#else // NOT HAVE_TPETRA_INST_INT_LONG
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INST_INT_LONG
            }
            else if( go == "long long int" ){
#ifdef HAVE_TPETRA_INT_LONG_LONG
              AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,long long int,DN);
#else // NOT HAVE_TPETRA_INT_LONG_LONG
              *fos << "GO=" << go << " was not enabled at configure time"
                   << std::endl;
#endif // HAVE_TPETRA_INT_LONG_LONG
            }
          }
          else if( lo == "long int" ){
            *fos << "Trilinos does not currently support LO=" << lo << std::endl;
          }
          else if( lo == "long long int" ){
            *fos << "Trilinos does not currently support LO=" << lo << std::endl;
          }
#else // NOT HAVE_TPETRA_INST_COMPLEX_DOUBLE
          *fos << "Scalar=complex<double> was not enabled at configure time"
               << std::endl;
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE
        } // end mag == "double"

        // mfh 15 Jan 2019: The code used to test
        // Scalar=std::complex<dd_real> and std::complex<qd_real>
        // here.  The C++ Standard does not permit portable use of
        // std::complex<T> for T != float, double, or long double.
        // Thus, I deleted those cases, which almost certainly were
        // never tested, even if they ever built.

      } // scalar == "complex"

      if( !test_done && verbosity > 1 ){
        *fos << "type parameters not recognized or enabled" << std::endl;
      }
    } // else there are no runs given; no default cases for tpetra
  }

  return( success );
}

template<typename Scalar,
         typename LocalOrdinal,
         typename Node,
         typename TpetraScalar>
bool do_kokkos_test_with_types(const string& mm_file,
                               const string& solver_name,
                               ParameterList solve_params)
{
  // Here I am using Tpetra as a helper to load and make solutions.
  typedef DefaultNode TpetraNode;

  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Teuchos::Comm;
  using Teuchos::ScalarTraits;
  using std::endl;
  using std::flush;

  typedef Kokkos::Device<typename Node::execution_space, typename Node::memory_space> device_t;
  typedef KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,device_t> MAT;
  typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft, device_t> view_t;

  const size_t numVecs = 5;     // arbitrary number
  const size_t numRHS = 5;      // also arbitrary

  bool transpose = solve_params.get<bool>("Transpose", false);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  if (verbosity > 2) {
    *fos << endl << "      Reading matrix from " << mm_file << " ... " << flush;
  }
  std::string path = filedir + mm_file;

  // not sure about the loading schemes for kokkos - just use Tpetra right now
  // and load it into the Kokkos CrsMatrix
  typedef Tpetra::Map<>::global_ordinal_type TpetraGO;
  typedef Tpetra::CrsMatrix<TpetraScalar,LocalOrdinal,TpetraGO,TpetraNode> tpetra_crsmatrix_t;
  RCP<tpetra_crsmatrix_t> tpetraM =
    Tpetra::MatrixMarket::Reader<tpetra_crsmatrix_t>::readSparseFile (path, comm);

  Teuchos::ArrayRCP<const size_t> rowPointers;
  Teuchos::ArrayRCP<const LocalOrdinal> columnIndices;
  Teuchos::ArrayRCP<const TpetraScalar> values;
  tpetraM->getAllValues(rowPointers, columnIndices, values);

  // convert Tpetra size_t row ptrs to kokkos crsmatrix LO type
  Teuchos::ArrayRCP<LocalOrdinal> kokkosRowPointers = Teuchos::arcp(new LocalOrdinal[rowPointers.size()], 0, rowPointers.size());
  for(int n = 0; n < rowPointers.size(); ++n) {
    kokkosRowPointers[n] = Teuchos::as<LocalOrdinal>(rowPointers[n]);
  }

  // convert Tpetra values to kokkos - for complex this will be converting
  // std::complex to Kokkos::complex
  Teuchos::ArrayRCP<Scalar> kokkosValues = Teuchos::arcp(new Scalar[values.size()], 0, values.size());
  for(int n = 0; n < values.size(); ++n) {
    kokkosValues[n] = Teuchos::as<Scalar>(values[n]);
  }

  auto num_rows = tpetraM->getNodeNumRows();
  auto num_cols = tpetraM->getNodeNumCols();

  // Kokkos CrsMatrix builds with non const ptrs while Tpetra loads into const
  Teuchos::ArrayRCP<Scalar> non_const_values = Teuchos::arcp_const_cast<Scalar>(kokkosValues);
  Teuchos::ArrayRCP<LocalOrdinal> non_const_kokkosRowPointers = Teuchos::arcp_const_cast<LocalOrdinal>(kokkosRowPointers);
  Teuchos::ArrayRCP<LocalOrdinal> non_const_columnIndices = Teuchos::arcp_const_cast<LocalOrdinal>(columnIndices);
  RCP<MAT> A = rcp(new MAT("Kokkos CrsMatrix",
             num_rows,
             num_cols,
             tpetraM->getGlobalNumEntries(),
             non_const_values.getRawPtr(),
             non_const_kokkosRowPointers.getRawPtr(),
             non_const_columnIndices.getRawPtr()));

  if (verbosity > 2) {
    *fos << endl << "      Creating right-hand side and solution vectors" << endl;
  }

  ETransp trans = transpose ? CONJ_TRANS : NO_TRANS;
  typedef Tpetra::MultiVector<TpetraScalar,LocalOrdinal,TpetraGO,TpetraNode> MV;
  RCP<const Tpetra::Map<LocalOrdinal,TpetraGO,TpetraNode> > dmnmap = tpetraM->getDomainMap();
  RCP<const Tpetra::Map<LocalOrdinal,TpetraGO,TpetraNode> > rngmap = tpetraM->getRangeMap();
  Array<RCP<MV> > xMV(numRHS);
  Array<RCP<MV> > bMV(numRHS);
  for( size_t i = 0; i < numRHS; ++i ){
    if( transpose ){
      xMV[i] = rcp(new MV(dmnmap,numVecs));
      bMV[i] = rcp(new MV(rngmap,numVecs));
    } else {
      xMV[i] = rcp(new MV(rngmap,numVecs));
      bMV[i] = rcp(new MV(dmnmap,numVecs));
    }
    std::ostringstream xlabel, blabel;
    xlabel << "x[" << i << "]";
    blabel << "b[" << i << "]";
    xMV[i]->setObjectLabel(xlabel.str());
    bMV[i]->setObjectLabel(blabel.str());

    xMV[i]->randomize();
    tpetraM->apply(*xMV[i], *bMV[i], trans);
  }

  Array<RCP<view_t>> x(numRHS);
  Array<RCP<view_t>> b(numRHS);
  for( size_t i = 0; i < numRHS; ++i ){
    std::ostringstream xlabel, blabel;
    xlabel << "x[" << i << "]";
    blabel << "b[" << i << "]";
    if( transpose ){
      x[i] = rcp(new view_t(xlabel.str(), num_cols, numRHS));
      b[i] = rcp(new view_t(blabel.str(), num_cols, numRHS));
    } else {
      x[i] = rcp(new view_t(xlabel.str(), num_rows, numRHS));
      b[i] = rcp(new view_t(blabel.str(), num_rows, numRHS));
    }

    // MDM Right now I'm employing the Tpetra version to generate the random
    // values. But probably should make a pure kokkos version though I'd like
    // them all to be the same.
    RCP<MV> xMV, bMV;
    if( transpose ){
      xMV = rcp(new MV(rngmap,numVecs));
      bMV = rcp(new MV(dmnmap,numVecs));
    }
    else {
      xMV = rcp(new MV(dmnmap,numVecs));
      bMV = rcp(new MV(rngmap,numVecs));
    }
    xMV->randomize();
    tpetraM->apply(*xMV, *bMV, trans);

    Kokkos::deep_copy(*x[i], xMV->getLocalViewHost());
    Kokkos::deep_copy(*b[i], bMV->getLocalViewHost());
  }

  RCP<tpetra_crsmatrix_t> temp_tpetraM =
    Tpetra::MatrixMarket::Reader<tpetra_crsmatrix_t>::readSparseFile (path, comm);

  RCP<MAT> A2;
  RCP<view_t> Xhat, x2, b2;

  if (transpose) {
    Xhat = rcp(new view_t("Xhat", num_rows, numVecs));
    if (refactor) {
      x2 = rcp(new view_t("x2", num_rows, numVecs));
      b2 = rcp(new view_t("b2", num_cols, numVecs));
    }
  } else {
    Xhat = rcp(new view_t("Xhat", num_cols, numVecs));
    if (refactor) {
      x2 = rcp(new view_t("x2", num_cols, numVecs));
      b2 = rcp(new view_t ("b2", num_rows, numVecs));
    }
  }

  if (refactor) {
    if (verbosity > 2) {
      *fos << endl << "      Creating near-copy of matrix for refactor test" << endl;
    }

    RCP<tpetra_crsmatrix_t> tpetraM2 =
      Tpetra::MatrixMarket::Reader<tpetra_crsmatrix_t>::readSparseFile (path, comm);

    // perturb the values just a bit (element-wise square of first row)
    size_t l_fst_row_nnz = tpetraM2->getNumEntriesInLocalRow(0);
    Array<LocalOrdinal> indices(l_fst_row_nnz);
    Array<TpetraScalar> values(l_fst_row_nnz);
    tpetraM2->getLocalRowCopy(0, indices, values, l_fst_row_nnz);
    for( size_t i = 0; i < l_fst_row_nnz; ++i ){
      values[i] = values[i] * values[i];
    }
    tpetraM2->resumeFill ();
    tpetraM2->replaceLocalValues (0, indices, values);
    tpetraM2->fillComplete (tpetraM->getDomainMap (), tpetraM->getRangeMap ());

    // Get Tpetra ahain
    Teuchos::ArrayRCP<const size_t> rowPointers2;
    Teuchos::ArrayRCP<const LocalOrdinal> columnIndices2;
    Teuchos::ArrayRCP<const TpetraScalar> values2;
    tpetraM2->getAllValues(rowPointers2, columnIndices2, values2);

    Teuchos::ArrayRCP<LocalOrdinal> kokkosRowPointers2 = Teuchos::arcp(new LocalOrdinal[rowPointers2.size()], 0, rowPointers2.size());
    for(int n = 0; n < rowPointers2.size(); ++n) {
      kokkosRowPointers2[n] = Teuchos::as<LocalOrdinal>(rowPointers2[n]);
    }

    auto num_rows2 = tpetraM2->getNodeNumRows();
    auto num_cols2 = tpetraM2->getNodeNumCols();

    A2 = rcp(new MAT("Kokkos CrsMatrix 2",
             num_rows2,
             num_cols2,
             tpetraM2->getGlobalNumEntries(),
             (Scalar*)values2.getRawPtr(),
             (LocalOrdinal*)kokkosRowPointers2.getRawPtr(),
             (LocalOrdinal*)columnIndices2.getRawPtr()));

    RCP<MV> x2MV, b2MV;
    if (transpose) {
      if (refactor) {
        x2MV = rcp(new MV(dmnmap,numVecs));
        b2MV = rcp(new MV(rngmap,numVecs));
      }
    } else {
      if (refactor) {
        x2MV = rcp(new MV(rngmap,numVecs));
        b2MV = rcp(new MV(dmnmap,numVecs));
      }
    }

    tpetraM2->apply (*x2MV, *b2MV, trans);

    Kokkos::deep_copy(*x2, x2MV->getLocalViewHost());
    Kokkos::deep_copy(*b2, b2MV->getLocalViewHost());
  } // else A2 is never read

  return do_solve_routine<MAT,view_t>(solver_name, A, A2,
                                  Xhat, x, b, x2, b2,
                                  numRHS, solve_params);
}

bool test_kokkos(const string& mm_file,
                 const string& solver_name,
                 const ParameterList& kokkos_runs,
                 ParameterList solve_params)
{
  bool success = true;
  ParameterList::ConstIterator run_it;
  for( run_it = kokkos_runs.begin(); run_it != kokkos_runs.end(); ++run_it ){
    if( kokkos_runs.entry(run_it).isList() ){
      ParameterList run_list = Teuchos::getValue<ParameterList>(kokkos_runs.entry(run_it));
      ParameterList solve_params_copy(solve_params);
      if( run_list.isSublist("solver_run_params") ){
        solve_params_copy.sublist(solver_name).setParameters( run_list.sublist("solver_run_params") );
      }
      if( run_list.isSublist("amesos2_params") ){
        solve_params_copy.setParameters( run_list.sublist("amesos2_params") );
      }

      string scalar, mag, lo, node;
      if( run_list.isParameter("Scalar") ){
        scalar = run_list.get<string>("Scalar");
        if( scalar == "complex" ){
#ifdef HAVE_TEUCHOS_COMPLEX
          // get the magnitude parameter
          if( run_list.isParameter("Magnitude") ){
            mag = run_list.get<string>("Magnitude");
          } else {
            *fos << "    Must provide a type for `Magnitude' when Scalar='complex', aborting run `"
                 << run_list.name() << "'..." << std::endl;
            continue;
          }
#else
          if( verbosity > 1 ){
            *fos << "    Complex support not enabled, skipping run `"
                 << run_list.name() << "'..." << std::endl;
          }
          continue;
#endif  // HAVE_TEUCHOS_COMPLEX
        }
      } else {
        *fos << "    Must provide a type for `Scalar', aborting run `"
             << run_list.name() << "'..." << std::endl;
        continue;
      }
      if( run_list.isParameter("LocalOrdinal") ){
        lo = run_list.get<string>("LocalOrdinal");
      } else {
        *fos << "    Must provide a type for `LocalOrdinal', aborting run `"
             << run_list.name() << "'..." << std::endl;
        continue;
      }
      if( run_list.isParameter("GlobalOrdinal") ){
        *fos << "    GlobalOrdinal setting is currently ignored for Kokkos adapter. Intended for Serial only currently." << std::endl;
      }
      if( run_list.isParameter("Node") ){
        node = run_list.get<string>("Node");
      } else {
        node = "default";
      }

      string run_name = kokkos_runs.name(run_it);
      if( verbosity > 1 ){
        *fos << "    Doing kokkos test run `"
             << run_name << "' with"
             << " s=" << scalar;
        if( scalar == "complex" ){
          *fos << "(" << mag << ")";
        }
        *fos << " lo=" << lo
             << " node=" << node
             << " ... " << std::endl << std::flush;
      }

      string timer_name = mm_file + "_" + scalar + "_" + lo + "_" + node;
      RCP<Time> timer = TimeMonitor::getNewTimer(timer_name);
      TimeMonitor LocalTimer(*timer);

      bool test_done = false;

#define AMESOS2_SOLVER_KOKKOS_TEST(S,LO,N,TpetraScalar)                 \
      test_done = true;                                                 \
      if (verbosity > 1) {                                              \
        *fos << std::endl                                               \
             << "Running kokkos test with types "                       \
             << "S=" << Teuchos::TypeNameTraits<S>::name ()             \
             << ", LO=" << Teuchos::TypeNameTraits<LO>::name ()         \
             << ", N=" << Teuchos::TypeNameTraits<N>::name ()           \
             << std::endl;                                              \
      }                                                                 \
      bool run_success =                                                \
        do_kokkos_test_with_types<S,LO,N,TpetraScalar>                  \
        (mm_file,solver_name, solve_params_copy);                       \
      if (verbosity > 1) {                                              \
        if (!run_success)                                               \
          *fos << "failed!" << std::endl;                               \
        else                                                            \
          *fos << "passed" << std::endl;                                \
      }                                                                 \
      success &= run_success

      // create an option to test UVM off inputs
      #ifdef KOKKOS_ENABLE_CUDA
      typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>  uvm_off_node_t;
      #endif

      if( lo != "int" ) {
        *fos << "Trilinos does not currently support LO=" << lo << std::endl;
      }
      else {
        if( scalar == "float" ) {
          #ifdef HAVE_TPETRA_INST_FLOAT// Because of Tpetra maps this is currently needed for Kokkos adapter
          if( node == "default" ) {
            AMESOS2_SOLVER_KOKKOS_TEST(float,int,DefaultNode,float);
          }
          else if( node == "serial" ) {
            #ifdef KOKKOS_ENABLE_SERIAL
            *fos << "KokkosSerialWrapperNode float ";
            AMESOS2_SOLVER_KOKKOS_TEST(float,int,Kokkos::Serial,float);
            #else
            *fos << "node=serial was not enabled at configure time" << std::endl;
            #endif
          }
          else if( node == "cuda" ) {
            #ifdef KOKKOS_ENABLE_CUDA
            *fos << "KokkosCudaWrapperNode float ";
            AMESOS2_SOLVER_KOKKOS_TEST(float,int,Kokkos::Cuda,float);
            #else
            *fos << "node=cuda was not enabled at configure time" << std::endl;
            #endif
          }
          else if( node == "cudauvmoff" ) {
            #ifdef KOKKOS_ENABLE_CUDA
            *fos << "KokkosCudaUVMOffWrapperNode float ";
            AMESOS2_SOLVER_KOKKOS_TEST(float,int,uvm_off_node_t,float);
            #else
            *fos << "node=cudauvmoff was not enabled at configure time" << std::endl;
            #endif
          }
          #else
          *fos << "scalar=float was not enabled at configure time" << std::endl;
          #endif
        }
        else if( scalar == "double" ) {
          #ifdef HAVE_TPETRA_INST_DOUBLE // Because of Tpetra maps this is currently needed for Kokkos adapter
          if( node == "default" ) {
            AMESOS2_SOLVER_KOKKOS_TEST(double,int,DefaultNode,double);
          }
          else if( node == "serial" ) {
            #ifdef KOKKOS_ENABLE_SERIAL
            *fos << "KokkosSerialWrapperNode double ";
            AMESOS2_SOLVER_KOKKOS_TEST(double,int,Kokkos::Serial,double);
            #else
            *fos << "node=serial was not enabled at configure time" << std::endl;
            #endif
          }
          else if( node == "cuda" ) {
            #ifdef KOKKOS_ENABLE_CUDA
            *fos << "KokkosCudaWrapperNode double ";
            AMESOS2_SOLVER_KOKKOS_TEST(double,int,Kokkos::Cuda,double);
            #else
            *fos << "node=cuda was not enabled at configure time" << std::endl;
            #endif
          }
          else if( node == "cudauvmoff" ) {
            #ifdef KOKKOS_ENABLE_CUDA
            *fos << "KokkosCudaUVMOffWrapperNode double ";
            AMESOS2_SOLVER_KOKKOS_TEST(double,int,uvm_off_node_t,double);
            #else
            *fos << "node=cudauvmoff was not enabled at configure time" << std::endl;
            #endif
          }
          #else
          *fos << "scalar=double was not enabled at configure time" << std::endl;
          #endif
        }
        else if( scalar == "double double" ) {
          *fos << "scalar=double double not implemented yet for kokkos adapter." << std::endl;
        }
        else if( scalar == "quad" || scalar == "quad double" ) {
          *fos << "scalar=qd real not implemented yet for kokkos adapter." << std::endl;
        }
        else if( scalar == "complex" ){
          if( mag == "float" ){
  #ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
            typedef Kokkos::complex<float> cmplx_float;
            typedef std::complex<float> tpetra_cmplx_float;
            if( lo == "int" ){
              if( node == "default" ) {
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_float,int,DefaultNode,tpetra_cmplx_float);
              }
              else if( node == "serial" ) {
                #ifdef KOKKOS_ENABLE_SERIAL
                *fos << "KokkosSerialWrapperNode complex<float> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_float,int,Kokkos::Serial,tpetra_cmplx_float);
                #else
                *fos << "node=serial was not enabled at configure time" << std::endl;
                #endif
              }
              else if( node == "cuda" ) {
                #ifdef KOKKOS_ENABLE_CUDA
                *fos << "KokkosCudaWrapperNode complex<float> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_float,int,Kokkos::Cuda,tpetra_cmplx_float);
                #else
                *fos << "node=cuda was not enabled at configure time" << std::endl;
                #endif
              }
              else if( node == "cudauvmoff" ) {
                #ifdef KOKKOS_ENABLE_CUDA
                *fos << "KokkosCudaUVMOffWrapperNode complex<float> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_float,int,uvm_off_node_t,tpetra_cmplx_float);
                #else
                *fos << "node=cudauvmoff was not enabled at configure time" << std::endl;
                #endif
              }
            }
            else if( lo == "long int" ){
              *fos << "Trilinos does not currently support LO=" << lo << std::endl;
            }
            else if( lo == "long long int" ){
              *fos << "Trilinos does not currently support LO=" << lo << std::endl;
            }
  #else // NOT HAVE_TPETRA_INST_COMPLEX_FLOAT
            *fos << "Scalar=complex<float> was not enabled at configure time"
                 << std::endl;
  #endif // HAVE_TPETRA_INST_COMPLEX_FLOAT
          }

          if( mag == "double" ){
  #ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
            typedef Kokkos::complex<double> cmplx_double;
            typedef std::complex<double> tpetra_cmplx_double;
            if( lo == "int" ){
              if( node == "default" ) {
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_double,int,DefaultNode,tpetra_cmplx_double);
              }
              else if( node == "serial" ) {
                #ifdef KOKKOS_ENABLE_SERIAL
                *fos << "KokkosSerialWrapperNode complex<double> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_double,int,Kokkos::Serial,tpetra_cmplx_double);
                #else
                *fos << "node=serial was not enabled at configure time" << std::endl;
                #endif
              }
              else if( node == "cuda" ) {
                #ifdef KOKKOS_ENABLE_CUDA
                *fos << "KokkosCudaWrapperNode complex<double> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_double,int,Kokkos::Cuda,tpetra_cmplx_double);
                #else
                *fos << "node=cuda was not enabled at configure time" << std::endl;
                #endif
              }
              else if( node == "cudauvmoff" ) {
                #ifdef KOKKOS_ENABLE_CUDA
                *fos << "KokkosCudaUVMOffWrapperNode complex<double> ";
                AMESOS2_SOLVER_KOKKOS_TEST(cmplx_double,int,uvm_off_node_t,tpetra_cmplx_double);
                #else
                *fos << "node=cudauvmoff was not enabled at configure time" << std::endl;
                #endif
              }
            }
            else if( lo == "long int" ){
              *fos << "Trilinos does not currently support LO=" << lo << std::endl;
            }
            else if( lo == "long long int" ){
              *fos << "Trilinos does not currently support LO=" << lo << std::endl;
            }
  #else // NOT HAVE_TPETRA_INST_COMPLEX_DOUBLE
            *fos << "Scalar=complex<double> was not enabled at configure time"
                 << std::endl;
  #endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE
          } // end mag == "double"

          // mfh 15 Jan 2019: The code used to test
          // Scalar=std::complex<dd_real> and std::complex<qd_real>
          // here.  The C++ Standard does not permit portable use of
          // std::complex<T> for T != float, double, or long double.
          // Thus, I deleted those cases, which almost certainly were
          // never tested, even if they ever built.

        } // scalar == "complex"
      }

      if( !test_done && verbosity > 1 ){
        *fos << "type parameters not recognized or enabled" << std::endl;
      }
    } // else there are no runs given; no default cases for kokkos
  }

  return( success );
}
