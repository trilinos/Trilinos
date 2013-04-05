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

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp> // For reading matrix-market files

#include "Amesos2.hpp"          // includes everything from Amesos2

#ifdef HAVE_AMESOS2_EPETRAEXT
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
#endif	// HAVE_AMESOS2_EPETRAEXT


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

#ifdef HAVE_AMESOS2_EPETRAEXT
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


typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Platform::NodeType DefaultNode;

int main(int argc, char*argv[])
{
  Teuchos::GlobalMPISession mpisession(&argc,&argv,&std::cout);

  TimeMonitor TotalTimer(*total_timer);

  Tpetra::DefaultPlatform::DefaultPlatformType& platform
    = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();

  int root = 0;

  string xml_file("solvers_test.xml"); // default xml file
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
  } catch (Teuchos::CommandLineProcessor::HelpPrinted hp) {
    return EXIT_SUCCESS;	// help was printed, exit gracefully.
  }

  // set up output streams based on command-line parameters
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if( !allprint ) fos->setOutputToRootOnly( root );

  if( verbosity > 3 ){
    compare_fos = fos;
  } else {
    Teuchos::oblackholestream blackhole;
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

bool test_mat_with_solver(const string& mm_file,
			  const string& solver_name,
			  const ParameterList& test_params,
			  ParameterList solve_params)
{
  bool success = true;

  if( test_params.isSublist("solver_params") ){
    solve_params.sublist(solver_name) = test_params.sublist("solver_params");
  }

  ParameterList::ConstIterator object_it;
  for( object_it = test_params.begin(); object_it != test_params.end(); ++object_it ){
    if( ! test_params.entry(object_it).isUsed() ){
      string object_name = test_params.name(object_it);

      // There may be a more flexible way (e.g. a `query' function) to do this check
      if( object_name == "epetra" ){
#ifdef HAVE_AMESOS2_EPETRAEXT
	const ParameterList epetra_runs = Teuchos::getValue<ParameterList>(test_params.entry(object_it));
	success &= test_epetra(mm_file, solver_name, epetra_runs, solve_params);
#else
	if( verbosity > 1 ){
	  *fos << "    EpetraExt must be enabled for testing of epetra objects.  Skipping this test."
	       << std::endl;
	}
#endif
      } else if( object_name == "tpetra" ){
	const ParameterList tpetra_runs = Teuchos::getValue<ParameterList>(test_params.entry(object_it));
	success &= test_tpetra(mm_file, solver_name, tpetra_runs, solve_params);
      } else {
	*fos << "    Linear algebra objects of " << object_name << " are not supported." << std::endl;
      }
    }
  }

  return( success );
}


template <class Vector>
struct solution_checker {
  bool operator()(RCP<Vector> true_solution, RCP<Vector> solution)
  {
    return false;
  }
};

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

#ifdef HAVE_AMESOS2_EPETRAEXT
template <>
struct solution_checker<Epetra_MultiVector> {
  bool operator()(RCP<Epetra_MultiVector> true_solution, RCP<Epetra_MultiVector> given_solution)
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
  typedef typename ArrayView<const RCP<Vector> >::iterator rhs_it_t;
  bool success = true;		// prove me wrong!

  solution_checker<Vector> checker;
  RCP<Amesos2::Solver<Matrix,Vector> > solver;

  int phase = Amesos2::CLEAN;	// start with no computation

  while( phase != Amesos2::SOLVE ) {
    // and X and B will never be needed, so we create the solver without them
    solver = Amesos2::create<Matrix,Vector>(solver_name, A1);
    switch( phase ){
    case Amesos2::CLEAN: break;
    case Amesos2::PREORDERING:
      solver->preOrdering(); break;
    case Amesos2::SYMBFACT:
      solver->symbolicFactorization(); break;
    case Amesos2::NUMFACT:
      solver->numericFactorization();
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

  while( style <= SOLVE_SHORT ){
    rhs_it_t rhs_it = b.begin();
    rhs_it_t x_it = x.begin();

    // Create our solver according to the current style
    switch( style ){
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
    if( !success ) return success; // bail out early if necessary

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


#ifdef HAVE_AMESOS2_EPETRAEXT

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
  const size_t numRHS  = 5;	// also quite arbitrary

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
    int ret = EpetraExt::MatrixMarketFileToCrsMatrix(path.c_str(), comm,
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
#endif	// HAVE_AMESOS2_EPETRAEXT

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

  typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MAT;
  typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  const size_t numVecs = 5;     // arbitrary number
  const size_t numRHS = 5;	  // also arbitrary

  bool transpose = solve_params.get<bool>("Transpose", false);

  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Comm<int> > comm = platform.getComm();
  RCP<Node>             node = platform.getNode();

  if( verbosity > 2 ){
    *fos << std::endl << "      Reading matrix from " << mm_file << " ... " << std::flush;
  }
  std::string path = filedir + mm_file;
  RCP<MAT> A =
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path, comm, node);

  if( verbosity > 2 ){
    *fos << "done" << std::endl;
    switch( verbosity ){
    case 6:
      A->describe(*fos, Teuchos::VERB_EXTREME); break;
    case 5:
      A->describe(*fos, Teuchos::VERB_HIGH); break;
    case 4:
      A->describe(*fos, Teuchos::VERB_LOW); break;
    }
  }

  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > dmnmap = A->getDomainMap();
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rngmap = A->getRangeMap();

  ETransp trans = transpose ? CONJ_TRANS : NO_TRANS;
    
  RCP<MAT> A2;
  RCP<MV> Xhat, x2, b2;
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
    
  if( refactor ){
    // There isn't a really nice way to get a deep copy of an entire
    // CrsMatrix, so we just read the file again.
    A2 = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path, comm, node);

    // perturb the values just a bit (element-wise square of first row)
    size_t l_fst_row_nnz = A2->getNumEntriesInLocalRow(0);
    Array<LocalOrdinal> indices(l_fst_row_nnz);
    Array<Scalar> values(l_fst_row_nnz);
    A2->getLocalRowCopy(0, indices, values, l_fst_row_nnz);
    for( size_t i = 0; i < l_fst_row_nnz; ++i ){
      values[i] = values[i] * values[i];
    }
    A2->resumeFill();
    A2->replaceLocalValues(0, indices, values);
    A2->fillComplete();

    A2->apply(*x2, *b2, trans);
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
#endif	// HAVE_TEUCHOS_COMPLEX
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

#define AMESOS2_SOLVER_TPETRA_TEST(S,LO,GO,N)				\
      test_done = true;							\
      bool run_success = do_tpetra_test_with_types<S,LO,GO,N>(mm_file,solver_name, \
							      solve_params_copy); \
      if( verbosity > 1 ){						\
	if (!run_success)						\
	  *fos << "failed!" << std::endl;				\
	else								\
	  *fos << "passed" << std::endl;				\
      }									\
      success &= run_success

	
      // I can't think of any better way to do the types as
      // specified in the parameter list at runtime but to branch
      // out on all possibilities, sorry...  Note: we're going to
      // ignore the `node' parameter for now
#if !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) || ((defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) && (defined HAVE_TPETRA_INST_FLOAT))
      if( scalar == "float" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,int,int,DN);
	  }
	  else if( go == "int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,int,int,DN);
	  }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,int,long long int,DN);
	  }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	}
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	else if( lo == "long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,long int,long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,long int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,long int,long long int,DN);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,long long int,long long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(float,long long int,long long int,DN);
	  }
	}
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
      } // end scalar == "float"
#endif	// HAVE_TPETRA_INST_FLOAT
#if !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) || ((defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) && (defined HAVE_TPETRA_INST_DOUBLE))
      if( scalar == "double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,int,int,DN);
	  }
	  else if( go == "int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,int,int,DN);
	  }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,int,long long int,DN);
	  }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	}
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	else if( lo == "long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,long int,long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,long int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,long int,long long int,DN);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,long long int,long long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(double,long long int,long long int,DN);
	  }
	}
#endif
#endif	  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
      } // end scalar == "double"
#endif	  // HAVE_TPETRA_INST_DOUBLE
#if (defined HAVE_TEUCHOS_QD) && !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION)
      if( scalar == "double double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,int,DN);
	  }
	  else if( go == "int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,int,DN);
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,int,long long int,DN);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,long int,long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,long int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,long int,long long int,DN);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,long long int,long long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(dd_real,long long int,long long int,DN);
	  }
	}
#endif
      } // end scalar == "double double"
      if( scalar == "quad" || scalar == "quad double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,int,DN);
	  }
	  else if( go == "int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,int,DN);
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,int,long long int,DN);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,long int,long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,long int,long int,DN);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,long int,long long int,DN);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,long long int,long long int,DN);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	  }
	  else if( go == "long long int" ){
	    AMESOS2_SOLVER_TPETRA_TEST(qd_real,long long int,long long int,DN);
	  }
	}
#endif
      } // end scalar == "quad double"
#endif    // HAVE_TEUCHOS_QD
#ifdef HAVE_TEUCHOS_COMPLEX
      if( scalar == "complex" ){
#if !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) || ((defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) && (defined HAVE_TPETRA_INST_COMPLEX_FLOAT))
	if( mag == "float" ){
	  typedef std::complex<float> cmplx;
	  if( lo == "int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,int,DN);
	    }
	    else if( go == "int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,int,DN);
	    }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,int,long long int,DN);
	    }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,long int,long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,long int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,long int,long long int,DN);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,long long int,long long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx,long long int,long long int,DN);
	    }
	  }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	}
#endif
#if !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) || ((defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION) && (defined HAVE_TPETRA_INST_COMPLEX_DOUBLE))
	if( mag == "double" ){
	  typedef std::complex<double> cmplx_double;
	  if( lo == "int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,int,DN);
	    }
	    else if( go == "int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,int,DN);
	    }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,int,long long int,DN);
	    }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  }
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,long int,long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,long int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,long int,long long int,DN);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,long long int,long long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_double,long long int,long long int,DN);
	    }
	  }
#endif
#endif	// HAVE_AMESOS2_EXPLICIT_INSTANTIATION
	} // end mag == "double"
#endif	     
#if (defined HAVE_TEUCHOS_QD) && !(defined HAVE_AMESOS2_EXPLICIT_INSTANTIATION)
	if( mag == "double double" ){
	  typedef std::complex<dd_real> cmplx_dd;
	  if( lo == "int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,int,int,DN);
	    }
	    else if( go == "int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,int,int,DN);
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,int,long long int,DN);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,long int,long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,long int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,long int,long long int,DN);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,long long int,long long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_dd,long long int,long long int,DN);
	    }
	  }
#endif
	} // end scalar == "double double"
	else if( mag == "quad" || mag == "quad double" ){
	  typedef std::complex<qd_real> cmplx_qd;
	  if( lo == "int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,int,int,DN);
	    }
	    else if( go == "int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,int,int,DN);
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,int,long long int,DN);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,long int,long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,long int,long int,DN);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,long int,long long int,DN);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,long long int,long long int,DN);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    }
	    else if( go == "long long int" ){
	      AMESOS2_SOLVER_TPETRA_TEST(cmplx_qd,long long int,long long int,DN);
	    }
	  }
#endif
	} // end scalar == "quad double"
#endif    // HAVE_TEUCHOS_QD
      } // end if( scalar == "complex" )
#endif  // HAVE_TEUCHOS_COMPLEX

      if( !test_done && verbosity > 1 ){
	*fos << "type parameters not recognized or enabled" << std::endl;
      }
    } // else there are no runs given; no default cases for tpetra
  }

  return( success );
}
