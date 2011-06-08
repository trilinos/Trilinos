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
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_DefaultPlatform.hpp>
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

using Teuchos::RCP;
using Teuchos::ETransp;
using Teuchos::TRANS;
using Teuchos::NO_TRANS;
using Teuchos::CONJ_TRANS;
using Teuchos::ParameterList;
using Teuchos::Time;
using Teuchos::TimeMonitor;

/*
 * An example input xml file can be found in the default: "solvers_test.xml"
 */

// TODO: flush out the timers
RCP<Time> total_timer = TimeMonitor::getNewTimer("total time");
RCP<Teuchos::FancyOStream> fos; // for global output
int verbosity = 0;
std::string filedir = "";

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
typedef typename Platform::NodeType DefaultNode;

int main(int argc, char*argv[])
{
  Teuchos::GlobalMPISession mpisession(&argc,&argv,&std::cout);

  TimeMonitor TotalTimer(*total_timer);

  Tpetra::DefaultPlatform::DefaultPlatformType& platform
    = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  int rank = comm->getRank();
  int root = 0;

  string xml_file("solvers_test.xml"); // default xml file
  bool allprint = false;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_params", &xml_file, "XML Parameters file");
  cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  cmdp.setOption("filedir", &filedir, "Directory to search for matrix files");
  cmdp.setOption("verbosity", &verbosity, "Set verbosity level of output");
  Teuchos::CommandLineProcessor::EParseCommandLineReturn ret = cmdp.parse(argc,argv);
  if( ret != Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ){
    if( ret != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
      throw std::runtime_error("Error parsing command-line.");
    } else if( ret == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION ){
      throw std::invalid_argument("Unrecognized option given");
    }
  }

  // set up output streams based on command-line parameters
  Teuchos::oblackholestream blackhole;
  std::ostream& out = ( (allprint || (rank == root)) ? std::cout : blackhole );
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));


  //Read the contents of the xml file into a ParameterList.
  if( verbosity > 0 ){
    *fos << "Every proc reading parameters from xml_file: "
	 << xml_file << std::endl;
  }
  ParameterList test_params =
    Teuchos::ParameterXMLFileReader(xml_file).getParameters();

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
  *fos << "End Result: ";
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
    *fos << "Test matrix " << mm_file << " ... ";
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

  ParameterList solve_params;
  if( parameters.isSublist("all_solver_params") ){
    solve_params = parameters.sublist("all_solver_params");
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
	/* The <int,int> template parameters for the Factory is a kludge;
	 * we shouldn't have to provide template parameters if we're just
	 * asking for support of a solver.  This will be fixed when we
	 * move to nonmember constructors
	 * (i.e. Amesos::create<MAT,MV>("solver"); )
	 */
	if( Amesos::query(solver_name) ){
	  // then we have support for this solver

	  if( verbosity > 1) *fos << "  - with " << solver_name << " : " << std::endl;

	  ParameterList test_params = Teuchos::getValue<ParameterList>(parameters.entry(solver_it));
	  bool solver_success = test_mat_with_solver(mm_file, solver_name, test_params, solve_params);

	  if( verbosity > 1 ){
	    *fos << "  - Testing with " << solver_name << " ";
	    if( solver_success ) *fos << "passed";
	    else *fos << "failed";
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
    solve_params.setParameters( test_params.sublist("solver_params") );
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
	  *fos << "EpetraExt must be enabled for testing of epetra objects.  Skipping this test."
	       << std::endl;
	}
#endif
      } else if( object_name == "tpetra" ){
	const ParameterList tpetra_runs = Teuchos::getValue<ParameterList>(test_params.entry(object_it));
	success &= test_tpetra(mm_file, solver_name, tpetra_runs, solve_params);
      } else {
	*fos << "Linear algebra objects of " << object_name << " are not supported." << std::endl;
      }
    }
  }

  return( success );
}

#ifdef HAVE_AMESOS2_EPETRAEXT
bool do_epetra_test(const string& mm_file,
		    const string& solver_name,
		    ParameterList solve_params)
{
  using Teuchos::ScalarTraits;

  if( verbosity > 1 ){
    *fos << "    Doing epetra test... " << std::endl;
  }

  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;
  typedef ScalarTraits<double> ST;
  typedef typename ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const size_t numVecs = 5;     // arbitrary number

  bool transpose = false;
  if( solve_params.isParameter("transpose") ){
    if( solve_params.get<bool>("transpose") ){
      transpose = true;
    }
  }

#ifdef HAVE_MPI
  const Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  const Epetra_SerialComm comm;
#endif

  if( verbosity > 2 ){
    *fos << "      Reading matrix from " << mm_file << std::endl;
  }
  std::string path = filedir + mm_file;
  MAT* A;
  int ret = EpetraExt::MatrixMarketFileToCrsMatrix(path.c_str(), comm, A, false, false);
  if( ret == -1 ){
    *fos << "error reading file from disk, aborting run." << std::endl;
    return( false );
  }

  const Epetra_Map dmnmap = A->DomainMap();
  const Epetra_Map rngmap = A->RangeMap();

  ETransp trans;
  RCP<MV> X, B, Xhat;
  if( transpose ){
    trans = CONJ_TRANS;
    X = rcp(new MV(dmnmap,numVecs));
    B = rcp(new MV(rngmap,numVecs));
    Xhat = rcp(new MV(dmnmap,numVecs));
  } else {
    trans = NO_TRANS;
    X = rcp(new MV(rngmap,numVecs));
    B = rcp(new MV(dmnmap,numVecs));
    Xhat = rcp(new MV(rngmap,numVecs));
  }
  X->SetLabel("X");
  B->SetLabel("B");
  Xhat->SetLabel("Xhat");
  X->Random();
  A->Multiply(transpose, *X, *B);

  RCP<MAT> A_rcp(A);
  RCP<Amesos::SolverBase> solver
    = Amesos::create<MAT,MV>(solver_name, A_rcp, Xhat, B );

  solver->setParameters( rcpFromRef(solve_params) );
  try {
    solver->symbolicFactorization().numericFactorization().solve();
  } catch ( std::exception e ){
    *fos << "Exception encountered during solution." << std::endl;
    return( false );
  }
  if( verbosity > 2 ){
    *fos << "    Solution achieved. Checking solution ... ";
  }

  Teuchos::Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
  Xhat->Norm2(xhatnorms.getRawPtr());
  X->Norm2(xnorms.getRawPtr());

  // Set up an appropriate output stream for the comparison function
  RCP<Teuchos::FancyOStream> compare_fos; // for global output
  if( verbosity > 3 ){
    compare_fos = fos;
  } else {
    Teuchos::oblackholestream blackhole;
    compare_fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(blackhole));
  }

  const bool result = Teuchos::compareFloatingArrays(xhatnorms, "xhatnorms",
						     xnorms, "xnorms",
						     0.005, *compare_fos);
  if (!result) {
    if( verbosity > 2 ){
      *fos << "failed!" << std::endl;
    }
    return( false );
  } else {
    if( verbosity > 2 ){
      *fos << "passed!" << std::endl;
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
      if( run_list.isSublist("run_params") ){
	solve_params.setParameters( run_list.sublist("run_params") );
	success &= do_epetra_test(mm_file, solver_name, solve_params);
      } else {
	do_default = true;
      }
    } else {
      do_default = true;
    }
  }

  // only do one default run
  if( do_default ){
    success &= do_epetra_test(mm_file, solver_name, solve_params);
  }

  return( success );
}
#endif	// HAVE_AMESOS2_EPETRAEXT

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
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const size_t numVecs = 5;     // arbitrary number

  bool transpose = false;
  if( solve_params.isParameter("transpose") ){
    if( solve_params.get<bool>("transpose") ){
      transpose = true;
    }
  }

  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Comm<int> > comm = platform.getComm();
  RCP<Node>             node = platform.getNode();

  if( verbosity > 2 ){
    *fos << "      Reading matrix from " << mm_file << std::endl;
  }
  std::string path = filedir + mm_file;
  RCP<MAT> A =
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path, comm, node);

  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > dmnmap = A->getDomainMap();
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rngmap = A->getRangeMap();

  ETransp trans;
  RCP<MV> X, B, Xhat;
  if( transpose ){
    trans = CONJ_TRANS;
    X = rcp(new MV(dmnmap,numVecs));
    B = rcp(new MV(rngmap,numVecs));
    Xhat = rcp(new MV(dmnmap,numVecs));
  } else {
    trans = NO_TRANS;
    X = rcp(new MV(rngmap,numVecs));
    B = rcp(new MV(dmnmap,numVecs));
    Xhat = rcp(new MV(rngmap,numVecs));
  }
  X->setObjectLabel("X");
  B->setObjectLabel("B");
  Xhat->setObjectLabel("Xhat");
  X->randomize();
  A->apply(*X,*B,trans);

  RCP<Amesos::SolverBase> solver
    = Amesos::create<MAT,MV>(solver_name, A, Xhat, B );

  solver->setParameters( rcpFromRef(solve_params) );
  try {
    solver->symbolicFactorization().numericFactorization().solve();
  } catch ( std::exception e ){
    *fos << "Exception encountered during solution." << std::endl;
    return( false );
  }
  if( verbosity > 2 ){
    *fos << "    Solution achieved. Checking solution ... ";
  }

  Teuchos::Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
  Xhat->norm2(xhatnorms());
  X->norm2(xnorms());

  // Set up an appropriate output stream for the comparison function
  RCP<Teuchos::FancyOStream> compare_fos; // for global output
  if( verbosity > 3 ){
    compare_fos = fos;
  } else {
    Teuchos::oblackholestream blackhole;
    compare_fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(blackhole));
  }

  const bool result = Teuchos::compareFloatingArrays(xhatnorms, "xhatnorms",
						     xnorms, "xnorms",
						     0.005, *compare_fos);
  if (!result) {
    if( verbosity > 2 ){
      *fos << "failed!" << std::endl;
    }
    return( false );
  } else {
    if( verbosity > 2 ){
      *fos << "passed" << std::endl;
    }
    return( true );
  }
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
      if( run_list.isSublist("run_params") ){
	solve_params.setParameters( run_list.sublist("run_params") );
      }

      string scalar, mag, lo, go, node;
      if( run_list.isParameter("Scalar") ){
	scalar = run_list.get<string>("Scalar");
	if( scalar == "complex" ){
	  // get the magnitude parameter
	  if( run_list.isParameter("Magnitude") ){
	    mag = run_list.get<string>("Magnitude");
	  } else {
	    *fos << "Must provide a type for `Magnitude' when Scalar='complex', aborting run `"
		 << run_list.name() << "'..." << std::endl;
	  }
	}
      } else {
	*fos << "Must provide a type for `Scalar', aborting run `"
	     << run_list.name() << "'..." << std::endl;
      }
      if( run_list.isParameter("LocalOrdinal") ){
	lo = run_list.get<string>("LocalOrdinal");
      } else {
	*fos << "Must provide a type for `LocalOrdinal', aborting run `"
	     << run_list.name() << "'..." << std::endl;
      }
      if( run_list.isParameter("GlobalOrdinal") ){
	go = run_list.get<string>("GlobalOrdinal");
      } else {
	go = "default";
      }
      if( run_list.isParameter("Scalar") ){
	scalar = run_list.get<string>("Scalar");
      } else {
	node = "default";
      }

      if( verbosity > 1 ){
	*fos << "    Doing tpetra test with"
	     << " scalar=" << scalar;
	if( scalar == "complex" ){
	  *fos << "(" << mag << ")";
	}
	*fos << " lo=" << lo
	     << " go=" << go
	     << std::endl;
      }

      string timer_name = mm_file + "_" + scalar + "_" + lo + "_" + go;
      RCP<Time> timer = TimeMonitor::getNewTimer(timer_name);
      TimeMonitor LocalTimer(*timer);

      // I can't think of any better way to do the types as
      // specified in the parameter list at runtime but to branch
      // out on all possibilities, sorry...
      // Note: we're going to ignore the `node' parameter for now
      if( scalar == "float" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<float,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    success &= do_tpetra_test_with_types<float,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<float,int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<float,int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<float,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<float,long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<float,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<float,long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<float,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<float,long long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<float,long long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<float,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	}
#endif
      }
      else if( scalar == "double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<double,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    success &= do_tpetra_test_with_types<double,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<double,int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<double,int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<double,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<double,long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<double,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<double,long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<double,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<double,long long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<double,long long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<double,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	}
#endif
      } // end scalar == "double"
#ifdef HAVE_TEUCHOS_QD
      else if( scalar == "double double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<dd_real,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    success &= do_tpetra_test_with_types<dd_real,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<dd_real,int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<dd_real,int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<dd_real,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<dd_real,long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<dd_real,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<dd_real,long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<dd_real,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<dd_real,long long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<dd_real,long long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<dd_real,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	}
#endif
      } // end scalar == "double double"
      else if( scalar == "quad" || scalar == "quad double" ){
	if( lo == "int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<qd_real,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    success &= do_tpetra_test_with_types<qd_real,int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<qd_real,int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<qd_real,int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
	else if( lo == "long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<qd_real,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<qd_real,long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    success &= do_tpetra_test_with_types<qd_real,long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<qd_real,long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
#endif
	}
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	else if( lo == "long long int" ){
	  if( go == "default" ){
	    success &= do_tpetra_test_with_types<qd_real,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<qd_real,long long int,int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long int" ){
	    *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	    // success &= do_tpetra_test_with_types<qd_real,long long int,long int,DN>(mm_file,solver_name,solve_params);
	  }
	  else if( go == "long long int" ){
	    success &= do_tpetra_test_with_types<qd_real,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	  }
	}
#endif
      } // end scalar == "quad double"
#endif    // HAVE_TEUCHOS_QD
#ifdef HAVE_TEUCHOS_COMPLEX
      else if( scalar == "complex" ){
	if( mag == "float" ){
	  typedef std::complex<float> cmplx;
	  if( lo == "int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      success &= do_tpetra_test_with_types<cmplx,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx,int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx,int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx,long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx,long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx,long long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx,long long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	  }
#endif
	}
	else if( mag == "double" ){
	  typedef std::complex<double> cmplx_double;
	  if( lo == "int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_double,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_double,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_double,long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_double,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_double,long long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_double,long long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_double,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	  }
#endif
	} // end scalar == "double"
#ifdef HAVE_TEUCHOS_QD
	else if( mag == "double double" ){
	  typdef std::complex<dd_real> cmplx_dd;
	  if( lo == "int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_dd,long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_dd,long long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_dd,long long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_dd,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	  }
#endif
	} // end scalar == "double double"
	else if( mag == "quad" || mag == "quad double" ){
	  typdef std::complex<qd_real> cmplx_qd;
	  if( lo == "int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
	  else if( lo == "long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_qd,long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
#endif
	  }
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
	  else if( lo == "long long int" ){
	    if( go == "default" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_qd,long long int,int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long int" ){
	      *fos << "May not have global ordinal with size smaller than local ordinal" << std::endl;
	      // success &= do_tpetra_test_with_types<cmplx_qd,long long int,long int,DN>(mm_file,solver_name,solve_params);
	    }
	    else if( go == "long long int" ){
	      success &= do_tpetra_test_with_types<cmplx_qd,long long int,long long int,DN>(mm_file,solver_name,solve_params);
	    }
	  }
#endif
	} // end scalar == "quad double"
#endif    // HAVE_TEUCHOS_QD
      }
#endif  // HAVE_TEUCHOS_COMPLEX
    } // else : no default cases for tpetra
  }

  return( success );
}

