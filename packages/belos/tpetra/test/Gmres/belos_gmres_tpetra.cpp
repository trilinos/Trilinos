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

/// \file belos_gmres_tpetra.cpp
/// \brief Test GmresSolMgr with Tpetra objects.
///
#include <BelosGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace {

  /// \fn getNode
  /// \brief Return an RCP to a Kokkos Node
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode (Teuchos::RCP<Teuchos::ParameterList> params) {
    throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
  }

  // Specialization of getNode for SerialNode
  template<>
  Teuchos::RCP<Kokkos::SerialNode>
  getNode (Teuchos::RCP<Teuchos::ParameterList> params) {
    // "Num Threads" specifies the number of threads.  Defaults to an
    // automatically chosen value.
    if (params.is_null())
      params = Teuchos::parameterList ();

    return Teuchos::rcp (new Kokkos::SerialNode (*params));
  }

#if defined(HAVE_KOKKOSCLASSIC_TBB)
  // Specialization of getNode for TBBNode
  template<>
  Teuchos::RCP<Kokkos::TBBNode>
  getNode (Teuchos::RCP<Teuchos::ParameterList> params) {
    // "Num Threads" specifies the number of threads.  Defaults to an
    // automatically chosen value.
    if (params.is_null())
      params = Teuchos::parameterList ();

    return Teuchos::rcp (new Kokkos::TBBNode (*params));
  }
#endif // defined(HAVE_KOKKOSCLASSIC_TBB)

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  // Specialization of getNode for TPINode
  template<>
  Teuchos::RCP<Kokkos::TPINode>
  getNode (Teuchos::RCP<Teuchos::ParameterList> params) {
    using Teuchos::isParameterType;

    // "Num Threads" (defaults to 0) specifies the number of threads,
    // and "Verbose" specifies verbosity (defaults to 0, but we set 1
    // as the default, so that you can see how many threads are being
    // used if you don't set a specific number.
    if (params.is_null())
      {
	params = Teuchos::parameterList ();
	int verbosity = 1;
	params->set ("Verbose", verbosity);
      }
    else if (isParameterType<int>(*params, "Num Threads") && params->get<int>("Num Threads") == -1)
      params->set ("Num Threads", static_cast<int>(0));

    return Teuchos::rcp (new Kokkos::TPINode (*params));
  }
#endif // defined(HAVE_KOKKOSCLASSIC_THREADPOOL)

  /// \brief Show MsgType as comma-delimited list of names.
  ///
  /// Belos::MsgType claims to be an enum, but is really a C-style bit
  /// set (where you bitwise OR together different names to get a
  /// combination of values).  This function returns a string
  /// representing the given MsgType (represented here as an int,
  /// mainly because Teuchos::ParameterList seems to prefer storing
  /// these kind of C-style bit sets as int rather than MsgType) as a
  /// comma-delimited, human-readable list of names.  This is useful
  /// for debugging.
  /// 
  std::string 
  msgTypeToString (const int msgType)
  {
    using std::ostringstream;
    using std::vector;
    typedef vector<int>::size_type size_type;

    // Wouldn't it be nice if C++ enums had introspection and could
    // be enumerated?
    const size_type numValidTypes = 8;
    const int validTypes[] = {
      Belos::Errors, 
      Belos::Warnings, 
      Belos::IterationDetails,
      Belos::OrthoDetails,
      Belos::FinalSummary,
      Belos::TimingDetails,
      Belos::StatusTestDetails,
      Belos::Debug
    };
    const char* typeNames[] = {
      "Errors", 
      "Warnings", 
      "IterationDetails",
      "OrthoDetails",
      "FinalSummary",
      "TimingDetails",
      "StatusTestDetails",
      "Debug"
    };

    // We first generate a list, and only then build a single string.
    // This helps us decide where to put the commas.  The list just
    // uses the indices of the valid names, rather than the valid
    // names themselves, in order to save space and time.  We use
    // size_type for the indices to avoid signed/unsigned comparisons.
    vector<size_type> theList;
    for (size_type nameIndex = 0; nameIndex < numValidTypes; ++nameIndex)
      {
	if (msgType & validTypes[nameIndex])
	  theList.push_back (nameIndex);
      }
    ostringstream os;
    for (size_type k = 0; k < theList.size(); ++k)
      {
	const size_type nameIndex = theList[k];
	os << typeNames[nameIndex];
	if (nameIndex < theList.size() - 1)
	  os << ",";
      }
    return os.str();
  }

  /// \class ProblemMaker
  /// \brief Make a Tpetra sparse linear problem to solve.
  template<class SparseMatrixType>
  class ProblemMaker {
  public:
    typedef SparseMatrixType sparse_matrix_type;
    typedef typename SparseMatrixType::scalar_type scalar_type;
    typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename SparseMatrixType::node_type node_type;
    typedef Tpetra::MultiVector<scalar_type, 
				local_ordinal_type, 
				global_ordinal_type, 
				node_type> multivector_type;

  private:
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    Teuchos::RCP<node_type> node_;
    const bool tolerant_;
    const bool debug_;

    /// Read in a Tpetra::CrsMatrix from a Matrix Market file.
    ///
    /// \param matrixFilename [in] Name of the Matrix Market format sparse
    ///   matrix file to read (on MPI Rank 0 only)
    /// \param pComm [in] Communicator, over whose MPI ranks to
    ///   distribute the returned Tpetra::CrsMatrix.
    /// \param tolerant [in] Whether or not to parse the file 
    ///   tolerantly
    ///
    /// \return Tpetra::CrsMatrix read in from the file 
    Teuchos::RCP<sparse_matrix_type>
    readSparseFile (const std::string& matrixFilename)
    {
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
      const bool callFillComplete = true;
      return reader_type::readSparseFile (matrixFilename, comm_, node_,
					  callFillComplete, tolerant_, debug_);
    }

    /// \brief Generate, distribute, and return a test sparse matrix.
    ///
    /// \param globalNumRows [in] Global (over all MPI process(es))
    ///   number of rows in the sparse matrix.
    ///
    /// \param symmetric [in] Whether to generate a symmetric test problem.
    ///   Storage is nonsymmetric regardless; symmetry here only applies to 
    ///   the entries' locations and values.
    ///
    /// \return The sparse matrix (global, distributed)
    Teuchos::RCP<sparse_matrix_type>
    generateTestMatrix (const global_ordinal_type globalNumRows,
			const bool symmetric)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::tuple;
      using Tpetra::createUniformContigMapWithNode;
      using std::endl;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef node_type NT;

      // For a square matrix, we only need a Map for the range of the matrix.
      RCP<const map_type> pRangeMap = 
	createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm_, node_);
      // The sparse matrix object to fill.
      RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

      const int myRank = comm_->getRank();
      if (myRank == 0) {
	const scalar_type leftVal = -STS::one();
	const scalar_type rightVal = symmetric ? -STS::one() : +2*STS::one();
	// Boost the diagonal if nonsymmetric, hopefully to make
	// convergence a little faster.
	const scalar_type centerVal = symmetric ? 2*STS::one() : 8*STS::one(); 

	for (GO curRow = 0; curRow < globalNumRows; ++curRow) {
	  if (curRow > 0) {
	    pMat->insertGlobalValues (curRow, tuple(curRow-1), tuple(leftVal));
	  }
	  pMat->insertGlobalValues (curRow, tuple(curRow), tuple(centerVal));
	  if (curRow < globalNumRows-1) {
	    pMat->insertGlobalValues (curRow, tuple(curRow+1), tuple(rightVal));
	  }
	}
      }
      // Make sure Rank 0 is done filling in the matrix.
      Teuchos::barrier (*comm_);
      // fillComplete() doesn't need any arguments if the domain and
      // range maps are the same.
      pMat->fillComplete ();
      return pMat;
    }

    Teuchos::RCP<sparse_matrix_type>
    generateGmresUnfriendlyMatrix (const global_ordinal_type globalNumRows)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::tuple;
      using Tpetra::createUniformContigMapWithNode;
      using std::endl;

      typedef typename SparseMatrixType::local_ordinal_type LO;
      typedef typename SparseMatrixType::global_ordinal_type GO;
      typedef typename SparseMatrixType::node_type NT;

      // For a square matrix, we only need a Map for the range of the matrix.
      RCP<const map_type> pRangeMap = 
	createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm_, node_);
      // The sparse matrix object to fill.
      RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

      const int myRank = comm_->getRank();
      if (myRank == 0) {
	const scalar_type val = STS::one();
	for (GO curRow = 0; curRow < globalNumRows; ++curRow) {
	  const GO curCol = (curRow == 0) ? (globalNumRows-1) : (curRow-1);
	  pMat->insertGlobalValues (curRow, tuple(curCol), tuple(val));
	}
      }
      // Make sure Rank 0 is done filling in the matrix.
      Teuchos::barrier (*comm_);
      // fillComplete() doesn't need any arguments if the domain and
      // range maps are the same.
      pMat->fillComplete ();
      return pMat;
    }

  public:
    ProblemMaker (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		  const Teuchos::RCP<node_type>& node,
		  const bool tolerant,
		  const bool debug) :
      comm_ (comm), node_ (node), tolerant_ (tolerant), debug_ (debug)
    {}

    Teuchos::RCP<sparse_matrix_type>
    makeMatrixFromFile (const std::string& inMatrixFilename)
    {
      using Teuchos::RCP;
      using std::endl;

      const int myRank = comm_->getRank ();
      Teuchos::oblackholestream blackHole;
      std::ostream& err = (debug_ && myRank == 0) ? std::cerr : blackHole; 

      // Read the sparse matrix A from the file.
      err << "Reading sparse matrix A from file \""
	  << inMatrixFilename << "\"...";
      RCP<sparse_matrix_type> A = readSparseFile (inMatrixFilename);
      err << "done." << endl;
      return A;
    }

    Teuchos::RCP<sparse_matrix_type>
    makeMatrix (const std::string& inMatrixFilename,
		const std::string& outMatrixFilename,
		const int globalNumRows,
		const bool generated,
		const bool symmetric,
		const bool gmresUnfriendly)
    {
      Teuchos::RCP<sparse_matrix_type> A;
      if (inMatrixFilename != "") {
	A = makeMatrixFromFile (inMatrixFilename);
      } else if (generated) {
	if (gmresUnfriendly) {
	  A = generateGmresUnfriendlyMatrix (globalNumRows);
	} else {
	  A = generateTestMatrix (globalNumRows, symmetric);
	}
      } else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
      }

      if (outMatrixFilename != "") {
	typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
	writer_type::writeSparseFile (outMatrixFilename, A, "", "", debug_);
      }
      return A;
    }

    void
    makeVectors (const Teuchos::RCP<const sparse_matrix_type>& A,
		 Teuchos::RCP<multivector_type>& X_guess,
		 Teuchos::RCP<multivector_type>& X_exact,
		 Teuchos::RCP<multivector_type>& B,
		 const std::string& inRhsFilename,
		 const std::string& outRhsFilename,
		 const bool gmresUnfriendly)
    {
      using std::endl;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef node_type NT;
      typedef multivector_type MV;

      const int myRank = comm_->getRank ();
      Teuchos::oblackholestream blackHole;
      std::ostream& err = (debug_ && myRank == 0) ? std::cerr : blackHole; 

      // Construct X_guess (the initial guess for the solution of
      // AX=B) from the domain of the matrix A, and fill it with
      // zeros.
      err << "Constructing initial guess vector X_guess...";
      X_guess = rcp (new MV (A->getDomainMap(), 1));
      X_guess->putScalar (STS::zero());
      err << "done." << endl;

      if (inRhsFilename != "") {
	// If reading the right-hand side(s) from a file, don't set
	// the exact solution(s).  Later we could try to read those
	// from a file too.
	err << "Reading B from Matrix Market file...";
	typedef Tpetra::MatrixMarket::Reader<SparseMatrixType> reader_type;
	RCP<const map_type> map = A->getRangeMap();
	B = reader_type::readDenseFile (inRhsFilename, comm_, A->getNode(), 
					map, tolerant_, debug_);
	err << "...done." << endl;
      } else {
	// Our choice of exact solution and right-hand side depend on
	// the test problem.  If we generated the GMRES-unfriendly
	// example, we need B = e_{globalNumRows} and therefore
	// X_exact = e_{globalNumRows} as well.  Otherwise, we pick
	// X_exact first and compute B via SpMV: B = A * X_exact.
	X_exact = rcp (new MV (A->getDomainMap(), 1));

	// Construct the right-hand side B from the range of the
	// matrix A.  Don't just clone X_guess, since the range may
	// differ from the domain.
	B = rcp (new MV (A->getRangeMap(), 1));

	if (gmresUnfriendly) {
	  err << "Constructing B and X_exact for canonical \"GMRES-"
	    "unfriendly\" example...";
	  X_exact->putScalar (STS::zero());
	  X_exact->replaceGlobalValue (A->getGlobalNumRows()-1, 0, STS::one());
	  B->putScalar (STS::zero());
	  B->replaceGlobalValue (A->getGlobalNumRows()-1, 0, STS::one());
	  err << "done." << endl;
	} else {
	  // Construct the exact solution vector and fill it with all
	  // ones.  Tacky, but deterministic.  Not so good if we
	  // expect A to be singular with rigid body modes.
	  err << "Setting X_exact = [1; ...; 1]...";
	  X_exact->putScalar (STS::one());
	  err << "done." << endl;
      
	  // Compute the right-hand side B := A*X_exact.
	  err << "Computing B := A*X_exact...";
	  A->apply (*X_exact, *B);
	  err << "done." << endl;
	}
      }
    }
  };


} // (anonymous namespace)

/// \fn main
/// \brief Test driver for Belos' new GMRES implementations
int 
main (int argc, char *argv[]) 
{
  using Belos::GmresSolMgr;
  using Belos::LinearProblem;
  using Belos::OrthoManagerFactory;
  using Teuchos::CommandLineProcessor;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cerr;
  using std::cout;
  using std::endl;

  ////////////////////////////////////////////////////////////////////
  // typedefs make main() as generic as possible
  ////////////////////////////////////////////////////////////////////

  typedef double scalar_type;
  typedef int local_ordinal_type;

  // mfh 22 Feb 2011: Explicit instantiation of templates in the checkin
  // test script makes my life hard sometimes...
#if 0
#  ifdef HAVE_TEUCHOS_LONG_LONG_INT
  // Half the point of Tpetra is the capability of using global
  // ordinals that are bigger than local ordinals.  We should test
  // this whenever possible, even in code that is not specifically
  // testing Tpetra.
  typedef long long int global_ordinal_type;
#  else // not HAVE_TEUCHOS_LONG_LONG_INT
  typedef long int global_ordinal_type;
#  endif // HAVE_TEUCHOS_LONG_LONG_INT
#else // not 0
  typedef int global_ordinal_type;
#endif // 0

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
  typedef Kokkos::TPINode node_type;
#else
#  if defined(HAVE_KOKKOSCLASSIC_TBB)
  typedef Kokkos::TBBNode node_type;
#  else
  typedef Kokkos::SerialNode node_type;
#  endif // HAVE_KOKKOSCLASSIC_TBB
#endif // HAVE_KOKKOSCLASSIC_THREADPOOL

  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, 
    global_ordinal_type, node_type> MV;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, 
    global_ordinal_type, node_type> OP;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, 
    global_ordinal_type, node_type> sparse_matrix_type;
  typedef Belos::MultiVecTraits<scalar_type, MV> MVT;
  typedef Belos::OperatorTraits<scalar_type, MV, OP> OPT;

  ////////////////////////////////////////////////////////////////////
  // main() begins execution here
  ////////////////////////////////////////////////////////////////////

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // Name of the Matrix Market - format sparse matrix file from which
  // to read the matrix A.
  std::string matrixFilename;

  // Name of the Matrix Market - format dense matrix file from which
  // to read the right-hand side(s) B.
  std::string rhsFilename;

  // Whether to parse Matrix Market files tolerantly.
  bool tolerant = false;

  // Name of XML file containing parameters for the GMRES solver
  // manager.
  std::string xmlInFile;

  // Name of XML file to which to write the parameters used by the
  // GMRES solver manager.
  std::string xmlOutFile;

  // Whether, instead of reading the matrix from a file, to generate a
  // test matrix.
  bool generated = false;

  // If generating the test problem: the global number of rows (and
  // columns) in the generated matrix.  The global number of rows
  // should be >= max(3, Teuchos::size(*comm)), so that each MPI
  // process has at least one row, and the matrix is nonsingular.
  // Ignored if a sparse matrix file is given.
  int globalNumRows = std::max (3, Teuchos::size (*comm));

  // Whether the generated test problem is symmetric.
  // This option is ignored if generated is false.
  bool symmetric = false;

  // Whether to generate a "GMRES-unfriendly" test problem (the
  // classical example of a permutation matrix).
  // This option is ignored if generated is false.
  bool gmresUnfriendly = false;

  // Name of the file to which to dump the matrix.  If "", don't dump
  // the matrix to a file.
  std::string outMatrixFilename;

  // Name of the file to which to dump the right-hand side(s).  If "",
  // don't dump the right-hand side(s) to a file.
  std::string outRhsFilename;

  // Whether to print verbose output.
  bool verbose = false;

  // Whether to print debugging information.
  bool debug = false;

  // Frequency of intermediate status output.  Defaults to the invalid
  // value -2, which will be reset to -1 if not set to something else.
  // (We use this to detect whether the frequency has indeed been set
  // on the command line.)
  int frequency = -2;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("matrixFilename", &matrixFilename, "Name of a Matrix Market - "
		  "format sparse matrix file from which to read the matrix A."
		  "Do not provide if the \"--generated\" option is set.");
  cmdp.setOption ("rhsFilename", &rhsFilename, "Name of a Matrix Market - "
		  "format dense matrix file from which to read the right-hand "
		  "side(s) B.  A right-hand side will be generated if not "
		  "provided.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse Matrix Market files tolerantly.");
  cmdp.setOption ("generated", "nongenerated", &generated, 
		  "Whether to generate a test matrix.  Do not provide if you "
		  "give the \"--matrixFilename=<file>\" option.");
  cmdp.setOption ("globalNumRows", &globalNumRows, "Global number of rows in "
		  "the generated test matrix.  Ignored if the sparse matrix is "
		  "to be read from a file.");
  cmdp.setOption ("symmetric", "nonsymmetric", &symmetric, "Whether the "
		  "generated test problem is symmetric.  This option is "
		  "ignored if the \"--generated\" option is not provided.");
  cmdp.setOption ("gmresUnfriendly", "notGmresUnfriendly", &gmresUnfriendly, 
		  "Whether to generate a \"GMRES-unfriendly\" test problem "
		  "(the classical permutation matrix example).  This option is "
		  "ignored if the \"--generated\" option is not provided.");
  cmdp.setOption ("outMatrixFilename", &outMatrixFilename, "Name of the file to "
		  "which to dump the sparse matrix (in Matrix Market format).  "
		  "Ignored if the filename is empty or not provided.");
  cmdp.setOption ("outRhsFilename", &outRhsFilename, "Name of the file to "
		  "which to dump the right-hand side(s) (in Matrix Market "
		  "format).  Ignored if the filename is empty or not provided.");
  cmdp.setOption ("paramFile", &xmlInFile, "Name of an XML file containing "
		  "parameters (in Teuchos::ParameterList form) for the solve.");
  cmdp.setOption ("outParamFile", &xmlOutFile, "Name of an XML file to which to "
		  "dump the parameters used by the GMRES solver manager.  These "
		  "may be different than those given to the solver manager.");
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  cmdp.setOption ("frequency", &frequency, "Frequency of the iterative solver's "
		  "intermediate status output, in terms of number of iterations."
		  "  1 means every iteration; -1 means no intermediate status "
		  "output.  Defaults to -1.  If set, overrides the parameters "
		  "in the XML file.");
  //
  // Parse the command-line arguments.
  //
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // specify which matrix to use (from a file, or generated), we let
    // the "test" pass trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED || 
	(matrixFilename == "" && ! generated))
      {
	out << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
    if (generated)
      {
	TEUCHOS_TEST_FOR_EXCEPTION(matrixFilename != "",
			   std::invalid_argument,
			   "The --generated and \"--matrixFilename=<file>\" "
			   "options may not both be used.");
	TEUCHOS_TEST_FOR_EXCEPTION(globalNumRows < std::max(3, Teuchos::size(*comm)),
			   std::invalid_argument,
			   "The number of rows in the test matrix to generate "
			   "must be at least max(3, # MPI processes), in order "
			   "for the matrix to be nonsingular and for each MPI "
			   "process to have at least one row of the matrix.");
      }
  }

  // Now we can use the 'debug' Boolean parameter to set the debug
  // output stream.
  std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

  // We can tell if frequency was set via the "--frequency=<freq>"
  // command-line option, since the default value of frequency (-2) is
  // invalid.  Replace the default value with the valid value of -1
  // (which means, never display intermediate status output).
  const bool setFrequencyAtCommandLine = (frequency == -2);
  if (setFrequencyAtCommandLine)
    frequency = -1;

  // Default verbosity: only print warnings and errors.  This is not
  // directly a command-line option, but is affected by the --verbose
  // and --debug Boolean command-line options.
  int verbosityLevel = Belos::Errors | Belos::Warnings;

  // If the name of an XML file of parameters was provided, read it on
  // Proc 0, and broadcast from Proc 0 to the other procs.
  //
  // Start with an empty ParameterList on all procs.  If a matrixFilename
  // was provided, fill in the list on Proc 0, and broadcast.
  RCP<ParameterList> params = Teuchos::parameterList();
  const bool readParamsFromFile = (xmlInFile != "");
  if (readParamsFromFile)
    {
      err << "Reading solve parameters from XML file \"" 
	  << xmlInFile << "\"...";
      using Teuchos::updateParametersFromXmlFileAndBroadcast;
      updateParametersFromXmlFileAndBroadcast (xmlInFile, params.getRawPtr(), 
					       *comm);
      err << "done." << endl;

      // Set the verbosity level if it was specified in the XML file.
      try {
	const int verbLevel = params->get<int>("Verbosity");
	verbosityLevel = verbLevel;
      } catch (Teuchos::Exceptions::InvalidParameter&) {
	// Do nothing; leave verbosityLevel at its default value
      }

      // Did the XML file set the "Output Frequency" parameter?  If
      // so, and if the user didn't override at the command line, read
      // it in.  If provided at the command line, it overrides
      // whatever was in the XML file.
      if (! setFrequencyAtCommandLine) {
	try {
	  const int newFreq = params->get<int>("Output Frequency");
	  frequency = newFreq;
	} catch (Teuchos::Exceptions::InvalidParameter&) {
	  // Do nothing; leave frequency at its (new) default value
	}
      }
    }

  // Override the XML file's "Output Frequency" parameter in the
  // parameter list with the corresponding command-line argument, if
  // one was provided.
  if (setFrequencyAtCommandLine) {
    params->set ("Output Frequency", frequency);
  }

  // In debug mode, set verbosity level to its maximum, including
  // Belos::Debug.  In verbose mode, include everything but
  // Belos::Debug.
  if (verbose || debug)
    { // Change verbosity level from its default.
      if (verbose) {
	verbosityLevel = Belos::IterationDetails | 
	  Belos::OrthoDetails |
	  Belos::FinalSummary |
	  Belos::TimingDetails |
	  Belos::StatusTestDetails |
	  Belos::Warnings | 
	  Belos::Errors;
      } else if (debug) {
	verbosityLevel = Belos::Debug |
	  Belos::Warnings | 
	  Belos::Errors;
      }
      err << "Setting \"Verbosity\" to " << msgTypeToString(verbosityLevel) 
	  << endl;
      params->set ("Verbosity", verbosityLevel);
    }
  //
  // Construct a node, since we have parameters now.
  //
  RCP<node_type> node;
  {
    // FIXME (mfh 11 Oct 2011) Read Node parameters from XML file?
    // For now, we use default parameters.
    RCP<ParameterList> nodeParams = parameterList ("Node Parameters");
    node = getNode<node_type> (nodeParams);
    TEUCHOS_TEST_FOR_EXCEPTION(node.is_null(), std::logic_error, 
		       "Failed to initialize Kokkos Node.");
  }
  //
  // Make a sparse matrix A.  Either read it from the given file, or
  // generate it, depending on the command-line options.
  //
  ProblemMaker<sparse_matrix_type> maker (comm, node, tolerant, debug);
  // This dumps the matrix to a file, if outMatrixFilename != "".
  RCP<sparse_matrix_type> A = 
    maker.makeMatrix (matrixFilename, outMatrixFilename, globalNumRows, 
		      generated, symmetric, gmresUnfriendly);
  //
  // Construct the initial guess and the right-hand side.  If the
  // right-hand side was read from a file, we don't get X_exact.
  //
  RCP<MV> X_guess, X_exact, B;
  // This dumps B to a file, if outRhsFilename != "".
  maker.makeVectors (A, X_guess, X_exact, B, rhsFilename, 
		     outRhsFilename, gmresUnfriendly);
  // If B was read in from a file, we won't get X_exact.
  const bool haveExactSolution = ! X_exact.is_null();

  //
  // In debug mode, compute and print the initial residual
  // independently of the solver framework.  The solver may also
  // compute residuals, but the code below verifies the solver's
  // computation of residuals.
  //
  if (debug) {
    RCP<MV> R = MVT::Clone (*B, MVT::GetNumberVecs (*B));
    // R := A * X_guess
    OPT::Apply (*A, *X_guess, *R);
    // R := B - R
    MVT::MvAddMv (STS::one(), *B, -STS::one(), *R, *R);
    std::vector<magnitude_type> theNorm (MVT::GetNumberVecs (*R));
    MVT::MvNorm (*R, theNorm);

    err << "Initial residual norm: ||B - A*X_guess||_2 = " 
	<< theNorm[0] << endl;

    MVT::MvNorm (*B, theNorm);
    err << "||B||_2 = " << theNorm[0] << endl;

    if (haveExactSolution) {
      MVT::MvNorm (*X_exact, theNorm);
      const magnitude_type X_exact_norm = theNorm[0];
      err << "||X_exact||_2 = " << X_exact_norm << endl;

      RCP<MV> X_diff = MVT::CloneCopy (*X_exact);
      MVT::MvAddMv (STS::one(), *X_guess, -STS::one(), *X_diff, *X_diff);
      MVT::MvNorm (*X_diff, theNorm);

      if (X_exact_norm == STM::zero()) {
	// Don't compute a relative norm if ||X_exact|| is zero.
	err << "||X_guess - X_exact||_2 = " << theNorm[0] << endl;
      } else {
	err << "||X_guess - X_exact||_2 / ||X_exact||_2 = " 
	    << theNorm[0] / X_exact_norm << endl;
      }
    }

    MVT::MvNorm (*X_guess, theNorm);
    err << "||X_guess||_2 = " << theNorm[0] << endl;
  }

  // Wrap the linear problem to solve in a Belos::LinearProblem
  // object.  The "X" argument of the LinearProblem constructor is
  // only copied shallowly and will be overwritten by the solve, so we
  // make a deep copy here.  That way we can compare the result
  // against the original X_guess.
  RCP<MV> X = MVT::CloneCopy (*X_guess);
  typedef LinearProblem<scalar_type, MV, OP> prob_type;
  RCP<prob_type> problem (new prob_type (A, X, B));

  // Create a GMRES solver manager to solve the problem.
  typedef GmresSolMgr<scalar_type, MV, OP> solver_type;
  err << "Constructing the solver manager object...";
  solver_type solver (problem, params, debug);
  err << "done." << endl;

  // If specified, dump the parameters that the GMRES solver manager
  // is using to an XML file.  These parameters may be different than
  // those specified on construction of the solver manager, since the
  // latter may fill in unspecified default values and / or silently
  // correct invalid values.
  if (myRank == 0 && xmlOutFile != "") {
    err << "Dumping parameters for GMRES solve to XML file \"" 
	<< xmlOutFile << "\"...";
    RCP<const ParameterList> curParams = solver.getCurrentParameters();
    Teuchos::writeParameterListToXmlFile (*curParams, xmlOutFile);
    err << "done." << endl;
  }

  // Attempt to solve the linear system(s).
  // If in verbose mode, report the result.
  Belos::ReturnType result = solver.solve();
  if (verbose) {
    const int numEquations = MVT::GetNumberVecs(*B);
    out << "Total number of iterations: " << solver.getNumIters() 
	<< endl
	<< "Relative residual norm" << (numEquations != 1 ? "s" : "") 
	<< ":" << endl;

    // Norm(s) of the right-hand side(s) are useful for computing
    // relative residuals.
    std::vector<magnitude_type> rhsNorms (numEquations);
    MVT::MvNorm (*B, rhsNorms);

    // Compute relative (or absolute, if RHS has norm zero) residual
    // norms for the initial guess(es) for each of the equation(s).
    std::vector<magnitude_type> initResNorms (numEquations);
    RCP<MV> R = MVT::Clone (*B, numEquations);
    OPT::Apply (*A, *X_guess, *R);
    MVT::MvAddMv (-STS::one(), *R, STS::one(), *B, *R);
    MVT::MvNorm (*R, initResNorms);

    // Compute relative (or absolute, if RHS has norm zero) residual
    // norms for the approximate solution(s) for each of the
    // equation(s).  X should be the same multivector as returned by
    // problem->getLHS().
    std::vector<magnitude_type> finalResNorms (numEquations);
    OPT::Apply (*A, *X, *R);
    MVT::MvAddMv (-STS::one(), *R, STS::one(), *B, *R);
    MVT::MvNorm (*R, finalResNorms);

    // Compute absolute solution error(s) ||X_exact - X_guess||_2.
    std::vector<magnitude_type> absSolNorms (numEquations);
    {
      // In general, R may not necessarily be in the same vector
      // space as X, if we're using a left preconditioner.  With no
      // preconditioning, X and R must be in the same space, but we
      // prefer a more general approach; hence, we clone X_diff from
      // X, rather than recycling R for X_diff.  (MVT doesn't have a
      // notion of vector space, so we can't check whether R and X
      // are in the same vector space.)
      RCP<MV> X_diff = MVT::Clone (*X, numEquations);
      MVT::MvAddMv (STS::one(), *X, -STS::one(), *X_exact, *X_diff);
      MVT::MvNorm (*X_diff, absSolNorms);
    }

    // Display resulting residual norm(s) on Rank 0.
    for (std::vector<magnitude_type>::size_type k = 0;
	 k < rhsNorms.size(); ++k)
      {
	out << "For problem " << k+1 << " of " << numEquations << ": " 
	    << endl
	    << "* ||b||_2 = " << rhsNorms[k] << endl
	    << "* ||A x_guess - b||_2 ";
	if (rhsNorms[k] == STM::zero()) {
	  out << "= " << initResNorms[k] << endl;
	} else {
	  out << "/ ||b||_2 = "
	      << initResNorms[k] / rhsNorms[k] << endl;
	}
	out << "* ||A x - b||_2 ";
	if (rhsNorms[k] == STM::zero()) {
	  out << "= " << finalResNorms[k] << endl;
	} else {
	  out << "/ ||b||_2 = "
	      << finalResNorms[k] / rhsNorms[k] << endl;
	}
	out << "* ||x - x_exact||_2 = " << absSolNorms[k] << endl;
      }

    // Did the solution manager return a reasonable result
    // (Converged or Unconverged)?  (We allow the iterations not to
    // converge, since the test parameters may not have allowed
    // sufficient iterations for convergence.)  Broadcast the result
    // from Rank 0.  The result should be the same on all processes,
    // but we want to make sure.  Teuchos communication doesn't have
    // a reduction for Boolean values, so we use an int instead,
    // with the usual Boolean interpretation as in C (0 means false,
    // 1 means true).
    //
    // Note: "badness" doesn't mean "Unconverged", it means "neither
    // Converged nor Unconverged".
    int badness = 0;
    if (result == Belos::Converged) {
      out << "Result: Converged." << endl;
    } else if (result == Belos::Unconverged) {
      cout << "Result: Unconverged." << endl;	
    } else {
      badness = 1;
    }

    // Make sure that all processes throw the exception, to make the
    // test fail quickly.
    Teuchos::broadcast (*comm, 0, 1, &badness);
    if (badness) {
      err << "GmresSolMgr::solve() returned unknown result value " 
	  << result << "." << endl;
      Teuchos::barrier (*comm);
      return EXIT_FAILURE;
    }
  }

  // Report final result (PASSED or FAILED) to stdout on Rank 0.  This
  // tells the Trilinos test framework whether the test passed or
  // failed.  The framework must see "End Result: TEST PASSED" (or
  // something like that -- it's a regular expression test) in order
  // for the test to be considered passed.
  //
  // Later, numFailed might signify the number of test(s) that failed.
  // For now, it's always zero.
  int numFailed = 0;
  if (numFailed != 0)
    {
      if (verbose)
	{
	  err << "There were " << numFailed << " error" 
	      << (numFailed != 1 ? "s." : ".") << endl;
	}
      out << "End Result: TEST FAILED" << endl;
      return EXIT_FAILURE;
    }
  else 
    {
      out << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
}



