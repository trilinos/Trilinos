//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
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
#include <Teuchos_XMLParameterListHelpers.hpp>

//
// These typedefs make main() as generic as possible.
//
typedef double scalar_type;
typedef int local_ordinal_type;

// mfh 22 Feb 2011: Explicit instantiation of templates in the checkin
// test script makes my life hard sometimes...
#if 0
#  ifdef HAVE_TEUCHOS_LONG_LONG_INT
// Half the point of Tpetra is the capability of using global ordinals
// that are bigger than local ordinals.  We should test this whenever
// possible, even in code that is not specifically testing Tpetra.
typedef long long int global_ordinal_type;
#  else // not HAVE_TEUCHOS_LONG_LONG_INT
typedef long int global_ordinal_type;
#  endif // HAVE_TEUCHOS_LONG_LONG_INT
#else // not 0
typedef int global_ordinal_type;
#endif // 0

#ifdef HAVE_KOKKOS_TBB
typedef Kokkos::TBBNode node_type;
#else
typedef Kokkos::SerialNode node_type;
#endif // HAVE_KOKKOS_TBB

typedef Teuchos::ScalarTraits<scalar_type> STS;
typedef STS::magnitudeType magnitude_type;
typedef Teuchos::ScalarTraits<magnitude_type> STM;

typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> OP;
typedef Belos::MultiVecTraits<scalar_type, MV> MVT;
typedef Belos::OperatorTraits<scalar_type, MV, OP> OPT;
typedef Teuchos::SerialDenseMatrix<int, scalar_type> serial_matrix_type;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> sparse_matrix_type;


namespace {
  /// \fn getNode
  /// \brief Return an RCP to a Kokkos Node
  ///
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode() {
    throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
  }
  
  template<>
  Teuchos::RCP<Kokkos::SerialNode>
  getNode() {
    Teuchos::ParameterList defaultParams;
    return Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
  }
  
#if defined(HAVE_KOKKOS_TBB)
  template<>
  Teuchos::RCP<Kokkos::TBBNode>
  getNode() {
    // "Num Threads" specifies the number of threads.  Defaults to an
    // automatically chosen value.
    Teuchos::ParameterList defaultParams;
    return Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
  }
#endif // defined(HAVE_KOKKOS_TBB)
  
  /// Test Tpetra::MatrixMarket::Reader::readFile()
  ///
  /// \param filename [in] Name of the Matrix Market format sparse
  ///   matrix file to read (on MPI Rank 0 only)
  /// \param pComm [in] Communicator, over whose MPI ranks to
  ///   distribute the returned Tpetra::CrsMatrix.
  /// \param tolerant [in] Whether or not to parse the file 
  ///   tolerantly
  ///
  /// \return Tpetra::CrsMatrix read in from the file 
  template<class SparseMatrixType>
  Teuchos::RCP<SparseMatrixType>
  readFile (const std::string& filename,
	    const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
	    const bool tolerant,
	    const bool debug)
  {
    typedef typename SparseMatrixType::node_type node_type;
    Teuchos::RCP<node_type> pNode = getNode<node_type>();
    
    typedef Tpetra::MatrixMarket::Reader<SparseMatrixType> reader_type;
    return reader_type::readFile (filename, pComm, pNode, tolerant, debug);
  }

  template<class Ordinal, class Scalar>
  void
  writeRow (std::ostream& out, 
	    const Ordinal rowInd, 
	    const Ordinal colInd, 
	    const Scalar& val)
  {
    using std::endl;
    typedef typename Teuchos::ScalarTraits<Scalar> STS;

    out << (rowInd+1) << " " << (colInd+1) << " " << STS::real(val);
    if (STS::isComplex)
      out << " " << STS::imag(val);
    out << endl;
  }

  /// \fn generateTestProblem
  /// \brief Generate, distribute, and return a test sparse matrix.
  ///
  /// SparseMatrixType, the type of sparse matrix returned, must be a
  /// Tpetra::CrsMatrix.
  ///
  /// \param pComm [in] Communicator over which the sparse matrix is
  ///   to be distributed.
  ///
  /// \param globalNumRows [in] Global (over all MPI process(es))
  ///   number of rows in the sparse matrix.
  ///
  /// \param symmetric [in] Whether to generate a symmetric test problem.
  ///   Storage is nonsymmetric regardless; symmetry here only applies to 
  ///   the entries' locations and values.
  ///
  /// \param dumpToFile [in] Whether the generated sparse matrix
  ///   should be dumped to a Matrix Market file by Rank 0.
  ///
  /// \param outFilename [in] If my rank is 0 and dumpToFile is true,
  ///   dump the generated sparse matrix (in Matrix Market format,
  ///   nonsymmetric storage even if the matrix is symmetric) to the
  ///   file outFileName.  Otherwise, outFileName is ignored.
  ///
  /// \return The sparse matrix (global, distributed)
  template<class SparseMatrixType>
  Teuchos::RCP<SparseMatrixType>
  generateTestProblem (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
		       const typename SparseMatrixType::global_ordinal_type globalNumRows,
		       const bool symmetric,
		       const bool dumpToFile,
		       const std::string& outFileName)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::tuple;
    using Tpetra::createUniformContigMapWithNode;
    using std::endl;

    typedef typename SparseMatrixType::local_ordinal_type LO;
    typedef typename SparseMatrixType::global_ordinal_type GO;
    typedef typename SparseMatrixType::node_type NT;
    typedef Tpetra::Map<LO, GO, NT> map_type;

    RCP<node_type> pNode = getNode<node_type>();
    // For a square matrix, we only need a Map for the range of the matrix.
    RCP<const map_type> pRangeMap = 
      createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, pComm, pNode);
    // The sparse matrix object to fill.
    RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

    const int myRank = Teuchos::rank (*pComm);
    if (myRank == 0)
      {
	typedef typename SparseMatrixType::scalar_type Scalar;
	typedef Teuchos::ScalarTraits<Scalar> STS;
	const Scalar leftVal = -STS::one();
	const Scalar rightVal = symmetric ? -STS::one() : +2*STS::one();
	// Boost the diagonal if nonsymmetric, hopefully to make
	// convergence a little faster.
	const Scalar centerVal = symmetric ? 2*STS::one() : 8*STS::one(); 

	// Output file to which to dump the generated sparse matrix,
	// if the caller specified that this should happen.
	std::ofstream outFile;
	if (dumpToFile)
	  {
	    outFile.open (outFileName.c_str());
	    // "Banner" line for the Matrix Market file.
	    if (STS::isComplex)
	      outFile << "%%MatrixMarket matrix coordinate complex general" << endl;
	    else
	      outFile << "%%MatrixMarket matrix coordinate real general" << endl;
	    // Number of nonzeros.  Every row has three nonzeros,
	    // except for the first and last row, which each have two.
	    // Special cases when there are less than three rows
	    // (which we don't allow anyway, but why not be general?).
	    GO numNonzeros = 0;
	    if (globalNumRows < 2)
	      numNonzeros = 1;
	    else if (globalNumRows < 3)
	      numNonzeros = 4;
	    else // if (globalNumRows >= 3)
	      numNonzeros = (globalNumRows - 2) * 3 + 4;
	    // Write "<numRows> <numCols> <numNonzeros>" line to output file.
	    outFile << globalNumRows << " " << globalNumRows << " "
		    << numNonzeros << endl;
	  }
	for (GO curRow = 0; curRow < globalNumRows; ++curRow)
	  {
	    if (curRow > 0) 
	      {
		pMat->insertGlobalValues (curRow, tuple(curRow-1), tuple(leftVal));
		if (dumpToFile)
		  writeRow<GO, Scalar> (outFile, curRow, curRow-1, leftVal);
	      }
	    pMat->insertGlobalValues (curRow, tuple(curRow), tuple(centerVal));
	    if (dumpToFile)
	      writeRow<GO, Scalar> (outFile, curRow, curRow, centerVal);
	    if (curRow < globalNumRows-1)
	      {
		pMat->insertGlobalValues (curRow, tuple(curRow+1), tuple(rightVal));
		if (dumpToFile)
		  writeRow<GO, Scalar> (outFile, curRow, curRow+1, rightVal);
	      }
	  }
	if (dumpToFile)
	  outFile.close ();
      }
    // Make sure Rank 0 is done filling in the matrix.
    Teuchos::barrier (*pComm);
    // fillComplete() doesn't need any arguments if the domain and
    // range maps are the same.
    pMat->fillComplete ();
    return pMat;
  }

  template<class SparseMatrixType>
  Teuchos::RCP<SparseMatrixType>
  generateGmresUnfriendlyMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
				 const typename SparseMatrixType::global_ordinal_type globalNumRows,
				 const bool dumpToFile,
				 const std::string& outFileName)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::tuple;
    using Tpetra::createUniformContigMapWithNode;
    using std::endl;

    typedef typename SparseMatrixType::local_ordinal_type LO;
    typedef typename SparseMatrixType::global_ordinal_type GO;
    typedef typename SparseMatrixType::node_type NT;
    typedef Tpetra::Map<LO, GO, NT> map_type;

    RCP<node_type> pNode = getNode<node_type>();
    // For a square matrix, we only need a Map for the range of the matrix.
    RCP<const map_type> pRangeMap = 
      createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, pComm, pNode);
    // The sparse matrix object to fill.
    RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

    const int myRank = Teuchos::rank (*pComm);
    if (myRank == 0)
      {
	typedef typename SparseMatrixType::scalar_type Scalar;
	typedef Teuchos::ScalarTraits<Scalar> STS;
	const Scalar val = STS::one();

	// Output file to which to dump the generated sparse matrix,
	// if the caller specified that this should happen.
	std::ofstream outFile;
	if (dumpToFile)
	  {
	    outFile.open (outFileName.c_str());
	    // "Banner" line for the Matrix Market file.
	    if (STS::isComplex)
	      outFile << "%%MatrixMarket matrix coordinate complex general" << endl;
	    else
	      outFile << "%%MatrixMarket matrix coordinate real general" << endl;
	    // Number of nonzeros.  Every row has exactly one nonzero.
	    const GO numNonzeros = globalNumRows;
	    // Write "<# rows> <# cols> <# nonzeros>" line to output file.
	    outFile << globalNumRows << " " << globalNumRows << " "
		    << numNonzeros << endl;
	  }
	for (GO curRow = 0; curRow < globalNumRows; ++curRow)
	  {
	    const GO curCol = (curRow == 0) ? (globalNumRows-1) : (curRow-1);
	    pMat->insertGlobalValues (curRow, tuple(curCol), tuple(val));
	    if (dumpToFile)
	      writeRow<GO, Scalar> (outFile, curRow, curCol, val);
	  }
	if (dumpToFile)
	  outFile.close ();
      }
    // Make sure Rank 0 is done filling in the matrix.
    Teuchos::barrier (*pComm);
    // fillComplete() doesn't need any arguments if the domain and
    // range maps are the same.
    pMat->fillComplete ();
    return pMat;
  }

  /// Print out the given Belos::MsgType to a comma-delimited list of
  /// names.
  std::string 
  msgTypeToString (const int msgType)
  {
    // Wouldn't it be nice if C++ enums had introspection and could
    // be enumerated?
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
    const int numValidTypes = 8;
    std::ostringstream os;
    for (int k = 0; k < numValidTypes; ++k)
      {
	if (msgType & validTypes[k])
	  os << typeNames[k];
	if (k > 0 && k < numValidTypes - 1)
	  os << ", ";
      }
    return os.str();
  }
} // (anonymous namespace)

/// \fn main
/// \brief Test driver for (Mat)OrthoManager subclasses
int 
main (int argc, char *argv[]) 
{
  using Belos::GmresSolMgr;
  using Belos::LinearProblem;
  using Belos::OrthoManagerFactory;
  using Teuchos::CommandLineProcessor;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cerr;
  using std::cout;
  using std::endl;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Teuchos::Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = Teuchos::rank (*pComm);

  // Name of the Matrix Market - format sparse matrix file from which
  // to read the matrix A.
  std::string filename;
  // Optional argument: Name of XML file containing parameters for
  // the GMRES solver manager.
  std::string xmlInFile;
  // Optional argument: Name of XML file to which to write the
  // parameters used by the GMRES solver manager.
  std::string xmlOutFile;
  // Whether, instead of reading the matrix from a file, to use a
  // generated test problem.
  bool generated = false;
  // If generating the test problem: the global number of rows (and
  // columns) in the generated matrix.  The global number of rows
  // should be >= max(3, Teuchos::size(*pComm)), so that each MPI
  // process has at least one row, and the matrix is nonsingular.
  // Ignored if a sparse matrix file is given.
  int globalNumRows = std::max (3, Teuchos::size (*pComm));

  // Whether the generated test problem is symmetric.
  // This option is ignored if generated is false.
  bool symmetric = false;
  // Whether to generate a "GMRES-unfriendly" test problem (the
  // classical example of a permutation matrix).
  // This option is ignored if generated is false.
  bool gmresUnfriendly = false;
  // Whether to dump the generated test matrix to a Matrix Market -
  // format output file (named outFileName, which must be provided).
  bool dumpToFile = false;
  // Name of the file to which to dump the generated test matrix.
  std::string outFileName;

  // Optional argument: Whether to print verbose output.
  bool verbose = false;
  // Optional argument: Whether to print debugging information.
  bool debug = false;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("filename", &filename, "Name of a Matrix Market - format "
		  "sparse matrix file from which to read the matrix A in the "
		  "linear system Ax=b to solve.  Should not be provided if "
		  "the \"--generated\" option is set.");
  cmdp.setOption ("generated", "nongenerated", &generated, "Whether to use a "
		  "generated test problem.  Should not be provided if the "
		  "\"--filename=<file>\" option is given.");
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
  cmdp.setOption ("dumpToFile", "noDumpToFile", &dumpToFile, "Whether to dump "
		  "the generated test matrix to a Matrix Market - format output "
		  "file (named outFileName, which must be provided).  This "
		  "option is ignored if the \"--generated\" option is not "
		  "provided.");
  cmdp.setOption ("outFileName", &outFileName, "Name of the file to which to "
		  "dump the generated test matrix.  This option is ignored if "
		  "the \"--generated\" and \"--dumpToFile\" options are not "
		  "provided.");
  cmdp.setOption ("params", &xmlInFile, "Name of an XML file containing "
		  "parameters (in Teuchos::ParameterList form) for the solve.");
  cmdp.setOption ("dump", &xmlOutFile, "Name of an XML file to which to dump "
		  "the parameters used by the GMRES solver manager.  These may "
		  "be different than those given to the solver manager.");
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // specify which matrix to use (from a file, or generated), we let
    // the "test" pass trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED || 
	(filename == "" && ! generated))
      {
	if (myRank == 0)
	  cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
    if (generated)
      {
	TEST_FOR_EXCEPTION(filename != "",
			   std::invalid_argument,
			   "The --generated and \"--filename=<file>\" options may "
			   "not both be used.");
	TEST_FOR_EXCEPTION(globalNumRows < std::max(3, Teuchos::size(*pComm)),
			   std::invalid_argument,
			   "The number of rows in the test matrix to generate "
			   "must be at least max(3, # MPI processes), in order "
			   "for the matrix to be nonsingular and for each MPI "
			   "process to have at least one row of the matrix.");
	TEST_FOR_EXCEPTION(dumpToFile && outFileName == "", 
			   std::invalid_argument,
			   "If the \"--dumpToFile\" option is given, an output "
			   "file must be given via \"--outFileName=<filename>\""
			   ".");
      }
  }
  // If the name of an XML file of parameters was provided, read it on
  // Proc 0, and broadcast from Proc 0 to the other procs.
  //
  // Start with an empty ParameterList on all procs.  If a filename
  // was provided, fill in the list on Proc 0, and broadcast.
  RCP<ParameterList> params = Teuchos::parameterList();
  const bool readParamsFromFile = (xmlInFile != "");
  // Default verbosity: only print warnings and errors.
  int verbosityLevel = Belos::Errors | Belos::Warnings;
  if (readParamsFromFile)
    {
      if (debug && myRank == 0)
	cerr << "Reading parameters for GMRES solve from XML file \"" 
	     << xmlInFile << "\"...";
      using Teuchos::updateParametersFromXmlFileAndBroadcast;
      updateParametersFromXmlFileAndBroadcast (xmlInFile, params.getRawPtr(), 
					       *pComm);
      if (debug && myRank == 0)
	cerr << "done." << endl;
      try {
	verbosityLevel = params->get<int>("Verbosity");
      } catch (Teuchos::Exceptions::InvalidParameter&) {
	verbosityLevel = Belos::Errors | Belos::Warnings;
      }
    }
  if (verbose || debug)
    { // Change verbosity level from its default.
      if (verbose)
	{
	  verbosityLevel = verbosityLevel | 
	    Belos::IterationDetails | 
	    Belos::OrthoDetails |
	    Belos::FinalSummary |
	    Belos::TimingDetails |
	    Belos::StatusTestDetails;
	}
      if (debug)
	verbosityLevel = verbosityLevel | Belos::Debug;
      params->set ("Verbosity", verbosityLevel);
      // Change the output frequency to print output at every
      // iteration.
      params->set ("Output Frequency", static_cast<int>(1));
    }
  if (debug)
    cerr << "Verbosity: " << msgTypeToString(verbosityLevel) << endl;

  //
  // Let's make ourselves a sparse matrix A.  Either read it from the
  // given file, or generate it, depending on the command-line
  // options.
  //
  RCP<sparse_matrix_type> A;
  if (filename != "")
    { // Read the sparse matrix A from the file.
      const bool tolerant = false;
      if (debug && myRank == 0)
	cerr << "Reading sparse matrix A from file \"" << filename << "\"...";
      A = readFile<sparse_matrix_type> (filename, pComm, tolerant, debug);
      if (debug && myRank == 0)
	cerr << "done." << endl;
      // The matrix A must be square, and have a nonzero number of
      // rows and columns, in order to solve Ax=b with GMRES.
      TEST_FOR_EXCEPTION(A->getGlobalNumRows() != A->getGlobalNumRows(), 
			 std::invalid_argument,
			 "The sparse matrix A must be square in order to solve "
			 "Ax=b with GMRES.");
      TEST_FOR_EXCEPTION(A->getGlobalNumRows() == 0,
			 std::invalid_argument,
			 "The sparse matrix A must have a nonzero number of rows "
			 "and columns in order to solve Ax=b with GMRES.");
    }
  else if (generated)
    {
      if (gmresUnfriendly)
	A = generateGmresUnfriendlyMatrix<sparse_matrix_type> (pComm, globalNumRows, dumpToFile, outFileName);
      else
	A = generateTestProblem<sparse_matrix_type> (pComm, globalNumRows, symmetric, dumpToFile, outFileName);
    }
  else 
    throw std::logic_error("Should never get here!");

  // Construct X_guess (the initial guess for the solution of AX=B)
  // from the domain of the matrix A, and fill it with zeros.
  if (debug && myRank == 0)
    cerr << "Constructing initial guess vector X_guess...";
  RCP<MV> X_guess = rcp (new MV (A->getDomainMap(), 1));
  MVT::MvInit (*X_guess, STS::zero());
  if (debug && myRank == 0)
    cerr << "done." << endl;

  // Our choice of exact solution and right-hand side depend on the
  // test problem.  If we generated the GMRES-unfriendly example, we
  // need B = e_{globalNumRows} and therefore X_exact =
  // e_{globalNumRows} as well.  Otherwise, we pick X_exact first and
  // compute B via SpMV: B = A * X_exact.
  RCP<MV> X_exact = rcp (new MV (A->getDomainMap(), 1));
  // Construct the right-hand side B from the range of the matrix A.
  // Don't just clone X_guess, since the range may differ from the
  // domain.  (While we require that the matrix A be square, we don't
  // require that the domain and range have the same distribution.)
  RCP<MV> B = rcp (new MV (A->getRangeMap(), 1));

  if (generated && gmresUnfriendly)
    {
      if (debug && myRank == 0)
	cerr << "Constructing B and X_exact for canonical \"GMRES-unfriendly\""
	  " example...";
      // We know at this point that MV is a Tpetra::MultiVector
      // object, so we can escape the MVT interface in order to set an
      // individual value.  We might as well use other
      // Tpetra::MultiVector methods here as well.
      X_exact->putScalar (STS::zero());
      X_exact->replaceGlobalValue (globalNumRows-1, 0, STS::one());
      B->putScalar (STS::zero());
      B->replaceGlobalValue (globalNumRows-1, 0, STS::one());
      if (debug && myRank == 0)
	cerr << "done." << endl;
    }
  else
    {
      // Construct the exact solution vector and fill it with all ones.
      // Tacky, but deterministic.  Not so good if we expect A to be
      // singular with rigid body modes.
      if (debug && myRank == 0)
	cerr << "Setting X_exact = [1; ...; 1]...";
      MVT::MvInit (*X_exact, STS::one());
      if (debug && myRank == 0)
	cerr << "done." << endl;
      
      // Compute the right-hand side B := A*X_exact.
      if (debug && myRank == 0)
	cerr << "Computing B := A*X_exact...";
      OPT::Apply (*A, *X_exact, *B);
      if (debug && myRank == 0)
	cerr << "done." << endl;
    }

  // In debug mode, compute and print the initial residual
  // independently of the solver framework.
  if (debug)
    {
      RCP<MV> R = MVT::Clone (*B, MVT::GetNumberVecs (*B));
      // R := A * X_guess
      OPT::Apply (*A, *X_guess, *R);
      // R := B - R
      MVT::MvAddMv (STS::one(), *B, -STS::one(), *R, *R);
      std::vector<magnitude_type> theNorm (MVT::GetNumberVecs (*R));
      MVT::MvNorm (*R, theNorm);

      if (myRank == 0)
	cerr << "Initial residual norm: ||B - A*X_guess||_2 = " 
	     << theNorm[0] << endl;

      MVT::MvNorm (*B, theNorm);
      if (myRank == 0)
	cerr << "||B||_2 = "
	     << theNorm[0] << endl;

      MVT::MvNorm (*X_exact, theNorm);
      if (myRank == 0)
	cerr << "||X_exact||_2 = " << theNorm[0] << endl;

      MVT::MvNorm (*X_guess, theNorm);
      if (myRank == 0)
	cerr << "||X_guess||_2 = " << theNorm[0] << endl;
    }

  // Construct the linear problem to solve.  X_guess is only copied
  // shallowly and will likely be changed by the solve.
  typedef LinearProblem<scalar_type, MV, OP> prob_type;
  RCP<prob_type> problem (new prob_type (A, X_guess, B));

  // Create a GMRES solver manager to solve the problem.
  typedef GmresSolMgr<scalar_type, MV, OP> solver_type;
  if (debug && myRank == 0)
    cerr << "Constructing the solver manager object...";
  solver_type solver (problem, params, debug);
  if (debug && myRank == 0)
    cerr << "done." << endl;

  // If specified, dump the parameters that the GMRES solver manager
  // is using to an XML file.  These parameters may be different than
  // those specified on construction of the solver manager, since the
  // latter may fill in unspecified default values and / or silently
  // correct invalid values.
  if (myRank == 0 && xmlOutFile != "")
    {
      if (debug && myRank == 0)
	cerr << "Dumping parameters for GMRES solve to XML file \"" 
	     << xmlOutFile << "\"...";
      RCP<const ParameterList> curParams = solver.getCurrentParameters();
      Teuchos::writeParameterListToXmlFile (*curParams, xmlOutFile);
      if (debug && myRank == 0)
	cerr << "done." << endl;
    }

  // Attempt to solve the linear system, and report the result.
  Belos::ReturnType result = solver.solve();
  if (verbose)
    {
      // Teuchos communication doesn't have a reduction for bool, so
      // we use an int instead.
      int badness = 0;
      if (myRank == 0)
	{
	  if (result == Belos::Converged)
	    cout << "Result: Converged." << endl;
	  else if (result == Belos::Unconverged)
	    cout << "Result: Unconverged." << endl;	
	  else
	    badness = 1;
	}	
      // Make sure that all processes throw the exception, to make the
      // test fail quickly.
      Teuchos::broadcast (*pComm, 0, 1, &badness);
      if (badness)
	{
	  if (myRank == 0)
	    cerr << "GmresSolMgr::solve() returned unknown result value " 
		 << result << "." << endl;
	  Teuchos::barrier (*pComm);
	  return EXIT_FAILURE;
	}
    }

  // Later, numFailed might signify the number of test(s) that failed.
  // For now, it's always zero.
  int numFailed = 0;
  if (numFailed != 0)
    {
      if (myRank == 0)
	{
	  if (verbose)
	    {
	      cerr << "There were " << numFailed << " error" 
		   << (numFailed != 1 ? "s." : ".") << endl;
	    }
	  // The Trilinos test framework depends on seeing this message,
	  // so don't rely on the OutputManager to report it correctly.
	  cout << "End Result: TEST FAILED" << endl;
	}
      return EXIT_FAILURE;
    }
  else 
    {
      if (myRank == 0)
	cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
}



