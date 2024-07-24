// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This benchmark exercises different ways of finding the number of
// entries in each row of a sparse matrix.  It uses three different
// data structures: Epetra_CrsMatrix (in the Epetra package),
// Tpetra::CrsMatrix (in this, the TpetraCore subpackage of the Tpetra
// package), and KokkosSparse::CrsMatrix (in the TpetraKernels
// subpackage of the Tpetra package).
//
// Both Epetra and Tpetra are meant for all users.  KokkosSparse is
// meant for developers and expert users.  Epetra is the older sparse
// linear algebra library, which supports only MPI and a little bit of
// OpenMP parallelism.  Tpetra is newer and supports both MPI and
// thread parallelism.  Tpetra uses the Kokkos programming model for
// thread parallelism and data structures.  Furthermore,
// Tpetra::CrsMatrix uses KokkosSparse::CrsMatrix to store the sparse
// matrix, after it has been constructed.
//
// All three data structures use compressed sparse row (CSR) storage
// for the sparse matrix.  KokkosSparse::CrsMatrix is really just CSR,
// and is only a local data structure; it does not do MPI
// communication or know about MPI.  The Epetra and Tpetra data
// structures both do MPI communication.  They also support different
// ways to create and modify a sparse matrix, which means that they
// are not _just_ CSR; they have other data structures which they
// enable or disable as needed.  In particular, Tpetra::CrsMatrix
// _uses_ KokkosSparse::CrsMatrix inside.  Thus, it's not fair to
// compare KokkosSparse::CrsMatrix directly to the Epetra and Tpetra
// data structures; one should properly compare Epetra to Tpetra.

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI

#ifdef HAVE_TPETRACORE_EPETRA
#  include "Epetra_Comm.h"
#  include "Epetra_CrsGraph.h"
#  include "Epetra_CrsMatrix.h"
#  include "Epetra_Map.h"
#  include "Epetra_SerialComm.h"
#  ifdef HAVE_TEUCHOS_MPI
#    include "Epetra_MpiComm.h"
#  endif // HAVE_TEUCHOS_MPI
#endif // HAVE_TPETRACORE_EPETRA

namespace { // (anonymous)

  // Options to read in from the command line
  struct CmdLineOpts {
    // Do the benchmark this many times in a single timing loop, in
    // case the timer's granularity is too coarse to capture run time
    // to adequate precision.
    int numTrials;
    // Local (per MPI process) number of rows in the sparse matrix.
    int lclNumRows;
    // Number of entries per row in the sparse matrix.  This benchmark
    // doesn't care so much about the matrix's actual structure.
    int numEntPerRow;
    // Whether to test Epetra_CrsMatrix::ExtractMyRowView
    bool testEpetra;
    // Whether to test Epetra_CrsMatrix::NumMyEntries
    bool testEpetraLen;
    // Whether to test Tpetra::CrsMatrix::getLocalRowView
    bool testTpetra;
    // Whether to test Tpetra::CrsMatrix::getNumEntriesInLocalRow
    bool testTpetraLen;
    // Whether to test KokkosSparse::CrsMatrix methods
    bool testKokkos;
  };

  // Use a utility from the Teuchos package of Trilinos to set up
  // command-line options for reading, and set default values of
  // command-line options.  clp is an output argument containing the
  // set-up options.  It retains pointers to fields in 'opts'.
  // Reading the command-line options will update those fields in
  // place.
  void
  setCmdLineOpts (CmdLineOpts& opts,
                  Teuchos::CommandLineProcessor& clp)
  {
    // Set default values of command-line options.

    opts.numTrials = 250;
    opts.lclNumRows = 10000;
    opts.numEntPerRow = 10;
    // This example lives in Tpetra, and has the option to use Epetra,
    // but it does not require Epetra.  Users are not allowed to ask
    // for the Epetra benchmarks to run if Epetra is not available.
#ifdef HAVE_TPETRACORE_EPETRA
    const bool testEpetraDefault = true;
#else
    const bool testEpetraDefault = false;
#endif //HAVE_TPETRACORE_EPETRA
    opts.testEpetra = testEpetraDefault;
    opts.testEpetraLen = testEpetraDefault;
    opts.testTpetra = true;
    opts.testTpetraLen = true;
    opts.testKokkos = true;

    clp.setOption ("numTrials", &(opts.numTrials),
                   "Number of trials per timing loop (to increase timer precision)");
    clp.setOption ("lclNumRows", &(opts.lclNumRows), "Number of rows in the "
                   "matrix per MPI process");
    clp.setOption ("numEntPerRow", &(opts.numEntPerRow),
                   "Number of entries per row of the matrix");
    clp.setOption ("testEpetra", "skipEpetra", &(opts.testEpetra),
                   "Whether to test Epetra_CrsMatrix::ExtractMyRowView");
    clp.setOption ("testEpetraLen", "skipEpetraLen", &(opts.testEpetraLen),
                   "Whether to test Epetra_CrsMatrix::NumMyEntries");
    clp.setOption ("testTpetra", "skipTpetra", &(opts.testTpetra),
                   "Whether to test Tpetra::CrsMatrix::getLocalRowView");
    clp.setOption ("testTpetraLen", "skipTpetraLen", &(opts.testTpetraLen),
                   "Whether to test Tpetra::CrsMatrix::getNumEntriesInLocalRow");
    clp.setOption ("testKokkos", "skipKokkos", &(opts.testKokkos),
                   "Whether to test KokkosSparse::CrsMatrix");
  }

  // Actually read the command-line options from the command line,
  // using the argc and argv arguments to main().  Use the clp output
  // argument of setCmdLineOpts.  The options actually get written to
  // the same CmdLineOpts struct instance that was passed into
  // setCmdLineOpts above.
  //
  // Return 0 if successful, 1 if help printed due to the user
  // specifying the "--help" option (indicates that the application
  // shouldn't run), and -1 on error.
  int
  parseCmdLineOpts (Teuchos::CommandLineProcessor& clp, int argc, char* argv[])
  {
    auto result = clp.parse (argc, argv);

    switch (result) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return 1;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return -1;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      return 0;
    default:
      return -1;
    }
  }

  // Check the command-line options that were read in by
  // parseCmdLineOpts.  Return 0 if all correct, else return nonzero,
  // using the LAPACK error reporting convention of the negative of
  // the argument in its original order (starting with 1) as the error
  // code.  Print informative error messages to the given output
  // stream \c out.
  int
  checkCmdLineOpts (std::ostream& out, const CmdLineOpts& opts)
  {
    int err = 0;

    if (opts.numTrials < 0) {
      out << "numTrials = " << opts.numTrials << " < 0." << std::endl;
      err = -1; // LAPACK error reporting convention
    }
    if (opts.lclNumRows < 0) {
      out << "lclNumRows = " << opts.lclNumRows << " < 0." << std::endl;
      err = -2; // LAPACK error reporting convention
    }
    if (opts.numEntPerRow < 0) {
      out << "numEntPerRow = " << opts.numEntPerRow << " < 0." << std::endl;
      err = -3; // LAPACK error reporting convention
    }
#ifndef HAVE_TPETRACORE_EPETRA
    if (opts.testEpetra || opts.testEpetraLen) {
      out << "If you want to test Epetra, you must enable the Epetra package "
        "when configuring Trilinos, by setting Trilinos_ENABLE_Epetra=ON."
          << std::endl;
      err = -4; // LAPACK error reporting convention
    }
#endif // ! HAVE_TPETRACORE_EPETRA

    return err;
  }

  // Return a pointer (RCP is like std::shared_ptr) to an output
  // stream.  It prints on Process 0 of the given MPI communicator,
  // but ignores all output on other MPI processes.
  Teuchos::RCP<Teuchos::FancyOStream>
  getOutputStream (const Teuchos::Comm<int>& comm)
  {
    using Teuchos::getFancyOStream;

    const int myRank = comm.getRank ();
    if (myRank == 0) {
      // Process 0 of the given communicator prints to std::cout.
      return getFancyOStream (Teuchos::rcpFromRef (std::cout));
    }
    else {
      // A "black hole output stream" ignores all output directed to it.
      return getFancyOStream (Teuchos::rcp (new Teuchos::oblackholestream ()));
    }
  }

  // Get a Tpetra::CrsMatrix for use in benchmarks.  This method takes
  // parameters that come from the command-line options read in by
  // parseCmdLineOpts.
  Teuchos::RCP<Tpetra::CrsMatrix<> >
  getTpetraMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   const CmdLineOpts& opts)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Tpetra::Map<> map_type;
    typedef Tpetra::CrsGraph<> graph_type;
    typedef Tpetra::CrsMatrix<> matrix_type;
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef Tpetra::CrsMatrix<>::scalar_type SC;
    typedef Tpetra::CrsMatrix<>::mag_type MT;
    typedef Tpetra::global_size_t GST;

    const LO lclNumRows = opts.lclNumRows;
    const GST gblNumRows = static_cast<GST> (opts.lclNumRows) *
      static_cast<GST> (comm->getSize ());
    const GO indexBase = 0;

    // A Map describes a distribution of data over MPI processes.
    // This "row Map" will describe the distribution of rows of the
    // sparse matrix that we will create.
    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, static_cast<size_t> (lclNumRows),
                         indexBase, comm));
    const GO gblNumCols = static_cast<GO> (rowMap->getGlobalNumElements ());
    // Create the graph structure of the sparse matrix.
    RCP<graph_type> G =
      rcp (new graph_type (rowMap, opts.numEntPerRow));
    // Fill in the graph structure of the sparse matrix.
    Teuchos::Array<GO> gblColInds (opts.numEntPerRow);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) { // for each of my rows
      const GO gblInd = rowMap->getGlobalElement (lclRow);
      // Just put some entries in the graph.  The actual column
      // indices don't matter so much, as long as they make the
      // resulting matrix square and don't go out of bounds.
      for (LO k = 0; k < static_cast<LO> (opts.numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
      }
      G->insertGlobalIndices (gblInd, gblColInds ());
    }
    // Make the graph ready for use by a matrix.
    G->fillComplete ();

    // Create the sparse matrix.  Tpetra and Epetra have many ways to
    // create a sparse matrix.  The most efficient way and preferred
    // is to use a previously fill-complete graph.  This constrains
    // the matrix's structure never to change (the input graph is
    // const), which lets Tpetra and Epetra make assumptions that
    // speed up changing entries of the matrix.
    RCP<matrix_type> A = rcp (new matrix_type (G));

    // Fill the matrix with values.  Their values don't matter for
    // this benchmark.
    Teuchos::Array<SC> vals (opts.numEntPerRow);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblInd = rowMap->getGlobalElement (lclRow);
      for (LO k = 0; k < static_cast<LO> (opts.numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
        // Cast first to MT, then SC, to avoid issues like
        // std::complex<double>'s constructor not taking int.
        vals[k] = static_cast<SC> (static_cast<MT> (gblColInds[k]));
      }
      A->replaceGlobalValues (gblInd, gblColInds (), vals ());
    }
    // We're done changing entries of the sparse matrix.
    A->fillComplete ();

    return A; // return is a shallow copy (RCP is like std::shared_ptr)
  }

#ifdef HAVE_TPETRACORE_EPETRA
  // Convert from a Teuchos::Comm MPI communicator wrapper (used by
  // Tpetra classes) to an Epetra communicator wrapper (used by Epetra
  // classes).
  Teuchos::RCP<const Epetra_Comm>
  makeEpetraComm (const Teuchos::Comm<int>& comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

    // Trilinos uses MPI, but does not require it.  We want to be able
    // to build and run this benchmark whether or not MPI is enabled.
#ifdef HAVE_TEUCHOS_MPI
    const Teuchos::MpiComm<int>* mpiComm =
      dynamic_cast<const Teuchos::MpiComm<int>* > (&comm);
    if (mpiComm == NULL) {
      const Teuchos::SerialComm<int>* serialComm =
        dynamic_cast<const Teuchos::SerialComm<int>* > (&comm);
      TEUCHOS_TEST_FOR_EXCEPTION
        (serialComm == NULL, std::logic_error, "Input Teuchos::Comm "
         "is neither a Teuchos::MpiComm nor a Teuchos::SerialComm.");
      RCP<const Epetra_SerialComm> epetraComm = rcp (new Epetra_SerialComm ());
      return rcp_implicit_cast<const Epetra_Comm> (epetraComm);
    }
    else {
      auto rawMpiCommWrapped = mpiComm->getRawMpiComm ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (rawMpiCommWrapped.is_null (), std::logic_error, "The input "
         "Teuchos::Comm is a Teuchos::MpiComm, but its getRawMpiComm() "
         "method return null.");
      MPI_Comm rawMpiComm = *rawMpiCommWrapped;
      RCP<const Epetra_MpiComm> epetraComm =
        rcp (new Epetra_MpiComm (rawMpiComm));
      return rcp_implicit_cast<const Epetra_Comm> (epetraComm);
    }
#else // ! HAVE_TEUCHOS_MPI
    RCP<const Epetra_SerialComm> epetraComm = rcp (new Epetra_SerialComm ());
    return rcp_implicit_cast<const Epetra_Comm> (epetraComm);
#endif // HAVE_TEUCHOS_MPI
  }

  // Get an Epetra_CrsMatrix for use in benchmarks.  This method takes
  // parameters that come from the command-line options read in by
  // parseCmdLineOpts.  Epetra and Tpetra are similar enough that
  // reading the comments in getTpetraMatrix() above should help you
  // understand what this function does inside.
  Teuchos::RCP<Epetra_CrsMatrix>
  getEpetraMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   const CmdLineOpts& opts)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::outArg;
    typedef Epetra_Map map_type;
    typedef Epetra_CrsGraph graph_type;
    typedef Epetra_CrsMatrix matrix_type;
    typedef double SC;
    typedef int LO;
    typedef int GO;
    typedef double MT;

    const LO lclNumRows = opts.lclNumRows;
    const GO gblNumRows = static_cast<GO> (opts.lclNumRows) *
      static_cast<GO> (comm->getSize ());
    const GO indexBase = 0;

    RCP<const Epetra_Comm> epetraComm = makeEpetraComm (*comm);
    const map_type rowMap (gblNumRows, lclNumRows, indexBase, *epetraComm);
    const GO gblNumCols = gblNumRows;
    const bool staticProfile = true;
    graph_type G (Copy, rowMap, opts.numEntPerRow, staticProfile);

    Teuchos::Array<GO> gblColInds (opts.numEntPerRow);
    int lclSuccess = 1;
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblInd = rowMap.GID (lclRow);
      for (LO k = 0; k < static_cast<LO> (opts.numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
      }
      const int err = G.InsertGlobalIndices (gblInd, opts.numEntPerRow,
                                             gblColInds.getRawPtr ());
      if (err != 0) {
        lclSuccess = 0;
      }
    }

    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEUCHOS_TEST_FOR_EXCEPTION
      (gblSuccess != 1, std::logic_error, "While filling Epetra_CrsGraph, "
       "error on one or more MPI processes.");

    G.FillComplete ();

    RCP<matrix_type> A = rcp (new matrix_type (Copy, G));

    Teuchos::Array<SC> vals (opts.numEntPerRow);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblInd = rowMap.GID (lclRow);
      for (LO k = 0; k < static_cast<LO> (opts.numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
        // Cast first to MT, then SC, to avoid issues like
        // std::complex<double>'s constructor not taking int.
        vals[k] = static_cast<SC> (static_cast<MT> (gblColInds[k]));
      }
      const int err = A->ReplaceGlobalValues (gblInd, opts.numEntPerRow,
                                              vals.getRawPtr (),
                                              gblColInds.getRawPtr ());
      if (err != 0) {
        lclSuccess = 0;
      }
    }
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEUCHOS_TEST_FOR_EXCEPTION
      (gblSuccess != 1, std::logic_error, "While filling Epetra_CrsMatrix, "
       "error on one or more MPI processes.");

    A->FillComplete ();
    return A;
  }
#endif // HAVE_TPETRACORE_EPETRA

} // namespace (anonymous)

#if defined(CUDA_VERSION) && (CUDA_VERSION < 8000)
// mfh 20 Aug 2017: I've had troubles in the past with Kokkos functors
// in anonymous namespaces, so I've gotten in the habit of putting
// them in named namespaces.

namespace TpetraExample {
  template<class KokkosCrsMatrixType,
           class LocalOrdinalType>
  class ParallelReduceLoopBody {
  public:
    ParallelReduceLoopBody (const KokkosCrsMatrixType& A_lcl) :
      A_lcl_ (A_lcl)
    {}

    KOKKOS_FUNCTION void
    operator () (const LocalOrdinalType& lclRow, size_t& count) const
    {
      auto rowView = A_lcl_.row (lclRow);
      auto length  = rowView.length;
      count += static_cast<size_t> (length);
    }

    KokkosCrsMatrixType A_lcl_;
  };
} // namespace TpetraExample

#endif // defined(CUDA_VERSION) && (CUDA_VERSION < 8000)


int
main (int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::TimeMonitor;
  using std::endl;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();

    // Output stream 'out' will ignore output not from Process 0.
    RCP<Teuchos::FancyOStream> pOut = getOutputStream (*comm);
    Teuchos::FancyOStream& out = *pOut;

    // Read command-line options into the 'opts' struct.
    CmdLineOpts opts;
    {
      Teuchos::CommandLineProcessor clp;
      setCmdLineOpts (opts, clp);
      int result = parseCmdLineOpts (clp, argc, argv);
      if (result == 1) { // help printed
	return EXIT_SUCCESS;
      }
      else if (result == -1) { // parse error
	return EXIT_FAILURE;
      }
      result = checkCmdLineOpts (out, opts);
      if (result != 0) {
	return EXIT_FAILURE;
      }
    }

    out << "Command-line options:" << endl;
    {
      Teuchos::OSTab tab1 (out); // push one tab in this scope
      out << "numTrials: " << opts.numTrials << endl
	  << "lclNumRows: " << opts.lclNumRows << endl
	  << "numEntPerRow: " << opts.numEntPerRow << endl
	  << "testEpetra: " << opts.testEpetra << endl
	  << "testEpetraLen: " << opts.testEpetraLen << endl
	  << "testTpetra: " << opts.testTpetra << endl
	  << "testTpetraLen: " << opts.testTpetraLen << endl
	  << "testKokkos: " << opts.testKokkos << endl
	  << endl;
    }

    // The benchmark is supposed to be self-checking.
    const size_t expectedTotalLclNumEnt =
      static_cast<size_t> (opts.numTrials) *
      static_cast<size_t> (opts.lclNumRows) *
      static_cast<size_t> (opts.numEntPerRow);

    size_t totalLclNumEnt = 0;
    int lclSuccess = 1;
    int gblSuccess = 0; // not proven successful yet

    totalLclNumEnt = 0;
    if (opts.testEpetra) {
#ifdef HAVE_TPETRACORE_EPETRA
      typedef int LO;

      auto timer = TimeMonitor::getNewCounter ("Epetra ExtractMyRowView");
      RCP<const Epetra_CrsMatrix> A = getEpetraMatrix (comm, opts);
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	const LO lclNumRows = opts.lclNumRows;
	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    int numEnt;
	    double* val;
	    int* ind;
	    A->ExtractMyRowView (lclRow, numEnt, val, ind);
	    totalLclNumEnt += numEnt;
	  }
	}
      }
#else
      // We've already checked this case when checking the command-line arguments.
      TEUCHOS_TEST_FOR_EXCEPTION
	(true, std::logic_error, "Epetra not enabled; should never get here!");
#endif // HAVE_TPETRACORE_EPETRA
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Epetra ExtractMyRowView validation FAILED.  On my process, "
	"totalLclNumEnt = " << totalLclNumEnt << " != expectedTotalLclNumEnt = "
	  << expectedTotalLclNumEnt << "." << endl;
    }

    totalLclNumEnt = 0;
    if (opts.testEpetraLen) {
#ifdef HAVE_TPETRACORE_EPETRA
      typedef int LO;

      auto timer = TimeMonitor::getNewCounter ("Epetra NumMyEntries");
      RCP<const Epetra_CrsMatrix> A_ptr = getEpetraMatrix (comm, opts);
      TEUCHOS_TEST_FOR_EXCEPTION
	(A_ptr.is_null (), std::logic_error, "getEpetraMatrix returned null!  "
	 "This should never happen.");
      const Epetra_CrsMatrix& A = *A_ptr;
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	const LO lclNumRows = opts.lclNumRows;
	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    const size_t len = static_cast<size_t> (A.NumMyEntries (lclRow));
	    totalLclNumEnt += len;
	  }
	}
      }
#else
      // We've already checked this case when checking the command-line arguments.
      TEUCHOS_TEST_FOR_EXCEPTION
	(true, std::logic_error, "Epetra not enabled; should never get here!");
#endif // HAVE_TPETRACORE_EPETRA
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Epetra NumMyEntries validation FAILED.  On my process, "
	"totalLclNumEnt = " << totalLclNumEnt << " != expectedTotalLclNumEnt = "
	  << expectedTotalLclNumEnt << "." << endl;
    }

    totalLclNumEnt = 0;
    if (opts.testEpetraLen) {
#ifdef HAVE_TPETRACORE_EPETRA
      typedef int LO;

      auto timer = TimeMonitor::getNewCounter ("Epetra NumMyRowEntries");
      RCP<const Epetra_CrsMatrix> A_ptr = getEpetraMatrix (comm, opts);
      TEUCHOS_TEST_FOR_EXCEPTION
	(A_ptr.is_null (), std::logic_error, "getEpetraMatrix returned null!  "
	 "This should never happen.");
      const Epetra_CrsMatrix& A = *A_ptr;
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	const LO lclNumRows = opts.lclNumRows;
	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    int numEnt;
	    (void) A.NumMyRowEntries (lclRow, numEnt); // ignore error code
	    totalLclNumEnt += static_cast<size_t> (numEnt);
	  }
	}
      }
#else
      // We've already checked this case when checking the command-line arguments.
      TEUCHOS_TEST_FOR_EXCEPTION
	(true, std::logic_error, "Epetra not enabled; should never get here!");
#endif // HAVE_TPETRACORE_EPETRA
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Epetra NumMyRowEntries validation FAILED.  On my process, "
	"totalLclNumEnt = " << totalLclNumEnt << " != expectedTotalLclNumEnt = "
	  << expectedTotalLclNumEnt << "." << endl;
    }

    totalLclNumEnt = 0;
    if (opts.testTpetra) {
      typedef Tpetra::CrsMatrix<>::local_ordinal_type LO;

      auto timer = TimeMonitor::getNewCounter ("Tpetra getLocalRowView");
      RCP<const Tpetra::CrsMatrix<> > A_ptr = getTpetraMatrix (comm, opts);
      TEUCHOS_TEST_FOR_EXCEPTION
	(A_ptr.is_null (), std::logic_error, "getTpetraMatrix returned null!  "
	 "This should never happen.");
      const Tpetra::CrsMatrix<>& A = *A_ptr;
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  const LO lclNumRows = opts.lclNumRows;
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    typename Tpetra::CrsMatrix<>::local_inds_host_view_type ind;
	    typename Tpetra::CrsMatrix<>::values_host_view_type val;
	    A.getLocalRowView (lclRow, ind, val);
	    const size_t len = static_cast<size_t> (ind.size ());
	    totalLclNumEnt += len;
	  }
	}
      }
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Tpetra::CrsMatrix::getLocalRowView validation FAILED.  "
	"On my process, totalLclNumEnt = "
	<< totalLclNumEnt << " != expectedTotalLclNumEnt = "
	<< expectedTotalLclNumEnt << "." << endl;
    }

    totalLclNumEnt = 0;
    if (opts.testTpetraLen) {
      using LO = Tpetra::CrsMatrix<>::local_ordinal_type;
      auto timer =
	TimeMonitor::getNewCounter ("Tpetra getNumEntriesInLocalRow");
      auto A_ptr = getTpetraMatrix (comm, opts);
      TEUCHOS_TEST_FOR_EXCEPTION
	(A_ptr.is_null (), std::logic_error, "getTpetraMatrix "
	 "returned null!  "
	 "This should never happen.");
      const Tpetra::CrsMatrix<>& A = *A_ptr;
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  const LO lclNumRows = opts.lclNumRows;
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    const size_t len = A.getNumEntriesInLocalRow (lclRow);
	    totalLclNumEnt += len;
	  }
	}
      }
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
			 outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Tpetra::CrsMatrix::getNumEntriesInLocalRow validation "
	"FAILED.  On my process, totalLclNumEnt = "
	<< totalLclNumEnt << " != expectedTotalLclNumEnt = "
	<< expectedTotalLclNumEnt << "." << endl;
    }

    totalLclNumEnt = 0;
    if (opts.testKokkos) {
      using LO = Tpetra::CrsMatrix<>::local_ordinal_type;

      auto timer = TimeMonitor::getNewCounter ("Kokkos sequential");
      auto A = getTpetraMatrix (comm, opts);
      auto A_lcl = A->getLocalMatrixDevice ();
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  const LO lclNumRows = opts.lclNumRows;
	  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
	    auto rowView = A_lcl.row (lclRow);
	    auto len = rowView.length;

	    (void) rowView;
	    totalLclNumEnt += static_cast<size_t> (len);
	  }
	}
      }
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
			 outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Kokkos sequential loop validation FAILED.  "
	"On my process, totalLclNumEnt = " << totalLclNumEnt << " != "
	"expectedTotalLclNumEnt = " << expectedTotalLclNumEnt << "."
	<< endl;
    }

    totalLclNumEnt = 0;
    if (opts.testKokkos) {
      using Kokkos::parallel_reduce;
      using LO = Tpetra::CrsMatrix<>::local_ordinal_type;
      using DT = Tpetra::CrsMatrix<>::device_type;
      using host_execution_space =
	Kokkos::View<LO*, DT>::HostMirror::execution_space;
      using policy_type = Kokkos::RangePolicy<host_execution_space, LO>;

      auto timer = TimeMonitor::getNewCounter ("Kokkos parallel");
      auto A = getTpetraMatrix (comm, opts);
      auto A_lcl = A->getLocalMatrixDevice ();
      { // Start timing after matrix creation
	TimeMonitor timeMon (*timer);

	for (int trial = 0; trial < opts.numTrials; ++trial) {
	  policy_type range (0, static_cast<LO> (opts.lclNumRows));

	  size_t curTotalLclNumEnt = 0;
#if ! defined(CUDA_VERSION) || (CUDA_VERSION >= 8000)
	  parallel_reduce ("loop", range,
			   KOKKOS_LAMBDA (const LO& lclRow, size_t& count) {
			     auto rowView = A_lcl.row (lclRow);
			     auto length  = rowView.length;
			     count += static_cast<size_t> (length);
			   }, curTotalLclNumEnt);
#else
	  using TpetraExample::ParallelReduceLoopBody;
	  typedef ParallelReduceLoopBody<decltype (A_lcl), LO> functor_type;
	  parallel_reduce ("loop", range,
			   functor_type (A_lcl),
			   curTotalLclNumEnt);
#endif // ! defined(CUDA_VERSION) || (CUDA_VERSION >= 8000)
	  totalLclNumEnt += curTotalLclNumEnt;
	}
      }
    }
    lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
			 outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Kokkos host parallel loop validation FAILED.  "
	"On my process, totalLclNumEnt = "
	<< totalLclNumEnt << " != expectedTotalLclNumEnt = "
	<< expectedTotalLclNumEnt << "." << endl;
    }

    TimeMonitor::report (comm.ptr (), out);
  }
  return EXIT_SUCCESS;
}

