// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This benchmark exercises Tpetra::BlockCrsMatrix's
// sparse matrix-vector multiply (the apply() method).  Tpetra
// implements sparse matrix and dense vector data structures and
// computational kernels for users and other Trilinos data structures.
// Tpetra uses MPI (Message Passing Interface) for distributed-memory
// parallelism, and Kokkos for shared-memory parallelism within an MPI
// process.

#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_BlockView.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

#include "Kokkos_Random.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"

namespace { // (anonymous)

template<class Scalar, class LO, class GO, class Node>
void
localApplyBlockNoTrans (Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>& A,
                        Tpetra::BlockMultiVector<Scalar, LO, GO, Node>& X,
                        Tpetra::BlockMultiVector<Scalar, LO, GO, Node>& Y,
                        const Scalar& alpha,
                        const Scalar& beta)
{
  using Tpetra::COPY;
  using Tpetra::FILL;
  using Tpetra::GEMV;
  using Tpetra::SCAL;
  typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>
    block_crs_matrix_type;
  typedef typename block_crs_matrix_type::impl_scalar_type IST;
  typedef Kokkos::ArithTraits<IST> KAT;
  typedef typename block_crs_matrix_type::little_vec_type little_vec_type;
  typedef typename block_crs_matrix_type::little_block_type little_blk_type;

  const auto G = A.getCrsGraph ();

  const IST zero = KAT::zero ();
  const IST one = KAT::one ();
  const LO numLocalMeshRows =
    static_cast<LO> (G.getRowMap ()->getLocalNumElements ());
  const LO numVecs = static_cast<LO> (X.getNumVectors ());
  const LO blockSize = A.getBlockSize ();

  // NOTE (mfh 01 Jun 2016) This is a host code, so we have to sync
  // all the objects to host.  We sync them back to device after we're
  // done.  That's why all objects come in by nonconst reference.
  // A.sync_host ();
  // X.sync_host ();
  // Y.sync_host ();
  // Y.modify_host (); // only Y gets modified here

  // Get the matrix values.  Blocks are stored contiguously, each
  // block in row-major order (Kokkos::LayoutRight).
  auto val = A.getValuesHostNonConst ();

  auto gblGraph = A.getCrsGraph ();
  auto lclGraph = G.getLocalGraphHost ();
  auto ptrHost = lclGraph.row_map;
  auto indHost = lclGraph.entries;
  Teuchos::Array<IST> localMem (blockSize);
  little_vec_type Y_lcl (localMem.getRawPtr (), blockSize, 1);

  for (LO j = 0; j < numVecs; ++j) {
    for (LO lclRow = 0; lclRow < numLocalMeshRows; ++lclRow) {
      auto Y_cur = Y.getLocalBlockHost (lclRow, j, Tpetra::Access::ReadWrite);
      if (beta == zero) {
        FILL (Y_lcl, zero);
      } else if (beta == one) {
        COPY (Y_cur, Y_lcl);
      } else {
        COPY (Y_cur, Y_lcl);
        SCAL (beta, Y_lcl);
      }

      const size_t meshBeg = ptrHost[lclRow];
      const size_t meshEnd = ptrHost[lclRow+1];
      for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
        const LO meshCol = indHost[absBlkOff];

        auto A_cur_1d = Kokkos::subview (val, absBlkOff * blockSize * blockSize);
        little_blk_type A_cur (A_cur_1d.data (), blockSize, blockSize);
        auto X_cur = X.getLocalBlockHost (meshCol, j, Tpetra::Access::ReadOnly);

        GEMV (alpha, A_cur, X_cur, Y_lcl); // Y_lcl += alpha*A_cur*X_cur
      } // for each entry in the current local row of the matrix

      COPY (Y_lcl, Y_cur);
    } // for each local row of the matrix
  } // for each column j of the input / output block multivector

  // Sync everything back to device when we're done.  This only
  // actually copies Y back to device, but it ensures that all the
  // modified flags are right.
  // A.template sync<device_memory_space> ();
  // X.template sync<device_memory_space> ();
  // Y.template sync<device_memory_space> ();
}


template<class Scalar, class LO, class GO, class Node>
bool
compareLocalMatVec (Teuchos::FancyOStream& out,
                    Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>& A,
                    Tpetra::MultiVector<Scalar, LO, GO, Node>& X_mv,
                    Tpetra::MultiVector<Scalar, LO, GO, Node>& Y_mv,
                    const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::one (),
                    const Scalar& beta = Teuchos::ScalarTraits<Scalar>::zero ())
{
  using std::endl;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
  typedef typename MV::mag_type mag_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Teuchos::ScalarTraits<mag_type> STM;

  Teuchos::OSTab tab0 (out);
  out << "Test BlockCrsMatrix::apply by comparing against hand-rolled "
    "sequential code" << endl;
  Teuchos::OSTab tab1 (out);

  const LO numVecs = Y_mv.getNumVectors ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (numVecs != static_cast<LO> (X_mv.getNumVectors ()),
     std::invalid_argument,
     "X_mv and Y_mv must have the same number of columns.");

  const auto G = A.getCrsGraph ();
  const size_t lclNumMeshRows = G.getRowMap ()->getLocalNumElements ();
  const LO blockSize = A.getBlockSize ();
  const size_t maxNumTermsInRowSum =
    static_cast<size_t> (G.getLocalMaxNumRowEntries ()) *
    static_cast<size_t> (blockSize);
  const mag_type tol =
    STM::squareroot (static_cast<mag_type> (maxNumTermsInRowSum)) *
    STS::eps ();

  out << "Number of mesh rows in A: " << lclNumMeshRows << endl
      << "Block size: " << blockSize << endl
      << "Number of columns in X and Y: " << numVecs << endl
      << "Test tolerance: " << tol << endl;

  auto meshRangeMap = G.getRangeMap ();
  auto meshDomainMap = G.getDomainMap ();

  BMV X (X_mv, *meshDomainMap, blockSize);
  // Make a copy of X too, just in case the code has a bug that makes
  // it overwrite X.
  MV X_mv_copy (X_mv, Teuchos::Copy);
  BMV X_copy (X_mv_copy, *meshDomainMap, blockSize);

  MV Y_mv_copy (Y_mv, Teuchos::Copy);
  BMV Y_copy (Y_mv_copy, *meshRangeMap, blockSize);

  Teuchos::Array<mag_type> norms (numVecs);
  X_mv.normInf (norms);
  {
    out << "Make sure copies worked" << endl;
    Teuchos::OSTab tab2 (out);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||X_mv(:," << j << ")||_inf = " << norms[j] << endl;
    }
    X_mv_copy.normInf (norms);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||X_mv_copy(:," << j << ")||_inf = " << norms[j] << endl;
    }
  }

  out << "Call A.apply(...)" << endl;
  {
    Teuchos::OSTab tab2 (out);
    A.apply (X_mv, Y_mv, Teuchos::NO_TRANS, alpha, beta);
    Y_mv.normInf (norms);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||Y_mv(:," << j << ")||_inf = " << norms[j] << endl;
    }
    Y_mv.norm2 (norms);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||Y_mv(:," << j << ")||_2 = " << norms[j] << endl;
    }
  }

  out << "Call the hand-rolled reference implementation" << endl;
  {
    Teuchos::OSTab tab2 (out);
    localApplyBlockNoTrans (A, X_copy, Y_copy, alpha, beta);
    Y_mv_copy.normInf (norms);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||Y_mv_copy(:," << j << ")||_inf = " << norms[j] << endl;
    }
    Y_mv_copy.norm2 (norms);
    for (LO j = 0; j < numVecs; ++j) {
      out << "||Y_mv_copy(:," << j << ")||_2 = " << norms[j] << endl;
    }
  }

  out << "Compare results" << endl;
  bool success = true;
  {
    Teuchos::OSTab tab2 (out);

    Y_mv_copy.update (STS::one (), Y_mv, -STS::one ());
    Y_mv_copy.normInf (norms);
    for (LO j = 0; j < numVecs; ++j) {
      if (norms[j] > tol) {
        out << "||Y(:," << j << ") - Y_copy(:," << j << ")||_inf = "
            << norms[j] << " > tol = " << tol << std::endl;
        success = false;
      }
    }

    Y_mv_copy.norm2 (norms);
    for (LO j = 0; j < numVecs; ++j) {
      if (norms[j] > tol) {
        out << "||Y(:," << j << ") - Y_copy(:," << j << ")||_2 = "
            << norms[j] << " > tol = " << tol << std::endl;
      }
    }
  }

  if (success) {
    out << "SUCCESS" << endl;
  }
  else {
    out << "FAILURE" << endl;
  }
  return success;
}

// Options to read in from the command line
struct CmdLineOpts {
  // Do the benchmark this many times in a single timing loop, in case
  // the timer's granularity is too coarse to capture run time to
  // adequate precision.
  int numTrials;
  // Number of rows per MPI process (hence "local") in the graph;
  // number of block rows per MPI process in the BlockCrsMatrix.
  int lclNumRows;
  // Number of entries per row in the sparses graph; thus, number of
  // blocks per block row of the BlockCrsMatrix.
  int numEntPerRow;
  // Block size (number of rows / columns per block).  Applications
  // hardly ever want powers of two.  Only optimizing your library for
  // powers of two block sizes is tacky.
  int blockSize;
  // Instead of doing the benchmark, run a test: Compare results of
  // BlockCrsMatrix::apply to those of a hand-rolled non-threaded
  // code.
  bool runTest;
};

// Use a utility from the Teuchos package of Trilinos to set up
// command-line options for reading, and set default values of
// command-line options.  clp is an output argument containing the
// set-up options.  It retains pointers to fields in 'opts'.  Reading
// the command-line options will update those fields in place.
void
setCmdLineOpts (CmdLineOpts& opts,
                Teuchos::CommandLineProcessor& clp)
{
  // Set default values of command-line options.

  opts.numTrials = 200;
  opts.lclNumRows = 10000;
  opts.numEntPerRow = 10;
  opts.blockSize = 7;
  opts.runTest = false;

  clp.setOption ("numTrials", &(opts.numTrials), "Number of trials per "
                 "timing loop (to increase timer precision)");
  clp.setOption ("lclNumRows", &(opts.lclNumRows), "Number of rows per MPI "
                 "process in the sparse graph; that is, the number of "
                 "block rows in the BlockCrsMatrix");
  clp.setOption ("numEntPerRow", &(opts.numEntPerRow), "Number of entries "
                 "per row in the sparse graph; that is, number of blocks "
                 "per block row of the BlockCrsMatrix");
  clp.setOption ("blockSize", &(opts.blockSize), "Block size; number of rows "
                 "/ columns per block");
  clp.setOption ("runTest", "runBenchmark", &(opts.runTest), "Instead of "
                 "benchmarking, test BlockCrsMatrix::apply");
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
// using the LAPACK error reporting convention of the negative of the
// argument in its original order (starting with 1) as the error code.
// Print informative error messages to the given output stream \c out.
int
checkCmdLineOpts (std::ostream& out,
                  const Teuchos::Comm<int>& comm,
                  const CmdLineOpts& opts)
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
  if (opts.blockSize < 0) {
    out << "blockSize = " << opts.blockSize << " < 0." << std::endl;
    err = -4; // LAPACK error reporting convention
  }
  if (opts.runTest && comm.getSize () > 1) {
    out << "runTest is true, but the number of processes in the communicator "
      "is > 1.  Currently, this test only works if the number of processes is"
      " 1." << std::endl;
    err = -5; // LAPACK error reporting convention
  }

  return err;
}

// Print values of the command-line options, as read in by
// parseCmdLineOpts, to the given output stream.
void
printCmdLineOpts (Teuchos::FancyOStream& out,
                  const CmdLineOpts& opts)
{
  using std::endl;
  // Convention for FancyOStream is to push one tab before printing in
  // a scope.  OSTab pops the tab when leaving the scope.
  Teuchos::OSTab tab1 (out);
  out << "numTrials: " << opts.numTrials << endl
      << "lclNumRows: " << opts.lclNumRows << endl
      << "numEntPerRow: " << opts.numEntPerRow << endl
      << "blockSize: " << opts.blockSize << endl
      << "runTest: " << (opts.runTest ? "true" : "false") << endl
      << endl;
}

// Return a pointer (RCP is like std::shared_ptr) to an output stream.
// It prints on Process 0 of the given MPI communicator, but ignores
// all output on other MPI processes.
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

// Get a Tpetra::CrsGraph for use in benchmarks.  This method takes
// parameters that come from the command-line options read in by
// parseCmdLineOpts.
Teuchos::RCP<Tpetra::CrsGraph<> >
getTpetraGraph (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                const CmdLineOpts& opts)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::CrsGraph<> graph_type;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  const LO lclNumRows = opts.lclNumRows;
  const GST gblNumRows = static_cast<GST> (opts.lclNumRows) *
    static_cast<GST> (comm->getSize ());
  const GO indexBase = 0;

  // A Map describes a distribution of data over MPI processes.
  // This "row Map" will describe the distribution of rows of the
  // sparse graph that we will create.
  RCP<const map_type> rowMap =
    rcp (new map_type (gblNumRows, static_cast<size_t> (lclNumRows),
                       indexBase, comm));
  const GO gblNumCols = static_cast<GO> (rowMap->getGlobalNumElements ());
  // Create the graph structure of the sparse matrix.
  RCP<graph_type> G =
    rcp (new graph_type (rowMap, opts.numEntPerRow));
  // Fill in the sparse graph.
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
  // Make the graph ready for use by (Block)CrsMatrix.
  G->fillComplete ();
  return G;
}

// Get a Tpetra::BlockCrsMatrix for use in benchmarks.
// This method takes the result of getTpetraGraph() (above) and
// parameters that come from the command-line options read in by
// parseCmdLineOpts.
Teuchos::RCP<Tpetra::BlockCrsMatrix<> >
getTpetraBlockCrsMatrix (Teuchos::FancyOStream& out,
                         const Teuchos::RCP<const Tpetra::CrsGraph<> >& graph,
                         const CmdLineOpts& opts)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::BlockCrsMatrix<> matrix_type;
  typedef matrix_type::impl_scalar_type SC;
  typedef Kokkos::ArithTraits<SC> KAT;
  typedef Tpetra::Map<>::local_ordinal_type LO;

  // mfh 02 Jun 2016: Prefer Kokkos::Serial to the HostMirror as the
  // pseudorandom number generator's execution space.  This is
  // because, for CudaUVMSpace, the HostMirror is the same as the
  // original.  This causes segfaults in the pseudorandom number
  // generator, due to CUDA code trying to access host memory.
#ifdef KOKKOS_ENABLE_SERIAL
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> host_device_type;
#else
  typedef Kokkos::View<SC**, Kokkos::LayoutRight, device_type>::host_mirror_space host_device_type;
#endif // KOKKOS_ENABLE_SERIAL
  typedef host_device_type::execution_space host_execution_space;
  typedef Kokkos::View<SC**, Kokkos::LayoutRight, host_device_type> block_type;

  // We're filling on the host, so generate random numbers on the host.
  typedef Kokkos::Random_XorShift64_Pool<host_execution_space> pool_type;

  Teuchos::OSTab tab0 (out);
  out << "Create BlockCrsMatrix for "
      << (opts.runTest ? "test" : "benchmark") << endl;
  Teuchos::OSTab tab1 (out);

  const auto meshRowMap = * (graph->getRowMap ());
  // Contrary to expectations, asking for the graph's number of
  // columns, or asking the column Map for the number of entries,
  // won't give the correct number of columns in the graph.
  // const GO gblNumCols = graph->getDomainMap ()->getGlobalNumElements ();
  const LO lclNumRows = meshRowMap.getLocalNumElements ();
  const LO blkSize = opts.blockSize;

  RCP<matrix_type> A = rcp (new matrix_type (*graph, blkSize));

  // We're filling on the host.
  // A->sync_host ();
  // A->modify_host ();

  // This only matters if filling with random values.  We only do that
  // if opts.runTest is true (that is, if we're testing correctness of
  // BlockCrsMatrix::apply, instead of benchmarking it).
  const uint64_t myRank =
    static_cast<uint64_t> (graph->getMap ()->getComm ()->getRank ());
  const uint64_t seed64 = static_cast<uint64_t> (std::rand ()) + myRank + 17311uLL;
  const unsigned int seed = static_cast<unsigned int> (seed64&0xffffffff);
  pool_type randPool (seed);

  // Create a "prototype block" of values to use when filling the
  // block sparse matrix.
  block_type curBlk ("curBlk", blkSize, blkSize);
  // We only use this if filling with random values.
  Kokkos::View<SC*, Kokkos::LayoutRight, host_device_type,
    Kokkos::MemoryUnmanaged> curBlk_1d (curBlk.data (),
                                        blkSize * blkSize);
  if (! opts.runTest) {
    // For benchmarks, we don't care so much about the values; we just
    // want them not to be Inf or NaN, in case the processor makes the
    // unfortunate choice to handle arithmetic with those via traps.
    for (LO j = 0; j < blkSize; ++j) {
      for (LO i = 0; i < blkSize; ++i) {
        curBlk(i,j) = 1.0;
      }
    }
  }

  // Fill in the block sparse matrix.
  out << "Fill the BlockCrsMatrix" << endl;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) { // for each of my rows
    Tpetra::CrsGraph<>::local_inds_host_view_type lclColInds;
    graph->getLocalRowView (lclRow, lclColInds);

    // Put some entries in the matrix.
    for (LO k = 0; k < static_cast<LO> (lclColInds.size ()); ++k) {
      if (opts.runTest) {
        // Fill the current block with random values between -1 and 1.
        Kokkos::fill_random (curBlk_1d, randPool, -KAT::one (), KAT::one ());
      }
      const LO lclColInd = lclColInds[k];
      const LO err =
        A->replaceLocalValues (lclRow, &lclColInd, curBlk.data (), 1);
      TEUCHOS_TEST_FOR_EXCEPTION(err != 1, std::logic_error, "Bug");
    }
  }

  // We're done filling on the host, so sync to device.
  //A->sync<device_type::memory_space> ();
  return A;
}

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::TimeMonitor;
  using std::endl;
  typedef Tpetra::Vector<>::scalar_type SC;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  bool success = true;
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
      result = checkCmdLineOpts (out, *comm, opts);
      if (result != 0) {
        return EXIT_FAILURE;
      }
    }

    out << "Command-line options:" << endl;
    printCmdLineOpts (out, opts);

    auto G = getTpetraGraph (comm, opts);
    auto A = getTpetraBlockCrsMatrix (out, G, opts);
    Tpetra::Vector<> X (A->getDomainMap ());
    Tpetra::Vector<> Y (A->getRangeMap ());

    // Fill X with values that don't increase the max-norm of results.
    // That way, repeated mat-vecs won't overflow.  This matters
    // because some processors do a silly thing and handle Inf or NaN
    // (or even denorms) via traps.  This is very expensive, so if the
    // norms increase or decrease a lot, that might trigger the slow
    // case.
    const SC X_val = static_cast<SC> (1.0) /
      static_cast<SC> (opts.numEntPerRow * opts.blockSize);
    X.putScalar (X_val);
    Y.putScalar (0.0);

    if (opts.runTest) {
      const bool lclSuccess = compareLocalMatVec (out, *A, X, Y);
      success = success && lclSuccess;
    }
    else {
      auto timer =
        TimeMonitor::getNewCounter ("Tpetra BlockCrsMatrix apply (mat-vec)");
      {
        TimeMonitor timeMon (*timer);
        for (int trial = 0; trial < opts.numTrials; ++trial) {
          A->apply (X, Y);
        }
      }
    }

    if (! opts.runTest) {
      TimeMonitor::report (comm.ptr (), out);
    }
  }

  if (success) {
    return EXIT_SUCCESS;
  }
  else {
    return EXIT_FAILURE;
  }
}

