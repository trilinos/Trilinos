
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_TEUCHOS_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI

#ifdef HAVE_TPETRACORE_EPETRA
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#endif // HAVE_TPETRACORE_EPETRA

#include "Epetra_SerialComm.h"
#ifdef HAVE_TEUCHOS_MPI
#include "Epetra_MpiComm.h"
#endif // HAVE_TEUCHOS_MPI


namespace { // (anonymous)

  struct CmdLineOpts {
    int numTrials;
    int lclNumRows;
    int numEntPerRow;
    bool testEpetra;
    bool testEpetraLen;
    bool testTpetra;
    bool testTpetraLen;
    bool testKokkos;
  };

  void
  setCmdLineOpts (CmdLineOpts& opts,
                  Teuchos::CommandLineProcessor& clp)
  {
    opts.numTrials = 10;
    opts.lclNumRows = 10000;
    opts.numEntPerRow = 10;
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
                   "Number of trials per timing loop");
    clp.setOption ("lclNumRows", &(opts.lclNumRows), "Number of rows in the "
                   "matrix per MPI process");
    clp.setOption ("numEntPerRow", &(opts.numEntPerRow),
                   "Number of entries per row");
    clp.setOption ("testEpetra", "skipEpetra", &(opts.testEpetra),
                   "Whether to test Epetra_CrsMatrix::ExtractMyRowView");
    clp.setOption ("testEpetraLen", "skipEpetraLen", &(opts.testEpetraLen),
                   "Whether to test Epetra_CrsMatrix::NumMyEntries");
    clp.setOption ("testTpetra", "skipTpetra", &(opts.testTpetra),
                   "Whether to test Tpetra::CrsMatrix::getLocalRowView");
    clp.setOption ("testTpetraLen", "skipTpetraLen", &(opts.testTpetraLen),
                   "Whether to test Tpetra::CrsMatrix::getNumEntriesInLocalRow");
  }

  // Return 0 if successful, 1 if help printed (indicates that the
  // application shouldn't run), and -1 on error.
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

  // Return 0 if all correct, else return nonzero.  Print informative
  // error messages to the given output stream \c out.
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

  Teuchos::RCP<Teuchos::FancyOStream>
  getOutputStream (const Teuchos::Comm<int>& comm)
  {
    using Teuchos::getFancyOStream;

    const int myRank = comm.getRank ();
    if (myRank == 0) {
      return getFancyOStream (Teuchos::rcpFromRef (std::cout));
    }
    else {
      return getFancyOStream (Teuchos::rcp (new Teuchos::oblackholestream ()));
    }
  }

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

    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, static_cast<size_t> (lclNumRows),
                         indexBase, comm));
    const GO gblNumCols = static_cast<GO> (rowMap->getGlobalNumElements ());
    RCP<graph_type> G =
      rcp (new graph_type (rowMap, opts.numEntPerRow,
                           Tpetra::StaticProfile));

    Teuchos::Array<GO> gblColInds (opts.numEntPerRow);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblInd = rowMap->getGlobalElement (lclRow);
      for (LO k = 0; k < static_cast<LO> (opts.numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
      }
      G->insertGlobalIndices (gblInd, gblColInds ());
    }
    G->fillComplete ();

    RCP<matrix_type> A = rcp (new matrix_type (G));

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
    A->fillComplete ();

    return A;
  }


  Teuchos::RCP<const Epetra_Comm>
  makeEpetraComm (const Teuchos::Comm<int>& comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

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

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::TimeMonitor;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  RCP<Teuchos::FancyOStream> pOut = getOutputStream (*comm);
  Teuchos::FancyOStream& out = *pOut;

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

  const size_t expectedTotalLclNumEnt =
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
    out << "Epetra ExtractMyRowView validation FAILED" << std::endl;
  }

  totalLclNumEnt = 0;
  if (opts.testEpetraLen) {
#ifdef HAVE_TPETRACORE_EPETRA
    typedef int LO;

    auto timer = TimeMonitor::getNewCounter ("Epetra NumMyEntries");
    RCP<const Epetra_CrsMatrix> A = getEpetraMatrix (comm, opts);
    { // Start timing after matrix creation
      TimeMonitor timeMon (*timer);

      const LO lclNumRows = opts.lclNumRows;
      for (int trial = 0; trial < opts.numTrials; ++trial) {
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const size_t len = static_cast<size_t> (A->NumMyEntries (lclRow));
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
    out << "Epetra NumMyEntries validation FAILED" << std::endl;
  }

  totalLclNumEnt = 0;
  if (opts.testTpetra) {
    typedef Tpetra::CrsMatrix<>::scalar_type SC;
    typedef Tpetra::CrsMatrix<>::local_ordinal_type LO;

    auto timer = TimeMonitor::getNewCounter ("Tpetra getLocalRowView");
    RCP<const Tpetra::CrsMatrix<> > A = getTpetraMatrix (comm, opts);
    { // Start timing after matrix creation
      TimeMonitor timeMon (*timer);

      for (int trial = 0; trial < opts.numTrials; ++trial) {
        const LO lclNumRows = opts.lclNumRows;
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          Teuchos::ArrayView<const LO> ind;
          Teuchos::ArrayView<const SC> val;
          A->getLocalRowView (lclRow, ind, val);
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
    out << "Tpetra validation FAILED" << std::endl;
  }

  totalLclNumEnt = 0;
  if (opts.testTpetraLen) {
    typedef Tpetra::CrsMatrix<>::local_ordinal_type LO;

    auto timer = TimeMonitor::getNewCounter ("Tpetra getNumEntriesInLocalRow");
    RCP<const Tpetra::CrsMatrix<> > A = getTpetraMatrix (comm, opts);
    { // Start timing after matrix creation
      TimeMonitor timeMon (*timer);

      for (int trial = 0; trial < opts.numTrials; ++trial) {
        const LO lclNumRows = opts.lclNumRows;
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const size_t len = A->getNumEntriesInLocalRow (lclRow);
          totalLclNumEnt += len;
        }
      }
    }
  }
  lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    out << "Tpetra validation FAILED" << std::endl;
  }

  totalLclNumEnt = 0;
  if (opts.testKokkos) {
    typedef Tpetra::CrsMatrix<>::local_ordinal_type LO;

    auto timer = TimeMonitor::getNewCounter ("Kokkos sequential");
    auto A = getTpetraMatrix (comm, opts);
    auto A_lcl = A->getLocalMatrix ();
    { // Start timing after matrix creation
      TimeMonitor timeMon (*timer);

      for (int trial = 0; trial < opts.numTrials; ++trial) {
        const LO lclNumRows = opts.lclNumRows;
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          auto rowView = A_lcl.template row<LO> (lclRow);
          auto len = rowView.length;

          (void) rowView;
          totalLclNumEnt += static_cast<size_t> (len);
        }
      }
    }
  }
  lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    out << "Kokkos sequential loop validation FAILED" << std::endl;
  }

  totalLclNumEnt = 0;
  if (opts.testKokkos) {
    using Kokkos::parallel_reduce;
    typedef Tpetra::CrsMatrix<>::local_ordinal_type LO;
    typedef Tpetra::CrsMatrix<>::device_type DT;
    typedef Kokkos::View<LO*, DT>::HostMirror::execution_space host_execution_space;
    typedef Kokkos::RangePolicy<host_execution_space, LO> policy_type;

    auto timer = TimeMonitor::getNewCounter ("Kokkos parallel");
    auto A = getTpetraMatrix (comm, opts);
    auto A_lcl = A->getLocalMatrix ();
    { // Start timing after matrix creation
      TimeMonitor timeMon (*timer);

      for (int trial = 0; trial < opts.numTrials; ++trial) {
        policy_type range (0, static_cast<LO> (opts.lclNumRows));
        parallel_reduce ("loop", range,
                         KOKKOS_LAMBDA (const LO& lclRow, size_t& count) {
            auto rowView = A_lcl.template row<LO> (lclRow);
            auto length  = rowView.length;
            count += static_cast<size_t> (length);
          });
      }
    }
  }
  lclSuccess = (totalLclNumEnt == expectedTotalLclNumEnt) ? 1 : 0;
  gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    out << "Kokkos parallel loop validation FAILED" << std::endl;
  }

  TimeMonitor::report (comm.ptr (), out);
  return EXIT_SUCCESS;
}

