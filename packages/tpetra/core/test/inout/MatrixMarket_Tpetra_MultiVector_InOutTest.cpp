// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <algorithm>

namespace {
  // Ensure that X and Y have the same dimensions, and that the
  // corresponding columns of X and Y have 2-norms that differ by no
  // more than a prespecified tolerance (that accounts for rounding
  // errors in computing the 2-norm).
  template<class MV>
  void
  assertMultiVectorsEqual (const Teuchos::RCP<const MV>& X,
                           const Teuchos::RCP<const MV>& Y)
  {
    using Teuchos::as;
    using Teuchos::RCP;
    using std::cerr;
    using std::cout;
    using std::endl;

    typedef typename MV::scalar_type scalar_type;

    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    TEUCHOS_TEST_FOR_EXCEPTION(X->getGlobalLength() != Y->getGlobalLength(),
      std::logic_error, "Y has a different number of rows than X.");
    TEUCHOS_TEST_FOR_EXCEPTION(X->getNumVectors() != Y->getNumVectors(),
      std::logic_error, "Y has a different number of columns than X.");

    Tpetra::global_size_t numRows = X->getGlobalLength();
    const size_t numVecs = X->getNumVectors();
    Teuchos::Array<magnitude_type> X_norm2 (numVecs);
    Teuchos::Array<magnitude_type> Y_norm2 (numVecs);
    X->norm2 (X_norm2 ());
    Y->norm2 (Y_norm2 ());

    // For the relative tolerance, I'm using the typical heuristic: a
    // (fudge factor) times (machine epsilon) times sqrt(the number of
    // floating-point numbers involved in computing the norm for one
    // column).  Our output routine is careful to use enough digits,
    // so the input matrix shouldn't be that much different.
    const magnitude_type tol = as<magnitude_type> (10) *
      STS::magnitude (STS::eps ()) *
      STM::squareroot (as<magnitude_type> (numRows));
    std::vector<size_t> badColumns;
    for (size_t j = 0; j < numVecs; ++j) {
      // If the norm of the current column of X is zero, use the
      // absolute difference; otherwise use the relative difference.
      if ((X_norm2[j] == STM::zero() && STS::magnitude (Y_norm2[j]) > tol) ||
          STS::magnitude (X_norm2[j] - Y_norm2[j]) > tol) {
        badColumns.push_back (j);
      }
    }

    if (badColumns.size() > 0) {
      const size_t numBad = badColumns.size();
      std::ostringstream os;
      std::copy (badColumns.begin(), badColumns.end(),
                 std::ostream_iterator<size_t> (os, " "));
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Column"
        << (numBad != 1 ? "s" : "") << " [" << os.str() << "] of X and Y have "
        "norms that differ relatively by more than " << tol << ".");
    }
  }

  // Test Tpetra::MatrixMarket::Reader::readDenseFile(), using no input Map.
  //
  // \tparam MV Tpetra::MultiVector specialization
  //
  // \param inputFilename [in] Name of the Matrix Market format dense
  //   matrix file to read (on Proc 0 of comm only).
  // \param outputFilename [in] Filename to which we can validly write
  //   (on Proc 0 of comm only); used as a sanity test.  If empty
  //   (""), the file is not written.
  // \param comm [in] Communicator, over whose process(es) to
  //   distribute the returned Tpetra::MultiVector.
  // \param tolerant [in] Whether or not to parse the file tolerantly.
  // \param verbose [in] Whether to print verbose output.
  // \param debug [in] Whether to print debugging output.
  template<class MV>
  Teuchos::RCP<MV>
  testReadDenseFile (const std::string& inputFilename,
                     const std::string& outputFilename,
                     const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const bool tolerant,
                     const bool verbose,
                     const bool debug)
  {
    using Teuchos::RCP;
    using std::cerr;
    using std::cout;
    using std::endl;

    using scalar_type = typename MV::scalar_type;
    using local_ordinal_type = typename MV::local_ordinal_type;
    using global_ordinal_type = typename MV::global_ordinal_type;
    using node_type = typename MV::node_type;
    using map_type =
      Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;

    // The reader and writer classes are templated on the
    // Tpetra::CrsMatrix specialization, from which the
    // Tpetra::MultiVector specialization is derived.
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> sparse_matrix_type;

    const int myRank = comm->getRank ();
    if (verbose && myRank == 0) {
      cout << "testReadDenseFile: reading Matrix Market file \""
           << inputFilename << "\":" << endl;
    }

    // Map describing multivector distribution; starts out null and is
    // an output argument of readDenseFile().
    RCP<const map_type> map = Teuchos::null;

    // Read the dense matrix from the given Matrix Market file.
    // This routine acts like an MPI barrier.
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
    RCP<MV> X = reader_type::readDenseFile (inputFilename, comm,
                                            map, tolerant, debug);
    TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::runtime_error,
      "The Tpetra::MultiVector returned from readDenseFile() is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(map.is_null(), std::runtime_error,
      "The Tpetra::Map returned from readDenseFile() is null.");
    if (! X.is_null() && verbose && myRank == 0) {
      cout << "Finished reading file." << endl;
    }

    // Sanity check: write the multivector X to a file, read it back
    // in as a different multivector Y, then make sure it has the same
    // column norms as X.  This assumes that writing a multivector
    // works.  Only do this if the caller supplied a nonempty output
    // file name.
    if (outputFilename != "") {
      typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
      writer_type::writeDenseFile (outputFilename, X);
      // For the sanity test, we also check that reading in the
      // multivector produces a Map which is compatible with X's Map.
      // It does NOT have to be the same as X's Map, and it won't be
      // if X's Map is noncontiguous, has an indexBase other than
      // zero, etc.
      //
      // We're not testing here the functionality of readDenseFile()
      // to reuse an existing Map.
      RCP<const map_type> testMap;
      RCP<MV> Y = reader_type::readDenseFile (inputFilename, comm,
                                              testMap, tolerant, debug);
      TEUCHOS_TEST_FOR_EXCEPTION(! map->isCompatible (*testMap), std::logic_error,
        "Writing the read-in Tpetra::MultiVector X (with no input Map provided) "
        "to a file and reading it back in resulted in a multivector Y with an "
        "incompatible Map.  Please report this bug to the Tpetra developers.");
      assertMultiVectorsEqual<MV> (X, Y);
    }
    return X;
  }

  // Test Tpetra::MatrixMarket::Reader::readDenseFile(), using a custom input Map.
  //
  // \tparam MV Tpetra::MultiVector specialization
  //
  // \param inputFilename [in] Name of the Matrix Market format dense
  //   matrix file to read (on Proc 0 of map's communicator only).
  // \param map [in] Map over which to distribute the returned
  //   Tpetra::MultiVector.  Must be nonnull.  Use testReadDenseFile()
  //   to test the null input map case.
  // \param tolerant [in] Whether or not to parse the file tolerantly.
  // \param verbose [in] Whether to print verbose output.
  // \param debug [in] Whether to print debugging output.
  template<class MV>
  Teuchos::RCP<MV>
  testReadDenseFileWithInputMap (const std::string& inputFilename,
                                 const std::string& outputFilename,
                                 const Teuchos::RCP<const Tpetra::Map<typename MV::local_ordinal_type, typename MV::global_ordinal_type, typename MV::node_type> >& map,
                                 const bool tolerant,
                                 const bool verbose,
                                 const bool debug)
  {
    using Teuchos::RCP;
    using std::cerr;
    using std::cout;
    using std::endl;

    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>
      map_type;

    // The reader and writer classes are templated on the
    // Tpetra::CrsMatrix specialization, from which the
    // Tpetra::MultiVector specialization is derived.
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> sparse_matrix_type;

    TEUCHOS_TEST_FOR_EXCEPTION(map.is_null(), std::invalid_argument,
      "testReadDenseFileWithInputMap requires a nonnull input map.");

    RCP<const Teuchos::Comm<int> > comm = map->getComm();
    const int myRank = comm->getRank ();
    if (verbose && myRank == 0) {
      cout << "testReadDenseFileWithInputMap:" << endl;
    }
    if (verbose && myRank == 0) {
      cout << "Reading Matrix Market file \"" << inputFilename << "\":" << endl;
    }

    // readDenseFile()'s interface reserves the right to modify the
    // input map pointer, as it takes a nonconst reference to an RCP
    // of a const Tpetra::Map.  Tpetra::Map doesn't have a copy or
    // clone routine (as of 07 Feb 2012), so we just copy the pointer
    // and hope for the best.  Ideally, we should copy the map, pass
    // readDenseFile() the copy, and ensure that the resulting output
    // map isSameAs() the input map.
    RCP<const map_type> theMap = map;

    TEUCHOS_TEST_FOR_EXCEPTION(theMap.is_null(), std::logic_error,
      "testReadDenseFileWithInputMap: theMap (a shallow copy of map) should not "
      "be null, but it is.  Please report this bug to the Tpetra developers.");

    // Read the dense matrix from the given Matrix Market file.
    // This routine acts like an MPI barrier.
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
    RCP<MV> X = reader_type::readDenseFile (inputFilename, comm,
                                            theMap, tolerant, debug);
    TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::runtime_error,
      "The Tpetra::MultiVector returned from readDenseFile() is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(theMap.is_null(), std::runtime_error,
      "The Tpetra::Map returned from readDenseFile() is null.");
    if (! X.is_null() && verbose && myRank == 0) {
      cout << "Finished reading file." << endl;
    }

    // Sanity check: write the multivector X to a file, read it back
    // in as a different multivector Y, then make sure it has the same
    // column norms as X.  This assumes that writing a multivector
    // works.  Only do this if the caller supplied a nonempty output
    // file name.
    if (outputFilename != "") {
      typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
      writer_type::writeDenseFile (outputFilename, X);

      // The input Map need not necessarily be compatible with Y's
      // Map.  This is because we are supplying a null input Map to
      // readDenseFile(), so readDenseFile() will construct a
      // contiguous, uniformly distributed Map.  The input Map need
      // not necessarily be contiguous or uniformly distributed.
      RCP<const map_type> testMap;
      RCP<MV> Y = reader_type::readDenseFile (inputFilename, comm,
                                              testMap, tolerant, debug);
      assertMultiVectorsEqual<MV> (X, Y);
    }
    return X;
  }


  // Test Tpetra::MatrixMarket::Reader::writeDenseFile()
  //
  // \tparam MV Tpetra::MultiVector specialization
  //
  // \param outputFilename [in] Name of the Matrix Market format dense
  //   matrix file to write (on Proc 0 in X's communicator only).
  // \param X [in] The nonnull Tpetra::MultiVector instance to write
  //   to the given file.
  // \param echo [in] Whether or not to echo the output to stdout (on
  //   Proc 0 in X's communicator only).
  // \param verbose [in] Whether to print verbose output.
  // \param debug [in] Whether to print debugging output.
  template<class MV>
  void
  testWriteDenseFile (const std::string& outputFilename,
                      const Teuchos::RCP<const MV>& X,
                      const bool echo,
                      const bool verbose,
                      const bool debug)
  {
    using Teuchos::RCP;
    using std::cerr;
    using std::cout;
    using std::endl;

    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>
      map_type;

    // The reader and writer classes are templated on the
    // Tpetra::CrsMatrix specialization, from which the
    // Tpetra::MultiVector specialization is derived.
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> sparse_matrix_type;

    TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::invalid_argument,
      "testWriteDenseFile: The input Tpetra::MultiVector instance X is null.");

    RCP<const map_type> map = X->getMap();
    RCP<const Teuchos::Comm<int> > comm = map->getComm();
    const int myRank = comm->getRank ();

    if (verbose && myRank == 0) {
      cout << "testWriteDenseFile: writing Tpetra::MultiVector instance to "
        "file \"" << outputFilename << "\"" << endl;
    }
    // If specified, write the read-in sparse matrix to a file and/or
    // echo it to stdout.
    typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
    if (outputFilename != "") {
      writer_type::writeDenseFile (outputFilename, X);
    }
    if (verbose && myRank == 0) {
      cout << "Finished writing file." << endl;
    }
    if (echo) {
      if (verbose && myRank == 0) {
        cout << "Echoing output to stdout." << endl;
      }
      writer_type::writeDense (cout, X);
    }

    // Sanity check: read in the multivector again, and make sure that
    // it has the same column norms as the input multivector.  This
    // assumes that reading a multivector works.
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
    const bool tolerant = false; // Our own output shouldn't need tolerant mode.
    RCP<MV> Y = reader_type::readDenseFile (outputFilename, comm,
                                            map, tolerant, debug);
    assertMultiVectorsEqual<MV> (X, Y);
  }

} // namespace (anonymous)


int
main (int argc, char *argv[])
{
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MAX;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::cout;
  using std::endl;

  typedef Tpetra::MultiVector<>::scalar_type scalar_type;
  typedef Tpetra::Map<>::local_ordinal_type local_ordinal_type;
  typedef Tpetra::Map<>::global_ordinal_type global_ordinal_type;
  typedef Tpetra::Map<>::node_type node_type;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  std::string inputFilename;  // Matrix Market file to read
  std::string temporaryFilename; // Matrix Market file to write (if applicable)
  std::string outputFilename; // Matrix Market file to write (if applicable)

  // Number of a specific test to run.  If nonzero, only run that
  // test.  We always run Test #1 since its result is needed for
  // subsequent tests.
  int testToRun = 0;

  // FIXME (mfh 07 Feb 2012) Currently, all tests with a different
  // index base FAIL.  Reading in the multivector appears to be
  // correct, but writing it results in a multivector of all zeros (on
  // _all_ processes).
  bool testDifferentIndexBase = false;
  bool testContiguousInputMap = true;
  bool testNoncontiguousInputMap = false;

  bool testWrite = true; // Test Matrix Market output?
  bool tolerant = false; // Parse the file tolerantly?
  bool echo = false;     // Echo the read-in matrix back?
  bool verbose = false;  // Verbosity of output
  bool debug = false;    // Print debugging info?
  // If true, stop after a single test failure.  Intended for
  // interactive use, so that you can examine a test's output file.
  // Not intended for batch or ctest use.
  bool stopAfterFailure = false;

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("inputFilename", &inputFilename,
                  "Name of the Matrix Market dense matrix file to read.");
  cmdp.setOption ("temporaryFilename", &temporaryFilename,
                  "If --testWrite is true, then use this file as temporary "
                  "storage on (MPI) Proc 0.  Otherwise, this argument is "
                  "ignored.");
  cmdp.setOption ("outputFilename", &outputFilename,
                  "If --testWrite is true, then write the read-in matrix to "
                  "this file in Matrix Market format on (MPI) Proc 0.  "
                  "Otherwise, this argument is ignored.  Note that the output "
                  "file may not be identical to the input file.");
  cmdp.setOption ("testToRun", &testToRun, "Number of a specific test to run.  "
                  "If nonzero, only run that test.  We always run Test #1 since"
                  " its result is needed for subsequent tests.");
  cmdp.setOption ("testDifferentIndexBase", "dontTestDifferentIndexBase",
                  &testDifferentIndexBase, "Whether to test input and output "
                  "for Maps with different index bases.");
  cmdp.setOption ("testContiguousInputMap", "dontTestContiguousInputMap",
                  &testContiguousInputMap,
                  "Whether to test input and output for nonnull contiguous "
                  "input Maps.");
  cmdp.setOption ("testNoncontiguousInputMap", "dontTestNoncontiguousInputMap",
                  &testNoncontiguousInputMap,
                  "Whether to test input and output for nonnull noncontiguous "
                  "input Maps.");
  cmdp.setOption ("testWrite", "noTestWrite", &testWrite,
                  "Whether to test Matrix Market file output.  Ignored if no "
                  "--outputFilename value is given.");
  cmdp.setOption ("tolerant", "strict", &tolerant,
                  "Whether to parse the input Matrix Market file tolerantly.");
  cmdp.setOption ("echo", "noecho", &echo,
                  "Whether to echo the read-in matrix back to stdout on Rank 0 "
                  "in Matrix Market format.  Note that the echoed matrix may "
                  "not be identical to the input file.");
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, "Print debugging information.");
  cmdp.setOption ("stopAfterFailure", "dontStopAfterFailure", &stopAfterFailure,
                  "Whether to stop after a single test failure.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      if (myRank == 0) {
        cout << "End Result: TEST PASSED" << endl;
      }
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::invalid_argument, "Failed to parse command-line arguments.");
  }

  // List of numbers of failed tests.
  std::vector<int> failedTests;
  // Current test number.  Increment before starting each test.  If a
  // test is only run conditionally, increment before evaluating the
  // condition.  This ensures that each test has the same number each
  // time, whether or not a particular test is run.
  int testNum = 0;

  // Run all the tests.  If no input filename was specified, we don't
  // invoke the test and we report a "TEST PASSED" message.
  if (inputFilename != "") {
    // Convenient abbreviations
    typedef scalar_type ST;
    typedef local_ordinal_type LO;
    typedef global_ordinal_type GO;
    typedef node_type NT;
    typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;
    typedef Tpetra::Map<LO, GO, NT> MT;

    // If not testing writes, don't do the sanity check that tests
    // input by comparing against output.
    std::string outFilename = testWrite ? outputFilename : "";
    std::string tmpFilename = testWrite ? temporaryFilename : "";

    // Test 1: null input Map.
    ++testNum;
    if (verbose && myRank == 0) {
      cout << "Test " << testNum << ": Null Map on input to readDenseFile()" << endl;
    }
    RCP<MV> X;
    try {
      X = testReadDenseFile<MV> (inputFilename, tmpFilename, comm,
                                 tolerant, verbose, debug);
      if (outFilename != "") {
        testWriteDenseFile<MV> (outFilename, X, echo, verbose, debug);
      }
    } catch (std::exception& e) {
      failedTests.push_back (testNum);
      // If Test 1 fails, the other tests shouldn't even run, since
      // they depend on the result of Test 1 (the multivector X).
      throw e;
    }

    // Test 2: nonnull contiguous input Map with the same index base
    // as X's Map.  This Map may or may not necessarily be the same as
    // (in the sense of isSameAs()) or even compatible with X's Map.
    ++testNum;
    if ((testToRun == 0 && testContiguousInputMap) ||
        (testToRun != 0 && testToRun == testNum)) {
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull contiguous Map (same index "
          "base) on input to readDenseFile()" << endl;
      }
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();
      const GO indexBase = X->getMap()->getIndexBase();
      // Create the Map.
      RCP<const MT> map =
        rcp (new Tpetra::Map<LO, GO, NT> (globalNumRows, indexBase, comm,
                                          Tpetra::GloballyDistributed));
      try {
        RCP<MV> X2 =
          testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                             map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X2, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    }

    // Test 3: nonnull contiguous input Map, with a different index
    // base than X's Map.  In this case, the index base is X's Map's
    // index base plus a small number (3).  For sufficiently long
    // vectors, this tests the case where the GID sets overlap.
    ++testNum;
    if ((testToRun == 0 && testContiguousInputMap && testDifferentIndexBase) ||
        (testToRun != 0 && testToRun == testNum)) {
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull contiguous Map (different "
          "index base) on input to readDenseFile()" << endl;
      }
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();
      const GO indexBase = X->getMap()->getIndexBase() + as<GO> (3);

      // Make sure that the index base is the same on all processes.
      // It definitely should be, since the Map's getMaxAllGlobalIndex()
      // method should return the same value on all processes.
      GO minIndexBase = indexBase;
      reduceAll (*comm, REDUCE_MIN, indexBase, ptr (&minIndexBase));
      GO maxIndexBase = indexBase;
      reduceAll (*comm, REDUCE_MAX, indexBase, ptr (&maxIndexBase));
      TEUCHOS_TEST_FOR_EXCEPTION(minIndexBase != maxIndexBase || minIndexBase != indexBase,
        std::logic_error, "Index base values do not match on all processes.  "
        "Min value is " << minIndexBase << " and max value is " << maxIndexBase
        << ".");

      // Create the Map.
      RCP<const MT> map =
        rcp (new Tpetra::Map<LO, GO, NT> (globalNumRows, indexBase, comm,
                                          Tpetra::GloballyDistributed));
      try {
        RCP<MV> X3 =
          testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                             map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X3, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    }

    // Test 4: nonnull contiguous input Map, with a different index
    // base than X's Map.  In this case, the new index base is chosen
    // so that the new GID set does not overlap with X's Map's GID
    // set.
    ++testNum;
    if ((testToRun == 0 && testContiguousInputMap && testDifferentIndexBase) ||
        (testToRun != 0 && testToRun == testNum)) {
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull contiguous Map (different "
          "index base) on input to readDenseFile()" << endl;
      }
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();
      // Choose the Map's index base so that the global ordinal sets
      // of X->getMap() and map don't overlap.  This will ensure that
      // we test something nontrivial.
      const GO indexBase = X->getMap()->getMaxAllGlobalIndex() + 1;

      // Make sure that the index base is the same on all processes.
      // It definitely should be, since the Map's getMaxAllGlobalIndex()
      // method should return the same value on all processes.
      GO minIndexBase = indexBase;
      reduceAll (*comm, REDUCE_MIN, indexBase, ptr (&minIndexBase));
      GO maxIndexBase = indexBase;
      reduceAll (*comm, REDUCE_MAX, indexBase, ptr (&maxIndexBase));
      TEUCHOS_TEST_FOR_EXCEPTION(minIndexBase != maxIndexBase || minIndexBase != indexBase,
        std::logic_error, "Index base values do not match on all processes.  "
        "Min value is " << minIndexBase << " and max value is " << maxIndexBase
        << ".");

      // Create the Map.
      RCP<const MT> map =
        rcp (new Tpetra::Map<LO, GO, NT> (globalNumRows, indexBase, comm,
                                          Tpetra::GloballyDistributed));
      try {
        RCP<MV> X3 =
          testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                             map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X3, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    }

    ++testNum;
    if ((testToRun == 0 && testNoncontiguousInputMap) ||
        (testToRun != 0 && testToRun == testNum)) {
      // Test 5: nonnull input Map with the same index base as X's
      // Map, and a "noncontiguous" distribution (in the sense that
      // the Map is constructed using the constructor that takes an
      // arbitrary list of GIDs; that doesn't necessarily mean that
      // the GIDs themselves are noncontiguous).
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull noncontiguous Map (same index "
          "base) on input to readDenseFile()" << endl;
      }
      const GO indexBase = X->getMap()->getIndexBase();
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();

      // Compute number of GIDs owned by each process.  We're
      // replicating Tpetra functionality here because we want to
      // trick Tpetra into thinking we have a noncontiguous
      // distribution.  This is the most general case and the most
      // likely to uncover bugs.
      const size_t quotient = globalNumRows / numProcs;
      const size_t remainder = globalNumRows - quotient * numProcs;
      const size_t localNumRows = (as<size_t> (myRank) < remainder) ?
        (quotient + 1) : quotient;

      // Build the list of GIDs owned by this process.
      Array<GO> elementList (localNumRows);
      GO myStartGID;
      if (as<size_t> (myRank) < remainder) {
        myStartGID = indexBase + as<GO> (myRank) * as<GO> (quotient + 1);
      }
      else {
        // This branch does _not_ assume that GO is a signed type.
        myStartGID = indexBase + as<GO> (remainder) * as<GO> (quotient + 1) +
          (as<GO> (myRank) - as<GO> (remainder)) * as<GO> (quotient);
      }
      for (GO i = 0; i < as<GO> (localNumRows); ++i) {
        elementList[i] = myStartGID + i;
      }

      if (debug) {
        for (int p = 0; p < numProcs; ++p) {
          if (p == myRank) {
            if (elementList.size() > 0) {
              const GO minGID = *std::min_element (elementList.begin(), elementList.end());
              const GO maxGID = *std::max_element (elementList.begin(), elementList.end());
              cerr << "On Proc " << p << ": min,max GID = " << minGID << "," << maxGID << endl;
            }
            else {
              cerr << "On Proc " << p << ": elementList is empty" << endl;
            }
            cerr << std::flush;
          }
          comm->barrier ();
          comm->barrier ();
          comm->barrier ();
        }
      }

      // Create the Map.
      using Tpetra::createNonContigMapWithNode;
      RCP<const MT> map =
        createNonContigMapWithNode<LO, GO, NT> (elementList(), comm);
      try {
        RCP<MV> X4 = testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                                        map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X4, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    } // if test noncontiguous input Map

    ++testNum;
    if ((testToRun == 0 && testNoncontiguousInputMap && testDifferentIndexBase) ||
        (testToRun != 0 && testToRun == testNum)) {
      // Test 6: nonnull input Map with a different index base than
      // X's Map, and a "noncontiguous" distribution (in the sense
      // that the Map is constructed using the constructor that takes
      // an arbitrary list of GIDs; that doesn't necessarily mean that
      // the GIDs themselves are noncontiguous).
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull noncontiguous Map (different "
          "index base) on input to readDenseFile()" << endl;
      }
      // Make sure that the global ordinal sets of X->getMap() and
      // map don't overlap.
      GO indexBase = X->getMap()->getMaxAllGlobalIndex() + 1;
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();

      // Compute number of GIDs owned by each process.  We're
      // replicating Tpetra functionality here because we want to
      // trick Tpetra into thinking we have a noncontiguous
      // distribution.  This is the most general case and the most
      // likely to uncover bugs.
      const size_t quotient = globalNumRows / numProcs;
      const size_t remainder = globalNumRows - quotient * numProcs;
      const size_t localNumRows = (as<size_t> (myRank) < remainder) ?
        (quotient + 1) : quotient;

      // Build the list of GIDs owned by this process.
      Array<GO> elementList (localNumRows);
      GO myStartGID;
      if (as<size_t> (myRank) < remainder) {
        myStartGID = indexBase + as<GO> (myRank) * as<GO> (quotient + 1);
      }
      else {
        // This branch does _not_ assume that GO is a signed type.
        myStartGID = indexBase + as<GO> (remainder) * as<GO> (quotient + 1) +
          (as<GO> (remainder) - as<GO> (myRank)) * as<GO> (quotient);
      }
      for (GO i = 0; i < as<GO> (localNumRows); ++i) {
        elementList[i] = myStartGID + i;
      }

      // Create the Map.
      using Tpetra::createNonContigMapWithNode;
      RCP<const MT> map =
        createNonContigMapWithNode<LO, GO, NT> (elementList(), comm);
      try {
        RCP<MV> X5 = testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                                        map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X5, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    } // if test noncontiguous input Map

    ++testNum;
    if ((testToRun == 0 && testNoncontiguousInputMap) ||
        (testToRun != 0 && testToRun == testNum)) {
      // Test 7: nonnull input Map with the same index base as X's
      // Map, and a "noncontiguous" distribution with GIDs that start
      // at 3.  This lets us easily observe any missing entries after
      // writing X and reading it back in again.
      if (verbose && myRank == 0) {
        cout << "Test " << testNum << ": Nonnull noncontiguous Map (same index "
          "base, GIDs not in 0 .. N-1) on input to readDenseFile()" << endl;
      }
      const Tpetra::global_size_t globalNumRows = X->getMap()->getGlobalNumElements();
      const GO globalStartGID = as<GO> (3);

      // Compute number of GIDs owned by each process.  We're
      // replicating Tpetra functionality here because we want to
      // trick Tpetra into thinking we have a noncontiguous
      // distribution.  This is the most general case and the most
      // likely to uncover bugs.
      const size_t quotient = globalNumRows / numProcs;
      const size_t remainder = globalNumRows - quotient * numProcs;
      const size_t localNumRows = (as<size_t> (myRank) < remainder) ?
        (quotient + 1) : quotient;

      // Build the list of GIDs owned by this process.
      Array<GO> elementList (localNumRows);
      GO myStartGID;
      if (as<size_t> (myRank) < remainder) {
        myStartGID = globalStartGID + as<GO> (myRank) * as<GO> (quotient + 1);
      }
      else {
        // This branch does _not_ assume that GO is a signed type.
        myStartGID = globalStartGID + as<GO> (remainder) * as<GO> (quotient + 1) +
          (as<GO> (myRank) - as<GO> (remainder)) * as<GO> (quotient);
      }
      for (GO i = 0; i < as<GO> (localNumRows); ++i) {
        elementList[i] = myStartGID + i;
      }

      if (debug) {
        for (int p = 0; p < numProcs; ++p) {
          if (p == myRank) {
            if (elementList.size() > 0) {
              const GO minGID = *std::min_element (elementList.begin(), elementList.end());
              const GO maxGID = *std::max_element (elementList.begin(), elementList.end());
              cerr << "On Proc " << p << ": min,max GID = " << minGID << "," << maxGID << endl;
            }
            else {
              cerr << "On Proc " << p << ": elementList is empty" << endl;
            }
            cerr << std::flush;
          }
          comm->barrier ();
          comm->barrier ();
          comm->barrier ();
        }
      }

      // Create the Map.
      using Tpetra::createNonContigMapWithNode;
      RCP<const MT> map =
        createNonContigMapWithNode<LO, GO, NT> (elementList(), comm);
      try {
        RCP<MV> X7 = testReadDenseFileWithInputMap<MV> (inputFilename, tmpFilename,
                                                        map, tolerant, verbose, debug);
        if (outFilename != "") {
          testWriteDenseFile<MV> (outFilename, X7, echo, verbose, debug);
        }
      } catch (std::exception& e) {
        failedTests.push_back (testNum);
        if (myRank == 0) {
          cerr << "Test " << testNum << " failed: " << e.what() << endl;
        }

        if (stopAfterFailure) {
          if (failedTests.size() > 0) {
            if (myRank == 0) {
              cout << "End Result: TEST FAILED" << endl;
            }
            return EXIT_FAILURE;
          }
          else {
            if (myRank == 0) {
              cout << "End Result: TEST PASSED" << endl;
            }
            return EXIT_SUCCESS;
          }
        } // if stop after failure
      }
    } // if test noncontiguous input Map
  }

  if (failedTests.size() > 0) {
    if (myRank == 0) {
      cout << "End Result: TEST FAILED" << endl;
    }
    return EXIT_FAILURE;
  }
  else {
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
    return EXIT_SUCCESS;
  }
}





