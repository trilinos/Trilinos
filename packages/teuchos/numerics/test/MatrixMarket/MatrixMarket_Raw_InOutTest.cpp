// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_MatrixMarket_Raw_Checker.hpp>
#include <Teuchos_MatrixMarket_Raw_Reader.hpp>
#include <Teuchos_MatrixMarket_Raw_Writer.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <algorithm>

using std::endl;

namespace {
  // Sample Matrix Market sparse matrix file.  We include this so we
  // can test without needing to read in a file.  Notice that all the
  // decimal floating-point values in this example can be represented
  // exactly in binary floating point.  This example has correct
  // syntax, so you won't need to use tolerant mode to parse it.
  const char sampleMatrixMarketFile[] =
    "%%MatrixMarket matrix coordinate real general\n"
    "5 5 10\n"
    "5 5 55.0\n"
    "4 4 44.0\n"
    "3 3 33.0\n"
    "2 2 22.0\n"
    "1 1 11.0\n"
    "4 5 45.0\n"
    "3 4 34.0\n"
    "2 3 23.0\n"
    "1 2 12.0\n"
    "1 5 15.0\n";

  // Sample Matrix Market sparse matrix file for testing symmetric
  // storage.  Matrix Market format for symmetric, skew-symemtric,
  // etc. specifies that only the lower triangle should be stored.
  const char symmetricMatrixMarketFile[] =
    "%%MatrixMarket matrix coordinate real symmetric\n"
    "5 5 10\n"
    "5 5 55.0\n"
    "4 4 44.0\n"
    "3 3 33.0\n"
    "2 2 22.0\n"
    "1 1 11.0\n"
    "5 4 54.0\n"
    "4 3 43.0\n"
    "3 2 32.0\n"
    "2 1 21.0\n"
    "5 1 51.0\n";

} // namespace (anonymous)

// Benchmark driver
  int
main (int argc, char *argv[])
{
  using Teuchos::MatrixMarket::Raw::Checker;
  using Teuchos::MatrixMarket::Raw::Reader;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::SerialComm;
  using std::cout;
  using std::cerr;
  typedef double scalar_type;
  typedef int ordinal_type;

  bool success = false;
  // Verbosity of output
  bool verbose = true;
  try {
    // Name of the Matrix Market sparse matrix file to read.  If empty,
    // use the Matrix Market example embedded as a string in this file.
    std::string filename;
    // If true, just check the sparse matrix file.  Otherwise,
    // do a full conversion to CSR (compressed sparse row) format.
    bool checkOnly = false;
    // Whether to echo the sparse matrix to stdout after reading it
    // successfully.
    bool echo = false;
    // Whether to parse the Matrix Market file tolerantly.
    bool tolerant = false;
    // Whether to print debugging-level output
    bool debug = false;

    CommandLineProcessor cmdp (false, true);
    cmdp.setOption ("filename", &filename,
        "Name of the Matrix Market sparse matrix file to read.");
    cmdp.setOption ("checkOnly", "fullTest", &checkOnly,
        "If true, just check the syntax of the input file.  "
        "Otherwise, do a full test.");
    cmdp.setOption ("echo", "noecho", &echo,
        "Whether to echo the sparse matrix contents to stdout "
        "after reading it successfully.");
    cmdp.setOption ("tolerant", "strict", &tolerant,
        "Whether to tolerate syntax errors in the Matrix Market file.");
    cmdp.setOption ("verbose", "quiet", &verbose,
        "Print status output to stdout.");
    cmdp.setOption ("debug", "nodebug", &debug,
        "Print possibly copious debugging output to stderr.");
    // Parse the command-line arguments.
    {
      const CommandLineProcessor::EParseCommandLineReturn parseResult =
        cmdp.parse (argc,argv);
      // If the caller asks us to print the documentation, or does not
      // explicitly say to run the benchmark, we let this "test" pass
      // trivially.
      if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
        std::cout << "End Result: TEST PASSED" << endl;
        return EXIT_SUCCESS;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
          parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
          std::invalid_argument, "Failed to parse command-line arguments.");
    }

    // Test reading in the sparse matrix.  If no filename or an empty
    // filename is specified, the test passes trivially.
    success = true;
    {
      // The following tests check reading in different banners. A bug was found
      // in the banner reader wherein banners with multiple consecutive spaces
      // were not read correctly. These tests assure that banners with or
      // without multiple consecutive spaces/tabs are read correctly.
      bool xs;
      if (verbose) {
        cout << "Checking MatrixMarket banner parsing\n";
      }
      {
        // Well formatted banner, passes trivially
        if (verbose) cout << "Reading first banner\n";
        typedef Checker<scalar_type, ordinal_type> checker_type;
        checker_type checker (echo, false, false);
        RCP<const Comm<int> > comm = rcp (new SerialComm<int>);
        std::stringstream in;
        in.str("%%MatrixMarket matrix coordinate real symmetric\n0 0 0\n");
        RCP<std::istream> inStream = rcpFromRef(in);
        xs = checker.read (*comm, inStream);
        if (verbose) {
          cout << "Banner read " << (!xs ? "un" : "") << "successfully\n";
        }
        success = success && xs;
      }

      {
        // Banner with multiple adjacent/consecutive spaces/tabs
        if (verbose) cout << "Reading second banner\n";
        typedef Checker<scalar_type, ordinal_type> checker_type;
        checker_type checker (echo, false, false);
        RCP<const Comm<int> > comm = rcp (new SerialComm<int>);
        std::stringstream in;
        in.str("%%MatrixMarket\tmatrix\t\tcoordinate   real  symmetric\n0 0 0\n");
        RCP<std::istream> inStream = rcpFromRef(in);
        xs = checker.read (*comm, inStream);
        if (verbose) {
          cout << "Banner read " << (!xs ? "un" : "") << "successfully\n";
        }
        success = success && xs;
      }

      {
        // Bad value in banner.  Should throw std::runtime_error
        if (verbose) cout << "Reading third banner\n";
        typedef Checker<scalar_type, ordinal_type> checker_type;
        checker_type checker (echo, false, false);
        RCP<const Comm<int> > comm = rcp (new SerialComm<int>);
        std::stringstream in;
        try {
          in.str("%%MatrixMarket matrix coordinate real xyz\n0 0 0\n");
          RCP<std::istream> inStream = rcpFromRef(in);
          checker.read (*comm, inStream);
          // The call to read *should* raise an error and the following line
          // should not be encountered
          xs = false;
        }
        catch (const std::runtime_error& e) {
          // The error message will include that "xyz" is a bad value. Check
          // that the string "xyz" is in the error mesage.
          std::string es(e.what());
          xs = es.find("xyz") != std::string::npos;
        }
        if (verbose) {
          cout << "Banner read " << (!xs ? "un" : "") << "successfully\n";
        }
        success = success && xs;
      }

      if (verbose) {
        if (success) {
          cout << "Banners parsed successfully\n";
        } else {
          cout << "Banners not parsed successfully\n";
        }
      }
    }

    if (checkOnly) {
      typedef Checker<scalar_type, ordinal_type> checker_type;
      checker_type checker (echo, tolerant, debug);

      RCP<const Comm<int> > comm = rcp (new SerialComm<int>);
      if (filename != "") {
        if (verbose) {
          cout << "Checking syntax of the Matrix Market file \"" << filename
            << "\"" << endl;
        }
        success = success && checker.readFile (*comm, filename);
        if (verbose) {
          if (success) {
            cout << "The given file is a valid Matrix Market file." << endl;
          }
          else {
            cout << "The given file has syntax errors." << endl;
          }
        }
      }
      else {
        if (verbose) {
          cout << "Checking syntax of the first built-in Matrix Market example" << endl
            << std::flush;// for debug output next
        }
        if (debug) {
          cerr << "First built-in Matrix Market example: " << endl
            << sampleMatrixMarketFile << endl;
        }
        std::istringstream in (sampleMatrixMarketFile);
        RCP<std::istream> inStream = rcpFromRef (in);
        success = success && checker.read (*comm, inStream);
        if (verbose) {
          if (success) {
            cout << "The example has valid Matrix Market syntax." << endl;
          }
          else {
            cout << "The example has syntax errors." << endl;
          }
        }
      }
    }
    else {
      typedef Reader<scalar_type, ordinal_type> reader_type;
      reader_type reader (tolerant, debug);
      ArrayRCP<ordinal_type> ptr, ind;
      ArrayRCP<scalar_type> val;
      ordinal_type numRows, numCols;
      //
      // Read the Matrix Market data, either from a file or from a
      // built-in string.
      //
      if (filename != "") {
        if (verbose) {
          cout << "Reading the Matrix Market file \"" << filename << "\"" << endl;
        }
        success = success && reader.readFile (ptr, ind, val,
            numRows, numCols, filename);
      }
      else {
        if (verbose) {
          cout << "Reading the first built-in Matrix Market example" << endl;
        }
        if (debug) {
          cerr << "First built-in Matrix Market example:" << endl
            << sampleMatrixMarketFile << endl;
        }
        std::istringstream inStr (sampleMatrixMarketFile);
        success = success && reader.read (ptr, ind, val, numRows, numCols, inStr);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error, "Matrix Market "
          "reader failed to read the given file or input stream.");
      if (success && verbose) {
        cout << "Returned from reading the Matrix Market data" << endl
          << std::flush; // for following debug output
      }
      if (debug) {
        cerr << "CSR output info:" << endl
          << "  ptr.size() = " << ptr.size()
          << ", ind.size() = " << ind.size()
          << ", val.size() = " << val.size()
          << ", numRows = " << numRows
          << ", numCols = " << numCols
          << endl;
      }

      // Here's the fun part.  Output the CSR data to an output stream.
      // Then read in the output stream.  The resulting matrix should be
      // exactly the same (unless the original file had elements at the
      // same location that were added together with rounding error).
      // This is a test for both Writer and Reader.
      std::ostringstream outStr;
      if (success && verbose) {
        cout << "Printing the CSR arrays to a Matrix Market output stream"
          << endl << std::flush;
      }
      Teuchos::MatrixMarket::Raw::Writer<scalar_type, ordinal_type> writer;
      writer.write (outStr, ptr (), ind (), val (), numRows, numCols);

      if (debug && echo) {
        cerr << "CSR data:" << endl
          << "- ptr = [";
        for (ordinal_type i = 0; i < ptr.size(); ++i) {
          cerr << ptr[i];
          if (i+1 != ptr.size()) { // don't subtract from zero if unsigned
            cerr << ", ";
          }
        }
        cerr << "]" << endl
          << "- ind = [";
        for (ordinal_type i = 0; i < ind.size(); ++i) {
          cerr << ind[i];
          if (i+1 != ind.size()) { // don't subtract from zero if unsigned
            cerr << ", ";
          }
        }
        cerr << "]" << endl
          << "- val = [";
        for (ordinal_type i = 0; i < val.size(); ++i) {
          cerr << val[i];
          if (i+1 != val.size()) { // don't subtract from zero if unsigned
            cerr << ", ";
          }
        }
        cerr << "]" << endl;

        cerr << "CSR data, converted back to Matrix Market format" << endl;
        writer.write (cerr, ptr (), ind (), val (), numRows, numCols);
        cerr << endl;
      }

      ArrayRCP<ordinal_type> newptr, newind;
      ArrayRCP<scalar_type> newval;
      ordinal_type newNumRows, newNumCols;
      if (success && verbose) {
        cout << "Reading the Matrix Market output back into CSR arrays" << endl;
      }
      {
        std::istringstream inStr (outStr.str ());
        success = success && reader.read (newptr, newind, newval,
            newNumRows, newNumCols, inStr);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(! success, std::logic_error, "Matrix Market "
          "reader failed to read the output back into CSR arrays.");
      if (success && verbose) {
        cout << "Successfully read the Matrix Market output back into CSR arrays"
          << endl << std::flush;
      }
      if (debug) {
        cerr << "CSR output info:" << endl
          << "  newptr.size() = " << newptr.size()
          << ", newind.size() = " << newind.size()
          << ", newval.size() = " << newval.size()
          << ", newNumRows = " << newNumRows
          << ", newNumCols = " << newNumCols
          << endl;
      }

      // The old arrays should equal the new arrays.
      TEUCHOS_TEST_FOR_EXCEPTION(ptr.size () != newptr.size (), std::logic_error,
          "New ptr array has a different length than old ptr array");
      TEUCHOS_TEST_FOR_EXCEPTION(ind.size () != newind.size (), std::logic_error,
          "New ind array has a different length than old ind array");
      TEUCHOS_TEST_FOR_EXCEPTION(val.size () != newval.size (), std::logic_error,
          "New val array has a different length than old val array");
      TEUCHOS_TEST_FOR_EXCEPTION(newNumRows != numRows || newNumCols != numCols,
          std::logic_error, "New dimensions differ from old dimensions");
      TEUCHOS_TEST_FOR_EXCEPTION(ptr.size () != numRows+1, std::logic_error,
          "ptr.size() != numRows+1");
      TEUCHOS_TEST_FOR_EXCEPTION(newptr.size () != newNumRows+1, std::logic_error,
          "newptr.size() != newNumRows+1");

      for (ordinal_type rowIndex = 0; rowIndex < numRows; ++rowIndex) {
        TEUCHOS_TEST_FOR_EXCEPTION(ptr[rowIndex] != newptr[rowIndex],
            std::logic_error, "At row index " << rowIndex << ", ptr[rowIndex] = "
            << ptr[rowIndex] << " != newptr[rowIndex] = " << newptr[rowIndex]
            << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(ptr[rowIndex+1] != newptr[rowIndex+1],
            std::logic_error, "At row index " << rowIndex << ", ptr[rowIndex+1] = "
            << ptr[rowIndex+1] << " != newptr[rowIndex+1] = " << newptr[rowIndex+1]
            << ".");
        for (ordinal_type k = ptr[rowIndex]; k < ptr[rowIndex+1]; ++k) {
          TEUCHOS_TEST_FOR_EXCEPTION(ind[k] != newind[k], std::logic_error,
              "At row index " << rowIndex << ", ind[k=" << k << "] = "
              << ind[k] << " != newind[k] = " << newind[k] << ".");
          // You may want to relax this inequality if the original
          // Matrix Market file had multiple entries at the same
          // location and if adding them together resulted in rounding
          // error.
          TEUCHOS_TEST_FOR_EXCEPTION(val[k] != newval[k], std::logic_error,
              "At row index " << rowIndex << ", val[k=" << k << "] = "
              << val[k] << " != newval[k] = " << newval[k] << ".");
        }
      }

      // Now test reading symmetric data, if no filename was specified.
      if (filename == "") {
        std::istringstream inStr (symmetricMatrixMarketFile);
        success = success && reader.read (ptr, ind, val, numRows, numCols, inStr);
        TEUCHOS_TEST_FOR_EXCEPTION(! success, std::logic_error,
            "Matrix Market reader failed to read the given example string.");
        if (success && verbose) {
          cout << "Returned from reading the Matrix Market data" << endl
            << std::flush; // for following debug output
        }
        if (debug) {
          cerr << "CSR output info:" << endl
            << "  ptr.size() = " << ptr.size()
            << ", ind.size() = " << ind.size()
            << ", val.size() = " << val.size()
            << ", numRows = " << numRows
            << ", numCols = " << numCols
            << endl;
        }

        // This is a bit of a hack, since we know the contents of the
        // example.  Since we "symmetrize" when reading in symmetric
        // data, there should be 15 entries in the resulting matrix.
        const ordinal_type correctNumEntries = 15;
        TEUCHOS_TEST_FOR_EXCEPTION(
            val.size() != correctNumEntries,
            std::logic_error,
            "Incorrect number of entries after symmetrization: There should be "
            << correctNumEntries << ", but there are " << val.size() << " entries "
            "instead.");
      }
    } // end of the file / string Reader tests

    if (success)
      std::cout << "End Result: TEST PASSED" << endl;
    else
      std::cout << "End Result: TEST FAILED" << endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
