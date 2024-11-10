// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Checker_hpp
#define __Teuchos_MatrixMarket_Raw_Checker_hpp

#include "Teuchos_MatrixMarket_Raw_Adder.hpp"
#include "Teuchos_MatrixMarket_SymmetrizingAdder.hpp"
#include "Teuchos_MatrixMarket_CoordDataReader.hpp"


namespace Teuchos {
  namespace MatrixMarket {
    namespace Raw {
      /// \class Checker
      /// \brief Tool for debugging the syntax of a Matrix Market file
      ///   containing a sparse matrix.
      ///
      /// \tparam Scalar The type of entries of the sparse matrix.
      /// \tparam Ordinal The type of indices of the sparse matrix.
      ///
      /// This class is useful for checking the integrity of a Matrix
      /// Market sparse matrix file and printing its contents.  For
      /// reading a sparse matrix from a Matrix Market file into raw
      /// compressed sparse row (CSR) arrays on a single (MPI)
      /// process, use the Reader class.  For reading in a
      /// Tpetra::CrsMatrix, use the Tpetra::MatrixMarket::Reader
      /// class.
      template<class Scalar, class Ordinal>
      class Checker {
      public:
        /// \brief Constructor that takes Boolean parameters.
        ///
        /// \param echo [in] Whether to echo the sparse matrix to
        ///   stdout on MPI Process 0, after reading the sparse matrix.
        /// \param tolerant [in] Whether to parse the Matrix Market
        ///   files tolerantly.
        /// \param debug [in] Whether to print debugging output to
        ///   stderr.  This will happen on all MPI processes, so it
        ///   could be a lot of output.
        Checker (const bool echo, const bool tolerant, const bool debug) :
          echo_ (echo), tolerant_ (tolerant), debug_ (debug)
        {}

        //! Constructor that sets default Boolean parameters.
        Checker () :
          echo_ (false), tolerant_ (false), debug_ (false)
        {}

        /// \brief Constructor that takes a ParameterList of parameters.
        ///
        /// Parameters (all of them have type bool, all default to false):
        /// - "Echo to stdout": Whether to echo the sparse matrix to
        ///   stdout on MPI Process 0, after reading the sparse matrix.
        /// - "Parse tolerantly": Whether to parse Matrix Market files
        ///   tolerantly.
        /// - "Debug mode" Whether to print debugging output to
        ///   stderr.  This will happen on all MPI processes, so it
        ///   could be a lot of output.
        Checker (const RCP<ParameterList>& params) :
          echo_ (false), tolerant_ (false), debug_ (false)
        {
          setParameters (params);
        }

        /// \brief Set parameters from the given ParameterList.
        ///
        /// See constructor documentation for the accepted parameters.
        void
        setParameters (const RCP<ParameterList>& params)
        {
          // Default parameter values.
          bool echo = false;
          bool tolerant = false;
          bool debug = false;

          // Read parameters.
          echo = params->get ("Echo to stdout", echo);
          tolerant = params->get ("Parse tolerantly", tolerant);
          debug = params->get ("Debug mode", debug);

          // No side effects on the class until ParameterList
          // processing is complete.
          echo_ = echo;
          tolerant_ = tolerant;
          debug_ = debug;
        }

        RCP<const ParameterList>
        getValidParameters () const
        {
          // Default parameter values.
          const bool echo = false;
          const bool tolerant = false;
          const bool debug = false;

          // Set default parameters with documentation.
          RCP<ParameterList> params = parameterList ("Matrix Market Checker");
          params->set ("Echo to stdout", echo, "Whether to echo the sparse "
                       "matrix to stdout after reading it");
          params->set ("Parse tolerantly", tolerant, "Whether to tolerate "
                       "syntax errors when parsing the Matrix Market file");
          params->set ("Debug mode", debug, "Whether to print debugging output "
                       "to stderr, on all participating MPI processes");

          return rcp_const_cast<const ParameterList> (params);
        }

        /// \brief Read the sparse matrix from the given file.
        ///
        /// This is a collective operation.  Only MPI Process 0 opens
        /// the file and reads data from it, but all ranks participate
        /// and wait for the final result.
        ///
        /// \note This whole "raw" reader is meant for debugging and
        ///   diagnostics of syntax errors in the Matrix Market file;
        ///   it's not performance-oriented.  That's why we do all the
        ///   broadcasts of and checks for "success".
        bool
        readFile (const Teuchos::Comm<int>& comm,
                  const std::string& filename)
        {
          using std::cerr;
          using std::endl;

          const int myRank = comm.getRank ();
          // Teuchos::broadcast doesn't accept a bool; we use an int
          // instead, with the usual 1->true, 0->false Boolean
          // interpretation.
          int didReadFile = 0;
          RCP<std::ifstream> in; // only valid on Rank 0
          if (myRank == 0) {
            if (debug_) {
              cerr << "Attempting to open file \"" << filename
                   << "\" on Rank 0...";
            }
            in = rcp (new std::ifstream (filename.c_str()));
            if (! *in) {
              didReadFile = 0;
              if (debug_) {
                cerr << "failed." << endl;
              }
            }
            else {
              didReadFile = 1;
              if (debug_) {
                cerr << "succeeded." << endl;
              }
            }
          }
          Teuchos::broadcast (comm, 0, &didReadFile);
          // All MPI processes should throw at the same time, or none.
          TEUCHOS_TEST_FOR_EXCEPTION(! didReadFile, std::runtime_error,
            "Failed to open input file \"" + filename + "\".");
          // Only Rank 0 will try to dereference "in".
          return read (comm, in);
        }

        /// \brief Read the sparse matrix from the given input stream.
        ///
        /// This is a collective operation.  Only MPI Process 0 reads
        /// from the given input stream, but all MPI processes
        /// participate and wait for the final result.
        ///
        /// \note This whole "raw" reader is meant for debugging and
        ///   diagnostics of syntax errors in the Matrix Market file;
        ///   it's not performance-oriented.  That's why we do all the
        ///   broadcasts of and checks for "success".
        bool
        read (const Teuchos::Comm<int>& comm,
              const RCP<std::istream>& in)
        {
          using std::cerr;
          using std::endl;

          const int myRank = comm.getRank ();
          std::pair<bool, std::string> result;
          int msgSize = 0; // Size of error message (if any)
          if (myRank == 0) {
            if (in.is_null()) {
              result.first = false;
              result.second = "Input stream is null on Rank 0";
            }
            else {
              if (debug_) {
                cerr << "About to read from input stream on Rank 0" << endl;
              }
              result = readOnRank0 (*in);
              if (debug_) {
                if (result.first) {
                  cerr << "Successfully read sparse matrix from "
                    "input stream on Rank 0" << endl;
                }
                else {
                  cerr << "Failed to read sparse matrix from input "
                    "stream on Rank 0" << endl;
                }
              }
            }
            if (result.first) {
              msgSize = 0;
            }
            else {
              msgSize = result.second.size();
            }
          }
          int success = result.first ? 1 : 0;
          Teuchos::broadcast (comm, 0, &success);
          if (! success) {
            if (! tolerant_) {
              // Tell all ranks how long the error message is, so
              // they can make space for it in order to receive
              // the broadcast of the error message.
              Teuchos::broadcast (comm, 0, &msgSize);

              if (msgSize > 0) {
                std::string errMsg (msgSize, ' ');
                if (myRank == 0) {
                  std::copy (result.second.begin(), result.second.end(),
                             errMsg.begin());
                }
                Teuchos::broadcast (comm, 0, static_cast<int> (msgSize), &errMsg[0]);
                TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error, errMsg);
              }
              else {
                TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error,
                  "Unknown error when reading Matrix Market sparse matrix file; "
                  "the error is \"unknown\" because the error message has length 0.");
              }
            }
            else if (myRank == 0) {
              using std::cerr;
              using std::endl;
              cerr << "The following error occurred when reading the "
                "sparse matrix: " << result.second << endl;
            }
          }
          return success;
        }

      private:
        //! Whether to echo the sparse matrix to stdout after reading it.
        bool echo_;
        //! Whether to parse the Matrix Market file tolerantly.
        bool tolerant_;
        //! Whether to print debugging output to stderr.
        bool debug_;

        /// \brief Read in the Banner line from the given input stream.
        ///
        /// \warning Only call this method on MPI Process 0.
        ///
        /// \param in [in/out] Input stream from which to read the
        ///   Banner line.  This must be valid on the calling process.
        ///
        /// \param lineNumber [in/out] On input: Current line number
        ///   of the input stream.  On output: if any line(s) were
        ///   successfully read from the input stream, this is
        ///   incremented by the number of line(s) read.  (This
        ///   includes comment lines.)
        ///
        /// \return The Matrix Market file's Banner line (never null).
        RCP<const Teuchos::MatrixMarket::Banner>
        readBanner (std::istream& in, size_t& lineNumber)
        {
          using std::cerr;
          using std::endl;
          std::string line; // The presumed banner line

          // The first line of the Matrix Market file should always be
          // the banner line.  In tolerant mode, we allow comment
          // lines before the banner line.  This complicates detection
          // of comment lines a bit.
          if (tolerant_) {
            // Keep reading lines until we get a noncomment line.
            const bool maybeBannerLine = true;
            size_t numLinesRead = 0;
            bool commentLine = false;
            do {
              // Try to read a line from the input stream.
              const bool readFailed = ! getline (in, line);
              TEUCHOS_TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
                "Failed to get Matrix Market banner line from input, after reading "
                << numLinesRead << "line" << (numLinesRead != 1 ? "s." : "."));
              // We read a line from the input stream.
              ++lineNumber;
              ++numLinesRead;
              size_t start, size; // Output args of checkCommentLine
              commentLine = checkCommentLine (line, start, size, lineNumber,
                                              tolerant_, maybeBannerLine);
            } while (commentLine); // Loop until we find a noncomment line.
          }
          else {
            const bool readFailed = ! getline (in, line);
            TEUCHOS_TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
              "Failed to get Matrix Market banner line from input.  This "
              "probably means that the file is empty (contains zero lines).");
          }

          if (debug_) {
            cerr << "Raw::Checker::readBanner: Here is the presumed banner line:"
                 << endl << line << endl;
          }

          // Assume that the noncomment line we found is the banner line.
          RCP<Banner> banner;
          try {
            banner = rcp (new Banner (line, tolerant_));
          } catch (std::exception& e) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
              "Matrix Market file's banner line contains syntax error(s): "
              << e.what ());
          }
          return rcp_const_cast<const Banner> (banner);
        }

        /// \brief Read the sparse matrix on MPI Rank 0.
        ///
        /// \warning To be called only on MPI Rank 0.
        ///
        /// \param in [in/out] The input stream from which to read.
        ///   This must be valid on the calling process.
        ///
        /// \return First value: Whether reading was successful.  If
        ///   false, the second value is the error message.
        std::pair<bool, std::string>
        readOnRank0 (std::istream& in)
        {
          using std::cerr;
          using std::cout;
          using std::endl;
          typedef ScalarTraits<Scalar> STS;

          // This "Adder" knows how to add sparse matrix entries,
          // given a line of data from the file.  It also stores the
          // entries and can sort them.
          typedef Adder<Scalar, Ordinal> raw_adder_type;
          // SymmetrizingAdder "advices" (yes, I'm using that as a verb)
          // the original Adder, so that additional entries are filled
          // in symmetrically, if the Matrix Market banner line
          // specified a symmetry type other than "general".
          typedef SymmetrizingAdder<raw_adder_type> adder_type;

          // Current line number of the input stream.
          size_t lineNumber = 1;

          // Construct the "Banner" (matrix metadata, including type
          // and symmetry information, but not dimensions).
          std::ostringstream err;
          RCP<const Banner> pBanner;
          try {
            pBanner = readBanner (in, lineNumber);
          }
          catch (std::exception& e) {
            err << "Failed to read Matrix Market file's Banner: " << e.what();
            return std::make_pair (false, err.str());
          }
          //
          // Validate the metadata in the Banner.
          //
          if (pBanner->matrixType () != "coordinate") {
            err << "Matrix Market input file must contain a \"coordinate\"-"
              "format sparse matrix in order to create a sparse matrix object "
              "from it.";
            return std::make_pair (false, err.str ());
          }
          else if (! STS::isComplex && pBanner->dataType () == "complex") {
            err << "The Matrix Market sparse matrix file contains complex-"
              "valued data, but you are try to read the data into a sparse "
              "matrix containing real values (your matrix's Scalar type is "
              "real).";
            return std::make_pair (false, err.str ());
          }
          else if (pBanner->dataType () != "real" &&
                   pBanner->dataType () != "complex") {
            err << "Only real or complex data types (no pattern or integer "
              "matrices) are currently supported.";
            return std::make_pair (false, err.str ());
          }
          if (debug_) {
            cerr << "Banner line:" << endl << *pBanner << endl;
          }

          // The reader will invoke the adder (see below) once for
          // each matrix entry it reads from the input stream.
          typedef CoordDataReader<adder_type, Ordinal, Scalar,
            STS::isComplex> reader_type;
          // We will set the adder below, after calling readDimensions().
          reader_type reader;

          // Read in the dimensions of the sparse matrix: (# rows, #
          // columns, # matrix entries (counting duplicates as
          // separate entries)).  The second element of the pair tells
          // us whether the values were gotten successfully.
          std::pair<Tuple<Ordinal, 3>, bool> dims =
            reader.readDimensions (in, lineNumber, tolerant_);
          if (! dims.second) {
            err << "Error reading Matrix Market sparse matrix "
              "file: failed to read coordinate dimensions.";
            return std::make_pair (false, err.str ());
          }
          // These are "expected" values read from the input stream's
          // metadata.  The actual matrix entries read from the input
          // stream might not conform to their constraints.  We allow
          // such nonconformity only in "tolerant" mode; otherwise, we
          // throw an exception.
          const Ordinal numRows = dims.first[0];
          const Ordinal numCols = dims.first[1];
          const Ordinal numEntries = dims.first[2];
          if (debug_) {
            cerr << "Reported dimensions: " << numRows << " x " << numCols
                 << ", with " << numEntries << " entries (counting possible "
                 << "duplicates)." << endl;
          }

          // The "raw" adder knows about the expected matrix
          // dimensions, but doesn't know about symmetry.
          RCP<raw_adder_type> rawAdder =
            rcp (new raw_adder_type (numRows, numCols, numEntries,
                                     tolerant_, debug_));
          // The symmetrizing adder knows about symmetry.
          RCP<adder_type> adder =
            rcp (new adder_type (rawAdder, pBanner->symmType ()));

          // Give the adder to the reader.
          reader.setAdder (adder);

          // Read the sparse matrix entries.  "results" just tells us if
          // and where there were any bad lines of input.  The actual
          // sparse matrix entries are stored in the (raw) Adder object.
          std::pair<bool, std::vector<size_t> > results =
            reader.read (in, lineNumber, tolerant_, debug_);
          if (debug_) {
            if (results.first) {
              cerr << "Matrix Market file successfully read" << endl;
            }
            else {
              cerr << "Failed to read Matrix Market file" << endl;
            }
          }

          // Report any bad line number(s).
          if (! results.first) {
            if (! tolerant_) {
              err << "The Matrix Market input stream had syntax error(s)."
                "  Here is the error report." << endl;
              reportBadness (err, results);
              err << endl;
              return std::make_pair (false, err.str ());
            }
            else {
              if (debug_) {
                reportBadness (cerr, results);
              }
            }
          }
          // We're done reading in the sparse matrix.  If we're in
          // "echo" mode, print out the matrix entries to stdout.  The
          // entries will have been symmetrized if applicable.
          if (echo_) {
            const bool doMerge = false;
            const bool replace = false;
            rawAdder->print (cout, doMerge, replace);
            cout << endl;
          }
          return std::make_pair (true, err.str());
        }

        //! To be called only on MPI Rank 0.
        void
        reportBadness (std::ostream& out,
                       const std::pair<bool, std::vector<size_t> >& results)
        {
          using std::endl;
          const size_t numErrors = results.second.size();
          const size_t maxNumErrorsToReport = 20;
          out << numErrors << " errors when reading Matrix Market sparse "
            "matrix file." << endl;
          if (numErrors > maxNumErrorsToReport) {
            out << "-- We do not report individual errors when there "
              "are more than " << maxNumErrorsToReport << ".";
          }
          else if (numErrors == 1) {
            out << "Error on line " << results.second[0] << endl;
          }
          else if (numErrors > 1) {
            out << "Errors on lines {";
            for (size_t k = 0; k < numErrors-1; ++k) {
              out << results.second[k] << ", ";
            }
            out << results.second[numErrors-1] << "}" << endl;
          }
        }
      }; // end of class Checker
    } // namespace Raw
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_Raw_Checker_hpp
