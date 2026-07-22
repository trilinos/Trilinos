// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Reader_hpp
#define __Teuchos_MatrixMarket_Raw_Reader_hpp

#include "Teuchos_MatrixMarket_Raw_Adder.hpp"
#include "Teuchos_MatrixMarket_SymmetrizingAdder.hpp"
#include "Teuchos_MatrixMarket_CoordDataReader.hpp"


namespace Teuchos {
  /// \namespace MatrixMarket
  /// \brief Matrix Market file utilities
  /// \author Mark Hoemmen
  ///
  /// The Matrix Market (see their <a
  /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
  /// details) defines a human-readable ASCII text file format
  /// ("Matrix Market format") for interchange of sparse and dense
  /// matrices.  This namespace defines utility classes for input and
  /// output of Matrix Market files or input streams.  Users will
  /// likely find Raw::Reader most useful; it reads from a Matrix
  /// Market file or input stream into raw compressed sparse row (CSR)
  /// arrays.  Other classes and functions will probably be more
  /// useful for Trilinos developers.
  ///
  /// Matrix Market files are designed for easy reading and writing of
  /// test matrices by both humans and computers.  They are not
  /// intended for high-performance or parallel file input and output.
  /// You should use a true parallel file format if you want to do
  /// parallel input and output of sparse or dense matrices.
  namespace MatrixMarket {
    /// \namespace Raw
    /// \brief "Raw" input of sparse matrices from Matrix Market files.
    ///
    /// "Raw" means serial (not MPI or otherwise distributed over
    /// parallel processes), with storage as a collection of matrix
    /// indices and values.  This is useful if you want to read the
    /// sparse matrix on one (MPI) process and store it in a custom
    /// format.
    ///
    /// For reading a sparse matrix from a Matrix Market file into raw
    /// compressed sparse row (CSR) arrays on a single (MPI) process,
    /// use the Reader class.  For reading in a Tpetra::CrsMatrix, use
    /// the Tpetra::MatrixMarket::Reader class.  Nearly everything
    /// else in this namespace is of interest only to Trilinos
    /// developers.
    namespace Raw {
      /// \class Reader
      /// \brief Read a sparse matrix from a Matrix Market file into
      ///   raw CSR (compressed sparse row) storage.
      ///
      /// \tparam Scalar The type of entries of the sparse matrix.
      /// \tparam Ordinal The type of indices of the sparse matrix.
      ///
      /// This class is useful for benchmarking local sparse kernels.
      /// It should only be called by one MPI process at a time and is
      /// not aware of parallel communication.  Use
      /// Tpetra::MatrixMarket::Reader if you want to read a
      /// Tpetra::CrsMatrix from a Matrix Market file.
      template<class Scalar, class Ordinal>
      class Reader {
      public:
        /// \brief Constructor that takes Boolean parameters.
        ///
        /// \param tolerant [in] Whether to parse the Matrix Market
        ///   files tolerantly.
        /// \param debug [in] Whether to print (possibly copious)
        ///   debugging output to stderr.
        Reader (const bool tolerant, const bool debug) :
          tolerant_ (tolerant), debug_ (debug)
        {
          init ();
        }

        //! Constructor that sets default Boolean parameters.
        Reader () :
          tolerant_ (false), debug_ (false)
        {
          init ();
        }

        /// \brief Constructor that takes a ParameterList of parameters.
        ///
        /// Parameters (all of them have type bool, all default to false):
        /// - "Parse tolerantly": Whether to parse Matrix Market files
        ///   tolerantly.
        /// - "Debug mode": Whether to print (possibly copious)
        ///   debugging output to stderr.
        Reader (const RCP<ParameterList>& params) :
          tolerant_ (false), debug_ (false)
        {
          setParameters (params);
          init ();
        }

        /// \brief Set parameters from the given ParameterList.
        ///
        /// See constructor documentation for the accepted parameters.
        void
        setParameters (const RCP<ParameterList>& params)
        {
          // Default parameter values.
          bool tolerant = false;
          bool debug = false;

          // Read parameters.
          tolerant = params->get ("Parse tolerantly", tolerant);
          debug = params->get ("Debug mode", debug);

          // No side effects on the class until ParameterList
          // processing is complete.
          tolerant_ = tolerant;
          debug_ = debug;
        }

        //! Get a list of valid default parameters, with documentation.
        RCP<const ParameterList>
        getValidParameters () const
        {
          // Default parameter values.
          const bool tolerant = false;
          const bool debug = false;

          // Set default parameters with documentation.
          RCP<ParameterList> params = parameterList ("Matrix Market Reader");
          params->set ("Parse tolerantly", tolerant, "Whether to tolerate "
                       "syntax errors when parsing the Matrix Market file");
          params->set ("Debug mode", debug, "Whether to print debugging output "
                       "to stderr, on all participating MPI processes");

          return rcp_const_cast<const ParameterList> (params);
        }

        /// \brief Read the sparse matrix from the given file into CSR storage.
        ///
        /// The outputs rowptr, colind, values together form the
        /// common three-arrays representation of compressed sparse
        /// row (CSR) storage.
        ///
        /// \param rowptr [out] Array of numRows+1 offsets, where
        ///   numRows is the number of rows in the sparse matrix.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        /// \param colind [out] Column indices of the matrix.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        /// \param values [out] Values stored in the matrix.
        /// \param numRows [out] Number of rows in the sparse matrix.
        ///   This is redundant, because rowptr.size() == numRows on
        ///   output.
        /// \param numCols [out] Number of columns in the sparse matrix.
        /// \param filename [in] Name of the Matrix Market file from
        ///   which to read the sparse matrix.
        ///
        /// \return If parsing tolerantly: false if the Matrix Market
        ///   file has any syntax errors, else true.  If not parsing
        ///   tolerantly, this should always return true.
        bool
        readFile (ArrayRCP<Ordinal>& rowptr,
                  ArrayRCP<Ordinal>& colind,
                  ArrayRCP<Scalar>& values,
                  Ordinal& numRows,
                  Ordinal& numCols,
                  const std::string& filename)
        {
          std::ifstream in (filename.c_str ());
          TEUCHOS_TEST_FOR_EXCEPTION(! in, std::runtime_error,
            "Failed to open file \"" << filename << "\" for reading.");
          return read (rowptr, colind, values, numRows, numCols, in);
        }

        /// \brief Read the sparse matrix from the given input stream
        ///   into CSR storage.
        ///
        /// The outputs rowptr, colind, values together form the
        /// common three-arrays representation of compressed sparse
        /// row (CSR) storage.
        ///
        /// \param rowptr [out] Array of numRows+1 offsets, where
        ///   numRows is the number of rows in the sparse matrix.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        /// \param colind [out] Column indices of the matrix.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        /// \param values [out] Values stored in the matrix.
        /// \param numRows [out] Number of rows in the sparse matrix.
        ///   This is redundant, because rowptr.size() == numRows on
        ///   output.
        /// \param numCols [out] Number of columns in the sparse matrix.
        /// \param in [in/out] Input stream from which to read the
        ///   sparse matrix.
        ///
        /// \return If parsing tolerantly: false if the Matrix Market
        ///   file has any syntax errors, else true.  If not parsing
        ///   tolerantly, this should always return true.
        bool
        read (ArrayRCP<Ordinal>& rowptr,
              ArrayRCP<Ordinal>& colind,
              ArrayRCP<Scalar>& values,
              Ordinal& numRows,
              Ordinal& numCols,
              std::istream& in)
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
          RCP<const Banner> banner;
          std::ostringstream err;
          try {
            banner = readBanner (in, lineNumber);
          }
          catch (std::exception& e) {
            err << "Failed to read Matrix Market input's Banner: " << e.what();
            if (tolerant_) {
              if (debug_) {
                cerr << err.str() << endl;
              }
              return false;
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err.str());
            }
          }

          //
          // Validate the metadata in the Banner.
          //
          bool ok = true;
          if (banner->matrixType () != "coordinate") {
            err << "Matrix Market input file must contain a \"coordinate\"-"
              "format sparse matrix in order to create a sparse matrix object "
              "from it.";
            ok = false;
          }
          else if (! STS::isComplex && banner->dataType () == "complex") {
            err << "The Matrix Market sparse matrix file contains complex-"
              "valued data, but you are try to read the data into a sparse "
              "matrix containing real values (your matrix's Scalar type is "
              "real).";
            ok = false;
          }
          else if (banner->dataType () != "real" &&
                   banner->dataType () != "complex") {
            err << "Only real or complex data types (no pattern or integer "
              "matrices) are currently supported.";
            ok = false;
          }
          if (! ok) {
            if (tolerant_) {
              if (debug_) {
                cerr << "Matrix Market banner is invalid: " << err.str () << endl;
                return false;
              }
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                "Matrix Market banner is invalid: " << err.str ());
            }
          }
          if (debug_) {
            cerr << "Matrix Market Banner line:" << endl << *banner << endl;
          }

          // The reader will invoke the adder (see below) once for
          // each matrix entry it reads from the input stream.
          typedef CoordDataReader<adder_type, Ordinal, Scalar, STS::isComplex> reader_type;
          // We will set the adder below, after calling readDimensions().
          reader_type reader;

          // Read in the dimensions of the sparse matrix: (# rows, #
          // columns, # matrix entries (counting duplicates as
          // separate entries)).  The second element of the pair tells
          // us whether the values were gotten successfully.
          std::pair<Tuple<Ordinal, 3>, bool> dims =
            reader.readDimensions (in, lineNumber, tolerant_);
          if (! dims.second) {
            err << "Error reading Matrix Market sparse matrix file: failed to "
              "read coordinate dimensions.";
            if (tolerant_) {
              if (debug_) {
                cerr << err.str () << endl;
              }
              return false;
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err.str ());
            }
          }

          // These are "expected" values read from the input stream's
          // metadata.  The actual matrix entries read from the input
          // stream might not conform to their constraints.  We allow
          // such nonconformity only in "tolerant" mode; otherwise, we
          // throw an exception.
          numRows = dims.first[0];
          numCols = dims.first[1];
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
          // The symmetrizing adder knows about symmetry.  It mediates
          // adding entries to the "raw" adder.  We'll use the raw
          // adder to compute the CSR arrays.
          RCP<adder_type> adder =
            rcp (new adder_type (rawAdder, banner->symmType ()));

          // Give the adder to the reader.
          reader.setAdder (adder);

          // Read the sparse matrix entries.  "results" just tells us if
          // and where there were any bad lines of input.  The actual
          // sparse matrix entries are stored in the (raw) Adder object.
          std::pair<bool, std::vector<size_t> > results =
            reader.read (in, lineNumber, tolerant_, debug_);

          // Report any bad line number(s).
          if (! results.first) {
            err << "The Matrix Market input stream had syntax error(s)."
              "  Here is the error report." << endl;
            reportBadness (err, results);
            if (tolerant_) {
              if (debug_) {
                cerr << err.str() << endl;
              }
              return false;
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err.str ());
            }
          }
          // Done reading the sparse matrix; now extract CSR arrays.
          size_t numUnique, numRemoved;
          ArrayRCP<Ordinal> ptr;
          ArrayRCP<Ordinal> ind;
          ArrayRCP<Scalar> val;
          try {
            rawAdder->mergeAndConvertToCSR (numUnique, numRemoved, ptr, ind, val);
          }
          catch (std::exception& e) {
            err << "Failed to convert sparse matrix data to CSR (compressed "
              "sparse row) format.  Reported error: " << e.what ();
            if (tolerant_) {
              if (debug_) {
                cerr << err.str () << endl;
              }
              return false;
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err.str ());
            }
          }
          rowptr = ptr;
          colind = ind;
          values = val;
          return true;
        }

      private:
        //! Whether to parse the Matrix Market file tolerantly.
        bool tolerant_;
        //! Whether to print debugging output to stderr.
        bool debug_;

        /// \brief "Initialize" the Reader.
        ///
        /// Right now, this means print debugging output on creation
        /// of the Reader, if the user selected that option.  We put
        /// this into a function to avoid duplicated debugging output
        /// code in the different constructors.
        void init () {
          using std::cerr;
          using std::endl;

          if (debug_) {
            cerr << "MatrixMarket::Raw::Reader:" << endl
                 << "- Tolerant mode: " << tolerant_ << endl
                 << "- Debug mode: " << debug_ << endl;
          }
        }

        /// \brief Read in the Banner line from the given input stream.
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
        RCP<const Banner>
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
            cerr << "Raw::Reader::readBanner: Here is the presumed banner line:"
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

        /// Report syntax errors in the input stream's sparse matrix data.
	///
	/// \param out [in/out] Output stream to which to report.
	/// \param results [in] Return value of CoordDataReader's read() method.
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
      }; // end of class Reader
    } // namespace Raw
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_Raw_Reader_hpp
