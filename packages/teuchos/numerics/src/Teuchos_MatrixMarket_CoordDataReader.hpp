// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef __Teuchos_MatrixMarket_CoordDataReader_hpp
#define __Teuchos_MatrixMarket_CoordDataReader_hpp

#include "Teuchos_MatrixMarket_generic.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Tuple.hpp"


namespace Teuchos {
  namespace MatrixMarket {
    /// \class CoordDataReaderBase
    /// \brief Common functionality of a coordinate-format sparse
    ///   matrix or graph data reader.
    ///
    /// This class provides common functionality for reading
    /// coordinate-format sparse matrix or graph data from a Matrix
    /// Market file.  In particular, this class does not depend on the
    /// type of matrix entries.  If you are writing a function to read
    /// a sparse matrix, use CoordDataReader.  If you are writing a
    /// function to read a sparse graph (a "pattern matrix" in Matrix
    /// Market terms), use CoordPatternReader.  If you just want to
    /// read a sparse matrix into raw compressed sparse row (CSR)
    /// arrays, use Raw::Reader.  This class is mainly for Teuchos
    /// developers.
    ///
    /// \tparam Callback The type of a callback (a.k.a. closure) that
    ///   knows how to add entries to the sparse graph or matrix.
    /// \param Ordinal The type of indices of the sparse graph or matrix.
    ///
    /// Instances of Callback must implement an operator() that adds a
    /// single entry to the sparse matrix.  If the matrix contains
    /// values (i.e., is not a "pattern" matrix), the Callback's
    /// operator() takes three arguments: the row index, the column
    /// index, and the value to add.  If the matrix is a "pattern"
    /// matrix, then the operator() takes two arguments: the row
    /// index, and the column index.  In either case, the Callback
    /// expects one-based row and column indices, as specified by the
    /// Matrix Market standard.  The Raw::Adder class provides a model
    /// for Callback with the three-argument version of operator().
    ///
    /// The Callback object defines for itself what happens if one
    /// attempts to insert two entries at the same location.  They
    /// must be "merged" in some way.  Usually the right thing to do
    /// is to add the values of the two entries, but you may want to
    /// define different behavior.  Your Callback implementation may
    /// defer merging entries at the same location for performance
    /// reasons.  (Raw::Adder does this.  Its merge() method merges
    /// the entries.)
    ///
    /// Note that this class does not have the notion of the type of
    /// entries of the matrix, only the type of the indices in the
    /// matrix (Ordinal).  The CoordDataReader subclass defines the
    /// type of matrix entries as its Scalar template parameter.  The
    /// CoordPatternReader subclass does not define the type of matrix
    /// entries, because it stores a graph (a "pattern matrix" in
    /// Matrix Market terms).
    template<class Callback, class Ordinal>
    class CoordDataReaderBase {
    protected:
      //! Closure that knows how to add entries to the sparse graph or matrix.
      Teuchos::RCP<Callback> adder_;

    public:
      /// \brief Constructor with "adder" argument.
      ///
      /// This is the favored way to construct an instance of this
      /// type.  Only use the no-argument constructor if you have a
      /// "chicken-and-egg" problem, where in order to create the
      /// Callback instance, you need the graph or matrix dimensions.
      ///
      /// \param adder [in/out] Closure (a.k.a. callback) whose
      ///   operator() adds an entry to the sparse graph or matrix on
      ///   each invocation.
      CoordDataReaderBase (const Teuchos::RCP<Callback>& adder) :
        adder_ (adder) {}

      /// \brief No-argument constructor.
      ///
      /// We offer this option in case the adder's constructor needs
      /// the graph or matrix dimensions, so that it's necessary to
      /// call readDimensions() first before constructing the adder.
      /// You should call setAdder() with a non-null argument before
      /// calling read() or readLine().
      CoordDataReaderBase () : adder_ (null) {}

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReaderBase () {}

      /// \brief Set the Adder object.
      ///
      /// Please don't call this after calling read() or readLine().
      /// The right time to call this is right after calling the
      /// no-argument constructor, if it's not possible to supply an
      /// Adder object before calling readDimensions().
      void setAdder (const Teuchos::RCP<Callback>& adder) {
        adder_ = adder;
      }

    protected:
      /// \brief Read in the data from a single line of the input stream.
      ///
      /// \param theLine [in] The line read in from the input stream.
      /// \param adder [in/out] The callback to invoke for adding an
      ///   entry to the sparse matrix.
      /// \param lineNumber [in] Current line number of the file.
      ///   We use this for generating informative exception messages.
      /// \param tolerant [in] Whether to parse tolerantly.
      ///
      /// \return In tolerant parsing mode (tolerant==true), then this
      ///   method returns true if parsing the current line succeeded,
      ///   else false.  Otherwise, this method throws an exception
      ///   (and does not invoke the adder) if parsing the current
      ///   line did not succeed.
      ///
      /// Subclasses must implement this method in order to read one
      /// entry of the sparse graph or matrix.  Implementations should
      /// use the callback (\c adder_) to add the entry.
      ///
      /// \note To implementers: We defer implementation of this
      ///   method to subclasses, because the callback for a graph
      ///   will take different arguments than the callback for a
      ///   matrix.  Abstracting around that using templates isn't
      ///   worth the trouble.  (Remember you're reading from a file
      ///   and parsing strings.  Virtual method call overhead isn't
      ///   significant by comparison.)
      virtual bool
      readLine (const std::string& theLine,
                const size_t lineNumber,
                const bool tolerant) = 0;

    public:

      /// \brief Read in all the data from the given input stream.
      ///
      /// \param in [in/out] The input stream from which to read
      ///
      /// \param startingLineNumber [in] The line number of the file
      ///   from which we begin reading.  (This is used for
      ///   informative error output, if an error is detected in the
      ///   file.)
      ///
      /// \param tolerant [in] If true, parse tolerantly.  The
      ///   resulting read-in data may be incorrect, but the parser
      ///   won't throw an exception if it can just let bad data
      ///   through and continue.
      ///
      /// \param debug [in] If true, print verbose debugging output.
      ///
      /// \return If tolerant==false, the returned pair is always true
      ///   and the empty vector.  If tolerant==true, the first entry
      ///   of the pair is whether all lines of data were read
      ///   successfully, and the second entry is the list of line
      ///   numbers containing bad data.  The first entry of the pair
      ///   is false if and only if the vector has nonzero size.
      ///
      /// \note This method is virtual in case derived classes want to
      ///   override the default behavior.  The specific example we
      ///   have in mind is a "pattern matrix" (i.e., a sparse graph),
      ///   which we would represent with Scalar=void.
      virtual std::pair<bool, std::vector<size_t> >
      read (std::istream& in,
            const size_t startingLineNumber,
            const bool tolerant,
            const bool debug = false)
      {
        (void) debug; // silence unused input argument warning
        TEUCHOS_TEST_FOR_EXCEPTION(! in, std::invalid_argument,
          "Input stream is invalid.");

        std::string line;
        size_t lineNumber = startingLineNumber;
        bool allSucceeded = true;
        std::vector<size_t> badLineNumbers;
        size_t validDataLines = 0;
        while (getline (in, line)) {
          size_t start, size;
          if (checkCommentLine (line, start, size, lineNumber, tolerant)) {
            ++lineNumber;
            continue; // it's a comment line
          }
          const std::string theLine = line.substr (start, size);

          const bool localSuccess = readLine (theLine, lineNumber, tolerant);
          ++lineNumber;
          allSucceeded = allSucceeded && localSuccess;
          if (! localSuccess) {
            badLineNumbers.push_back (lineNumber);
          }
          else {
            ++validDataLines;
          }
        }
        return std::make_pair (allSucceeded, badLineNumbers);
      }

      /// \brief Read (numRows, numCols, numNonzeros).
      ///
      /// Read one line from the given input stream, and parse it into
      /// the number of rows, the number of columns, and the number of
      /// nonzeros in the sparse matrix.  We assume that those data
      /// are in whitespace-delimited format and are read from a
      /// Matrix Market - format file.
      ///
      /// \param in [in/out] The input stream from which to attempt to
      ///   read one line.
      ///
      /// \param lineNumber The starting line number from which we
      ///   begin reading from the input stream.  Used for diagnostic
      ///   error output.
      ///
      /// \param tolerant [in] Whether to read "tolerantly" (setting
      ///   defaults and returning whether we were successful) or
      ///   "intolerantly" (throwing an exception on any deviation
      ///   from the expected format).
      ///
      /// \return ((numRows, numCols, numNonzeros), success).  In
      ///   tolerant mode, success may be false, meaning that the
      ///   read-in triple may not be valid.
      std::pair<Teuchos::Tuple<Ordinal, 3>, bool>
      readDimensions (std::istream& in,
                      size_t& lineNumber,
                      const bool tolerant = false)
      {
        Teuchos::Tuple<Ordinal, 3> dims;
        // Fill in (numRows, numCols, numNonzeros) with reasonable
        // defaults.  If we don't succeed in reading all the data
        // from the current line of the input stream, the remaining
        // values not read will be these default values.
        dims[0] = 0;
        dims[1] = 0;
        dims[2] = 0;

        // Keep reading lines from the input stream until we find a
        // non-comment line, or until we run out of lines.  The latter
        // is an error, since every "coordinate" format Matrix Market
        // file must have a dimensions line after the banner (even if
        // the matrix has zero rows or columns, or zero entries).
        std::string line;
        bool commentLine = true;
        while (commentLine) {
          // Is it even valid to read from the input stream?
          if (in.eof() || in.fail()) {
            if (tolerant) {
              return std::make_pair (dims, false);
            }
            else {
              std::ostringstream os;
              os << "Unable to get coordinate dimensions line (at all) "
                "from (line " << lineNumber << ") of input stream; the "
                "input stream claims that it is at \"end-of-file\" or has "
                "an otherwise \"fail\"ed state.";
              throw std::invalid_argument(os.str());
            }
          }
          // Try to get the next line from the input stream.
          if (getline(in, line)) {
            lineNumber++; // We did actually read a line
          }
          else {
            if (tolerant) {
              return std::make_pair (dims, false);
            }
            else {
              std::ostringstream os;
              os << "Failed to read coordinate dimensions line (at all) "
                "from (line " << lineNumber << " from input stream.  The "
                "line should contain the coordinate matrix dimensions in "
                 << " the form \"<numRows> <numCols> <numNonzeros>\".";
              throw std::invalid_argument (os.str());
            }
          }
          // Is the current line a comment line?  Ignore start and
          // size; they are only useful for reading the actual matrix
          // entries.  (We could use them here as an optimization, but
          // we've chosen not to.)
          size_t start = 0, size = 0;
          commentLine = checkCommentLine (line, start, size,
                                          lineNumber, tolerant);
        }
        //
        // Read in <numRows> <numCols> <numNonzeros> from input line
        //
        std::istringstream istr (line);
        // Does line contain anything at all?  Can we safely read from
        // the input stream wrapping the line?
        if (istr.eof() || istr.fail()) {
          if (tolerant) {
            return std::make_pair (dims, false);
          }
          std::ostringstream os;
          os << "Unable to read any data from line " << lineNumber
             << " of input; the line should contain the coordinate matrix "
             << "dimensions \"<numRows> <numCols> <numNonzeros>\".";
          throw std::invalid_argument(os.str());
        }
        // Read in <numRows>
        {
          Ordinal theNumRows = 0;
          istr >> theNumRows;
          if (istr.fail()) {
            if (tolerant) {
              return std::make_pair (dims, false);
            }
            std::ostringstream os;
            os << "Failed to get number of rows from line " << lineNumber
               << " of input; the line should contain the coordinate matrix "
               << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
            throw std::invalid_argument(os.str());
          }
          else { // Capture the validly read result before checking for eof.
            dims[0] = theNumRows;
          }
        }
        // There should be two more things to read.
        if (istr.eof()) {
          if (tolerant) {
            return std::make_pair (dims, false);
          }
          std::ostringstream os;
          os << "No more data after number of rows on line " << lineNumber
             << " of input; the line should contain the coordinate matrix "
             << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
          throw std::invalid_argument(os.str());
        }
        // Read in <numCols>
        {
          Ordinal theNumCols = 0;
          istr >> theNumCols;
          if (istr.fail()) {
            if (tolerant) {
              return std::make_pair (dims, false);
            }
            std::ostringstream os;
            os << "Failed to get number of columns from line " << lineNumber
               << " of input; the line should contain the coordinate matrix "
               << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
            throw std::invalid_argument(os.str());
          }
          else { // Capture the validly read result before checking for eof.
            dims[1] = theNumCols;
          }
        }
        // There should be one more thing to read.
        if (istr.eof()) {
          if (tolerant) {
            return std::make_pair (dims, false);
          }
          std::ostringstream os;
          os << "No more data after number of columns on line " << lineNumber
             << " of input; the line should contain the coordinate matrix "
             << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
          throw std::invalid_argument(os.str());
        }
        // Read in <numNonzeros>
        {
          Ordinal theNumNonzeros = 0;
          istr >> theNumNonzeros;
          if (istr.fail()) {
            if (tolerant) {
              return std::make_pair (dims, false);
            }
            std::ostringstream os;
            os << "Failed to get number of (structural) nonzeros from line "
               << lineNumber
               << " of input; the line should contain the coordinate matrix "
               << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
            throw std::invalid_argument(os.str());
          }
          else { // Capture the validly read result
            dims[2] = theNumNonzeros;
          }
        }
        // It would be nice to validate the read-in data further.  The
        // only thing we can do now is test if it's negative.  However,
        // we don't know syntactically whether Ordinal is a signed or
        // unsigned type, so we shouldn't even test for negativity.
        return std::make_pair (dims, true);
      }
    };

    /// \class CoordDataReader
    /// \brief Coordinate-format sparse matrix data reader.
    ///
    /// Use this class to read in sparse matrix data and add it to the
    /// sparse matrix using your Callback implementation (see the
    /// documentation of CoordDataReaderBase).
    ///
    /// \tparam Callback Same as in CoordDataReaderBase.  Its
    ///   operator() takes three arguments: the row index, the column
    ///   index, and the matrix value.  The indices are all one-based.
    ///
    /// \tparam Ordinal Same as in CoordDataReaderBase.  The type of
    ///   indices in the sparse matrix.
    ///
    /// \tparam Scalar The type of entries of the sparse matrix.  For
    ///   a real Scalar type (isComplex=false), any type for which a
    ///   Teuchos::ScalarTraits specialization exists is valid.  If
    ///   Scalar is complex (isComplex=true), then Scalar's operator=
    ///   must accept an std::complex<typename
    ///   Teuchos::ScalarTraits<Scalar>::magnitudeType> input.
    ///
    /// \tparam isComplex Whether Scalar is a complex-valued type or a
    ///   real-valued type.  The default value here is fine for most
    ///   users.
    ///
    /// This class completes the implementation of CoordDataReaderBase
    /// for sparse matrices.  We provide two partial specializations:
    /// one for real-valued data (isComplex=false), and the other for
    /// complex-valued data (isComplex=true).
    ///
    /// If you have an existing sparse matrix implementation to which
    /// you want to add entries, you'll either have to make it
    /// implement Callback's three-argument operator(), or write your
    /// own Callback implementation that wraps the sparse matrix.
    template<class Callback,
             class Ordinal,
             class Scalar,
             bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex>
    class CoordDataReader : public CoordDataReaderBase<Callback, Ordinal> {
    public:
      /// \brief Constructor with "adder" argument.
      ///
      /// This is the favored way to construct an instance of this
      /// type.  Only use the no-argument constructor if you have a
      /// "chicken-and-egg" problem, where in order to create the
      /// Callback instance, you need the matrix dimensions.  This is
      /// the case if your Callback wraps Tpetra::CrsMatrix, for
      /// example.
      ///
      /// \param adder [in/out] Closure (a.k.a. callback) whose
      ///   operator() adds an entry to the sparse graph or matrix on
      ///   each invocation.
      CoordDataReader (const Teuchos::RCP<Callback>& adder);

      /// \brief No-argument constructor.
      ///
      /// We offer this option in case the adder's constructor needs
      /// the matrix dimensions, so that it's necessary to call
      /// readDimensions() first before constructing the adder.  You
      /// should call setAdder() with a non-null argument before
      /// calling read() or readLine().
      CoordDataReader ();

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReader();

    protected:
      bool
      readLine (const std::string& theLine,
                const size_t lineNumber,
                const bool tolerant);
    };

#ifdef HAVE_TEUCHOS_COMPLEX
    // Partial specialization for complex Scalar types.
    template<class Callback, class Ordinal, class Scalar>
    class CoordDataReader<Callback, Ordinal, Scalar, true> :
      public CoordDataReaderBase<Callback, Ordinal> {
    public:
      CoordDataReader (const Teuchos::RCP<Callback>& adder) :
        CoordDataReaderBase<Callback, Ordinal> (adder)
      {}

      CoordDataReader() :
        CoordDataReaderBase<Callback, Ordinal> (null)
      {}

      virtual ~CoordDataReader() {};

    protected:
      bool
      readLine (const std::string& theLine,
                const size_t lineNumber,
                const bool tolerant)
      {
        typedef Teuchos::ScalarTraits<Scalar> STS;
        typedef typename STS::magnitudeType Real;

        Ordinal rowIndex;
        Ordinal colIndex;
        Scalar value;

        Real realPart, imagPart;
        const bool localSuccess =
          readComplexLine (theLine, rowIndex, colIndex, realPart, imagPart,
                           lineNumber, tolerant);
        if (localSuccess) {
          // Assume that assignment from std::complex<Real> to Scalar
          // (which itself is complex-valued) is valid.  We have to do
          // this, since the C++ compiler may not be smart enough to
          // assume here (when it instantiates the templates) that
          // Scalar is an std::complex<Real> -- even though it has to
          // be, if STS::isComplex is true (which as of 31 Jan 2011,
          // only holds for std::complex<T>).
          value = std::complex<Real> (realPart, imagPart);

          // Now that we've read in the (i, j, A_ij) triple
          // successfully, we can add the entry to the sparse matrix.
          (*(this->adder_)) (rowIndex, colIndex, value);
        }
        return localSuccess;
      }
    };
#endif // HAVE_TEUCHOS_COMPLEX

    // Partial specialization for real Scalar types.
    template<class Callback, class Ordinal, class Scalar>
    class CoordDataReader<Callback, Ordinal, Scalar, false> :
      public CoordDataReaderBase<Callback, Ordinal> {
    public:
      CoordDataReader (const Teuchos::RCP<Callback>& adder) :
        CoordDataReaderBase<Callback, Ordinal> (adder)
      {}

      CoordDataReader() :
        CoordDataReaderBase<Callback, Ordinal> (null)
      {}

      virtual ~CoordDataReader() {};

    protected:
      bool
      readLine (const std::string& theLine,
                const size_t lineNumber,
                const bool tolerant)
      {
        Ordinal rowIndex;
        Ordinal colIndex;
        Scalar value;
        const bool localSuccess = readRealLine (theLine, rowIndex, colIndex,
                                                value, lineNumber, tolerant);
        if (localSuccess) {
          // Now that we've read in the (i, j, A_ij) triple
          // successfully, we can add the entry to the sparse matrix.
          (*(this->adder_)) (rowIndex, colIndex, value);
        }
        return localSuccess;
      }
    };


    /// \class CoordPatternReader
    /// \brief Coordinate-format sparse graph data reader.
    ///
    /// Use this class to read in sparse graph data and add it to the
    /// sparse graph using your Callback implementation (see the
    /// documentation of CoordDataReaderBase).  CoordDataReaderBase is
    /// sufficiently general that we can extend it to read in a sparse
    /// graph.  "Pattern" refers to the Matrix Market keyword
    /// "pattern" that indicates a sparse graph.
    ///
    /// \tparam Callback Same as the \c Callback template parameter of
    ///   \c CoordDataReaderBase.  The type of the callback
    ///   (a.k.a. closure) for adding entries to the sparse graph.
    ///
    /// \tparam Ordinal Same as the \c Ordinal template parameter of
    ///   \c CoordDataReaderBase.  The type of indices of the sparse
    ///   graph.
    ///
    /// If you have an existing sparse graph implementation to which
    /// you want to add entries, you'll either have to make it
    /// implement Callback's two-argument operator(), or write your
    /// own Callback implementation that wraps the sparse graph.
    template<class Callback, class Ordinal>
    class CoordPatternReader : public CoordDataReaderBase<Callback, Ordinal> {
    public:
      /// \brief Constructor with "adder" argument.
      ///
      /// This is the favored way to construct an instance of this
      /// type.  Only use the no-argument constructor if you have a
      /// "chicken-and-egg" problem, where in order to create the
      /// Callback instance, you need the graph dimensions.  This is
      /// the case if your Callback wraps Tpetra::CrsGraph, for
      /// example.
      ///
      /// \param adder [in/out] Closure (a.k.a. callback) whose
      ///   operator() adds an entry to the sparse graph on each
      ///   invocation.
      CoordPatternReader (const Teuchos::RCP<Callback>& adder) :
        CoordDataReaderBase<Callback, Ordinal> (adder)
      {}

      /// \brief No-argument constructor.
      ///
      /// We offer this option in case the adder's constructor needs
      /// the graph dimensions, so that it's necessary to call
      /// readDimensions() first before constructing the adder.  You
      /// should call setAdder() with a non-null argument before
      /// calling read() or readLine().
      CoordPatternReader() :
        CoordDataReaderBase<Callback, Ordinal> (null)
      {}

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordPatternReader() {};

    protected:
      bool
      readLine (const std::string& theLine,
                const size_t lineNumber,
                const bool tolerant)
      {
        Ordinal rowIndex;
        Ordinal colIndex;
        const bool localSuccess =
          readPatternLine (theLine, rowIndex, colIndex, lineNumber, tolerant);
        if (localSuccess) {
          // Now that we've read in the (i, j) pair successfully, we
          // can add the entry to the sparse graph.
          (*(this->adder_)) (rowIndex, colIndex);
        }
        return localSuccess;
      }
    };

  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_CoordDataReader_hpp
