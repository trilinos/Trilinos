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

#ifndef __MatrixMarket_raw_hpp
#define __MatrixMarket_raw_hpp

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "MatrixMarket_Banner.hpp"
#include "MatrixMarket_CoordDataReader.hpp"
#include "MatrixMarket_util.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>


namespace Teuchos {
  namespace MatrixMarket {
    /// \namespace Raw
    /// \brief "Raw" input of sparse matrices from Matrix Market files.
    ///
    /// "Raw" means serial (not MPI or otherwise distributed over
    /// parallel processes), with storage as a collection of matrix
    /// indices and values.  This is useful if you want to read the
    /// sparse matrix on one (MPI) process and store it in a custom
    /// format.  For reading in a Tpetra::CrsMatrix, use the
    /// Tpetra::MatrixMarket::Reader class.
    namespace Raw {
      /// \class Element
      /// \author Mark Hoemmen
      /// \brief Stores one entry of a sparse matrix.
      ///
      /// \tparam Scalar The type of entries of the sparse matrix.
      /// \tparam Ordinal The type of indices of the sparse matrix.
      ///
      /// This class is mainly useful as an implementation detail of
      /// Adder.  We expose it to users only if they wish to convert
      /// the sparse matrix read in by Adder into a storage format
      /// other than CSR (compressed sparse row).
      ///
      /// An array of Elements implements the so-called "array of
      /// structs" representation of a coordinate format sparse
      /// matrix.  An Element has a row and column index (each of type
      /// Ordinal) and a value (of type Scalar).  Elements also have
      /// equality and ordering comparisons.  The equality comparison
      /// only tests the row and column index, and is intended to
      /// simplify merging matrix entries with the same row and column
      /// indices.  The ordering comparison means that std::sort of a
      /// sequence of Elements will put them in an order suitable for
      /// extracting the CSR (compressed sparse row) representation of
      /// the sparse matrix.
      template<class Scalar, class Ordinal>
      class Element {
      public:
        /// \brief Default constructor: an invalid entry of the matrix.
        ///
        /// FIXME (mfh 31 May 2012) This currently only works for
        /// signed Ordinal types to which -1 can be assigned.  It
        /// would be better to use
        /// Teuchos::OrdinalTraits<Ordinal>::invalid() as the default
        /// invalid value.
        Element () : rowIndex_ (-1), colIndex_ (-1), value_ (0) {}

        //! Create a sparse matrix entry at (i,j) with value Aij.
        Element (const Ordinal i, const Ordinal j, const Scalar& Aij) :
          rowIndex_ (i), colIndex_ (j), value_ (Aij) {}

        //! Ignore the matrix value for comparisons.
        bool operator== (const Element& rhs) {
          return rowIndex_ == rhs.rowIndex_ && colIndex_ == rhs.colIndex_;
        }

        //! Ignore the matrix value for comparisons.
        bool operator!= (const Element& rhs) {
          return ! (*this == rhs);
        }

        //! Lexicographic order first by row index, then by column index.
        bool operator< (const Element& rhs) const {
          if (rowIndex_ < rhs.rowIndex_)
            return true;
          else if (rowIndex_ > rhs.rowIndex_)
            return false;
          else { // equal
            return colIndex_ < rhs.colIndex_;
          }
        }

        /// \brief Merge rhs into this Element, using custom binary function.
        ///
        /// This replaces the current value Aij with f(rhs.value_,
        /// Aij).  The object f must be a binary function that takes
        /// two Scalar arguments and returns a Scalar.
        template<class BinaryFunction>
        void merge (const Element& rhs, const BinaryFunction& f) {
          if (rowIndex() != rhs.rowIndex() || colIndex() != rhs.colIndex()) {
            throw std::invalid_argument ("Can only merge elements at the same "
                                         "location in the sparse matrix");
          }
          else {
            value_ = f (rhs.value_, value_);
          }
        }

        /// \brief Merge rhs into this Element.
        ///
        /// "Replace" means replace this Element's value with that of
        /// rhs.  Otherwise, this Element's value is added to rhs's
        /// value.
        void merge (const Element& rhs, const bool replace=false) {
          if (rowIndex() != rhs.rowIndex() || colIndex() != rhs.colIndex())
            throw std::logic_error("Can only merge elements at the same "
                                   "location in the sparse matrix");
          else if (replace)
            value_ = rhs.value_;
          else
            value_ += rhs.value_;
        }

        //! Row index (zero-based) of this Element.
        Ordinal rowIndex() const { return rowIndex_; }

        //! Column index (zero-based) of this Element.
        Ordinal colIndex() const { return colIndex_; }

        //! Value (A(rowIndex(), colIndex()) of this Element.
        Scalar value() const { return value_; }

      private:
        Ordinal rowIndex_, colIndex_;
        Scalar value_;
      };

      /// \brief Print out an Element to the given output stream.
      ///
      /// This method is suitable for printing a sparse matrix to a
      /// Matrix Market file.  We try to print out floating-point
      /// values with enough digits to reproduce the results, but in
      /// general this routine does not promise that the matrix entry
      /// later read in from this printing will be bitwise identical
      /// to the matrix entry supplied as input.  There _are_ printing
      /// algorithms that make this guarantee; we just haven't
      /// implemented them yet here.
      template<class Scalar, class Ordinal>
      std::ostream&
      operator<< (std::ostream& out, const Element<Scalar, Ordinal>& elt)
      {
        typedef ScalarTraits<Scalar> STS;
        // Non-Ordinal types are floating-point types.  In order not to
        // lose information when we print a floating-point type, we have
        // to set the number of digits to print.  C++ standard behavior
        // in the default locale seems to be to print only five decimal
        // digits after the decimal point; this does not suffice for
        // double precision.  We solve the problem of how many digits to
        // print more generally below.  It's a rough solution so please
        // feel free to audit and revise it.
        //
        // FIXME (mfh 01 Feb 2011)
        // This really calls for the following approach:
        //
        // Guy L. Steele and Jon L. White, "How to print floating-point
        // numbers accurately", 20 Years of the ACM/SIGPLAN Conference
        // on Programming Language Design and Implementation
        // (1979-1999): A Selection, 2003.
        if (! STS::isOrdinal) {
          // std::scientific, std::fixed, and default are the three
          // output states for floating-point numbers.  A reasonable
          // user-defined floating-point type should respect these
          // flags; hopefully it does.
          out << std::scientific;

          // Decimal output is standard for Matrix Market format.
          out << std::setbase (10);

          // Compute the number of decimal digits required for expressing
          // a Scalar, by comparing with IEEE 754 double precision (16
          // decimal digits, 53 binary digits).  This would be easier if
          // Teuchos exposed std::numeric_limits<T>::digits10, alas.
          const double numDigitsAsDouble =
            16 * ((double) STS::t() / (double) ScalarTraits<double>::t());
          // Adding 0.5 and truncating is a portable "floor".
          const int numDigits = static_cast<int> (numDigitsAsDouble + 0.5);

          // Precision to which a Scalar should be written.  Add one
          // for good measure, since 17 is necessary for IEEE 754
          // doubles.
          out << std::setprecision (numDigits + 1);
        }
        out << elt.rowIndex () << " " << elt.colIndex () << " ";
        if (STS::isComplex) {
          out << STS::real (elt.value ()) << " " << STS::imag (elt.value ());
        }
        else {
          out << elt.value ();
        }
        return out;
      }

      /// \class Adder
      /// \brief To be used with Reader for "raw" sparse matrix input.
      ///
      /// \tparam Scalar The type of entries in the sparse matrix.
      /// \tparam Ordinal The type of indices in the sparse matrix.
      ///
      /// This class implements the interface required by the Callback
      /// template parameter of
      /// Teuchos::MatrixMarket::CoordDataReader.  It provides a
      /// simple implementation of this interface which is useful for
      /// things like printing out a sparse matrix's entries, or
      /// converting between storage formats.
      template<class Scalar, class Ordinal>
      class Adder {
      public:
        typedef Ordinal index_type;
        typedef Scalar value_type;
        typedef Element<Scalar, Ordinal> element_type;
        typedef typename std::vector<element_type>::size_type size_type;

        /// \brief Default constructor.
        ///
        /// If you call the default constructor, we assume that you
        /// want tolerant mode (in which the Adder tries to infer the
        /// matrix dimensions and number of entries from the actual
        /// matrix data, not from any metadata).  Tolerant mode is
        /// similar to what Matlab does if you give it an ASCII file
        /// of (i,j,Aij) triples.  It may get the matrix dimensions
        /// (m,n) wrong if the lower right entry of the matrix is zero
        /// and is not supplied explicitly by calling operator().
        Adder () :
          expectedNumRows_(0),
          expectedNumCols_(0),
          expectedNumEntries_(0),
          seenNumRows_(0),
          seenNumCols_(0),
          seenNumEntries_(0),
          tolerant_ (true),
          debug_ (false)
        {}

        /// \brief Standard constructor.
        ///
        /// \param expectedNumRows [in] Number of rows in the matrix,
        ///   as specified by the matrix metadata.
        ///
        /// \param expectedNumCols [in] Number of columns in the
        ///   matrix, as specified by the matrix metadata.
        ///
        /// \param expectedNumEntries [in] Number of entries in the
        ///   matrix, as specified by the matrix metadata.
        ///
        /// \param tolerant [in] Whether the "expected" metadata is
        ///   required to match what the read-in matrix entries tell
        ///   us.
        ///
        /// \param debug [in] If true, we may print copious status
        ///   output for debugging purposes.
        Adder (const Ordinal expectedNumRows,
               const Ordinal expectedNumCols,
               const Ordinal expectedNumEntries,
               const bool tolerant=false,
               const bool debug=false) :
          expectedNumRows_(expectedNumRows),
          expectedNumCols_(expectedNumCols),
          expectedNumEntries_(expectedNumEntries),
          seenNumRows_(0),
          seenNumCols_(0),
          seenNumEntries_(0),
          tolerant_ (tolerant),
          debug_ (debug)
        {}

        /// \brief Add an entry to the sparse matrix.
        ///
        /// If tolerant==false, this method will perform error
        /// checking to ensure that the matrix data matches the
        /// metadata.  For example, it will check that i and j are in
        /// bounds, and that we haven't added more than the expected
        /// number of matrix entries.  Regardless, this method will
        /// update the "actual" metadata.
        ///
        /// \param i [in] (1-based) row index
        /// \param j [in] (1-based) column index
        /// \param Aij [in] Value of the entry A(i,j)
        void
        operator() (const Ordinal i, const Ordinal j, const Scalar& Aij)
        {
          if (! tolerant_) {
            const bool indexPairOutOfRange = i < 1 || j < 1 ||
              i > expectedNumRows_ || j > expectedNumCols_;

            TEUCHOS_TEST_FOR_EXCEPTION(indexPairOutOfRange,
              std::invalid_argument, "Matrix is " << expectedNumRows_ << " x "
              << expectedNumCols_ << ", so entry A(" << i << "," << j << ") = "
              << Aij << " is out of range.");
            TEUCHOS_TEST_FOR_EXCEPTION(seenNumEntries_ >= expectedNumEntries_,
              std::invalid_argument, "Cannot add entry A(" << i << "," << j
              << ") = " << Aij << " to matrix; already have expected number "
              "of entries " << expectedNumEntries_ << ".");
          }
          // i and j are 1-based indices, but we store them as 0-based.
          elts_.push_back (element_type (i-1, j-1, Aij));

          // Keep track of the rightmost column containing a matrix
          // entry, and the bottommost row containing a matrix entry.
          // This gives us a lower bound for the dimensions of the
          // matrix, and a check for the reported dimensions of the
          // matrix in the Matrix Market file.
          seenNumRows_ = std::max (seenNumRows_, i);
          seenNumCols_ = std::max (seenNumCols_, j);
          seenNumEntries_++;
        }

        /// \brief Print the sparse matrix data.
        ///
        /// We always print the data sorted.  You may also merge
        /// duplicate entries if you prefer.
        ///
        /// \param out [out] Output stream to which to print
        /// \param doMerge [in] Whether to merge entries before printing
        /// \param replace [in] If merging, whether to replace duplicate
        ///   entries; otherwise their values are added together.
        void
        print (std::ostream& out, const bool doMerge, const bool replace=false)
        {
          if (doMerge) {
            merge (replace);
          } else {
            std::sort (elts_.begin(), elts_.end());
          }
          // Print out the results, delimited by newlines.
          typedef std::ostream_iterator<element_type> iter_type;
          std::copy (elts_.begin(), elts_.end(), iter_type (out, "\n"));
        }

        /// \brief Merge duplicate elements.
        ///
        /// Merge elements of the sparse matrix that have the same row
        /// and column indices ("duplicates").  Resize the array of
        /// elements to fit just the "unique" (not duplicate)
        /// elements.
        ///
        /// \param replace [in] If true, replace each duplicate
        ///   element with the next element sharing the same row and
        ///   column index.  This means that results will depend on
        ///   the order in which the duplicate elements were added.
        ///   Otherwise, duplicate elements have their values added
        ///   together; in that case, the result is independent (in
        ///   exact arithmetic, not in finite-precision arithmetic) of
        ///   their order.
        ///
        /// \return (# unique elements, # removed elements)
        ///
        /// \note This method does not change the "expected" or "seen"
        ///   numbers of entries, since both of those count entries
        ///   with the same row and column indices as separate
        ///   entries.
        std::pair<size_type, size_type>
        merge (const bool replace=false)
        {
          typedef typename std::vector<element_type>::iterator iter_type;

          // Start with a sorted container.  Element objects sort in
          // lexicographic order of their (row, column) indices, for
          // easy conversion to CSR format.  If you expect that the
          // elements will usually be sorted in the desired order, you
          // can check first whether they are already sorted.  We have
          // no such expectation, so we don't even bother to spend the
          // extra O(# entries) operations to check.
          std::sort (elts_.begin(), elts_.end());

          // Walk through the array of elements in place, merging
          // duplicates and pushing unique elements up to the front of
          // the array.  We can't use std::unique for this because it
          // doesn't let us merge duplicate elements; it only removes
          // them from the sequence.
          size_type numUnique = 0;
          iter_type cur = elts_.begin();
          if (cur == elts_.end()) { // No elements to merge
            return std::make_pair (numUnique, size_type (0));
          }
          else {
            iter_type next = cur;
            ++next; // There is one unique element
            ++numUnique;
            while (next != elts_.end()) {
              if (*cur == *next) {
                // Merge in the duplicated element *next
                cur->merge (*next, replace);
              } else {
                // *cur is already a unique element.  Move over one to
                // *make space for the new unique element.
                ++cur;
                *cur = *next; // Add the new unique element
                ++numUnique;
              }
              // Look at the "next" not-yet-considered element
              ++next;
            }
            // Remember how many elements we removed before resizing.
            const size_type numRemoved = elts_.size() - numUnique;
            elts_.resize (numUnique);
            return std::make_pair (numUnique, numRemoved);
          }
        }

        /// \brief Merge duplicate elements and convert to zero-based CSR.
        ///
        /// Merge elements of the sparse matrix that have the same row
        /// and column indices ("duplicates").  Resize the array of
        /// elements to fit just the "unique" (not duplicate)
        /// elements.  Return a CSR (compressed sparse row) version of
        /// the data, with zero-based indices.
        ///
        /// We combine merge and conversion to CSR because the latter
        /// requires the former.
        ///
        /// \param numUniqueElts [out] Same as the first return value
        ///   of merge().
        ///
        /// \param numRemovedElts [out] Same as the second return
        ///   value of merge().
        ///
        /// \param rowptr [out] Array of numRows+1 offsets, where
        ///   numRows is the number of rows in the sparse matrix.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        ///
        /// \param colind [out] Column indices of the matrix.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        ///
        /// \param values [out] Values stored in the matrix.
        ///
        /// \param replace [in] If true, replace each duplicate
        ///   element with the next element sharing the same row and
        ///   column index.  This means that results will depend on
        ///   the order in which the duplicate elements were added.
        ///   Otherwise, duplicate elements have their values added
        ///   together; in that case, the result is independent (in
        ///   exact arithmetic, not in finite-precision arithmetic) of
        ///   their order.
        ///
        /// \note This method does not change the "expected" or "seen"
        ///   numbers of entries, since both of those count entries
        ///   with the same row and column indices as separate
        ///   entries.
        void
        mergeAndConvertToCSR (size_type& numUniqueElts,
                              size_type& numRemovedElts,
                              Teuchos::ArrayRCP<Ordinal>& rowptr,
                              Teuchos::ArrayRCP<Ordinal>& colind,
                              Teuchos::ArrayRCP<Scalar>& values,
                              const bool replace=false)
        {
          using Teuchos::arcp;
          using Teuchos::ArrayRCP;

          std::pair<size_type, size_type> mergeResult = merge (replace);

          // At this point, elts_ is already in CSR order.
          // Now we can allocate and fill the ind and val arrays.
          ArrayRCP<Ordinal> ind = arcp<Ordinal> (elts_.size ());
          ArrayRCP<Scalar> val = arcp<Scalar> (elts_.size ());

          // Number of rows in the matrix.
          const Ordinal numRows = tolerant_ ? seenNumRows_ : expectedNumRows_;
          ArrayRCP<Ordinal> ptr = arcp<Ordinal> (numRows + 1);

          // Copy over the elements, and fill in the ptr array with
          // offsets.  Note that merge() sorted the entries by row
          // index, so we can assume the row indices are increasing in
          // the list of entries.
          Ordinal curRow = 0;
          Ordinal curInd = 0;
          typedef typename std::vector<element_type>::const_iterator iter_type;
          for (iter_type it = elts_.begin(); it != elts_.end(); ++it) {
            const Ordinal i = it->rowIndex ();
            const Ordinal j = it->colIndex ();
            const Scalar Aij = it->value ();

            TEUCHOS_TEST_FOR_EXCEPTION(i < curRow, std::logic_error, "The "
              "current matrix entry's row index " << i << " is less then what "
              "should be the current row index lower bound " << curRow << ".");
            for (Ordinal k = curRow; k < i; ++k) {
              ptr[k] = curInd;
            }
            curRow = i;

            TEUCHOS_TEST_FOR_EXCEPTION(curInd >= elts_.size (),
              std::logic_error, "The current index " << curInd << " into ind "
              "and val is >= the number of matrix entries " << elts_.size ()
              << ".");
            ind[curInd] = j;
            val[curInd] = Aij;
            ++curInd;
          }
          for (Ordinal k = curRow; k <= numRows; ++k) {
            ptr[k] = curInd;
          }

          // Assign to outputs here, to ensure the strong exception
          // guarantee (assuming that ArrayRCP's operator= doesn't
          // throw).
          rowptr = ptr;
          colind = ind;
          values = val;
          numUniqueElts = mergeResult.first;
          numRemovedElts = mergeResult.second;
        }

        //! A temporary const view of the entries of the matrix.
        const std::vector<element_type>& getEntries() const {
          return elts_;
        }

        //! Clear all the added matrix entries and reset metadata.
        void clear() {
          seenNumRows_ = 0;
          seenNumCols_ = 0;
          seenNumEntries_ = 0;
          elts_.resize (0);
        }

        //! Computed number of rows.
        const Ordinal numRows() const { return seenNumRows_; }

        //! Computed number of columns.
        const Ordinal numCols() const { return seenNumCols_; }

      private:
        Ordinal expectedNumRows_, expectedNumCols_, expectedNumEntries_;
        Ordinal seenNumRows_, seenNumCols_, seenNumEntries_;
        bool tolerant_;
        bool debug_;
        std::vector<element_type> elts_;
      };

      /// \class Reader
      /// \brief "Raw" reader for debugging a Matrix Market file.
      ///
      /// This class' methods are useful for examining the contents of
      /// a Matrix Market file, and checking the integrity of its
      /// data.  See MatrixMarket_Tpetra.hpp for a Matrix Market
      /// reader that constructs a Tpetra::CrsMatrix object.
      template<class Scalar, class Ordinal>
      class Reader {
      public:
        /// \brief Read the sparse matrix from the given file.
        ///
        /// This is a collective operation.  Only Rank 0 opens the
        /// file and reads data from it, but all ranks participate and
        /// wait for the final result.
        ///
        /// \note This whole "raw" reader is meant for debugging and
        ///   diagnostics of syntax errors in the Matrix Market file;
        ///   it's not performance-oriented.  That's why we do all the
        ///   broadcasts of and checks for "success".
        static bool
        readFile (const Teuchos::Comm<int>& comm,
                  const std::string& filename,
                  const bool echo,
                  const bool tolerant,
                  const bool debug=false)
        {
          using std::cerr;
          using std::endl;

          const int myRank = Teuchos::rank (comm);
          // Teuchos::broadcast doesn't accept a bool; we use an int
          // instead, with the usual 1->true, 0->false Boolean
          // interpretation.
          int readFile = 0;
          RCP<std::ifstream> in; // only valid on Rank 0
          if (myRank == 0) {
            if (debug) {
              cerr << "Attempting to open file \"" << filename
                   << "\" on Rank 0...";
            }
            in = rcp (new std::ifstream (filename.c_str()));
            if (! *in) {
              readFile = 0;
              if (debug)
                cerr << "failed." << endl;
            }
            else {
              readFile = 1;
              if (debug)
                cerr << "succeeded." << endl;
            }
          }
          Teuchos::broadcast (comm, 0, &readFile);
          TEUCHOS_TEST_FOR_EXCEPTION(! readFile, std::runtime_error,
                             "Failed to open input file \"" + filename + "\".");
          // Only Rank 0 will try to dereference "in".
          return read (comm, in, echo, tolerant, debug);
        }


        /// \brief Read the sparse matrix from the given input stream.
        ///
        /// This is a collective operation.  Only Rank 0 reads from
        /// the given input stream, but all ranks participate and wait
        /// for the final result.
        ///
        /// \note This whole "raw" reader is meant for debugging and
        ///   diagnostics of syntax errors in the Matrix Market file;
        ///   it's not performance-oriented.  That's why we do all the
        ///   broadcasts of and checks for "success".
        static bool
        read (const Teuchos::Comm<int>& comm,
              const RCP<std::istream>& in,
              const bool echo,
              const bool tolerant,
              const bool debug=false)
        {
          using std::cerr;
          using std::endl;

          const int myRank = Teuchos::rank (comm);
          std::pair<bool, std::string> result;
          int msgSize = 0; // Size of error message (if any)
          if (myRank == 0) {
            if (in.is_null()) {
              result.first = false;
              result.second = "Input stream is null on Rank 0";
            }
            else {
              if (debug) {
                cerr << "About to read from input stream on Rank 0" << endl;
              }
              result = readOnRank0 (*in, echo, tolerant, debug);
              if (debug) {
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
            if (! tolerant) {
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
                Teuchos::broadcast (comm, 0, static_cast<int>(msgSize), &errMsg[0]);
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
        /// \brief Read in the Banner line from the given input stream.
        ///
        /// Only call this method on Rank 0.
        ///
        /// \param in [in/out] Input stream from which to read the
        ///   Banner line.
        ///
        /// \param lineNumber [in/out] On input: Current line number
        ///   of the input stream.  On output: if any line(s) were
        ///   successfully read from the input stream, this is
        ///   incremented by the number of line(s) read.  (This
        ///   includes comment lines.)
        ///
        /// \return Banner [non-null]
        static RCP<const Teuchos::MatrixMarket::Banner>
        readBanner (std::istream& in,
                    size_t& lineNumber,
                    const bool tolerant=false,
                    const bool debug=false)
        {
          typedef ScalarTraits<Scalar> STS;

          // The pointer will be non-null on return only on MPI Rank 0.
          // Using a pointer lets the data persist outside the
          // "myRank==0" scopes.
          RCP<Banner> pBanner;

          // Keep reading lines until we get a noncomment line.
          std::string line;
          size_t numLinesRead = 0;
          bool commentLine = false;
          do {
            // Try to read a line from the input stream.
            const bool readFailed = ! getline(in, line);
            TEUCHOS_TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
                               "Failed to get Matrix Market banner line "
                               "from input, after reading " << numLinesRead
                               << "line" << (numLinesRead != 1 ? "s." : "."));
            // We read a line from the input stream.
            lineNumber++;
            numLinesRead++;
            size_t start, size; // Output args of checkCommentLine
            commentLine = checkCommentLine (line, start, size,
                                            lineNumber, tolerant);
          } while (commentLine); // Loop until we find a noncomment line.

          // Assume that the noncomment line we found is the banner line.
          try {
            pBanner = rcp (new Banner (line, tolerant));
          } catch (std::exception& e) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                               "Matrix Market banner line contains syntax "
                               "error(s): " << e.what());
          }
          return pBanner;
        }

        //! To be called only on MPI Rank 0.
        static std::pair<bool, std::string>
        readOnRank0 (std::istream& in,
                     const bool echo,
                     const bool tolerant,
                     const bool debug=false)
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
            pBanner = readBanner (in, lineNumber, tolerant, debug);
          }
          catch (std::exception& e) {
            err << "Failed to read Banner: " << e.what();
            return std::make_pair (false, err.str());
          }
          //
          // Validate the metadata in the Banner.
          //
          if (pBanner->matrixType() != "coordinate") {
            err << "Matrix Market input file must contain a "
              "\"coordinate\"-format sparse matrix in "
              "order to create a sparse matrix object "
              "from it.";
            return std::make_pair (false, err.str());
          }
          else if (! STS::isComplex && pBanner->dataType() == "complex") {
            err << "The Matrix Market sparse matrix file contains complex-"
              "valued data, but you are try to read the data into a sparse "
              "matrix containing real values (your matrix's Scalar type is "
              "real).";
            return std::make_pair (false, err.str());
          }
          else if (pBanner->dataType() != "real" &&
                   pBanner->dataType() != "complex") {
            err << "Only real or complex data types (no pattern or integer "
              "matrices) are currently supported.";
            return std::make_pair (false, err.str());
          }
          if (debug) {
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
            reader.readDimensions (in, lineNumber, tolerant);
          if (! dims.second) {
            err << "Error reading Matrix Market sparse matrix "
              "file: failed to read coordinate dimensions.";
            return std::make_pair (false, err.str());
          }
          // These are "expected" values read from the input stream's
          // metadata.  The actual matrix entries read from the input
          // stream might not conform to their constraints.  We allow
          // such nonconformity only in "tolerant" mode; otherwise, we
          // throw an exception.
          const Ordinal numRows = dims.first[0];
          const Ordinal numCols = dims.first[1];
          const Ordinal numEntries = dims.first[2];
          if (debug) {
            cerr << "Reported dimensions: " << numRows << " x " << numCols
                 << ", with " << numEntries << " entries (counting possible "
                 << "duplicates)." << endl;
          }

          // The "raw" adder knows about the expected matrix
          // dimensions, but doesn't know about symmetry.
          RCP<raw_adder_type> rawAdder =
            rcp (new raw_adder_type (numRows, numCols, numEntries,
                                     tolerant, debug));
          // The symmetrizing adder knows about symmetry.
          RCP<adder_type> adder =
            rcp (new adder_type (rawAdder, pBanner->symmType()));

          // Give the adder to the reader.
          reader.setAdder (adder);

          // Read the sparse matrix entries.  "results" just tells us if
          // and where there were any bad lines of input.  The actual
          // sparse matrix entries are stored in the (raw) Adder object.
          std::pair<bool, std::vector<size_t> > results =
            reader.read (in, lineNumber, tolerant, debug);
          if (debug) {
            if (results.first) {
              cerr << "Matrix Market file successfully read" << endl;
            }
            else {
              cerr << "Failed to read Matrix Market file" << endl;
            }
          }

          // Report any bad line number(s).
          if (! results.first) {
            if (! tolerant) {
              err << "The Matrix Market input stream had syntax error(s)."
                "  Here is the error report." << endl;
              reportBadness (err, results);
              err << endl;
              return std::make_pair (false, err.str());
            }
            else {
              if (debug) {
                reportBadness (cerr, results);
              }
            }
          }
          // We're done reading in the sparse matrix.  If we're in
          // "echo" mode, print out the matrix entries to stdout.  The
          // entries will have been symmetrized if applicable.
          if (echo) {
            const bool doMerge = false;
            const bool replace = false;
            rawAdder->print (cout, doMerge, replace);
            cout << endl;
          }
          return std::make_pair (true, err.str());
        }

        //! To be called only on MPI Rank 0.
        static void
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
      };

    } // namespace Raw
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __MatrixMarket_raw_hpp
