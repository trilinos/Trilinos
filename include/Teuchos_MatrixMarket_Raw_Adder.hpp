// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Adder_hpp
#define __Teuchos_MatrixMarket_Raw_Adder_hpp

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_MatrixMarket_Banner.hpp"
#include "Teuchos_MatrixMarket_CoordDataReader.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>


namespace Teuchos {
  namespace MatrixMarket {
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
        //! Default constructor: an invalid entry of the matrix.
        Element () :
          rowIndex_ (Teuchos::OrdinalTraits<Ordinal>::invalid ()),
          colIndex_ (Teuchos::OrdinalTraits<Ordinal>::invalid ()),
          value_ (Teuchos::ScalarTraits<Scalar>::zero ())
        {}

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
          TEUCHOS_TEST_FOR_EXCEPTION(
            rowIndex() != rhs.rowIndex() || colIndex() != rhs.colIndex(),
            std::invalid_argument,
            "Attempt to merge elements at different locations in the sparse "
            "matrix.  The current element is at (" << rowIndex() << ", "
            << colIndex() << ") and the element you asked me to merge with it "
            "is at (" << rhs.rowIndex() << ", " << rhs.colIndex() << ").  This "
            "probably indicates a bug in the sparse matrix reader.");

          value_ = f (rhs.value_, value_);
        }

        /// \brief Merge rhs into this Element, either by addition or replacement.
        ///
        /// \param rhs [in] Element to merge in.
        /// \param replace [in] If true, replace this Element's value
        ///   with that of rhs.  If false, add rhs to this Element's
        ///   value.
        void merge (const Element& rhs, const bool replace=false) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            rowIndex() != rhs.rowIndex() || colIndex() != rhs.colIndex(),
            std::invalid_argument,
            "Attempt to merge elements at different locations in the sparse "
            "matrix.  The current element is at (" << rowIndex() << ", "
            << colIndex() << ") and the element you asked me to merge with it "
            "is at (" << rhs.rowIndex() << ", " << rhs.colIndex() << ").  This "
            "probably indicates a bug in the sparse matrix reader.");

          if (replace) {
            value_ = rhs.value_;
          }
          else {
            value_ += rhs.value_;
          }
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
        std::ios::fmtflags f( out.flags() );
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
        // Restore flags
        out.flags( f );
        return out;
      }

      /// \class Adder
      /// \brief To be used with Checker for "raw" sparse matrix input.
      ///
      /// \tparam Scalar The type of entries in the sparse matrix.
      /// \tparam Ordinal The type of indices in the sparse matrix.
      ///
      /// This class implements the following interface, which is
      /// required by the Callback template parameter of
      /// Teuchos::MatrixMarket::CoordDataReader:
      /// \code
      /// class AdderType {
      /// public:
      ///   typedef ... index_type; // Ellipsis represents the actual type
      ///   typedef ... value_type; // Ellipsis represents the actual type
      ///   void operator() (const index_type, const index_type, const value_type&);
      /// };
      /// \endcode
      /// For Adder, the Scalar template parameter is value_type, and
      /// the Ordinal template parameter is index_type.  Adder
      /// provides a simple implementation of the above interface
      /// which is useful for things like printing out a sparse
      /// matrix's entries, or converting between storage formats.
      ///
      /// It is possible to nest classes that implement the above
      /// interface, in order to modify the definition of inserting
      /// values into a sparse matrix.  (If you are familiar with
      /// Emacs Lisp, this is called "advising" (the insertion
      /// function, in this case).  See the
      /// <a href="http://www.gnu.org/software/emacs/manual/html_node/elisp/Advising-Functions.html">Emacs Lisp Manual</a>
      /// for details.)  If you are building a chain of classes, each
      /// of which implements the above interface by calling the next
      /// lower class' operator() beneath it, this class is a good
      /// start, since it implements insertion of values directly.
      /// SymmetrizingAdder is an example; it "advises" Adder by
      /// inserting an entry at (j,i) whenever Adder inserts an entry
      /// at (i,j) with i != j.
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
        /// bounds.  If countAgainstTotal is true, it will also check
        /// to make sure you haven't added more than the expected
        /// number of matrix entries.  Regardless, this method will
        /// update the "actual" metadata.
        ///
        /// \param i [in] (1-based) row index
        /// \param j [in] (1-based) column index
        /// \param Aij [in] Value of the entry A(i,j)
        /// \param countAgainstTotal [in] Whether to count the entry
        ///   to insert against the total expected number of entries.
        ///   The default is true.  Make this false if you are
        ///   inserting an entry that wasn't stored in the original
        ///   Matrix Market file, which you're adding in order to
        ///   preserve symmetry or some other related structural
        ///   property of the matrix.
        void
        operator() (const Ordinal i,
                    const Ordinal j,
                    const Scalar& Aij,
                    const bool countAgainstTotal=true)
        {
          if (! tolerant_) {
            const bool indexPairOutOfRange = i < 1 || j < 1 ||
              i > expectedNumRows_ || j > expectedNumCols_;

            TEUCHOS_TEST_FOR_EXCEPTION(indexPairOutOfRange,
              std::invalid_argument, "Matrix is " << expectedNumRows_ << " x "
              << expectedNumCols_ << ", so entry A(" << i << "," << j << ") = "
              << Aij << " is out of range.");
            if (countAgainstTotal) {
              TEUCHOS_TEST_FOR_EXCEPTION(seenNumEntries_ >= expectedNumEntries_,
                std::invalid_argument, "Cannot add entry A(" << i << "," << j
                << ") = " << Aij << " to matrix; already have expected number "
                "of entries " << expectedNumEntries_ << ".");
            }
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
          if (countAgainstTotal) {
            ++seenNumEntries_;
          }
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
          const Ordinal nrows = tolerant_ ? seenNumRows_ : expectedNumRows_;
          ArrayRCP<Ordinal> ptr = arcp<Ordinal> (nrows + 1);

          // Copy over the elements, and fill in the ptr array with
          // offsets.  Note that merge() sorted the entries by row
          // index, so we can assume the row indices are increasing in
          // the list of entries.
          Ordinal curRow = 0;
          Ordinal curInd = 0;
          typedef typename std::vector<element_type>::const_iterator iter_type;
          ptr[0] = 0; // ptr always has at least one entry.
          for (iter_type it = elts_.begin(); it != elts_.end(); ++it) {
            const Ordinal i = it->rowIndex ();
            const Ordinal j = it->colIndex ();
            const Scalar Aij = it->value ();

            TEUCHOS_TEST_FOR_EXCEPTION(i < curRow, std::logic_error, "The "
              "current matrix entry's row index " << i << " is less then what "
              "should be the current row index lower bound " << curRow << ".");
            for (Ordinal k = curRow+1; k <= i; ++k) {
              ptr[k] = curInd;
            }
            curRow = i;

            TEUCHOS_TEST_FOR_EXCEPTION(
              static_cast<size_t> (curInd) >= elts_.size (),
              std::logic_error, "The current index " << curInd << " into ind "
              "and val is >= the number of matrix entries " << elts_.size ()
              << ".");
            ind[curInd] = j;
            val[curInd] = Aij;
            ++curInd;
          }
          for (Ordinal k = curRow+1; k <= nrows; ++k) {
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

        /// \brief Computed number of rows.
        ///
        /// "Computed" means "as seen from the matrix data."
        const Ordinal numRows() const { return seenNumRows_; }

        /// \brief Computed number of columns.
        ///
        /// "Computed" means "as seen from the matrix data."
        const Ordinal numCols() const { return seenNumCols_; }

        /// \brief Computed number of columns.
        ///
        /// "Computed" means "as seen from the matrix data."
        const Ordinal numEntries() const { return seenNumEntries_; }


      private:
        Ordinal expectedNumRows_, expectedNumCols_, expectedNumEntries_;
        Ordinal seenNumRows_, seenNumCols_, seenNumEntries_;
        bool tolerant_;
        bool debug_;

        //! The actual matrix entries, stored as an array of structs.
        std::vector<element_type> elts_;
      };
    } // namespace Raw
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // #ifndef __Teuchos_MatrixMarket_Raw_Adder_hpp
