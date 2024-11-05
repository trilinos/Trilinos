// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Graph_Adder_hpp
#define __Teuchos_MatrixMarket_Raw_Graph_Adder_hpp

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
      /// \class GraphElement
      /// \author Alicia Klinvex
      /// \brief Stores one entry of a sparse graph.
      ///
      /// \tparam Ordinal The type of indices of the sparse graph.
      ///
      /// This class is mainly useful as an implementation detail of
      /// GraphAdder.  We expose it to users only if they wish to convert
      /// the sparse graph read in by GraphAdder into a storage format
      /// other than CSR (compressed sparse row).
      ///
      /// An array of Elements implements the so-called "array of
      /// structs" representation of a coordinate format sparse
      /// graph.  A GraphElement has a row and column index (each of type
      /// Ordinal).  Elements also have
      /// equality and ordering comparisons.  The equality comparison
      /// only tests the row and column index, and is intended to
      /// simplify merging matrix entries with the same row and column
      /// indices.  The ordering comparison means that std::sort of a
      /// sequence of Elements will put them in an order suitable for
      /// extracting the CSR (compressed sparse row) representation of
      /// the sparse graph.
      template<class Ordinal>
      class GraphElement {
      public:
        //! Default constructor: an invalid entry of the graph.
        GraphElement () :
          rowIndex_ (Teuchos::OrdinalTraits<Ordinal>::invalid ()),
          colIndex_ (Teuchos::OrdinalTraits<Ordinal>::invalid ())
        {}

        //! Create a sparse graph entry at (i,j).
        GraphElement (const Ordinal i, const Ordinal j) :
          rowIndex_ (i), colIndex_ (j) {}

        //! Compare row and column indices.
        bool operator== (const GraphElement& rhs) {
          return rowIndex_ == rhs.rowIndex_ && colIndex_ == rhs.colIndex_;
        }

        //! Compare row and column indices.
        bool operator!= (const GraphElement& rhs) {
          return ! (*this == rhs);
        }

        //! Lexicographic order first by row index, then by column index.
        bool operator< (const GraphElement& rhs) const {
          if (rowIndex_ < rhs.rowIndex_)
            return true;
          else if (rowIndex_ > rhs.rowIndex_)
            return false;
          else { // equal
            return colIndex_ < rhs.colIndex_;
          }
        }

        //! Row index (zero-based) of this GraphElement.
        Ordinal rowIndex() const { return rowIndex_; }

        //! Column index (zero-based) of this GraphElement.
        Ordinal colIndex() const { return colIndex_; }

      private:
        Ordinal rowIndex_, colIndex_;
      };

      /// \brief Print out a GraphElement to the given output stream.
      ///
      /// This method is suitable for printing a sparse graph to a
      /// Matrix Market file.
      template<class Ordinal>
      std::ostream&
      operator<< (std::ostream& out, const GraphElement<Ordinal>& elt)
      {
        out << elt.rowIndex () << " " << elt.colIndex ();
        return out;
      }

      /// \class GraphAdder
      /// \brief To be used with Checker for "raw" sparse matrix input.
      ///
      /// \tparam Ordinal The type of indices in the sparse matrix.
      ///
      /// This class implements the following interface, which is
      /// required by the Callback template parameter of
      /// Teuchos::MatrixMarket::CoordPatternReader:
      /// \code
      /// class AdderType {
      /// public:
      ///   typedef ... index_type; // Ellipsis represents the actual type
      ///   void operator() (const index_type, const index_type, const value_type&);
      /// };
      /// \endcode
      /// For GraphAdder, the Ordinal template parameter is index_type.  GraphAdder
      /// provides a simple implementation of the above interface
      /// which is useful for things like printing out a sparse
      /// graph's entries, or converting between storage formats.
      template<class Ordinal>
      class GraphAdder {
      public:
        typedef Ordinal index_type;
        typedef GraphElement<Ordinal> element_type;
        typedef typename std::vector<element_type>::size_type size_type;

        /// \brief Default constructor.
        ///
        /// If you call the default constructor, we assume that you
        /// want tolerant mode (in which the GraphAdder tries to infer the
        /// graph dimensions and number of entries from the actual
        /// graph data, not from any metadata).  Tolerant mode is
        /// similar to what Matlab does if you give it an ASCII file
        /// of (i,j,Aij) triples.  It may get the graph dimensions
        /// (m,n) wrong if the lower right entry of the graph is zero
        /// and is not supplied explicitly by calling operator().
        GraphAdder () :
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
        /// \param expectedNumRows [in] Number of rows in the graph,
        ///   as specified by the matrix metadata.
        ///
        /// \param expectedNumCols [in] Number of columns in the
        ///   graph, as specified by the matrix metadata.
        ///
        /// \param expectedNumEntries [in] Number of entries in the
        ///   graph, as specified by the matrix metadata.
        ///
        /// \param tolerant [in] Whether the "expected" metadata is
        ///   required to match what the read-in graph entries tell
        ///   us.
        ///
        /// \param debug [in] If true, we may print copious status
        ///   output for debugging purposes.
        GraphAdder (const Ordinal expectedNumRows,
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

        /// \brief Add an entry to the sparse graph.
        ///
        /// If tolerant==false, this method will perform error
        /// checking to ensure that the graph data matches the
        /// metadata.  For example, it will check that i and j are in
        /// bounds.  If countAgainstTotal is true, it will also check
        /// to make sure you haven't added more than the expected
        /// number of graph entries.  Regardless, this method will
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
                    const bool countAgainstTotal=true)
        {
          if (! tolerant_) {
            const bool indexPairOutOfRange = i < 1 || j < 1 ||
              i > expectedNumRows_ || j > expectedNumCols_;

            TEUCHOS_TEST_FOR_EXCEPTION(indexPairOutOfRange,
              std::invalid_argument, "Graph is " << expectedNumRows_ << " x "
              << expectedNumCols_ << ", so entry A(" << i << "," << j
              << ") is out of range.");
            if (countAgainstTotal) {
              TEUCHOS_TEST_FOR_EXCEPTION(seenNumEntries_ >= expectedNumEntries_,
                std::invalid_argument, "Cannot add entry A(" << i << "," << j
                << ")  to graph; already have expected number "
                "of entries " << expectedNumEntries_ << ".");
            }
          }
          // i and j are 1-based indices, but we store them as 0-based.
          elts_.push_back (element_type (i-1, j-1));

          // Keep track of the rightmost column containing a matrix
          // entry, and the bottommost row containing a matrix entry.
          // This gives us a lower bound for the dimensions of the
          // graph, and a check for the reported dimensions of the
          // graph in the Matrix Market file.
          seenNumRows_ = std::max (seenNumRows_, i);
          seenNumCols_ = std::max (seenNumCols_, j);
          if (countAgainstTotal) {
            ++seenNumEntries_;
          }
        }

        /// \brief Print the sparse graph data.
        ///
        /// We always print the data sorted.  You may also merge
        /// duplicate entries if you prefer.
        ///
        /// \param out [out] Output stream to which to print
        ///
        /// \param doMerge [in] Whether to merge entries before printing
        ///
        /// \param replace [in] If merging, whether to replace
        ///   duplicate entries; otherwise their values are added
        ///   together.
        ///
        /// \warning It never makes sense for replace to be true.
        ///   Perhaps we should get rid of this argument at some
        ///   point.
        void
        print (std::ostream& out, const bool doMerge, const bool replace=false)
        {
          if (doMerge) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (replace, std::logic_error, "replace = true not implemented!");
            //merge (replace);
            merge ();
          } else {
            std::sort (elts_.begin(), elts_.end());
          }
          // Print out the results, delimited by newlines.
          typedef std::ostream_iterator<element_type> iter_type;
          std::copy (elts_.begin(), elts_.end(), iter_type (out, "\n"));
        }

        /// \brief Merge duplicate elements.
        ///
        /// Merge elements of the sparse graph that have the same row
        /// and column indices ("duplicates").  Resize the array of
        /// elements to fit just the "unique" (not duplicate)
        /// elements.
        ///
        /// \return (# unique elements, # removed elements)
        ///
        /// \note This method does not change the "expected" or "seen"
        ///   numbers of entries, since both of those count entries
        ///   with the same row and column indices as separate
        ///   entries.
        std::pair<size_type, size_type>
        merge ()
        {
          typedef typename std::vector<element_type>::iterator iter_type;

          // Start with a sorted container.  GraphElement objects sort in
          // lexicographic order of their (row, column) indices, for
          // easy conversion to CSR format.  If you expect that the
          // elements will usually be sorted in the desired order, you
          // can check first whether they are already sorted.  We have
          // no such expectation, so we don't even bother to spend the
          // extra O(# entries) operations to check.
          std::sort (elts_.begin(), elts_.end());

          // Remove duplicate elements from the sequence
          iter_type it;
          it = std::unique(elts_.begin(), elts_.end());
          size_type numUnique = std::distance(elts_.begin(),it);
          const size_type numRemoved = elts_.size() - numUnique;
          elts_.resize( std::distance(elts_.begin(),it) );
          elts_.resize (numUnique);
          return std::make_pair (numUnique, numRemoved);
        }

        /// \brief Merge duplicate elements and convert to zero-based CSR.
        ///
        /// Merge elements of the sparse graph that have the same row
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
        ///   numRows is the number of rows in the sparse graph.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        ///
        /// \param colind [out] Column indices of the graph.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        ///
        /// \note This method does not change the "expected" or "seen"
        ///   numbers of entries, since both of those count entries
        ///   with the same row and column indices as separate
        ///   entries.
        void
        mergeAndConvertToCSR (size_type& numUniqueElts,
                              size_type& numRemovedElts,
                              Teuchos::ArrayRCP<Ordinal>& rowptr,
                              Teuchos::ArrayRCP<Ordinal>& colind)
        {
          using Teuchos::arcp;
          using Teuchos::ArrayRCP;

          std::pair<size_type, size_type> mergeResult = merge();

          // At this point, elts_ is already in CSR order.
          // Now we can allocate and fill the ind array.
          ArrayRCP<Ordinal> ind = arcp<Ordinal> (elts_.size ());

          // Number of rows in the graph.
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

            TEUCHOS_TEST_FOR_EXCEPTION(i < curRow, std::logic_error, "The "
              "current graph entry's row index " << i << " is less then what "
              "should be the current row index lower bound " << curRow << ".");
            for (Ordinal k = curRow+1; k <= i; ++k) {
              ptr[k] = curInd;
            }
            curRow = i;

            TEUCHOS_TEST_FOR_EXCEPTION(
              static_cast<size_t> (curInd) >= elts_.size (),
              std::logic_error, "The current index " << curInd << " into ind "
              "is >= the number of matrix entries " << elts_.size ()
              << ".");
            ind[curInd] = j;
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
          numUniqueElts = mergeResult.first;
          numRemovedElts = mergeResult.second;
        }

        //! A temporary const view of the entries of the graph.
        const std::vector<element_type>& getEntries() const {
          return elts_;
        }

        //! Clear all the added graph entries and reset metadata.
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

#endif // #ifndef __Teuchos_MatrixMarket_Raw_Graph_Adder_hpp
