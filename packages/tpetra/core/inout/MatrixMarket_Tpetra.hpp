// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __MatrixMarket_Tpetra_hpp
#define __MatrixMarket_Tpetra_hpp

/// \file MatrixMarket_Tpetra.hpp
/// \brief Matrix Market file readers and writers for Tpetra objects.
/// \author Mark Hoemmen
///
/// This header file implements Matrix Market file readers and writers
/// for both sparse and dense matrices (as \c Tpetra::CrsMatrix
/// resp. \c Tpetra::MultiVector).  The Matrix Market (see their <a
/// href="http://math.nist.gov/MatrixMarket"> web site </a> for
/// details) defines a human-readable ASCII text file format ("Matrix
/// Market format") for interchange of sparse and dense matrices.
///
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_ComputeGatherMap.hpp"
#include "Teuchos_MatrixMarket_Raw_Adder.hpp"
#include "Teuchos_MatrixMarket_Raw_Graph_Adder.hpp"
#include "Teuchos_MatrixMarket_SymmetrizingAdder.hpp"
#include "Teuchos_MatrixMarket_SymmetrizingGraphAdder.hpp"
#include "Teuchos_MatrixMarket_assignScalar.hpp"
#include "Teuchos_MatrixMarket_Banner.hpp"
#include "Teuchos_MatrixMarket_CoordDataReader.hpp"
#include "Teuchos_SetScientific.hpp"
#include "Teuchos_TimeMonitor.hpp"

extern "C" {
#include "mmio_Tpetra.h"
}
#include "Tpetra_Distribution.hpp"


#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>
#include <numeric>

namespace Tpetra {
  /// \namespace MatrixMarket
  /// \brief Matrix Market file readers and writers for sparse and
  ///   dense matrices (as \c CrsMatrix resp. \c MultiVector).
  /// \author Mark Hoemmen
  ///
  /// The Matrix Market (see their <a
  /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
  /// details) defines a human-readable ASCII text file format
  /// ("Matrix Market format") for interchange of sparse and dense
  /// matrices.  This namespace defines classes for reading and
  /// writing sparse or dense matrices from a Matrix Market file or
  /// input stream.
  ///
  /// Matrix Market files are designed for easy reading and writing of
  /// test matrices by both humans and computers.  They are <i>not</i>
  /// intended for high-performance or parallel file input and output.
  /// You should use a true parallel file format if you want to do
  /// parallel input and output of sparse or dense matrices.
  ///
  /// Since the Matrix Market format is not optimized for performance
  /// or parallelism, some of our readers and writers assume that the
  /// entire matrix can fit in a single MPI process.  Nevertheless, we
  /// always do all of the file input or output on the MPI process
  /// with rank 0 in the given communicator ("Process 0").
  /// Distributed input matrices may be gathered from all MPI
  /// processes in the participating communicator, and distributed
  /// output matrices are broadcast from Process 0 to all MPI
  /// processes in the participating communicator.  In some cases, we
  /// have optimized this to save memory.
  namespace MatrixMarket {
    /// \class Reader
    /// \brief Matrix Market file reader for CrsMatrix and MultiVector.
    /// \author Mark Hoemmen
    ///
    /// The Matrix Market (see its <a
    /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
    /// details) defines a human-readable text file format for storing
    /// sparse and dense matrices.  This class defines methods for
    /// reading sparse and dense matrices from a Matrix Market file or
    /// input stream.  It represents sparse matrices as CrsMatrix and
    /// dense vectors and matrices as MultiVector.  Reader can also
    /// read a Map (in the format produced by Writer) from a file or
    /// input stream.
    ///
    /// All methods of this class only open files or read from input
    /// streams on MPI Process 0, with respect to the MPI communicator
    /// over which the given CrsMatrix or MultiVector is to be
    /// distributed.
    ///
    /// We define the MultiVector type returned by readDense() and
    /// readDenseFile() using the scalar_type, local_ordinal_type,
    /// global_ordinal_type, and node_type typedefs in
    /// SparseMatrixType.  This ensures that the multivectors returned
    /// by those methods have a type compatible with the CrsMatrix
    /// sparse matrices returned by readSparse() and readSparseFile().
    /// We do this because the typical use case of Matrix Market files
    /// in Trilinos is to test sparse matrix methods, which usually
    /// involves reading a sparse matrix A and perhaps also a dense
    /// right-hand side b.
    ///
    /// \tparam SparseMatrixType A specialization of CrsMatrix.
    ///
    /// Templating on the specialization of CrsMatrix means that the
    /// Reader expects matrix data of a type compatible with the
    /// CrsMatrix's scalar_type.  In general, Matrix Market files may
    /// contain data of integer, real, or complex type.  However, the
    /// reader methods have to return a CrsMatrix of a specific type,
    /// so we require that you declare a Reader with the CrsMatrix
    /// type that you want and that you expect the file(s) to contain.
    ///
    /// We didn't find any of the alternatives to this approach
    /// acceptable.  One possibility would have been to have the
    /// reader methods return a "container that can hold anything,"
    /// like a boost::any.  However, then you would have to know all
    /// five template arguments of the CrsMatrix in order to get the
    /// actual CrsMatrix object out.  C++ doesn't have algebraic data
    /// types (see the Wikipedia entry for a good definition) that are
    /// disjoint unions of different types.  Thus, we couldn't have
    /// had the readers return a CrsMatrix with scalar_type = "int or
    /// double or complex<double>."  While you can implement such a
    /// type in C++ (see e.g., boost::variant), it would not be
    /// interchangeable for its component types.  This is because it
    /// may not have the same memory layout (e.g., copying an array of
    /// boost::variant<int, double, complex<double> > bitwise into an
    /// array of int may not work).
    template<class SparseMatrixType>
    class Reader {
    public:
      //! This class' template parameter; a specialization of CrsMatrix.
      typedef SparseMatrixType sparse_matrix_type;
      typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;

      /// Type of the entries of the sparse matrix.
      /// The first template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// Type of the local indices of the sparse matrix.
      /// The second template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// Type of the global indices of the sparse matrix.
      ///
      /// The third template parameter of CrsMatrix and MultiVector.
      /// This is also the type of indices as read from the Matrix
      /// Market file.  Indices of the sparse matrix are read in as
      /// global ordinals, since Matrix Market files represent the
      /// whole matrix and don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type
        global_ordinal_type;
      //! The fourth template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::node_type node_type;

      //! The CrsGraph specialization associated with SparseMatrixType.
      typedef CrsGraph<local_ordinal_type,
                       global_ordinal_type,
                       node_type> sparse_graph_type;

      //! The MultiVector specialization associated with SparseMatrixType.
      typedef MultiVector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> multivector_type;

      //! The Vector specialization associated with SparseMatrixType.
      typedef Vector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> vector_type;

      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

      //! Type of the MPI communicator.
      using trcp_tcomm_t = Teuchos::RCP<const Teuchos::Comm<int>>;

    private:
      /// \typedef size_type
      /// \brief Handy typedef for entries of arrays such as rowPtr.
      ///
      /// This assumes that the ArrayRCP<T>::size_type typedef is the
      /// same for all T, which we know to be true.
      typedef Teuchos::ArrayRCP<int>::size_type size_type;

      /// \brief Compute initial range map.
      ///
      /// Range maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pComm [in] Global communicator.
      /// \param numRows [in] Global number of rows in the matrix.
      ///
      /// \return Range map to be used for constructing a CrsMatrix.
      static Teuchos::RCP<const map_type>
      makeRangeMap (const trcp_tcomm_t& pComm,
                    const global_ordinal_type numRows)
      {
        // Return a conventional, uniformly partitioned, contiguous map.
        return rcp (new map_type (static_cast<global_size_t> (numRows),
                                  static_cast<global_ordinal_type> (0),
                                  pComm, GloballyDistributed));
      }

      /// \brief Compute initial row map, or verify an existing one.
      ///
      /// The typical case when reading a sparse matrix from a file is
      /// for the reader itself to create a new row map, in particular
      /// a standard uniform contiguous one-to-one row map.  However,
      /// we also give the option to use an existing row map, if you
      /// are already using a particular distribution for (say) vector
      /// data and don't want to stop using it.  In the latter case
      /// (pRowMap is not null), we validate the communicator and node
      /// of the existing row map that you pass in.  In either case,
      /// you need to know the (global) number of rows in the matrix.
      ///
      /// \param pRowMap [in] If non-null, test pRowMap for validity,
      ///   and return it if valid.  Otherwise, if pRowMap is null,
      ///   initialize and return a (uniform contiguous one-to-one)
      ///   row map.  "Validity" here means that the map's
      ///   communicator and node are the same objects (pointerwise)
      ///   as the corresponding arguments.  (Note that the global
      ///   number of elements may not be the same as the number of
      ///   rows; a row map is not required to be one-to-one.)  The
      ///   typical case is to pass in null here, which is why we call
      ///   this routine "makeRowMap".
      /// \param pComm [in] Global communicator.
      /// \param numRows [in] Global number of rows in the matrix.  If
      ///   pRowMap is nonnull, used only for error checking.
      ///
      /// \return If pRowMap is null, a new row map, otherwise pRowMap.
      static Teuchos::RCP<const map_type>
      makeRowMap (const Teuchos::RCP<const map_type>& pRowMap,
                  const trcp_tcomm_t& pComm,
                  const global_ordinal_type numRows)
      {
        // If the caller didn't provide a map, return a conventional,
        // uniformly partitioned, contiguous map.
        if (pRowMap.is_null ()) {
          return rcp (new map_type (static_cast<global_size_t> (numRows),
                                    static_cast<global_ordinal_type> (0),
                                    pComm, GloballyDistributed));
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION
            (! pRowMap->isDistributed () && pComm->getSize () > 1,
             std::invalid_argument, "The specified row map is not distributed, "
             "but the given communicator includes more than one process (in "
             "fact, there are " << pComm->getSize () << " processes).");
          TEUCHOS_TEST_FOR_EXCEPTION
            (pRowMap->getComm () != pComm, std::invalid_argument,
             "The specified row Map's communicator (pRowMap->getComm()) "
             "differs from the given separately supplied communicator pComm.");
          return pRowMap;
        }
      }

      /// \brief Compute domain map.
      ///
      /// Domain maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pRangeMap [in] Valid range map of the matrix,
      ///   as returned by \c makeRangeMap().
      /// \param numRows [in] Global number of rows in the matrix.
      /// \param numCols [in] Global number of columns in the matrix.
      ///
      /// \return The domain map.  If numRows == numCols, this is
      ///   identical to the range map, otherwise we make a new map
      ///   for the domain.
      static Teuchos::RCP<const map_type>
      makeDomainMap (const Teuchos::RCP<const map_type>& pRangeMap,
                     const global_ordinal_type numRows,
                     const global_ordinal_type numCols)
      {
        // Abbreviations so that the map creation call isn't too long.
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;

        if (numRows == numCols) {
          return pRangeMap;
        } else {
          return createUniformContigMapWithNode<LO,GO,NT> (numCols,
                                                           pRangeMap->getComm ());
        }
      }

      /// \brief Distribute the sparse matrix entries.
      ///
      /// This is one of those routines that just has to be messy.  We try
      /// to comment everything so you can see what's going on.  All of the
      /// matrix data starts out on Rank 0, stored in CSR format (rowPtr,
      /// colInd, values).  numEntriesPerRow applies to this data; it can be
      /// computed from rowPtr, and has one less entry than rowPtr.  Rank 0
      /// then engages in a dialog with Ranks 1 .. numProcs-1 to determine
      /// which part of the matrix data is theirs, and to send them their
      /// data.
      ///
      /// \param myNumEntriesPerRow [out] For my row indices, the number of
      ///   entries per row.  This array has
      ///   pRowMap->getLocalElementList().size() entries, and is indexed in
      ///   the same way, so that myNumEntriesPerRow[k] is the number of
      ///   entries in row pRowMap->getLocalElementList()[k].
      ///
      /// \param myRowPtr [out] The row pointer array for the rows
      ///   that belong to me.  This array has one more entry than
      ///   myNumEntriesPerRow, and is the prefix sum (with
      ///   myRowPtr[0] = 0) of myNumEntriesPerRow.
      ///
      /// \param myColInd [out] My rows' column indices.  If myRows =
      ///   pRowMap->getLocalElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then
      ///   myColInd[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param myValues [out] My rows' stored matrix values.  If myRows =
      ///   pRowMap->getLocalElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then
      ///   myValues[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param pRowMap [in] Map describing the distribution of rows among
      ///   processors.
      ///
      /// \param numEntriesPerRow [in/out] For all row indices, the number
      ///   of entries per row.  You can construct this from the usual CSR
      ///   matrix "rowPtr" array by differences: numEntriesPerRow[k] =
      ///   rowPtr[k+1] - rowPtr[k].  As a result, numEntriesPerRow has one
      ///   less entry than rowPtr.  On output, the reference is invalidated
      ///   to save space.
      ///
      /// \param rowPtr [in/out] On input: the usual CSR matrix row pointer
      ///   array for the whole matrix.  On output, the reference is
      ///   invalidated to save space.
      ///
      /// \param colInd [in/out] On input: all rows' column indices.  If
      ///   start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then colInd[start .. end-1] are the column
      ///   indices for row k.  On output, the reference is invalidated to
      ///   save space.
      ///
      /// \param values [in/out] On input: all rows' stored matrix values.
      ///   If start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then values[start .. end-1] are the values
      ///   for row k.  On output, the reference is invalidated to save
      ///   space.
      ///
      /// \param debug [in] If true, print copious debugging output to
      ///   stderr on Rank 0.  This option is unlikely to be useful to
      ///   anyone but a Tpetra developer debugging this code.
      ///
      /// \note It was only after I wrote this routine that I learned
      ///   it was completely unnecessary.  All the data
      ///   redistribution can be done in about 10 lines of code by
      ///   using Tpetra::Map objects, and either Import or Export.
      ///   (For example, you could read the file into the sparse
      ///   matrix entirely owned by Proc 0, then export it using a
      ///   distributed Map.)  However, this routine works and I
      ///   haven't had time to rewrite it yet.  Just expect that the
      ///   implementation of \c readSparse() may become a lot shorter
      ///   in the future.
      static void
      distribute (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& pRowMap,
                  Teuchos::ArrayRCP<size_t>& numEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& rowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& colInd,
                  Teuchos::ArrayRCP<scalar_type>& values,
                  const bool debug=false)
      {
         using Teuchos::arcp;
         using Teuchos::ArrayRCP;
         using Teuchos::ArrayView;
         using Teuchos::as;
         using Teuchos::Comm;
         using Teuchos::CommRequest;
         using Teuchos::null;
         using Teuchos::RCP;
         using Teuchos::receive;
         using Teuchos::send;
         using std::cerr;
         using std::endl;

         const bool extraDebug = false;
         trcp_tcomm_t pComm = pRowMap->getComm ();
         const int numProcs = pComm->getSize ();
         const int myRank = pComm->getRank ();
         const int rootRank = 0;

         // Type abbreviations to make the code more concise.
         typedef global_ordinal_type GO;

         // List of the global indices of my rows.  They may or may
         // not be contiguous, and the row map need not be one-to-one.
         ArrayView<const GO> myRows = pRowMap->getLocalElementList();
         const size_type myNumRows = myRows.size();
         TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(myNumRows) !=
                            pRowMap->getLocalNumElements(),
                            std::logic_error,
                            "pRowMap->getLocalElementList().size() = "
                            << myNumRows
                            << " != pRowMap->getLocalNumElements() = "
                            << pRowMap->getLocalNumElements() << ".  "
                            "Please report this bug to the Tpetra developers.");
         TEUCHOS_TEST_FOR_EXCEPTION(myRank == 0 && numEntriesPerRow.size() < myNumRows,
                            std::logic_error,
                            "On Proc 0: numEntriesPerRow.size() = "
                            << numEntriesPerRow.size()
                            << " != pRowMap->getLocalElementList().size() = "
                            << myNumRows << ".  Please report this bug to the "
                            "Tpetra developers.");

         // Space for my proc's number of entries per row.  Will be
         // filled in below.
         myNumEntriesPerRow = arcp<size_t> (myNumRows);

         if (myRank != rootRank) {
           // Tell the root how many rows we have.  If we're sending
           // none, then we don't have anything else to send, nor does
           // the root have to receive anything else.
           send (*pComm, myNumRows, rootRank);
           if (myNumRows != 0) {
             // Now send my rows' global indices.  Hopefully the cast
             // to int doesn't overflow.  This is unlikely, since it
             // should fit in a LO, even though it is a GO.
             send (*pComm, static_cast<int> (myNumRows),
                   myRows.getRawPtr(), rootRank);

             // I (this proc) don't care if my global row indices are
             // contiguous, though the root proc does (since otherwise
             // it needs to pack noncontiguous data into contiguous
             // storage before sending).  That's why we don't check
             // for contiguousness here.

             // Ask the root process for my part of the array of the
             // number of entries per row.
             receive (*pComm, rootRank,
                      static_cast<int> (myNumRows),
                      myNumEntriesPerRow.getRawPtr());

             // Use the resulting array to figure out how many column
             // indices and values I should ask from the root process.
             const local_ordinal_type myNumEntries =
               std::accumulate (myNumEntriesPerRow.begin(),
                                myNumEntriesPerRow.end(), 0);

             // Make space for my entries of the sparse matrix.  Note
             // that they don't have to be sorted by row index.
             // Iterating through all my rows requires computing a
             // running sum over myNumEntriesPerRow.
             myColInd = arcp<GO> (myNumEntries);
             myValues = arcp<scalar_type> (myNumEntries);
             if (myNumEntries > 0) {
               // Ask for that many column indices and values, if
               // there are any.
               receive (*pComm, rootRank,
                        static_cast<int> (myNumEntries),
                        myColInd.getRawPtr());
               receive (*pComm, rootRank,
                        static_cast<int> (myNumEntries),
                        myValues.getRawPtr());
             }
           } // If I own at least one row
         } // If I am not the root processor
         else { // I _am_ the root processor
           if (debug) {
             cerr << "-- Proc 0: Copying my data from global arrays" << endl;
           }
           // Proc 0 still needs to (allocate, if not done already)
           // and fill its part of the matrix (my*).
           for (size_type k = 0; k < myNumRows; ++k) {
             const GO myCurRow = myRows[k];
             const local_ordinal_type numEntriesInThisRow = numEntriesPerRow[myCurRow];
             myNumEntriesPerRow[k] = numEntriesInThisRow;
           }
           if (extraDebug && debug) {
             cerr << "Proc " << pRowMap->getComm ()->getRank ()
                  << ": myNumEntriesPerRow[0.." << (myNumRows-1) << "] = [";
             for (size_type k = 0; k < myNumRows; ++k) {
               cerr << myNumEntriesPerRow[k];
               if (k < myNumRows-1) {
                 cerr << " ";
               }
             }
             cerr << "]" << endl;
           }
           // The total number of matrix entries that my proc owns.
           const local_ordinal_type myNumEntries =
             std::accumulate (myNumEntriesPerRow.begin(),
                              myNumEntriesPerRow.end(), 0);
           if (debug) {
             cerr << "-- Proc 0: I own " << myNumRows << " rows and "
                  << myNumEntries << " entries" << endl;
           }
           myColInd = arcp<GO> (myNumEntries);
           myValues = arcp<scalar_type> (myNumEntries);

           // Copy Proc 0's part of the matrix into the my* arrays.
           // It's important that myCurPos be updated _before_ k,
           // otherwise myCurPos will get the wrong number of entries
           // per row (it should be for the row in the just-completed
           // iteration, not for the next iteration's row).
           local_ordinal_type myCurPos = 0;
           for (size_type k = 0; k < myNumRows;
                myCurPos += myNumEntriesPerRow[k], ++k) {
             const local_ordinal_type curNumEntries = myNumEntriesPerRow[k];
             const GO myRow = myRows[k];
             const size_t curPos = rowPtr[myRow];
             // Only copy if there are entries to copy, in order not
             // to construct empty ranges for the ArrayRCP views.
             if (curNumEntries > 0) {
               ArrayView<GO> colIndView = colInd (curPos, curNumEntries);
               ArrayView<GO> myColIndView = myColInd (myCurPos, curNumEntries);
               std::copy (colIndView.begin(), colIndView.end(),
                          myColIndView.begin());

               ArrayView<scalar_type> valuesView =
                 values (curPos, curNumEntries);
               ArrayView<scalar_type> myValuesView =
                 myValues (myCurPos, curNumEntries);
               std::copy (valuesView.begin(), valuesView.end(),
                          myValuesView.begin());
             }
           }

           // Proc 0 processes each other proc p in turn.
           for (int p = 1; p < numProcs; ++p) {
             if (debug) {
               cerr << "-- Proc 0: Processing proc " << p << endl;
             }

             size_type theirNumRows = 0;
             // Ask Proc p how many rows it has.  If it doesn't
             // have any, we can move on to the next proc.  This
             // has to be a standard receive so that we can avoid
             // the degenerate case of sending zero data.
             receive (*pComm, p, &theirNumRows);
             if (debug) {
               cerr << "-- Proc 0: Proc " << p << " owns "
                    << theirNumRows << " rows" << endl;
             }
             if (theirNumRows != 0) {
               // Ask Proc p which rows it owns.  The resulting global
               // row indices are not guaranteed to be contiguous or
               // sorted.  Global row indices are themselves indices
               // into the numEntriesPerRow array.
               ArrayRCP<GO> theirRows = arcp<GO> (theirNumRows);
               receive (*pComm, p, as<int> (theirNumRows),
                        theirRows.getRawPtr ());
               // Extra test to make sure that the rows we received
               // are all sensible.  This is a good idea since we are
               // going to use the global row indices we've received
               // to index into the numEntriesPerRow array.  Better to
               // catch any bugs here and print a sensible error
               // message, rather than segfault and print a cryptic
               // error message.
               {
                 const global_size_t numRows = pRowMap->getGlobalNumElements ();
                 const GO indexBase = pRowMap->getIndexBase ();
                 bool theirRowsValid = true;
                 for (size_type k = 0; k < theirNumRows; ++k) {
                   if (theirRows[k] < indexBase ||
                       as<global_size_t> (theirRows[k] - indexBase) >= numRows) {
                     theirRowsValid = false;
                   }
                 }
                 if (! theirRowsValid) {
                   TEUCHOS_TEST_FOR_EXCEPTION(
                     ! theirRowsValid, std::logic_error,
                     "Proc " << p << " has at least one invalid row index.  "
                     "Here are all of them: " <<
                     Teuchos::toString (theirRows ()) << ".  Valid row index "
                     "range (zero-based): [0, " << (numRows - 1) << "].");
                 }
               }

               // Perhaps we could save a little work if we check
               // whether Proc p's row indices are contiguous.  That
               // would make lookups in the global data arrays
               // faster.  For now, we just implement the general
               // case and don't prematurely optimize.  (Remember
               // that you're making Proc 0 read the whole file, so
               // you've already lost scalability.)

               // Compute the number of entries in each of Proc p's
               // rows.  (Proc p will compute its row pointer array
               // on its own, after it gets the data from Proc 0.)
               ArrayRCP<size_t> theirNumEntriesPerRow;
               theirNumEntriesPerRow = arcp<size_t> (theirNumRows);
               for (size_type k = 0; k < theirNumRows; ++k) {
                 theirNumEntriesPerRow[k] = numEntriesPerRow[theirRows[k]];
               }

               // Tell Proc p the number of entries in each of its
               // rows.  Hopefully the cast to int doesn't overflow.
               // This is unlikely, since it should fit in a LO,
               // even though it is a GO.
               send (*pComm, static_cast<int> (theirNumRows),
                     theirNumEntriesPerRow.getRawPtr(), p);

               // Figure out how many entries Proc p owns.
               const local_ordinal_type theirNumEntries =
                 std::accumulate (theirNumEntriesPerRow.begin(),
                                  theirNumEntriesPerRow.end(), 0);

               if (debug) {
                 cerr << "-- Proc 0: Proc " << p << " owns "
                      << theirNumEntries << " entries" << endl;
               }

               // If there are no entries to send, then we're done
               // with Proc p.
               if (theirNumEntries == 0) {
                 continue;
               }

               // Construct (views of) proc p's column indices and
               // values.  Later, we might like to optimize for the
               // (common) contiguous case, for which we don't need to
               // copy data into separate "their*" arrays (we can just
               // use contiguous views of the global arrays).
               ArrayRCP<GO> theirColInd (theirNumEntries);
               ArrayRCP<scalar_type> theirValues (theirNumEntries);
               // Copy Proc p's part of the matrix into the their*
               // arrays.  It's important that theirCurPos be updated
               // _before_ k, otherwise theirCurPos will get the wrong
               // number of entries per row (it should be for the row
               // in the just-completed iteration, not for the next
               // iteration's row).
               local_ordinal_type theirCurPos = 0;
               for (size_type k = 0; k < theirNumRows;
                    theirCurPos += theirNumEntriesPerRow[k], k++) {
                 const local_ordinal_type curNumEntries = theirNumEntriesPerRow[k];
                 const GO theirRow = theirRows[k];
                 const local_ordinal_type curPos = rowPtr[theirRow];

                 // Only copy if there are entries to copy, in order
                 // not to construct empty ranges for the ArrayRCP
                 // views.
                 if (curNumEntries > 0) {
                   ArrayView<GO> colIndView =
                     colInd (curPos, curNumEntries);
                   ArrayView<GO> theirColIndView =
                     theirColInd (theirCurPos, curNumEntries);
                   std::copy (colIndView.begin(), colIndView.end(),
                              theirColIndView.begin());

                   ArrayView<scalar_type> valuesView =
                     values (curPos, curNumEntries);
                   ArrayView<scalar_type> theirValuesView =
                     theirValues (theirCurPos, curNumEntries);
                   std::copy (valuesView.begin(), valuesView.end(),
                              theirValuesView.begin());
                 }
               }
               // Send Proc p its column indices and values.
               // Hopefully the cast to int doesn't overflow.  This
               // is unlikely, since it should fit in a LO, even
               // though it is a GO.
               send (*pComm, static_cast<int> (theirNumEntries),
                     theirColInd.getRawPtr(), p);
               send (*pComm, static_cast<int> (theirNumEntries),
                     theirValues.getRawPtr(), p);

               if (debug) {
                 cerr << "-- Proc 0: Finished with proc " << p << endl;
               }
             } // If proc p owns at least one row
           } // For each proc p not the root proc 0
         } // If I'm (not) the root proc 0

         // Invalidate the input data to save space, since we don't
         // need it anymore.
         numEntriesPerRow = null;
         rowPtr = null;
         colInd = null;
         values = null;

         if (debug && myRank == 0) {
           cerr << "-- Proc 0: About to fill in myRowPtr" << endl;
         }

         // Allocate and fill in myRowPtr (the row pointer array for
         // my rank's rows).  We delay this until the end because we
         // don't need it to compute anything else in distribute().
         // Each proc can do this work for itself, since it only needs
         // myNumEntriesPerRow to do so.
         myRowPtr = arcp<size_t> (myNumRows+1);
         myRowPtr[0] = 0;
         for (size_type k = 1; k < myNumRows+1; ++k) {
           myRowPtr[k] = myRowPtr[k-1] + myNumEntriesPerRow[k-1];
         }
         if (extraDebug && debug) {
           cerr << "Proc " << Teuchos::rank (*(pRowMap->getComm()))
                << ": myRowPtr[0.." << myNumRows << "] = [";
           for (size_type k = 0; k < myNumRows+1; ++k) {
             cerr << myRowPtr[k];
             if (k < myNumRows) {
               cerr << " ";
             }
           }
           cerr << "]" << endl << endl;
         }

         if (debug && myRank == 0) {
           cerr << "-- Proc 0: Done with distribute" << endl;
         }
      }

      /// \brief Given my proc's data, return the completed sparse matrix.
      ///
      /// Each proc inserts its data into the sparse matrix, and then,
      /// if callFillComplete is true, all procs call fillComplete().
      /// (For whatever reason, you might not be done with the matrix
      /// yet, so you might want to call fillComplete() yourself.
      /// CrsMatrix::fillResume() doesn't currently work as you might
      /// expect when storage optimization is enabled; it fixes the
      /// graph of the matrix, so that you can't add new entries.)
      ///
      /// Column indices are zero-based on input.  This method will
      /// change them in place to match the index base of the input
      /// row Map (\c pRowMap).
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& pRowMap,
                  const Teuchos::RCP<const map_type>& pRangeMap,
                  const Teuchos::RCP<const map_type>& pDomainMap,
                  const bool callFillComplete = true)
      {
        using Teuchos::ArrayView;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::cerr;
        using std::endl;
        // Typedef to make certain type declarations shorter.
        typedef global_ordinal_type GO;

        // The row pointer array always has at least one entry, even
        // if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
        // and myValues would all be empty arrays in that degenerate
        // case, but the row and domain maps would still be nonnull
        // (though they would be trivial maps).
        TEUCHOS_TEST_FOR_EXCEPTION(myRowPtr.is_null(), std::logic_error,
          "makeMatrix: myRowPtr array is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pDomainMap.is_null(), std::logic_error,
          "makeMatrix: domain map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRangeMap.is_null(), std::logic_error,
          "makeMatrix: range map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRowMap.is_null(), std::logic_error,
          "makeMatrix: row map is null.  "
          "Please report this bug to the Tpetra developers.");

        // Construct the CrsMatrix, using the row map, with the
        // constructor specifying the number of nonzeros for each row.
        RCP<sparse_matrix_type> A =
          rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow ()));

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = pRowMap->getLocalElementList ();
        const size_type myNumRows = myRows.size ();

        // Add this processor's matrix entries to the CrsMatrix.
        const GO indexBase = pRowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const local_ordinal_type curNumEntries = myNumEntriesPerRow[i];
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          // Avoid constructing empty views of ArrayRCP objects.
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        if (callFillComplete) {
          A->fillComplete (pDomainMap, pRangeMap);
        }
        return A;
      }

      /// \brief Variant of makeMatrix() that takes parameters for
      ///   CrsMatrix's constructor and for fillComplete().
      ///
      /// Each process inserts its data into the sparse matrix, and
      /// then all processes call fillComplete().
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& pRowMap,
                  const Teuchos::RCP<const map_type>& pRangeMap,
                  const Teuchos::RCP<const map_type>& pDomainMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                  const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams)
      {
        using Teuchos::ArrayView;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::cerr;
        using std::endl;
        // Typedef to make certain type declarations shorter.
        typedef global_ordinal_type GO;

        // The row pointer array always has at least one entry, even
        // if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
        // and myValues would all be empty arrays in that degenerate
        // case, but the row and domain maps would still be nonnull
        // (though they would be trivial maps).
        TEUCHOS_TEST_FOR_EXCEPTION(
          myRowPtr.is_null(), std::logic_error,
          "makeMatrix: myRowPtr array is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pDomainMap.is_null(), std::logic_error,
          "makeMatrix: domain map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pRangeMap.is_null(), std::logic_error,
          "makeMatrix: range map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pRowMap.is_null(), std::logic_error,
          "makeMatrix: row map is null.  "
          "Please report this bug to the Tpetra developers.");

        // Construct the CrsMatrix, using the row map, with the
        // constructor specifying the number of nonzeros for each row.
        RCP<sparse_matrix_type> A =
          rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow(),
                                       constructorParams));

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = pRowMap->getLocalElementList();
        const size_type myNumRows = myRows.size();

        // Add this processor's matrix entries to the CrsMatrix.
        const GO indexBase = pRowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const local_ordinal_type curNumEntries = myNumEntriesPerRow[i];
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        A->fillComplete (pDomainMap, pRangeMap, fillCompleteParams);
        return A;
      }

      /// \brief Variant of makeMatrix() that takes an optional column Map.
      ///
      /// This method computes \c colMap only if it is null on input,
      /// and if \c callFillComplete is true.
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& rowMap,
                  Teuchos::RCP<const map_type>& colMap,
                  const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const bool callFillComplete = true)
      {
        using Teuchos::ArrayView;
        using Teuchos::as;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcp;
        typedef global_ordinal_type GO;

        // Construct the CrsMatrix.

        RCP<sparse_matrix_type> A; // the matrix to return.
        if (colMap.is_null ()) { // the user didn't provide a column Map
          A = rcp (new sparse_matrix_type (rowMap, myNumEntriesPerRow));
        } else { // the user provided a column Map
          A = rcp (new sparse_matrix_type (rowMap, colMap, myNumEntriesPerRow));
        }

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = rowMap->getLocalElementList ();
        const size_type myNumRows = myRows.size ();

        // Add this process' matrix entries to the CrsMatrix.
        const GO indexBase = rowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const size_type curNumEntries = as<size_type> (myNumEntriesPerRow[i]);
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        if (callFillComplete) {
          A->fillComplete (domainMap, rangeMap);
          if (colMap.is_null ()) {
            colMap = A->getColMap ();
          }
        }
        return A;
      }

    private:

      /// \brief Read in the Banner line from the given input stream.
      ///
      /// Only call this method on one (MPI communicator) process.
      ///
      /// \param in [in/out] Input stream from which to read the
      ///   Banner line.
      /// \param lineNumber [in/out] On input: Current line number of
      ///   the input stream.  On output: if any line(s) were
      ///   successfully read from the input stream, this is
      ///   incremented by the number of line(s) read.  (This includes
      ///   comment lines.)
      /// \param tolerant [in] Whether to parse tolerantly.
      /// \param debug [in] Whether to write debugging output to
      ///   stderr.
      ///
      /// \return Banner [non-null]
      static Teuchos::RCP<const Teuchos::MatrixMarket::Banner>
      readBanner (std::istream& in,
                  size_t& lineNumber,
                  const bool tolerant=false,
                  const bool /* debug */=false,
                  const bool isGraph=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        RCP<Banner> pBanner; // On output, if successful: the read-in Banner.
        std::string line; // If read from stream successful: the Banner line

        // Try to read a line from the input stream.
        const bool readFailed = ! getline(in, line);        
        TEUCHOS_TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
                                   "Failed to get Matrix Market banner line from input.");

        // We read a line from the input stream.
        lineNumber++;

        // Assume that the line we found is the Banner line.
        try {
          pBanner = rcp (new Banner (line, tolerant));
        } catch (std::exception& e) {

          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                     "Matrix Market banner line contains syntax error(s): "
                                     << e.what());
        }

        TEUCHOS_TEST_FOR_EXCEPTION(pBanner->objectType() != "matrix",
          std::invalid_argument, "The Matrix Market file does not contain "
          "matrix data.  Its Banner (first) line says that its object type is \""
          << pBanner->matrixType() << "\", rather than the required \"matrix\".");

        // Validate the data type of the matrix, with respect to the
        // Scalar type of the CrsMatrix entries.
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! STS::isComplex && pBanner->dataType() == "complex",
          std::invalid_argument,
          "The Matrix Market file contains complex-valued data, but you are "
          "trying to read it into a matrix containing entries of the real-"
          "valued Scalar type \""
          << Teuchos::TypeNameTraits<scalar_type>::name() << "\".");
        TEUCHOS_TEST_FOR_EXCEPTION(
          !isGraph &&
          pBanner->dataType() != "real" &&
          pBanner->dataType() != "complex" &&
          pBanner->dataType() != "integer",
          std::invalid_argument,
          "When reading Matrix Market data into a Tpetra::CrsMatrix, the "
          "Matrix Market file may not contain a \"pattern\" matrix.  A "
          "pattern matrix is really just a graph with no weights.  It "
          "should be stored in a CrsGraph, not a CrsMatrix.");

        TEUCHOS_TEST_FOR_EXCEPTION(
          isGraph &&
          pBanner->dataType() != "pattern",
          std::invalid_argument,
          "When reading Matrix Market data into a Tpetra::CrsGraph, the "
          "Matrix Market file must contain a \"pattern\" matrix.");

        return pBanner;
      }

      /// \brief Read sparse matrix dimensions on Rank 0, and broadcast.
      ///
      /// Call on all MPI ranks.  MPI Rank 0 attempts to read in the
      /// coordinate dimensions from the input stream.  If it
      /// succeeds, it broadcasts them to all the other MPI ranks.
      /// (All ranks need to know the matrix dimensions in order to
      /// create domain, range, and column Maps.)
      ///
      /// \param in [in/out, valid only on Rank 0] Input stream from
      ///   which to read the sparse matrix dimensions.
      /// \param lineNumber [in/out, valid only on Rank 0] On input:
      ///   Current line number of the input stream.  On output: if
      ///   any line(s) were successfully read from the input stream,
      ///   this is incremented by the number of line(s) read.  (This
      ///   includes comment lines.)
      /// \param pComm [in, global] Communicator over which the matrix
      ///   will (eventually -- not here) be distributed.
      /// \param tolerant [in] Whether to parse tolerantly.
      /// \param debug [in] Whether to write debugging output to
      ///   stderr on MPI Proc 0.
      ///
      /// \return (numRows, numCols, numNonzeros)
      static Teuchos::Tuple<global_ordinal_type, 3>
      readCoordDims (std::istream& in,
                     size_t& lineNumber,
                     const Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                     const trcp_tcomm_t& pComm,
                     const bool tolerant = false,
                     const bool /* debug */ = false)
      {
        using Teuchos::MatrixMarket::readCoordinateDimensions;
        using Teuchos::Tuple;

        // Packed coordinate matrix dimensions (numRows, numCols,
        // numNonzeros); computed on Rank 0 and broadcasted to all MPI
        // ranks.
        Tuple<global_ordinal_type, 3> dims;

        // Read in the coordinate matrix dimensions from the input
        // stream.  "success" tells us whether reading in the
        // coordinate matrix dimensions succeeded ("Guilty unless
        // proven innocent").
        bool success = false;
        if (pComm->getRank() == 0) {
          TEUCHOS_TEST_FOR_EXCEPTION(pBanner->matrixType() != "coordinate",
            std::invalid_argument, "The Tpetra::CrsMatrix Matrix Market reader "
            "only accepts \"coordinate\" (sparse) matrix data.");
          // Unpacked coordinate matrix dimensions
          global_ordinal_type numRows, numCols, numNonzeros;
          // Only MPI Rank 0 reads from the input stream
          success = readCoordinateDimensions (in, numRows, numCols,
                                              numNonzeros, lineNumber,
                                              tolerant);
          // Pack up the data into a Tuple so we can send them with
          // one broadcast instead of three.
          dims[0] = numRows;
          dims[1] = numCols;
          dims[2] = numNonzeros;
        }
        // Only Rank 0 did the reading, so it decides success.
        //
        // FIXME (mfh 02 Feb 2011) Teuchos::broadcast doesn't know how
        // to send bools.  For now, we convert to/from int instead,
        // using the usual "true is 1, false is 0" encoding.
        {
          int the_success = success ? 1 : 0; // only matters on MPI Rank 0
          Teuchos::broadcast (*pComm, 0, &the_success);
          success = (the_success == 1);
        }
        if (success) {
          // Broadcast (numRows, numCols, numNonzeros) from Rank 0
          // to all the other MPI ranks.
          Teuchos::broadcast (*pComm, 0, dims);
        }
        else {
          // Perhaps in tolerant mode, we could set all the
          // dimensions to zero for now, and deduce correct
          // dimensions by reading all of the file's entries and
          // computing the max(row index) and max(column index).
          // However, for now we just error out in that case.
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            "Error reading Matrix Market sparse matrix: failed to read "
            "coordinate matrix dimensions.");
        }
        return dims;
      }

      /// \typedef adder_type
      /// \brief Type of object that "adds" entries to the sparse matrix.
      ///
      /// "Adds" here means that it collects and makes note of matrix
      /// entries read in from the input stream.  This object doesn't
      /// call insertGlobalEntries() or fillComplete() on the
      /// CrsMatrix.  Depending on the Matrix Market banner line
      /// information, it may "symmetrize" the matrix by adding entry
      /// A(j,i) (with the appropriate value depending on the symmetry
      /// type) if entry A(i,j) is seen.
      typedef Teuchos::MatrixMarket::SymmetrizingAdder<Teuchos::MatrixMarket::Raw::Adder<scalar_type, global_ordinal_type> > adder_type;

      typedef Teuchos::MatrixMarket::SymmetrizingGraphAdder<Teuchos::MatrixMarket::Raw::GraphAdder<global_ordinal_type> > graph_adder_type;

      /// \brief Make an "adder" object for processing matrix data.
      ///
      /// \param pComm [in] Communicator (across whose MPI ranks
      ///   the sparse matrix will be distributed)
      ///
      /// \param banner [in, nonnull and valid on Rank 0 only]
      ///   Object describing the type and symmetry of matrix data.
      ///
      /// \param dims [in] (numRows, numCols, numEntries).  These are
      ///   the "expected" values as read from the top of the Matrix
      ///   Market input stream.  Whether they are the final values
      ///   depends on the "tolerant" parameter and the actual matrix
      ///   data read from the input stream.
      ///
      /// \param tolerant [in] Whether the adder should be "tolerant"
      ///   of syntax errors and missing/incorrect metadata.  (In
      ///   particular, this refers to the number of rows, columns,
      ///   and entries in the matrix.)
      /// \param debug [in] Whether to print verbose debug output
      ///   to stderr on Rank 0.
      ///
      /// \return An adder_type object [nonnull and valid on Rank 0
      ///   only] that optionally symmetrizes the entries of the
      ///   sparse matrix.
      ///
      static Teuchos::RCP<adder_type>
      makeAdder (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                 Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                 const Teuchos::Tuple<global_ordinal_type, 3>& dims,
                 const bool tolerant=false,
                 const bool debug=false)
      {
        if (pComm->getRank () == 0) {
          typedef Teuchos::MatrixMarket::Raw::Adder<scalar_type,
                                                    global_ordinal_type>
            raw_adder_type;
          Teuchos::RCP<raw_adder_type> pRaw =
            Teuchos::rcp (new raw_adder_type (dims[0], dims[1], dims[2],
                                              tolerant, debug));
          return Teuchos::rcp (new adder_type (pRaw, pBanner->symmType ()));
        }
        else {
          return Teuchos::null;
        }
      }

      /// \brief Make an "adder" object for processing graph data.
      ///
      /// \param pComm [in] Communicator (across whose MPI ranks
      ///   the sparse matrix will be distributed)
      ///
      /// \param banner [in, nonnull and valid on Rank 0 only]
      ///   Object describing the type and symmetry of matrix data.
      ///
      /// \param dims [in] (numRows, numCols, numEntries).  These are
      ///   the "expected" values as read from the top of the Matrix
      ///   Market input stream.  Whether they are the final values
      ///   depends on the "tolerant" parameter and the actual graph
      ///   data read from the input stream.
      ///
      /// \param tolerant [in] Whether the adder should be "tolerant"
      ///   of syntax errors and missing/incorrect metadata.  (In
      ///   particular, this refers to the number of rows, columns,
      ///   and entries in the graph.)
      /// \param debug [in] Whether to print verbose debug output
      ///   to stderr on Rank 0.
      ///
      /// \return An adder_type object [nonnull and valid on Rank 0
      ///   only] that optionally symmetrizes the entries of the
      ///   sparse matrix.
      ///
      static Teuchos::RCP<graph_adder_type>
      makeGraphAdder (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                      Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                      const Teuchos::Tuple<global_ordinal_type, 3>& dims,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        if (pComm->getRank () == 0) {
          typedef Teuchos::MatrixMarket::Raw::GraphAdder<global_ordinal_type> raw_adder_type;
          Teuchos::RCP<raw_adder_type> pRaw =
            Teuchos::rcp (new raw_adder_type (dims[0], dims[1], dims[2],
                                              tolerant, debug));
          return Teuchos::rcp (new graph_adder_type (pRaw, pBanner->symmType ()));
        }
        else {
          return Teuchos::null;
        }
      }

      /// \brief Read sparse graph from the given Matrix Market input stream.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraphHelper (std::istream& in,
                             const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                             const Teuchos::RCP<const map_type>& rowMap,
                             Teuchos::RCP<const map_type>& colMap,
                             const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                             const bool tolerant,
                             const bool debug)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::RCP;
        using Teuchos::ptr;
        using Teuchos::Tuple;
        using std::cerr;
        using std::endl;

        const int myRank = pComm->getRank ();
        const int rootRank = 0;

        // Current line number in the input stream.  Various calls
        // will modify this depending on the number of lines that are
        // read from the input stream.  Only Rank 0 modifies this.
        size_t lineNumber = 1;

        if (debug && myRank == rootRank) {
          cerr << "Matrix Market reader: readGraph:" << endl
               << "-- Reading banner line" << endl;
        }

        // The "Banner" tells you whether the input stream represents
        // a sparse matrix, the symmetry type of the matrix, and the
        // type of the data it contains.
        //
        // pBanner will only be nonnull on MPI Rank 0.  It will be
        // null on all other MPI processes.
        RCP<const Banner> pBanner;
        {
          // We read and validate the Banner on Proc 0, but broadcast
          // the validation result to all processes.
          // Teuchos::broadcast doesn't currently work with bool, so
          // we use int (true -> 1, false -> 0).
          int bannerIsCorrect = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            // Read the Banner line from the input stream.
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug, true);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                  "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }

            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse graph.
              // We validate on Proc 0, since it reads the Banner.

              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                    "\"coordinate\"-format sparse graph in order to create a "
                    "Tpetra::CrsGraph object from it, but the file's matrix "
                    "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                    "format sparse graph in order to create a Tpetra::CrsGraph "
                    "object from it, but the file's matrix type is \"array\" "
                    "instead.  That probably means the file contains dense matrix "
                    "data.";
              }
            }
          } // Proc 0: Done reading the Banner, hopefully successfully.

          // Broadcast from Proc 0 whether the Banner was read correctly.
          broadcast (*pComm, rootRank, ptr (&bannerIsCorrect));

          // If the Banner is invalid, all processes throw an
          // exception.  Only Proc 0 gets the exception message, but
          // that's OK, since the main point is to "stop the world"
          // (rather than throw an exception on one process and leave
          // the others hanging).
          TEUCHOS_TEST_FOR_EXCEPTION(bannerIsCorrect == 0,
              std::invalid_argument, errMsg.str ());
        } // Done reading the Banner line and broadcasting success.
        if (debug && myRank == rootRank) {
          cerr << "-- Reading dimensions line" << endl;
        }

        // Read the graph dimensions from the Matrix Market metadata.
        // dims = (numRows, numCols, numEntries).  Proc 0 does the
        // reading, but it broadcasts the results to all MPI
        // processes.  Thus, readCoordDims() is a collective
        // operation.  It does a collective check for correctness too.
        Tuple<global_ordinal_type, 3> dims =
          readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Making Adder for collecting graph data" << endl;
        }

        // "Adder" object for collecting all the sparse graph entries
        // from the input stream.  This is only nonnull on Proc 0.
        // The Adder internally converts the one-based indices (native
        // Matrix Market format) into zero-based indices.
        RCP<graph_adder_type> pAdder =
          makeGraphAdder (pComm, pBanner, dims, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Reading graph data" << endl;
        }
        //
        // Read the graph entries from the input stream on Proc 0.
        //
        {
          // We use readSuccess to broadcast the results of the read
          // (succeeded or not) to all MPI processes.  Since
          // Teuchos::broadcast doesn't currently know how to send
          // bools, we convert to int (true -> 1, false -> 0).
          int readSuccess = 1;
          std::ostringstream errMsg; // Exception message (only valid on Proc 0)
          if (myRank == rootRank) {
            try {
              // Reader for "coordinate" format sparse graph data.
              typedef Teuchos::MatrixMarket::CoordPatternReader<graph_adder_type,
                global_ordinal_type> reader_type;
              reader_type reader (pAdder);

              // Read the sparse graph entries.
              std::pair<bool, std::vector<size_t> > results =
                reader.read (in, lineNumber, tolerant, debug);
              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }
          }
          broadcast (*pComm, rootRank, ptr (&readSuccess));

          // It would be nice to add a "verbose" flag, so that in
          // tolerant mode, we could log any bad line number(s) on
          // Proc 0.  For now, we just throw if the read fails to
          // succeed.
          //
          // Question: If we're in tolerant mode, and if the read did
          // not succeed, should we attempt to call fillComplete()?
          TEUCHOS_TEST_FOR_EXCEPTION(readSuccess == 0, std::runtime_error,
            "Failed to read the Matrix Market sparse graph file: "
            << errMsg.str());
        } // Done reading the graph entries (stored on Proc 0 for now)

        if (debug && myRank == rootRank) {
          cerr << "-- Successfully read the Matrix Market data" << endl;
        }

        // In tolerant mode, we need to rebroadcast the graph
        // dimensions, since they may be different after reading the
        // actual graph data.  We only need to broadcast the number
        // of rows and columns.  Only Rank 0 needs to know the actual
        // global number of entries, since (a) we need to merge
        // duplicates on Rank 0 first anyway, and (b) when we
        // distribute the entries, each rank other than Rank 0 will
        // only need to know how many entries it owns, not the total
        // number of entries.
        if (tolerant) {
          if (debug && myRank == rootRank) {
            cerr << "-- Tolerant mode: rebroadcasting graph dimensions"
                 << endl
                 << "----- Dimensions before: "
                 << dims[0] << " x " << dims[1]
                 << endl;
          }
          // Packed coordinate graph dimensions (numRows, numCols).
          Tuple<global_ordinal_type, 2> updatedDims;
          if (myRank == rootRank) {
            // If one or more bottom rows of the graph contain no
            // entries, then the Adder will report that the number
            // of rows is less than that specified in the
            // metadata.  We allow this case, and favor the
            // metadata so that the zero row(s) will be included.
            updatedDims[0] =
              std::max (dims[0], pAdder->getAdder()->numRows());
            updatedDims[1] = pAdder->getAdder()->numCols();
          }
          broadcast (*pComm, rootRank, updatedDims);
          dims[0] = updatedDims[0];
          dims[1] = updatedDims[1];
          if (debug && myRank == rootRank) {
            cerr << "----- Dimensions after: " << dims[0] << " x "
                 << dims[1] << endl;
          }
        }
        else {
          // In strict mode, we require that the graph's metadata and
          // its actual data agree, at least somewhat.  In particular,
          // the number of rows must agree, since otherwise we cannot
          // distribute the graph correctly.

          // Teuchos::broadcast() doesn't know how to broadcast bools,
          // so we use an int with the standard 1 == true, 0 == false
          // encoding.
          int dimsMatch = 1;
          if (myRank == rootRank) {
            // If one or more bottom rows of the graph contain no
            // entries, then the Adder will report that the number of
            // rows is less than that specified in the metadata.  We
            // allow this case, and favor the metadata, but do not
            // allow the Adder to think there are more rows in the
            // graph than the metadata says.
            if (dims[0] < pAdder->getAdder ()->numRows ()) {
              dimsMatch = 0;
            }
          }
          broadcast (*pComm, 0, ptr (&dimsMatch));
          if (dimsMatch == 0) {
            // We're in an error state anyway, so we might as well
            // work a little harder to print an informative error
            // message.
            //
            // Broadcast the Adder's idea of the graph dimensions
            // from Proc 0 to all processes.
            Tuple<global_ordinal_type, 2> addersDims;
            if (myRank == rootRank) {
              addersDims[0] = pAdder->getAdder()->numRows();
              addersDims[1] = pAdder->getAdder()->numCols();
            }
            broadcast (*pComm, 0, addersDims);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dimsMatch == 0, std::runtime_error,
              "The graph metadata says that the graph is " << dims[0] << " x "
              << dims[1] << ", but the actual data says that the graph is "
              << addersDims[0] << " x " << addersDims[1] << ".  That means the "
              "data includes more rows than reported in the metadata.  This "
              "is not allowed when parsing in strict mode.  Parse the graph in "
              "tolerant mode to ignore the metadata when it disagrees with the "
              "data.");
          }
        } // Matrix dimensions (# rows, # cols, # entries) agree.

        // Create a map describing a distribution where the root owns EVERYTHING
        RCP<map_type> proc0Map;
        global_ordinal_type indexBase;
        if(Teuchos::is_null(rowMap)) {
          indexBase = 0;
        }
        else {
          indexBase = rowMap->getIndexBase();
        }
        if(myRank == rootRank) {
          proc0Map = rcp(new map_type(dims[0],dims[0],indexBase,pComm));
        }
        else {
          proc0Map = rcp(new map_type(dims[0],0,indexBase,pComm));
        }

        // Create the graph where the root owns EVERYTHING
        std::map<global_ordinal_type, size_t> numEntriesPerRow_map;
        if (myRank == rootRank) {
          const auto& entries = pAdder()->getAdder()->getEntries();
          // This will count duplicates, but it's better than dense.
          // An even better approach would use a classic algorithm,
          // likely in Saad's old textbook, for converting COO (entries)
          // to CSR (the local part of the sparse matrix data structure).
          for (const auto& entry : entries) {
            const global_ordinal_type gblRow = entry.rowIndex () + indexBase;
            ++numEntriesPerRow_map[gblRow];
          }
        }

        Teuchos::Array<size_t> numEntriesPerRow (proc0Map->getLocalNumElements ());
        for (const auto& ent : numEntriesPerRow_map) {
          const local_ordinal_type lclRow = proc0Map->getLocalElement (ent.first);
          numEntriesPerRow[lclRow] = ent.second;
        }
        // Free anything we don't need before allocating the graph.
        // Swapping with an empty data structure is the standard idiom
        // for freeing memory used by Standard Library containers.
        // (Just resizing to 0 doesn't promise to free memory.)
        {
          std::map<global_ordinal_type, size_t> empty_map;
          std::swap (numEntriesPerRow_map, empty_map);
        }

        RCP<sparse_graph_type> proc0Graph =
          rcp(new sparse_graph_type(proc0Map,numEntriesPerRow (),
                                    constructorParams));
        if(myRank == rootRank) {
          typedef Teuchos::MatrixMarket::Raw::GraphElement<global_ordinal_type> element_type;

          // Get the entries
          const std::vector<element_type>& entries =
            pAdder->getAdder()->getEntries();

          // Insert them one at a time
          for(size_t curPos=0; curPos<entries.size(); curPos++) {
            const element_type& curEntry = entries[curPos];
            const global_ordinal_type curRow = curEntry.rowIndex()+indexBase;
            const global_ordinal_type curCol = curEntry.colIndex()+indexBase;
            Teuchos::ArrayView<const global_ordinal_type> colView(&curCol,1);
            proc0Graph->insertGlobalIndices(curRow,colView);
          }
        }
        proc0Graph->fillComplete();

        RCP<sparse_graph_type> distGraph;
        if(Teuchos::is_null(rowMap))
        {
          // Create a map describing the distribution we actually want
          RCP<map_type> distMap =
              rcp(new map_type(dims[0],0,pComm,GloballyDistributed));

          // Create the graph with that distribution too
          distGraph = rcp(new sparse_graph_type(distMap,colMap,0,constructorParams));

          // Create an importer/exporter/vandelay to redistribute the graph
          typedef Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
          import_type importer (proc0Map, distMap);

          // Import the data
          distGraph->doImport(*proc0Graph,importer,INSERT);
        }
        else {
          distGraph = rcp(new sparse_graph_type(rowMap,colMap,0,constructorParams));

          // Create an importer/exporter/vandelay to redistribute the graph
          typedef Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
          import_type importer (proc0Map, rowMap);

          // Import the data
          distGraph->doImport(*proc0Graph,importer,INSERT);
        }

        return distGraph;
      }

    public:
      /// \brief Read sparse graph from the given Matrix Market file.
      ///
      /// Open the given file on MPI Rank 0 (with respect to the given
      /// communicator).  The file should contain Matrix Market
      /// "coordinate" format sparse matrix data.  Read that data on
      /// Rank 0, and distribute it to all processors.  Return the
      /// resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 opens the
      ///   file and reads data from it, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraphFile (const std::string& filename,
                           const trcp_tcomm_t& comm,
                           const bool callFillComplete=true,
                           const bool tolerant=false,
                           const bool debug=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readSparseGraph (in, comm,
                                callFillComplete,
                                tolerant, debug);
        // We can rely on the destructor of the input stream to close
        // the file on scope exit, even if readSparseGraph() throws an
        // exception.
      }


      /// \brief Read sparse graph from the given Matrix Market file.
      ///
      /// This is a variant of readSparseGraph() that lets you pass
      /// parameters to the CrsGraph's constructor and to its
      /// fillComplete() method.
      ///
      /// Open the given file on Process 0 in the given communicator.
      /// The file must contain Matrix Market "coordinate" format
      /// sparse matrix data.  Read that data on Process 0, and
      /// distribute it to all processes.  Return the resulting
      /// distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 opens
      ///   the file and reads data from it, but all processes
      ///   participate and wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all process(es)
      ///   over which the sparse matrix will be distributed.
      /// \param constructorParams [in/out] Parameters for the
      ///   CrsMatrix constructor.
      /// \param fillCompleteParams [in/out] Parameters for
      ///   CrsMatrix's fillComplete call.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraphFile (const std::string& filename,
                           const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                           const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                           const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                           const bool tolerant=false,
                           const bool debug=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(pComm, filename, true);

        return readSparseGraph (in, pComm,
                                constructorParams,
                                fillCompleteParams, tolerant, debug);
        // We can rely on the destructor of the input stream to close
        // the file on scope exit, even if readSparseGraph() throws an
        // exception.
      }


      /// \brief Read sparse graph from the given Matrix Market file,
      ///   with provided Maps.
      ///
      /// This version of readSparseGraph() requires you to provide a
      /// row Map, domain Map, and range Map.  You may, if you wish,
      /// provide a column Map as well, but this is not required.
      ///
      /// Open the given file on Process 0 (with respect to the given
      /// Maps' communicator).  The file should contain Matrix Market
      /// "coordinate" format sparse matrix data.  Read that data on
      /// Process 0, and distribute it to all processors.  Return the
      /// resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 opens
      ///   the file and reads data from it, but all ranks participate
      ///   and wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param rowMap [in] The Map over which to distribute rows
      ///   of the sparse matrix.  This must be nonnull.
      /// \param colMap [in/out] If nonnull: the Map over which to
      ///   distribute columns of the sparse matrix.  If null and if
      ///   callFillComplete is true, we create this for you.
      /// \param domainMap [in] The sparse matrix's domain Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param rangeMap [in] The sparse matrix's range Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsGraph, after adding all the entries
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraphFile (const std::string& filename,
                           const Teuchos::RCP<const map_type>& rowMap,
                           Teuchos::RCP<const map_type>& colMap,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const bool callFillComplete=true,
                           const bool tolerant=false,
                           const bool debug=false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (rowMap.is_null (), std::invalid_argument,
           "Input rowMap must be nonnull.");
        trcp_tcomm_t comm = rowMap->getComm ();
        if (comm.is_null ()) {
          // If the input communicator is null on some process, then
          // that process does not participate in the collective.
          return Teuchos::null;
        }

        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readSparseGraph (in, rowMap, colMap, domainMap, rangeMap,
                                callFillComplete, tolerant, debug);
      }

      /// \brief Read sparse graph from the given Matrix Market input stream.
      ///
      /// The given input stream need only be readable by MPI Rank 0
      /// (with respect to the given communicator).  The input stream
      /// should contain Matrix Market "coordinate" format sparse
      /// matrix data.  Read that data on Rank 0, and distribute it to
      /// all processors.  Return the resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 reads
      ///   data from the input stream, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param in [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraph (std::istream& in,
                       const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                       const bool callFillComplete=true,
                       const bool tolerant=false,
                       const bool debug=false)
      {
        Teuchos::RCP<const map_type> fakeRowMap;
        Teuchos::RCP<const map_type> fakeColMap;
        Teuchos::RCP<Teuchos::ParameterList> fakeCtorParams;

        Teuchos::RCP<sparse_graph_type> graph =
          readSparseGraphHelper (in, pComm,
                                 fakeRowMap, fakeColMap,
                                 fakeCtorParams, tolerant, debug);
        if (callFillComplete) {
          graph->fillComplete ();
        }
        return graph;
      }


      /// \brief Read sparse graph from the given Matrix Market input
      ///   stream.
      ///
      /// This is a variant of readSparse() that lets you pass
      /// parameters to the CrsMatrix's constructor and to its
      /// fillComplete() method.
      ///
      /// The given input stream need only be readable by Process 0 in
      /// the given communicator.  The input stream must contain
      /// Matrix Market "coordinate" format sparse matrix data.  Read
      /// that data on Process 0, and distribute it to all Processes.
      /// Return the resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Proces 0 reads
      ///   data from the input stream, but all processes participate
      ///   and wait for the final result.
      ///
      /// \param in [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all process(es)
      ///   over which the sparse matrix will be distributed.
      /// \param constructorParams [in/out] Parameters for the
      ///   CrsMatrix constructor.
      /// \param fillCompleteParams [in/out] Parameters for
      ///   CrsMatrix's fillComplete call.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraph (std::istream& in,
                       const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                       const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                       const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                       const bool tolerant=false,
                       const bool debug=false)
      {
        Teuchos::RCP<const map_type> fakeRowMap;
        Teuchos::RCP<const map_type> fakeColMap;
        Teuchos::RCP<sparse_graph_type> graph =
          readSparseGraphHelper (in, pComm,
                                 fakeRowMap, fakeColMap,
                                 constructorParams, tolerant, debug);
        graph->fillComplete (fillCompleteParams);
        return graph;
      }


      /// \brief Read sparse matrix from the given Matrix Market input
      ///   stream, with provided Maps.
      ///
      /// This version of readSparse() requires you to provide a row
      /// Map, domain Map, and range Map.  You may, if you wish,
      /// provide a column Map as well, but this is not required.
      ///
      /// The given input stream need only be readable by Process 0
      /// (with respect to the given Maps' communicator).  The input
      /// stream must contain Matrix Market "coordinate" format sparse
      /// matrix data.  Read that data on Process 0, and distribute it
      /// to all processors.  Return the resulting distributed
      /// CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 reads
      ///   data from the input stream, but all processes participate
      ///   and wait for the final result.
      ///
      /// \param in [in/out] The input stream from which to read.
      /// \param rowMap [in] The Map over which to distribute rows
      ///   of the sparse matrix.  This must be nonnull.
      /// \param colMap [in/out] If nonnull: the Map over which to
      ///   distribute columns of the sparse matrix.  If null and if
      ///   callFillComplete is true, we create this for you.
      /// \param domainMap [in] The sparse matrix's domain Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param rangeMap [in] The sparse matrix's range Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_graph_type>
      readSparseGraph (std::istream& in,
                       const Teuchos::RCP<const map_type>& rowMap,
                       Teuchos::RCP<const map_type>& colMap,
                       const Teuchos::RCP<const map_type>& domainMap,
                       const Teuchos::RCP<const map_type>& rangeMap,
                       const bool callFillComplete=true,
                       const bool tolerant=false,
                       const bool debug=false)
      {
        Teuchos::RCP<sparse_graph_type> graph =
          readSparseGraphHelper (in, rowMap->getComm (),
                                 rowMap, colMap, Teuchos::null, tolerant,
                                 debug);
        if (callFillComplete) {
          graph->fillComplete (domainMap, rangeMap);
        }
        return graph;
      }

#include "MatrixMarket_TpetraNew.hpp"

      /// \brief Read sparse matrix from the given Matrix Market file.
      ///
      /// Open the given file on MPI Rank 0 (with respect to the given
      /// communicator).  The file should contain Matrix Market
      /// "coordinate" format sparse matrix data.  Read that data on
      /// Rank 0, and distribute it to all processors.  Return the
      /// resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 opens the
      ///   file and reads data from it, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparseFile (const std::string& filename,
                      const trcp_tcomm_t& comm,
                      const bool callFillComplete=true,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readSparse (in, comm, callFillComplete, tolerant, debug);
        // We can rely on the destructor of the input stream to close
        // the file on scope exit, even if readSparse() throws an
        // exception.
      }


      /// \brief Read sparse matrix from the given Matrix Market file.
      ///
      /// This is a variant of readSparseFile() that lets you pass
      /// parameters to the CrsMatrix's constructor and to its
      /// fillComplete() method.
      ///
      /// Open the given file on Process 0 in the given communicator.
      /// The file must contain Matrix Market "coordinate" format
      /// sparse matrix data.  Read that data on Process 0, and
      /// distribute it to all processes.  Return the resulting
      /// distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 opens
      ///   the file and reads data from it, but all processes
      ///   participate and wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all process(es)
      ///   over which the sparse matrix will be distributed.
      /// \param constructorParams [in/out] Parameters for the
      ///   CrsMatrix constructor.
      /// \param fillCompleteParams [in/out] Parameters for
      ///   CrsMatrix's fillComplete call.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparseFile (const std::string& filename,
                      const trcp_tcomm_t& comm,
                      const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                      const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readSparse (in, comm, constructorParams,
                           fillCompleteParams, tolerant, debug);
      }


      /// \brief Read sparse matrix from the given Matrix Market file,
      ///   with provided Maps.
      ///
      /// This version of readSparseFile() requires you to provide a
      /// row Map, domain Map, and range Map.  You may, if you wish,
      /// provide a column Map as well, but this is not required.
      ///
      /// Open the given file on Process 0 (with respect to the given
      /// Maps' communicator).  The file should contain Matrix Market
      /// "coordinate" format sparse matrix data.  Read that data on
      /// Process 0, and distribute it to all processors.  Return the
      /// resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 opens
      ///   the file and reads data from it, but all ranks participate
      ///   and wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param rowMap [in] The Map over which to distribute rows
      ///   of the sparse matrix.  This must be nonnull.
      /// \param colMap [in/out] If nonnull: the Map over which to
      ///   distribute columns of the sparse matrix.  If null and if
      ///   callFillComplete is true, we create this for you.
      /// \param domainMap [in] The sparse matrix's domain Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param rangeMap [in] The sparse matrix's range Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparseFile (const std::string& filename,
                      const Teuchos::RCP<const map_type>& rowMap,
                      Teuchos::RCP<const map_type>& colMap,
                      const Teuchos::RCP<const map_type>& domainMap,
                      const Teuchos::RCP<const map_type>& rangeMap,
                      const bool callFillComplete=true,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          rowMap.is_null (), std::invalid_argument,
          "Row Map must be nonnull.");

        trcp_tcomm_t comm = rowMap->getComm ();

        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readSparse (in, rowMap, colMap, domainMap, rangeMap,
                           callFillComplete, tolerant, debug);
      }

      /// \brief Read sparse matrix from the given Matrix Market input stream.
      ///
      /// The given input stream need only be readable by MPI Rank 0
      /// (with respect to the given communicator).  The input stream
      /// should contain Matrix Market "coordinate" format sparse
      /// matrix data.  Read that data on Rank 0, and distribute it to
      /// all processors.  Return the resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 reads
      ///   data from the input stream, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param in [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparse (std::istream& in,
                  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                  const bool callFillComplete=true,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::arcp;
        using Teuchos::ArrayRCP;
        using Teuchos::broadcast;
        using Teuchos::null;
        using Teuchos::ptr;
        using Teuchos::RCP;
        using Teuchos::REDUCE_MAX;
        using Teuchos::reduceAll;
        using Teuchos::Tuple;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        const bool extraDebug = false;
        const int myRank = pComm->getRank ();
        const int rootRank = 0;

        // Current line number in the input stream.  Various calls
        // will modify this depending on the number of lines that are
        // read from the input stream.  Only Rank 0 modifies this.
        size_t lineNumber = 1;

        if (debug && myRank == rootRank) {
          cerr << "Matrix Market reader: readSparse:" << endl
               << "-- Reading banner line" << endl;
        }

        // The "Banner" tells you whether the input stream represents
        // a sparse matrix, the symmetry type of the matrix, and the
        // type of the data it contains.
        //
        // pBanner will only be nonnull on MPI Rank 0.  It will be
        // null on all other MPI processes.
        RCP<const Banner> pBanner;
        {
          // We read and validate the Banner on Proc 0, but broadcast
          // the validation result to all processes.
          // Teuchos::broadcast doesn't currently work with bool, so
          // we use int (true -> 1, false -> 0).
          int bannerIsCorrect = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            // Read the Banner line from the input stream.
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }

            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse matrix.
              // We validate on Proc 0, since it reads the Banner.

              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                  "\"coordinate\"-format sparse matrix in order to create a "
                  "Tpetra::CrsMatrix object from it, but the file's matrix "
                  "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                  "format sparse matrix in order to create a Tpetra::CrsMatrix "
                  "object from it, but the file's matrix type is \"array\" "
                  "instead.  That probably means the file contains dense matrix "
                  "data.";
              }
            }
          } // Proc 0: Done reading the Banner, hopefully successfully.

          // Broadcast from Proc 0 whether the Banner was read correctly.
          broadcast (*pComm, rootRank, ptr (&bannerIsCorrect));

          // If the Banner is invalid, all processes throw an
          // exception.  Only Proc 0 gets the exception message, but
          // that's OK, since the main point is to "stop the world"
          // (rather than throw an exception on one process and leave
          // the others hanging).
          TEUCHOS_TEST_FOR_EXCEPTION(bannerIsCorrect == 0,
            std::invalid_argument, errMsg.str ());
        } // Done reading the Banner line and broadcasting success.
        if (debug && myRank == rootRank) {
          cerr << "-- Reading dimensions line" << endl;
        }

        // Read the matrix dimensions from the Matrix Market metadata.
        // dims = (numRows, numCols, numEntries).  Proc 0 does the
        // reading, but it broadcasts the results to all MPI
        // processes.  Thus, readCoordDims() is a collective
        // operation.  It does a collective check for correctness too.
        Tuple<global_ordinal_type, 3> dims =
          readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Making Adder for collecting matrix data" << endl;
        }

        // "Adder" object for collecting all the sparse matrix entries
        // from the input stream.  This is only nonnull on Proc 0.
        RCP<adder_type> pAdder =
          makeAdder (pComm, pBanner, dims, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Reading matrix data" << endl;
        }
        //
        // Read the matrix entries from the input stream on Proc 0.
        //
        {
          // We use readSuccess to broadcast the results of the read
          // (succeeded or not) to all MPI processes.  Since
          // Teuchos::broadcast doesn't currently know how to send
          // bools, we convert to int (true -> 1, false -> 0).
          int readSuccess = 1;
          std::ostringstream errMsg; // Exception message (only valid on Proc 0)
          if (myRank == rootRank) {
            try {
              // Reader for "coordinate" format sparse matrix data.
              typedef Teuchos::MatrixMarket::CoordDataReader<adder_type,
                global_ordinal_type, scalar_type, STS::isComplex> reader_type;
              reader_type reader (pAdder);

              // Read the sparse matrix entries.
              std::pair<bool, std::vector<size_t> > results =
                reader.read (in, lineNumber, tolerant, debug);
              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }
          }
          broadcast (*pComm, rootRank, ptr (&readSuccess));

          // It would be nice to add a "verbose" flag, so that in
          // tolerant mode, we could log any bad line number(s) on
          // Proc 0.  For now, we just throw if the read fails to
          // succeed.
          //
          // Question: If we're in tolerant mode, and if the read did
          // not succeed, should we attempt to call fillComplete()?
          TEUCHOS_TEST_FOR_EXCEPTION(readSuccess == 0, std::runtime_error,
            "Failed to read the Matrix Market sparse matrix file: "
            << errMsg.str());
        } // Done reading the matrix entries (stored on Proc 0 for now)

        if (debug && myRank == rootRank) {
          cerr << "-- Successfully read the Matrix Market data" << endl;
        }

        // In tolerant mode, we need to rebroadcast the matrix
        // dimensions, since they may be different after reading the
        // actual matrix data.  We only need to broadcast the number
        // of rows and columns.  Only Rank 0 needs to know the actual
        // global number of entries, since (a) we need to merge
        // duplicates on Rank 0 first anyway, and (b) when we
        // distribute the entries, each rank other than Rank 0 will
        // only need to know how many entries it owns, not the total
        // number of entries.
        if (tolerant) {
          if (debug && myRank == rootRank) {
            cerr << "-- Tolerant mode: rebroadcasting matrix dimensions"
                 << endl
                 << "----- Dimensions before: "
                 << dims[0] << " x " << dims[1]
                 << endl;
          }
          // Packed coordinate matrix dimensions (numRows, numCols).
          Tuple<global_ordinal_type, 2> updatedDims;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number
            // of rows is less than that specified in the
            // metadata.  We allow this case, and favor the
            // metadata so that the zero row(s) will be included.
            updatedDims[0] =
              std::max (dims[0], pAdder->getAdder()->numRows());
            updatedDims[1] = pAdder->getAdder()->numCols();
          }
          broadcast (*pComm, rootRank, updatedDims);
          dims[0] = updatedDims[0];
          dims[1] = updatedDims[1];
          if (debug && myRank == rootRank) {
            cerr << "----- Dimensions after: " << dims[0] << " x "
                 << dims[1] << endl;
          }
        }
        else {
          // In strict mode, we require that the matrix's metadata and
          // its actual data agree, at least somewhat.  In particular,
          // the number of rows must agree, since otherwise we cannot
          // distribute the matrix correctly.

          // Teuchos::broadcast() doesn't know how to broadcast bools,
          // so we use an int with the standard 1 == true, 0 == false
          // encoding.
          int dimsMatch = 1;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number of
            // rows is less than that specified in the metadata.  We
            // allow this case, and favor the metadata, but do not
            // allow the Adder to think there are more rows in the
            // matrix than the metadata says.
            if (dims[0] < pAdder->getAdder ()->numRows ()) {
              dimsMatch = 0;
            }
          }
          broadcast (*pComm, 0, ptr (&dimsMatch));
          if (dimsMatch == 0) {
            // We're in an error state anyway, so we might as well
            // work a little harder to print an informative error
            // message.
            //
            // Broadcast the Adder's idea of the matrix dimensions
            // from Proc 0 to all processes.
            Tuple<global_ordinal_type, 2> addersDims;
            if (myRank == rootRank) {
              addersDims[0] = pAdder->getAdder()->numRows();
              addersDims[1] = pAdder->getAdder()->numCols();
            }
            broadcast (*pComm, 0, addersDims);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dimsMatch == 0, std::runtime_error,
              "The matrix metadata says that the matrix is " << dims[0] << " x "
              << dims[1] << ", but the actual data says that the matrix is "
              << addersDims[0] << " x " << addersDims[1] << ".  That means the "
              "data includes more rows than reported in the metadata.  This "
              "is not allowed when parsing in strict mode.  Parse the matrix in "
              "tolerant mode to ignore the metadata when it disagrees with the "
              "data.");
          }
        } // Matrix dimensions (# rows, # cols, # entries) agree.

        if (debug && myRank == rootRank) {
          cerr << "-- Converting matrix data into CSR format on Proc 0" << endl;
        }

        // Now that we've read in all the matrix entries from the
        // input stream into the adder on Proc 0, post-process them
        // into CSR format (still on Proc 0).  This will facilitate
        // distributing them to all the processors.
        //
        // These arrays represent the global matrix data as a CSR
        // matrix (with numEntriesPerRow as redundant but convenient
        // metadata, since it's computable from rowPtr and vice
        // versa).  They are valid only on Proc 0.
        ArrayRCP<size_t> numEntriesPerRow;
        ArrayRCP<size_t> rowPtr;
        ArrayRCP<global_ordinal_type> colInd;
        ArrayRCP<scalar_type> values;

        // Proc 0 first merges duplicate entries, and then converts
        // the coordinate-format matrix data to CSR.
        {
          int mergeAndConvertSucceeded = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            try {
              typedef Teuchos::MatrixMarket::Raw::Element<scalar_type,
                global_ordinal_type> element_type;

              // Number of rows in the matrix.  If we are in tolerant
              // mode, we've already synchronized dims with the actual
              // matrix data.  If in strict mode, we should use dims
              // (as read from the file's metadata) rather than the
              // matrix data to determine the dimensions.  (The matrix
              // data will claim fewer rows than the metadata, if one
              // or more rows have no entries stored in the file.)
              const size_type numRows = dims[0];

              // Additively merge duplicate matrix entries.
              pAdder->getAdder()->merge ();

              // Get a temporary const view of the merged matrix entries.
              const std::vector<element_type>& entries =
                pAdder->getAdder()->getEntries();

              // Number of matrix entries (after merging).
              const size_t numEntries = (size_t)entries.size();

              if (debug) {
                cerr << "----- Proc 0: Matrix has numRows=" << numRows
                     << " rows and numEntries=" << numEntries
                     << " entries." << endl;
              }

              // Make space for the CSR matrix data.  Converting to
              // CSR is easier if we fill numEntriesPerRow with zeros
              // at first.
              numEntriesPerRow = arcp<size_t> (numRows);
              std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 0);
              rowPtr = arcp<size_t> (numRows+1);
              std::fill (rowPtr.begin(), rowPtr.end(), 0);
              colInd = arcp<global_ordinal_type> (numEntries);
              values = arcp<scalar_type> (numEntries);

              // Convert from array-of-structs coordinate format to CSR
              // (compressed sparse row) format.
              global_ordinal_type prvRow = 0;
              size_t curPos = 0;
              rowPtr[0] = 0;
              for (curPos = 0; curPos < numEntries; ++curPos) {
                const element_type& curEntry = entries[curPos];
                const global_ordinal_type curRow = curEntry.rowIndex();
                TEUCHOS_TEST_FOR_EXCEPTION(
                  curRow < prvRow, std::logic_error,
                  "Row indices are out of order, even though they are supposed "
                  "to be sorted.  curRow = " << curRow << ", prvRow = "
                  << prvRow << ", at curPos = " << curPos << ".  Please report "
                  "this bug to the Tpetra developers.");
                if (curRow > prvRow) {
                  for (global_ordinal_type r = prvRow+1; r <= curRow; ++r) {
                    rowPtr[r] = curPos;
                  }
                  prvRow = curRow;
                }
                numEntriesPerRow[curRow]++;
                colInd[curPos] = curEntry.colIndex();
                values[curPos] = curEntry.value();
              }
              // rowPtr has one more entry than numEntriesPerRow.  The
              // last entry of rowPtr is the number of entries in
              // colInd and values.
              rowPtr[numRows] = numEntries;
            } // Finished conversion to CSR format
            catch (std::exception& e) {
              mergeAndConvertSucceeded = 0;
              errMsg << "Failed to merge sparse matrix entries and convert to "
                "CSR format: " << e.what();
            }

            if (debug && mergeAndConvertSucceeded) {
              // Number of rows in the matrix.
              const size_type numRows = dims[0];
              const size_type maxToDisplay = 100;

              cerr << "----- Proc 0: numEntriesPerRow[0.."
                   << (numEntriesPerRow.size()-1) << "] ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > 0) {
                if (numRows > maxToDisplay) {
                  for (size_type k = 0; k < 2; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                  cerr << "... ";
                  for (size_type k = numRows-2; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                else {
                  for (size_type k = 0; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                cerr << numEntriesPerRow[numRows-1];
              } // numRows > 0
              cerr << "]" << endl;

              cerr << "----- Proc 0: rowPtr ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > maxToDisplay) {
                for (size_type k = 0; k < 2; ++k) {
                  cerr << rowPtr[k] << " ";
                }
                cerr << "... ";
                for (size_type k = numRows-2; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              else {
                for (size_type k = 0; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              cerr << rowPtr[numRows] << "]" << endl;
            }
          } // if myRank == rootRank
        } // Done converting sparse matrix data to CSR format

        // Now we're done with the Adder, so we can release the
        // reference ("free" it) to save space.  This only actually
        // does anything on Rank 0, since pAdder is null on all the
        // other MPI processes.
        pAdder = null;

        if (debug && myRank == rootRank) {
          cerr << "-- Making range, domain, and row maps" << endl;
        }

        // Make the maps that describe the matrix's range and domain,
        // and the distribution of its rows.  Creating a Map is a
        // collective operation, so we don't have to do a broadcast of
        // a success Boolean.
        RCP<const map_type> pRangeMap = makeRangeMap (pComm, dims[0]);
        RCP<const map_type> pDomainMap =
          makeDomainMap (pRangeMap, dims[0], dims[1]);
        RCP<const map_type> pRowMap = makeRowMap (null, pComm, dims[0]);

        if (debug && myRank == rootRank) {
          cerr << "-- Distributing the matrix data" << endl;
        }

        // Distribute the matrix data.  Each processor has to add the
        // rows that it owns.  If you try to make Proc 0 call
        // insertGlobalValues() for _all_ the rows, not just those it
        // owns, then fillComplete() will compute the number of
        // columns incorrectly.  That's why Proc 0 has to distribute
        // the matrix data and why we make all the processors (not
        // just Proc 0) call insertGlobalValues() on their own data.
        //
        // These arrays represent each processor's part of the matrix
        // data, in "CSR" format (sort of, since the row indices might
        // not be contiguous).
        ArrayRCP<size_t> myNumEntriesPerRow;
        ArrayRCP<size_t> myRowPtr;
        ArrayRCP<global_ordinal_type> myColInd;
        ArrayRCP<scalar_type> myValues;
        // Distribute the matrix data.  This is a collective operation.
        distribute (myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap,
                    numEntriesPerRow, rowPtr, colInd, values, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Inserting matrix entries on each processor";
          if (callFillComplete) {
            cerr << " and calling fillComplete()";
          }
          cerr << endl;
        }
        // Each processor inserts its part of the matrix data, and
        // then they all call fillComplete().  This method invalidates
        // the my* distributed matrix data before calling
        // fillComplete(), in order to save space.  In general, we
        // never store more than two copies of the matrix's entries in
        // memory at once, which is no worse than what Tpetra
        // promises.
        RCP<sparse_matrix_type> pMatrix =
          makeMatrix (myNumEntriesPerRow, myRowPtr, myColInd, myValues,
                      pRowMap, pRangeMap, pDomainMap, callFillComplete);
        // Only use a reduce-all in debug mode to check if pMatrix is
        // null.  Otherwise, just throw an exception.  We never expect
        // a null pointer here, so we can save a communication.
        if (debug) {
          int localIsNull = pMatrix.is_null () ? 1 : 0;
          int globalIsNull = 0;
          reduceAll (*pComm, REDUCE_MAX, localIsNull, ptr (&globalIsNull));
          TEUCHOS_TEST_FOR_EXCEPTION(globalIsNull != 0, std::logic_error,
            "Reader::makeMatrix() returned a null pointer on at least one "
            "process.  Please report this bug to the Tpetra developers.");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(pMatrix.is_null(), std::logic_error,
            "Reader::makeMatrix() returned a null pointer.  "
            "Please report this bug to the Tpetra developers.");
        }

        // We can't get the dimensions of the matrix until after
        // fillComplete() is called.  Thus, we can't do the sanity
        // check (dimensions read from the Matrix Market data,
        // vs. dimensions reported by the CrsMatrix) unless the user
        // asked makeMatrix() to call fillComplete().
        //
        // Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
        // what one might think it does, so you have to ask the range
        // resp. domain map for the number of rows resp. columns.
        if (callFillComplete) {
          const int numProcs = pComm->getSize ();

          if (extraDebug && debug) {
            const global_size_t globalNumRows =
              pRangeMap->getGlobalNumElements ();
            const global_size_t globalNumCols =
              pDomainMap->getGlobalNumElements ();
            if (myRank == rootRank) {
              cerr << "-- Matrix is "
                   << globalNumRows << " x " << globalNumCols
                   << " with " << pMatrix->getGlobalNumEntries()
                   << " entries, and index base "
                   << pMatrix->getIndexBase() << "." << endl;
            }
            pComm->barrier ();
            for (int p = 0; p < numProcs; ++p) {
              if (myRank == p) {
                cerr << "-- Proc " << p << " owns "
                     << pMatrix->getLocalNumCols() << " columns, and "
                     << pMatrix->getLocalNumEntries() << " entries." << endl;
              }
              pComm->barrier ();
            }
          } // if (extraDebug && debug)
        } // if (callFillComplete)

        if (debug && myRank == rootRank) {
          cerr << "-- Done creating the CrsMatrix from the Matrix Market data"
               << endl;
        }
        return pMatrix;
      }


      /// \brief Read sparse matrix from the given Matrix Market input stream.
      ///
      /// This is a variant of readSparse() that lets you pass
      /// parameters to the CrsMatrix's constructor and to its
      /// fillComplete() method.
      ///
      /// The given input stream need only be readable by Process 0 in
      /// the given communicator.  The input stream must contain
      /// Matrix Market "coordinate" format sparse matrix data.  Read
      /// that data on Process 0, and distribute it to all Processes.
      /// Return the resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Proces 0 reads
      ///   data from the input stream, but all processes participate
      ///   and wait for the final result.
      ///
      /// \param in [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all process(es)
      ///   over which the sparse matrix will be distributed.
      /// \param constructorParams [in/out] Parameters for the
      ///   CrsMatrix constructor.
      /// \param fillCompleteParams [in/out] Parameters for
      ///   CrsMatrix's fillComplete call.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparse (std::istream& in,
                  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
                  const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                  const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::arcp;
        using Teuchos::ArrayRCP;
        using Teuchos::broadcast;
        using Teuchos::null;
        using Teuchos::ptr;
        using Teuchos::RCP;
        using Teuchos::reduceAll;
        using Teuchos::Tuple;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        const bool extraDebug = false;
        const int myRank = pComm->getRank ();
        const int rootRank = 0;

        // Current line number in the input stream.  Various calls
        // will modify this depending on the number of lines that are
        // read from the input stream.  Only Rank 0 modifies this.
        size_t lineNumber = 1;

        if (debug && myRank == rootRank) {
          cerr << "Matrix Market reader: readSparse:" << endl
               << "-- Reading banner line" << endl;
        }

        // The "Banner" tells you whether the input stream represents
        // a sparse matrix, the symmetry type of the matrix, and the
        // type of the data it contains.
        //
        // pBanner will only be nonnull on MPI Rank 0.  It will be
        // null on all other MPI processes.
        RCP<const Banner> pBanner;
        {
          // We read and validate the Banner on Proc 0, but broadcast
          // the validation result to all processes.
          // Teuchos::broadcast doesn't currently work with bool, so
          // we use int (true -> 1, false -> 0).
          int bannerIsCorrect = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            // Read the Banner line from the input stream.
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }

            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse matrix.
              // We validate on Proc 0, since it reads the Banner.

              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                  "\"coordinate\"-format sparse matrix in order to create a "
                  "Tpetra::CrsMatrix object from it, but the file's matrix "
                  "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                  "format sparse matrix in order to create a Tpetra::CrsMatrix "
                  "object from it, but the file's matrix type is \"array\" "
                  "instead.  That probably means the file contains dense matrix "
                  "data.";
              }
            }
          } // Proc 0: Done reading the Banner, hopefully successfully.

          // Broadcast from Proc 0 whether the Banner was read correctly.
          broadcast (*pComm, rootRank, ptr (&bannerIsCorrect));

          // If the Banner is invalid, all processes throw an
          // exception.  Only Proc 0 gets the exception message, but
          // that's OK, since the main point is to "stop the world"
          // (rather than throw an exception on one process and leave
          // the others hanging).
          TEUCHOS_TEST_FOR_EXCEPTION(bannerIsCorrect == 0,
            std::invalid_argument, errMsg.str ());
        } // Done reading the Banner line and broadcasting success.
        if (debug && myRank == rootRank) {
          cerr << "-- Reading dimensions line" << endl;
        }

        // Read the matrix dimensions from the Matrix Market metadata.
        // dims = (numRows, numCols, numEntries).  Proc 0 does the
        // reading, but it broadcasts the results to all MPI
        // processes.  Thus, readCoordDims() is a collective
        // operation.  It does a collective check for correctness too.
        Tuple<global_ordinal_type, 3> dims =
          readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Making Adder for collecting matrix data" << endl;
        }

        // "Adder" object for collecting all the sparse matrix entries
        // from the input stream.  This is only nonnull on Proc 0.
        RCP<adder_type> pAdder =
          makeAdder (pComm, pBanner, dims, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Reading matrix data" << endl;
        }
        //
        // Read the matrix entries from the input stream on Proc 0.
        //
        {
          // We use readSuccess to broadcast the results of the read
          // (succeeded or not) to all MPI processes.  Since
          // Teuchos::broadcast doesn't currently know how to send
          // bools, we convert to int (true -> 1, false -> 0).
          int readSuccess = 1;
          std::ostringstream errMsg; // Exception message (only valid on Proc 0)
          if (myRank == rootRank) {
            try {
              // Reader for "coordinate" format sparse matrix data.
              typedef Teuchos::MatrixMarket::CoordDataReader<adder_type,
                global_ordinal_type, scalar_type, STS::isComplex> reader_type;
              reader_type reader (pAdder);

              // Read the sparse matrix entries.
              std::pair<bool, std::vector<size_t> > results =
                reader.read (in, lineNumber, tolerant, debug);
              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }
          }
          broadcast (*pComm, rootRank, ptr (&readSuccess));

          // It would be nice to add a "verbose" flag, so that in
          // tolerant mode, we could log any bad line number(s) on
          // Proc 0.  For now, we just throw if the read fails to
          // succeed.
          //
          // Question: If we're in tolerant mode, and if the read did
          // not succeed, should we attempt to call fillComplete()?
          TEUCHOS_TEST_FOR_EXCEPTION(readSuccess == 0, std::runtime_error,
            "Failed to read the Matrix Market sparse matrix file: "
            << errMsg.str());
        } // Done reading the matrix entries (stored on Proc 0 for now)

        if (debug && myRank == rootRank) {
          cerr << "-- Successfully read the Matrix Market data" << endl;
        }

        // In tolerant mode, we need to rebroadcast the matrix
        // dimensions, since they may be different after reading the
        // actual matrix data.  We only need to broadcast the number
        // of rows and columns.  Only Rank 0 needs to know the actual
        // global number of entries, since (a) we need to merge
        // duplicates on Rank 0 first anyway, and (b) when we
        // distribute the entries, each rank other than Rank 0 will
        // only need to know how many entries it owns, not the total
        // number of entries.
        if (tolerant) {
          if (debug && myRank == rootRank) {
            cerr << "-- Tolerant mode: rebroadcasting matrix dimensions"
                 << endl
                 << "----- Dimensions before: "
                 << dims[0] << " x " << dims[1]
                 << endl;
          }
          // Packed coordinate matrix dimensions (numRows, numCols).
          Tuple<global_ordinal_type, 2> updatedDims;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number
            // of rows is less than that specified in the
            // metadata.  We allow this case, and favor the
            // metadata so that the zero row(s) will be included.
            updatedDims[0] =
              std::max (dims[0], pAdder->getAdder()->numRows());
            updatedDims[1] = pAdder->getAdder()->numCols();
          }
          broadcast (*pComm, rootRank, updatedDims);
          dims[0] = updatedDims[0];
          dims[1] = updatedDims[1];
          if (debug && myRank == rootRank) {
            cerr << "----- Dimensions after: " << dims[0] << " x "
                 << dims[1] << endl;
          }
        }
        else {
          // In strict mode, we require that the matrix's metadata and
          // its actual data agree, at least somewhat.  In particular,
          // the number of rows must agree, since otherwise we cannot
          // distribute the matrix correctly.

          // Teuchos::broadcast() doesn't know how to broadcast bools,
          // so we use an int with the standard 1 == true, 0 == false
          // encoding.
          int dimsMatch = 1;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number of
            // rows is less than that specified in the metadata.  We
            // allow this case, and favor the metadata, but do not
            // allow the Adder to think there are more rows in the
            // matrix than the metadata says.
            if (dims[0] < pAdder->getAdder ()->numRows ()) {
              dimsMatch = 0;
            }
          }
          broadcast (*pComm, 0, ptr (&dimsMatch));
          if (dimsMatch == 0) {
            // We're in an error state anyway, so we might as well
            // work a little harder to print an informative error
            // message.
            //
            // Broadcast the Adder's idea of the matrix dimensions
            // from Proc 0 to all processes.
            Tuple<global_ordinal_type, 2> addersDims;
            if (myRank == rootRank) {
              addersDims[0] = pAdder->getAdder()->numRows();
              addersDims[1] = pAdder->getAdder()->numCols();
            }
            broadcast (*pComm, 0, addersDims);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dimsMatch == 0, std::runtime_error,
              "The matrix metadata says that the matrix is " << dims[0] << " x "
              << dims[1] << ", but the actual data says that the matrix is "
              << addersDims[0] << " x " << addersDims[1] << ".  That means the "
              "data includes more rows than reported in the metadata.  This "
              "is not allowed when parsing in strict mode.  Parse the matrix in "
              "tolerant mode to ignore the metadata when it disagrees with the "
              "data.");
          }
        } // Matrix dimensions (# rows, # cols, # entries) agree.

        if (debug && myRank == rootRank) {
          cerr << "-- Converting matrix data into CSR format on Proc 0" << endl;
        }

        // Now that we've read in all the matrix entries from the
        // input stream into the adder on Proc 0, post-process them
        // into CSR format (still on Proc 0).  This will facilitate
        // distributing them to all the processors.
        //
        // These arrays represent the global matrix data as a CSR
        // matrix (with numEntriesPerRow as redundant but convenient
        // metadata, since it's computable from rowPtr and vice
        // versa).  They are valid only on Proc 0.
        ArrayRCP<size_t> numEntriesPerRow;
        ArrayRCP<size_t> rowPtr;
        ArrayRCP<global_ordinal_type> colInd;
        ArrayRCP<scalar_type> values;

        // Proc 0 first merges duplicate entries, and then converts
        // the coordinate-format matrix data to CSR.
        {
          int mergeAndConvertSucceeded = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            try {
              typedef Teuchos::MatrixMarket::Raw::Element<scalar_type,
                global_ordinal_type> element_type;

              // Number of rows in the matrix.  If we are in tolerant
              // mode, we've already synchronized dims with the actual
              // matrix data.  If in strict mode, we should use dims
              // (as read from the file's metadata) rather than the
              // matrix data to determine the dimensions.  (The matrix
              // data will claim fewer rows than the metadata, if one
              // or more rows have no entries stored in the file.)
              const size_type numRows = dims[0];

              // Additively merge duplicate matrix entries.
              pAdder->getAdder()->merge ();

              // Get a temporary const view of the merged matrix entries.
              const std::vector<element_type>& entries =
                pAdder->getAdder()->getEntries();

              // Number of matrix entries (after merging).
              const size_t numEntries = (size_t)entries.size();

              if (debug) {
                cerr << "----- Proc 0: Matrix has numRows=" << numRows
                     << " rows and numEntries=" << numEntries
                     << " entries." << endl;
              }

              // Make space for the CSR matrix data.  Converting to
              // CSR is easier if we fill numEntriesPerRow with zeros
              // at first.
              numEntriesPerRow = arcp<size_t> (numRows);
              std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 0);
              rowPtr = arcp<size_t> (numRows+1);
              std::fill (rowPtr.begin(), rowPtr.end(), 0);
              colInd = arcp<global_ordinal_type> (numEntries);
              values = arcp<scalar_type> (numEntries);

              // Convert from array-of-structs coordinate format to CSR
              // (compressed sparse row) format.
              global_ordinal_type prvRow = 0;
              size_t curPos = 0;
              rowPtr[0] = 0;
              for (curPos = 0; curPos < numEntries; ++curPos) {
                const element_type& curEntry = entries[curPos];
                const global_ordinal_type curRow = curEntry.rowIndex();
                TEUCHOS_TEST_FOR_EXCEPTION(
                  curRow < prvRow, std::logic_error,
                  "Row indices are out of order, even though they are supposed "
                  "to be sorted.  curRow = " << curRow << ", prvRow = "
                  << prvRow << ", at curPos = " << curPos << ".  Please report "
                  "this bug to the Tpetra developers.");
                if (curRow > prvRow) {
                  for (global_ordinal_type r = prvRow+1; r <= curRow; ++r) {
                    rowPtr[r] = curPos;
                  }
                  prvRow = curRow;
                }
                numEntriesPerRow[curRow]++;
                colInd[curPos] = curEntry.colIndex();
                values[curPos] = curEntry.value();
              }
              // rowPtr has one more entry than numEntriesPerRow.  The
              // last entry of rowPtr is the number of entries in
              // colInd and values.
              rowPtr[numRows] = numEntries;
            } // Finished conversion to CSR format
            catch (std::exception& e) {
              mergeAndConvertSucceeded = 0;
              errMsg << "Failed to merge sparse matrix entries and convert to "
                "CSR format: " << e.what();
            }

            if (debug && mergeAndConvertSucceeded) {
              // Number of rows in the matrix.
              const size_type numRows = dims[0];
              const size_type maxToDisplay = 100;

              cerr << "----- Proc 0: numEntriesPerRow[0.."
                   << (numEntriesPerRow.size()-1) << "] ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > 0) {
                if (numRows > maxToDisplay) {
                  for (size_type k = 0; k < 2; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                  cerr << "... ";
                  for (size_type k = numRows-2; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                else {
                  for (size_type k = 0; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                cerr << numEntriesPerRow[numRows-1];
              } // numRows > 0
              cerr << "]" << endl;

              cerr << "----- Proc 0: rowPtr ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > maxToDisplay) {
                for (size_type k = 0; k < 2; ++k) {
                  cerr << rowPtr[k] << " ";
                }
                cerr << "... ";
                for (size_type k = numRows-2; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              else {
                for (size_type k = 0; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              cerr << rowPtr[numRows] << "]" << endl;
            }
          } // if myRank == rootRank
        } // Done converting sparse matrix data to CSR format

        // Now we're done with the Adder, so we can release the
        // reference ("free" it) to save space.  This only actually
        // does anything on Rank 0, since pAdder is null on all the
        // other MPI processes.
        pAdder = null;

        if (debug && myRank == rootRank) {
          cerr << "-- Making range, domain, and row maps" << endl;
        }

        // Make the maps that describe the matrix's range and domain,
        // and the distribution of its rows.  Creating a Map is a
        // collective operation, so we don't have to do a broadcast of
        // a success Boolean.
        RCP<const map_type> pRangeMap = makeRangeMap (pComm, dims[0]);
        RCP<const map_type> pDomainMap =
          makeDomainMap (pRangeMap, dims[0], dims[1]);
        RCP<const map_type> pRowMap = makeRowMap (null, pComm, dims[0]);

        if (debug && myRank == rootRank) {
          cerr << "-- Distributing the matrix data" << endl;
        }

        // Distribute the matrix data.  Each processor has to add the
        // rows that it owns.  If you try to make Proc 0 call
        // insertGlobalValues() for _all_ the rows, not just those it
        // owns, then fillComplete() will compute the number of
        // columns incorrectly.  That's why Proc 0 has to distribute
        // the matrix data and why we make all the processors (not
        // just Proc 0) call insertGlobalValues() on their own data.
        //
        // These arrays represent each processor's part of the matrix
        // data, in "CSR" format (sort of, since the row indices might
        // not be contiguous).
        ArrayRCP<size_t> myNumEntriesPerRow;
        ArrayRCP<size_t> myRowPtr;
        ArrayRCP<global_ordinal_type> myColInd;
        ArrayRCP<scalar_type> myValues;
        // Distribute the matrix data.  This is a collective operation.
        distribute (myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap,
                    numEntriesPerRow, rowPtr, colInd, values, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Inserting matrix entries on each process "
            "and calling fillComplete()" << endl;
        }
        // Each processor inserts its part of the matrix data, and
        // then they all call fillComplete().  This method invalidates
        // the my* distributed matrix data before calling
        // fillComplete(), in order to save space.  In general, we
        // never store more than two copies of the matrix's entries in
        // memory at once, which is no worse than what Tpetra
        // promises.
        Teuchos::RCP<sparse_matrix_type> pMatrix =
          makeMatrix (myNumEntriesPerRow, myRowPtr, myColInd, myValues,
                      pRowMap, pRangeMap, pDomainMap, constructorParams,
                      fillCompleteParams);
        // Only use a reduce-all in debug mode to check if pMatrix is
        // null.  Otherwise, just throw an exception.  We never expect
        // a null pointer here, so we can save a communication.
        if (debug) {
          int localIsNull = pMatrix.is_null () ? 1 : 0;
          int globalIsNull = 0;
          reduceAll (*pComm, Teuchos::REDUCE_MAX, localIsNull, ptr (&globalIsNull));
          TEUCHOS_TEST_FOR_EXCEPTION(globalIsNull != 0, std::logic_error,
            "Reader::makeMatrix() returned a null pointer on at least one "
            "process.  Please report this bug to the Tpetra developers.");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(pMatrix.is_null(), std::logic_error,
            "Reader::makeMatrix() returned a null pointer.  "
            "Please report this bug to the Tpetra developers.");
        }

        // Sanity check for dimensions (read from the Matrix Market
        // data, vs. dimensions reported by the CrsMatrix).
        //
        // Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
        // what one might think it does, so you have to ask the range
        // resp. domain map for the number of rows resp. columns.
        if (extraDebug && debug) {
          const int numProcs = pComm->getSize ();
          const global_size_t globalNumRows =
            pRangeMap->getGlobalNumElements();
          const global_size_t globalNumCols =
            pDomainMap->getGlobalNumElements();
          if (myRank == rootRank) {
            cerr << "-- Matrix is "
                 << globalNumRows << " x " << globalNumCols
                 << " with " << pMatrix->getGlobalNumEntries()
                 << " entries, and index base "
                 << pMatrix->getIndexBase() << "." << endl;
          }
          pComm->barrier ();
          for (int p = 0; p < numProcs; ++p) {
            if (myRank == p) {
              cerr << "-- Proc " << p << " owns "
                   << pMatrix->getLocalNumCols() << " columns, and "
                   << pMatrix->getLocalNumEntries() << " entries." << endl;
            }
            pComm->barrier ();
          }
        } // if (extraDebug && debug)

        if (debug && myRank == rootRank) {
          cerr << "-- Done creating the CrsMatrix from the Matrix Market data"
               << endl;
        }
        return pMatrix;
      }


      /// \brief Read sparse matrix from the given Matrix Market input
      ///   stream, with provided Maps.
      ///
      /// This version of readSparse() requires you to provide a row
      /// Map, domain Map, and range Map.  You may, if you wish,
      /// provide a column Map as well, but this is not required.
      ///
      /// The given input stream need only be readable by Process 0
      /// (with respect to the given Maps' communicator).  The input
      /// stream must contain Matrix Market "coordinate" format sparse
      /// matrix data.  Read that data on Process 0, and distribute it
      /// to all processors.  Return the resulting distributed
      /// CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Process 0 reads
      ///   data from the input stream, but all processes participate
      ///   and wait for the final result.
      ///
      /// \param in [in/out] The input stream from which to read.
      /// \param rowMap [in] The Map over which to distribute rows
      ///   of the sparse matrix.  This must be nonnull.
      /// \param colMap [in/out] If nonnull: the Map over which to
      ///   distribute columns of the sparse matrix.  If null and if
      ///   callFillComplete is true, we create this for you.
      /// \param domainMap [in] The sparse matrix's domain Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param rangeMap [in] The sparse matrix's range Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparse (std::istream& in,
                  const Teuchos::RCP<const map_type>& rowMap,
                  Teuchos::RCP<const map_type>& colMap,
                  const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const bool callFillComplete=true,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::arcp;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::as;
        using Teuchos::broadcast;
        using Teuchos::Comm;
        using Teuchos::null;
        using Teuchos::ptr;
        using Teuchos::RCP;
        using Teuchos::reduceAll;
        using Teuchos::Tuple;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        RCP<const Comm<int> > pComm = rowMap->getComm ();
        const int myRank = pComm->getRank ();
        const int rootRank = 0;
        const bool extraDebug = false;

        // Fast checks for invalid input.  We can't check other
        // attributes of the Maps until we've read in the matrix
        // dimensions.
        TEUCHOS_TEST_FOR_EXCEPTION(
          rowMap.is_null (), std::invalid_argument,
          "Row Map must be nonnull.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          rangeMap.is_null (), std::invalid_argument,
          "Range Map must be nonnull.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          domainMap.is_null (), std::invalid_argument,
          "Domain Map must be nonnull.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          rowMap->getComm().getRawPtr() != pComm.getRawPtr(),
          std::invalid_argument,
          "The specified row Map's communicator (rowMap->getComm())"
          "differs from the given separately supplied communicator pComm.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          domainMap->getComm().getRawPtr() != pComm.getRawPtr(),
          std::invalid_argument,
          "The specified domain Map's communicator (domainMap->getComm())"
          "differs from the given separately supplied communicator pComm.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          rangeMap->getComm().getRawPtr() != pComm.getRawPtr(),
          std::invalid_argument,
          "The specified range Map's communicator (rangeMap->getComm())"
          "differs from the given separately supplied communicator pComm.");

        // Current line number in the input stream.  Various calls
        // will modify this depending on the number of lines that are
        // read from the input stream.  Only Rank 0 modifies this.
        size_t lineNumber = 1;

        if (debug && myRank == rootRank) {
          cerr << "Matrix Market reader: readSparse:" << endl
               << "-- Reading banner line" << endl;
        }

        // The "Banner" tells you whether the input stream represents
        // a sparse matrix, the symmetry type of the matrix, and the
        // type of the data it contains.
        //
        // pBanner will only be nonnull on MPI Rank 0.  It will be
        // null on all other MPI processes.
        RCP<const Banner> pBanner;
        {
          // We read and validate the Banner on Proc 0, but broadcast
          // the validation result to all processes.
          // Teuchos::broadcast doesn't currently work with bool, so
          // we use int (true -> 1, false -> 0).
          int bannerIsCorrect = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            // Read the Banner line from the input stream.
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }

            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse matrix.
              // We validate on Proc 0, since it reads the Banner.

              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                  "\"coordinate\"-format sparse matrix in order to create a "
                  "Tpetra::CrsMatrix object from it, but the file's matrix "
                  "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                  "format sparse matrix in order to create a Tpetra::CrsMatrix "
                  "object from it, but the file's matrix type is \"array\" "
                  "instead.  That probably means the file contains dense matrix "
                  "data.";
              }
            }
          } // Proc 0: Done reading the Banner, hopefully successfully.

          // Broadcast from Proc 0 whether the Banner was read correctly.
          broadcast (*pComm, rootRank, ptr (&bannerIsCorrect));

          // If the Banner is invalid, all processes throw an
          // exception.  Only Proc 0 gets the exception message, but
          // that's OK, since the main point is to "stop the world"
          // (rather than throw an exception on one process and leave
          // the others hanging).
          TEUCHOS_TEST_FOR_EXCEPTION(bannerIsCorrect == 0,
            std::invalid_argument, errMsg.str ());
        } // Done reading the Banner line and broadcasting success.
        if (debug && myRank == rootRank) {
          cerr << "-- Reading dimensions line" << endl;
        }

        // Read the matrix dimensions from the Matrix Market metadata.
        // dims = (numRows, numCols, numEntries).  Proc 0 does the
        // reading, but it broadcasts the results to all MPI
        // processes.  Thus, readCoordDims() is a collective
        // operation.  It does a collective check for correctness too.
        Tuple<global_ordinal_type, 3> dims =
          readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Making Adder for collecting matrix data" << endl;
        }

        // "Adder" object for collecting all the sparse matrix entries
        // from the input stream.  This is only nonnull on Proc 0.
        // The Adder internally converts the one-based indices (native
        // Matrix Market format) into zero-based indices.
        RCP<adder_type> pAdder =
          makeAdder (pComm, pBanner, dims, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Reading matrix data" << endl;
        }
        //
        // Read the matrix entries from the input stream on Proc 0.
        //
        {
          // We use readSuccess to broadcast the results of the read
          // (succeeded or not) to all MPI processes.  Since
          // Teuchos::broadcast doesn't currently know how to send
          // bools, we convert to int (true -> 1, false -> 0).
          int readSuccess = 1;
          std::ostringstream errMsg; // Exception message (only valid on Proc 0)
          if (myRank == rootRank) {
            try {
              // Reader for "coordinate" format sparse matrix data.
              typedef Teuchos::MatrixMarket::CoordDataReader<adder_type,
                global_ordinal_type, scalar_type, STS::isComplex> reader_type;
              reader_type reader (pAdder);

              // Read the sparse matrix entries.
              std::pair<bool, std::vector<size_t> > results =
                reader.read (in, lineNumber, tolerant, debug);
              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }
          }
          broadcast (*pComm, rootRank, ptr (&readSuccess));

          // It would be nice to add a "verbose" flag, so that in
          // tolerant mode, we could log any bad line number(s) on
          // Proc 0.  For now, we just throw if the read fails to
          // succeed.
          //
          // Question: If we're in tolerant mode, and if the read did
          // not succeed, should we attempt to call fillComplete()?
          TEUCHOS_TEST_FOR_EXCEPTION(readSuccess == 0, std::runtime_error,
            "Failed to read the Matrix Market sparse matrix file: "
            << errMsg.str());
        } // Done reading the matrix entries (stored on Proc 0 for now)

        if (debug && myRank == rootRank) {
          cerr << "-- Successfully read the Matrix Market data" << endl;
        }

        // In tolerant mode, we need to rebroadcast the matrix
        // dimensions, since they may be different after reading the
        // actual matrix data.  We only need to broadcast the number
        // of rows and columns.  Only Rank 0 needs to know the actual
        // global number of entries, since (a) we need to merge
        // duplicates on Rank 0 first anyway, and (b) when we
        // distribute the entries, each rank other than Rank 0 will
        // only need to know how many entries it owns, not the total
        // number of entries.
        if (tolerant) {
          if (debug && myRank == rootRank) {
            cerr << "-- Tolerant mode: rebroadcasting matrix dimensions"
                 << endl
                 << "----- Dimensions before: "
                 << dims[0] << " x " << dims[1]
                 << endl;
          }
          // Packed coordinate matrix dimensions (numRows, numCols).
          Tuple<global_ordinal_type, 2> updatedDims;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number
            // of rows is less than that specified in the
            // metadata.  We allow this case, and favor the
            // metadata so that the zero row(s) will be included.
            updatedDims[0] =
              std::max (dims[0], pAdder->getAdder()->numRows());
            updatedDims[1] = pAdder->getAdder()->numCols();
          }
          broadcast (*pComm, rootRank, updatedDims);
          dims[0] = updatedDims[0];
          dims[1] = updatedDims[1];
          if (debug && myRank == rootRank) {
            cerr << "----- Dimensions after: " << dims[0] << " x "
                 << dims[1] << endl;
          }
        }
        else {
          // In strict mode, we require that the matrix's metadata and
          // its actual data agree, at least somewhat.  In particular,
          // the number of rows must agree, since otherwise we cannot
          // distribute the matrix correctly.

          // Teuchos::broadcast() doesn't know how to broadcast bools,
          // so we use an int with the standard 1 == true, 0 == false
          // encoding.
          int dimsMatch = 1;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number of
            // rows is less than that specified in the metadata.  We
            // allow this case, and favor the metadata, but do not
            // allow the Adder to think there are more rows in the
            // matrix than the metadata says.
            if (dims[0] < pAdder->getAdder ()->numRows ()) {
              dimsMatch = 0;
            }
          }
          broadcast (*pComm, 0, ptr (&dimsMatch));
          if (dimsMatch == 0) {
            // We're in an error state anyway, so we might as well
            // work a little harder to print an informative error
            // message.
            //
            // Broadcast the Adder's idea of the matrix dimensions
            // from Proc 0 to all processes.
            Tuple<global_ordinal_type, 2> addersDims;
            if (myRank == rootRank) {
              addersDims[0] = pAdder->getAdder()->numRows();
              addersDims[1] = pAdder->getAdder()->numCols();
            }
            broadcast (*pComm, 0, addersDims);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dimsMatch == 0, std::runtime_error,
              "The matrix metadata says that the matrix is " << dims[0] << " x "
              << dims[1] << ", but the actual data says that the matrix is "
              << addersDims[0] << " x " << addersDims[1] << ".  That means the "
              "data includes more rows than reported in the metadata.  This "
              "is not allowed when parsing in strict mode.  Parse the matrix in "
              "tolerant mode to ignore the metadata when it disagrees with the "
              "data.");
          }
        } // Matrix dimensions (# rows, # cols, # entries) agree.

        if (debug && myRank == rootRank) {
          cerr << "-- Converting matrix data into CSR format on Proc 0" << endl;
        }

        // Now that we've read in all the matrix entries from the
        // input stream into the adder on Proc 0, post-process them
        // into CSR format (still on Proc 0).  This will facilitate
        // distributing them to all the processors.
        //
        // These arrays represent the global matrix data as a CSR
        // matrix (with numEntriesPerRow as redundant but convenient
        // metadata, since it's computable from rowPtr and vice
        // versa).  They are valid only on Proc 0.
        ArrayRCP<size_t> numEntriesPerRow;
        ArrayRCP<size_t> rowPtr;
        ArrayRCP<global_ordinal_type> colInd;
        ArrayRCP<scalar_type> values;

        // Proc 0 first merges duplicate entries, and then converts
        // the coordinate-format matrix data to CSR.
        {
          int mergeAndConvertSucceeded = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            try {
              typedef Teuchos::MatrixMarket::Raw::Element<scalar_type,
                global_ordinal_type> element_type;

              // Number of rows in the matrix.  If we are in tolerant
              // mode, we've already synchronized dims with the actual
              // matrix data.  If in strict mode, we should use dims
              // (as read from the file's metadata) rather than the
              // matrix data to determine the dimensions.  (The matrix
              // data will claim fewer rows than the metadata, if one
              // or more rows have no entries stored in the file.)
              const size_type numRows = dims[0];

              // Additively merge duplicate matrix entries.
              pAdder->getAdder()->merge ();

              // Get a temporary const view of the merged matrix entries.
              const std::vector<element_type>& entries =
                pAdder->getAdder()->getEntries();

              // Number of matrix entries (after merging).
              const size_t numEntries = (size_t)entries.size();

              if (debug) {
                cerr << "----- Proc 0: Matrix has numRows=" << numRows
                     << " rows and numEntries=" << numEntries
                     << " entries." << endl;
              }

              // Make space for the CSR matrix data.  Converting to
              // CSR is easier if we fill numEntriesPerRow with zeros
              // at first.
              numEntriesPerRow = arcp<size_t> (numRows);
              std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 0);
              rowPtr = arcp<size_t> (numRows+1);
              std::fill (rowPtr.begin(), rowPtr.end(), 0);
              colInd = arcp<global_ordinal_type> (numEntries);
              values = arcp<scalar_type> (numEntries);

              // Convert from array-of-structs coordinate format to CSR
              // (compressed sparse row) format.
              global_ordinal_type prvRow = 0;
              size_t curPos = 0;
              rowPtr[0] = 0;
              for (curPos = 0; curPos < numEntries; ++curPos) {
                const element_type& curEntry = entries[curPos];
                const global_ordinal_type curRow = curEntry.rowIndex();
                TEUCHOS_TEST_FOR_EXCEPTION(
                  curRow < prvRow, std::logic_error,
                  "Row indices are out of order, even though they are supposed "
                  "to be sorted.  curRow = " << curRow << ", prvRow = "
                  << prvRow << ", at curPos = " << curPos << ".  Please report "
                  "this bug to the Tpetra developers.");
                if (curRow > prvRow) {
                  prvRow = curRow;
                }
                numEntriesPerRow[curRow]++;
                colInd[curPos] = curEntry.colIndex();
                values[curPos] = curEntry.value();
              }

              rowPtr[0] = 0;
              for (global_ordinal_type row = 1; row <= numRows; ++row) {
                rowPtr[row] = numEntriesPerRow[row-1] + rowPtr[row-1];
              }
            } // Finished conversion to CSR format
            catch (std::exception& e) {
              mergeAndConvertSucceeded = 0;
              errMsg << "Failed to merge sparse matrix entries and convert to "
                "CSR format: " << e.what();
            }

            if (debug && mergeAndConvertSucceeded) {
              // Number of rows in the matrix.
              const size_type numRows = dims[0];
              const size_type maxToDisplay = 100;

              cerr << "----- Proc 0: numEntriesPerRow[0.."
                   << (numEntriesPerRow.size()-1) << "] ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > 0) {
                if (numRows > maxToDisplay) {
                  for (size_type k = 0; k < 2; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                  cerr << "... ";
                  for (size_type k = numRows-2; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                else {
                  for (size_type k = 0; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                cerr << numEntriesPerRow[numRows-1];
              } // numRows > 0
              cerr << "]" << endl;

              cerr << "----- Proc 0: rowPtr ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > maxToDisplay) {
                for (size_type k = 0; k < 2; ++k) {
                  cerr << rowPtr[k] << " ";
                }
                cerr << "... ";
                for (size_type k = numRows-2; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              else {
                for (size_type k = 0; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              cerr << rowPtr[numRows] << "]" << endl;

              cerr << "----- Proc 0: colInd = [";
              for (size_t k = 0; k < rowPtr[numRows]; ++k) {
                cerr << colInd[k] << " ";
              }
              cerr << "]" << endl;
            }
          } // if myRank == rootRank
        } // Done converting sparse matrix data to CSR format

        // Now we're done with the Adder, so we can release the
        // reference ("free" it) to save space.  This only actually
        // does anything on Rank 0, since pAdder is null on all the
        // other MPI processes.
        pAdder = null;

        // Verify details of the Maps.  Don't count the global number
        // of entries in the row Map, since that number doesn't
        // correctly count overlap.
        if (debug && myRank == rootRank) {
          cerr << "-- Verifying Maps" << endl;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          as<global_size_t> (dims[0]) != rangeMap->getGlobalNumElements(),
          std::invalid_argument,
          "The range Map has " << rangeMap->getGlobalNumElements ()
          << " entries, but the matrix has a global number of rows " << dims[0]
          << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(
          as<global_size_t> (dims[1]) != domainMap->getGlobalNumElements (),
          std::invalid_argument,
          "The domain Map has " << domainMap->getGlobalNumElements ()
          << " entries, but the matrix has a global number of columns "
          << dims[1] << ".");

        // Create a row Map which is entirely owned on Proc 0.
        RCP<Teuchos::FancyOStream> err = debug ?
          Teuchos::getFancyOStream (Teuchos::rcpFromRef (cerr)) : null;

        RCP<const map_type> gatherRowMap = Details::computeGatherMap (rowMap, err, debug);
        ArrayView<const global_ordinal_type> myRows =
            gatherRowMap->getLocalElementList ();
        const size_type myNumRows = myRows.size ();
        const global_ordinal_type indexBase = gatherRowMap->getIndexBase ();

        ArrayRCP<size_t> gatherNumEntriesPerRow = arcp<size_t>(myNumRows);
        for (size_type i_ = 0; i_ < myNumRows; i_++) {
          gatherNumEntriesPerRow[i_] = numEntriesPerRow[myRows[i_]-indexBase];
        }

        // Create a matrix using this Map, and fill in on Proc 0.  We
        // know how many entries there are in each row, so we can use
        // static profile.
        RCP<sparse_matrix_type> A_proc0 =
          rcp (new sparse_matrix_type (gatherRowMap, gatherNumEntriesPerRow()));
        if (myRank == rootRank) {
          if (debug) {
            cerr << "-- Proc 0: Filling gather matrix" << endl;
          }
          if (debug) {
            cerr << "---- Rows: " << Teuchos::toString (myRows) << endl;
          }

          // Add Proc 0's matrix entries to the CrsMatrix.
          for (size_type i_ = 0; i_ < myNumRows; ++i_) {
            size_type i = myRows[i_] - indexBase;

            const size_type curPos = as<size_type> (rowPtr[i]);
            const local_ordinal_type curNumEntries = numEntriesPerRow[i];
            ArrayView<global_ordinal_type> curColInd =
              colInd.view (curPos, curNumEntries);
            ArrayView<scalar_type> curValues =
              values.view (curPos, curNumEntries);

            // Modify the column indices in place to have the right index base.
            for (size_type k = 0; k < curNumEntries; ++k) {
              curColInd[k] += indexBase;
            }
            if (debug) {
              cerr << "------ Columns: " << Teuchos::toString (curColInd) << endl;
              cerr << "------ Values: " << Teuchos::toString (curValues) << endl;
            }
            // Avoid constructing empty views of ArrayRCP objects.
            if (curNumEntries > 0) {
              A_proc0->insertGlobalValues (myRows[i_], curColInd, curValues);
            }
          }
          // Now we can save space by deallocating numEntriesPerRow,
          // rowPtr, colInd, and values, since we've already put those
          // data in the matrix.
          numEntriesPerRow = null;
          rowPtr = null;
          colInd = null;
          values = null;
        } // if myRank == rootRank

        RCP<sparse_matrix_type> A;
        typedef Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
        export_type exp (gatherRowMap, rowMap);

        // Communicate the precise number of nonzeros per row, which was already
        // calculated above.
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef Tpetra::MultiVector<GO, LO, GO, node_type> mv_type_go;
        mv_type_go target_nnzPerRow(rowMap,1);
        mv_type_go source_nnzPerRow(gatherRowMap,1);
        Teuchos::ArrayRCP<GO> srcData = source_nnzPerRow.getDataNonConst(0);
        for (int i=0; i<myNumRows; i++)
          srcData[i] = gatherNumEntriesPerRow[i];
        srcData = Teuchos::null;
        target_nnzPerRow.doExport(source_nnzPerRow,exp,Tpetra::INSERT);
        Teuchos::ArrayRCP<GO> targetData = target_nnzPerRow.getDataNonConst(0);
        ArrayRCP<size_t> targetData_size_t = arcp<size_t>(targetData.size());
        for (int i=0; i<targetData.size(); i++)
          targetData_size_t[i] = targetData[i];

        if (colMap.is_null ()) {
          A = rcp (new sparse_matrix_type (rowMap, targetData_size_t()));
        } else {
          A = rcp (new sparse_matrix_type (rowMap, colMap, targetData_size_t()));
        }
        A->doExport (*A_proc0, exp, INSERT);
        if (callFillComplete) {
          A->fillComplete (domainMap, rangeMap);
        }

        // We can't get the dimensions of the matrix until after
        // fillComplete() is called.  Thus, we can't do the sanity
        // check (dimensions read from the Matrix Market data,
        // vs. dimensions reported by the CrsMatrix) unless the user
        // asked us to call fillComplete().
        //
        // Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
        // what one might think it does, so you have to ask the range
        // resp. domain map for the number of rows resp. columns.
        if (callFillComplete) {
          const int numProcs = pComm->getSize ();

          if (extraDebug && debug) {
            const global_size_t globalNumRows = rangeMap->getGlobalNumElements ();
            const global_size_t globalNumCols = domainMap->getGlobalNumElements ();
            if (myRank == rootRank) {
              cerr << "-- Matrix is "
                   << globalNumRows << " x " << globalNumCols
                   << " with " << A->getGlobalNumEntries()
                   << " entries, and index base "
                   << A->getIndexBase() << "." << endl;
            }
            pComm->barrier ();
            for (int p = 0; p < numProcs; ++p) {
              if (myRank == p) {
                cerr << "-- Proc " << p << " owns "
                     << A->getLocalNumCols() << " columns, and "
                     << A->getLocalNumEntries() << " entries." << endl;
              }
              pComm->barrier ();
            }
          } // if (extraDebug && debug)
        } // if (callFillComplete)

        if (debug && myRank == rootRank) {
          cerr << "-- Done creating the CrsMatrix from the Matrix Market data"
               << endl;
        }
        return A;
      }

      /// \brief Read dense matrix (as a MultiVector) from the given
      ///   Matrix Market file.
      ///
      /// Open the given file on MPI Process 0 (with respect to the
      /// given communicator).  The file should contain Matrix Market
      /// "array" format dense matrix data.  Read that data on Process
      /// 0, and distribute it to all processors.  Return the
      /// resulting distributed MultiVector.
      ///
      /// See documentation of readDense() for details.
      ///
      /// \param filename [in] Name of the Matrix Market file from
      ///   which to read.  Both the filename and the file itself are
      ///   only accessed on Rank 0 of the given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the dense matrix will be distributed.
      /// \param map [in/out] On input: if nonnull, the map describing
      ///   how to distribute the vector (not modified).  In this
      ///   case, the map's communicator and node must equal \c comm
      ///   resp. \c node.  If null on input, then on output, a
      ///   sensible (contiguous and uniformly distributed over the
      ///   given communicator) map describing the distribution of the
      ///   output multivector.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      /// \param debug [in] If true, read in binary mode.
      static Teuchos::RCP<multivector_type>
      readDenseFile (const std::string& filename,
                     const trcp_tcomm_t& comm,
                     Teuchos::RCP<const map_type>& map,
                     const bool tolerant=false,
                     const bool debug=false,
                     const bool binary=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true, binary ? std::ios::binary : std::ios::in);

        return readDense (in, comm, map, tolerant, debug, binary);
      }

      /**
       * @brief Open an input file stream safely on rank zero.
       *
       * Only open the file on rank zero process. Test carefully to make
       * sure that the file opened successfully and broadcast that
       * result to all processes to prevent a hang on exception
       * throw.
       *
       * @note On processes that are not the rank zero process, the stream is left uninitialized.
       */
      static std::ifstream openInFileOnRankZero(
        const trcp_tcomm_t& comm,
        const std::string& filename, const bool safe = true,
        std::ios_base::openmode mode = std::ios_base::in
      ){
        // Input stream.
        std::ifstream in;

        // Placeholder for broadcasting in-stream state.
        int all_should_stop = false;

        // Try to open the in-stream on root rank.
        if (comm->getRank() == 0) {
            in.open(filename, mode);
            all_should_stop = !in && safe;
        }

        // Broadcast state and possibly throw.
        if(comm) Teuchos::broadcast(*comm, 0, &all_should_stop);

        TEUCHOS_TEST_FOR_EXCEPTION(
          all_should_stop,
          std::runtime_error,
          "Could not open input file '" + filename + "' on root rank 0."
        );

        return in;
      }


      /// \brief Read a Vector from the given Matrix Market file.
      ///
      /// Same as readDenseFile() but will work with 1D data only.
      ///
      /// Open the given file on MPI Process 0 (with respect to the
      /// given communicator).  The file should contain Matrix Market
      /// "array" format dense matrix data with number of columns==1.
      /// Read that data on Process 0, and distribute it to all processors.
      /// Return the resulting distributed Vector.
      ///
      /// See documentation of readVector() for details.
      ///
      /// \param filename [in] Name of the Matrix Market file from
      ///   which to read.  Both the filename and the file itself are
      ///   only accessed on Rank 0 of the given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the dense matrix will be distributed.
      /// \param map [in/out] On input: if nonnull, the map describing
      ///   how to distribute the vector (not modified).  In this
      ///   case, the map's communicator and node must equal \c comm
      ///   resp. \c node.  If null on input, then on output, a
      ///   sensible (contiguous and uniformly distributed over the
      ///   given communicator) map describing the distribution of the
      ///   output multivector.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<vector_type>
      readVectorFile (const std::string& filename,
                      const trcp_tcomm_t& comm,
                      Teuchos::RCP<const map_type>& map,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true);

        return readVector (in, comm, map, tolerant, debug);
      }


      /// \brief Read dense matrix (as a MultiVector) from the given
      ///   Matrix Market input stream.
      ///
      /// The given input stream need only be readable by MPI Rank 0
      /// (with respect to the given communicator).  The input stream
      /// should contain Matrix Market "array" format dense matrix
      /// data.  Read that data on Process 0, and distribute it to all
      /// processes in the communicator.  Return the resulting
      /// distributed MultiVector.
      ///
      /// Unlike readSparse(), this method allows callers to supply a
      /// Map over which to distribute the resulting MultiVector.  The
      /// Map argument is optional; if null, we construct our own
      /// reasonable Map.  We let users supply their own Map, because
      /// a common case in Tpetra is to read in or construct a sparse
      /// matrix first, and then create dense (multi)vectors
      /// distributed with the sparse matrix's domain or range Map.
      ///
      /// \note This is a collective operation.  Only Process 0 in the
      ///   communicator opens the file and reads data from it, but
      ///   all processes in the communicator participate and wait for
      ///   the final result.
      ///
      /// \note "Tolerant" parsing mode means something different for
      ///   dense matrices than it does for sparse matrices.  Since
      ///   Matrix Market dense matrix files don't store indices with
      ///   each value, unlike sparse matrices, we can't determine the
      ///   matrix dimensions from the data alone.  Thus, we require
      ///   the metadata to include a valid number of rows and
      ///   columns.  "Tolerant" in the dense case refers to the data;
      ///   in tolerant mode, the number of stored matrix entries may
      ///   be more or less than the number reported by the metadata
      ///   (number of rows times number of columns).  If more, the
      ///   extra data are ignored; if less, the remainder is filled
      ///   in with zeros.
      ///
      /// \note On Map compatibility: Suppose that you write a
      ///   multivector X to a file.  Then, you read it back in as a
      ///   different multivector Y distributed over the same
      ///   communicator, but with a Map constructed by the input
      ///   routine (i.e., a null Map on input to readDenseFile() or
      ///   readDense()).  In that case, the only properties shared by
      ///   the maps of X and Y are that they have the same
      ///   communicator and the same number of GIDs.  The Maps need
      ///   not necessarily be compatible (in the sense of
      ///   Map::isCompatible()), and they certainly need not
      ///   necessarily be the same Map (in the sense of
      ///   Map::isSameAs()).
      ///
      /// \param in [in] The input stream from which to read.  The
      ///   stream is only accessed on Process 0 of the given
      ///   communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the dense matrix will be distributed.
      /// \param map [in/out] On input: if nonnull, the map describing
      ///   how to distribute the vector (not modified).  In this
      ///   case, the map's communicator and node must equal comm
      ///   resp. node.  If null on input, then on output, a sensible
      ///   (contiguous and uniformly distributed over the given
      ///   communicator) map describing the distribution of the
      ///   output multivector.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the stream.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      /// \param debug [in] If true, read in binary mode.
      
      static Teuchos::RCP<multivector_type>
      readDense (std::istream& in,
                 const trcp_tcomm_t& comm,
                 Teuchos::RCP<const map_type>& map,
                 const bool tolerant=false,
                 const bool debug=false,
                 const bool binary=false)
      {
        Teuchos::RCP<Teuchos::FancyOStream> err =
          Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
        return readDenseImpl<scalar_type> (in, comm, map, err, tolerant, debug, binary);
      }


      //! Read Vector from the given Matrix Market input stream.
      static Teuchos::RCP<vector_type>
      readVector (std::istream& in,
                  const trcp_tcomm_t& comm,
                  Teuchos::RCP<const map_type>& map,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        Teuchos::RCP<Teuchos::FancyOStream> err =
          Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
        return readVectorImpl<scalar_type> (in, comm, map, err, tolerant, debug);
      }


      /// \brief Read Map (as a MultiVector) from the given
      ///   Matrix Market file.
      ///
      /// Open the given file on MPI Process 0 (with respect to the
      /// given communicator).  The file should contain Matrix Market
      /// "array" format dense matrix data with two columns, as
      /// generated by Writer::writeMap() or Writer::writeMapFile().
      /// Read that data on Process 0, and distribute it to all
      /// processes.  Return the resulting Map.
      ///
      /// \param filename [in] Name of the Matrix Market file from
      ///   which to read.  Both the filename and the file itself are
      ///   only accessed on Process 0 of the given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the Map will be distributed.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      /// \param debug [in] If true, read in binary mode.
      static Teuchos::RCP<const map_type>
      readMapFile (const std::string& filename,
                   const trcp_tcomm_t& comm,
                   const bool tolerant=false,
                   const bool debug=false,
                   const bool binary=false)
      {
        std::ifstream in = Reader::openInFileOnRankZero(comm, filename, true, binary ? std::ios::binary : std::ios::in);

        return readMap (in, comm, tolerant, debug, binary);
      }


    private:
      template<class MultiVectorScalarType>
      static Teuchos::RCP<Tpetra::MultiVector<MultiVectorScalarType,
                                     local_ordinal_type,
                                     global_ordinal_type,
                                     node_type> >
      readDenseImpl (std::istream& in,
                     const trcp_tcomm_t& comm,
                     Teuchos::RCP<const map_type>& map,
                     const Teuchos::RCP<Teuchos::FancyOStream>& err,
                     const bool tolerant=false,
                     const bool debug=false,
                     const bool binary=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::MatrixMarket::checkCommentLine;
        using Teuchos::ArrayRCP;
        using Teuchos::as;
        using Teuchos::broadcast;
        using Teuchos::outArg;
        using Teuchos::RCP;
        using Teuchos::Tuple;
        using std::endl;
        typedef MultiVectorScalarType ST;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;
        typedef Teuchos::ScalarTraits<ST> STS;
        typedef typename STS::magnitudeType MT;
        typedef Teuchos::ScalarTraits<MT> STM;
        typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;

        // Rank 0 is the only (MPI) process allowed to read from the
        // input stream.
        const int myRank = comm->getRank ();

        if (! err.is_null ()) {
          err->pushTab ();
        }
        if (debug) {
          *err << myRank << ": readDenseImpl" << endl;
        }
        if (! err.is_null ()) {
          err->pushTab ();
        }

        // mfh 17 Feb 2013: It's not strictly necessary that the Comm
        // instances be identical and that the Node instances be
        // identical.  The essential condition is more complicated to
        // test and isn't the same for all Node types.  Thus, we just
        // leave it up to the user.

        // // If map is nonnull, check the precondition that its
        // // communicator resp. node equal comm resp. node.  Checking
        // // now avoids doing a lot of file reading before we detect the
        // // violated precondition.
        // TEUCHOS_TEST_FOR_EXCEPTION(
        //   ! map.is_null() && (map->getComm() != comm || map->getNode () != node,
        //   std::invalid_argument, "If you supply a nonnull Map, the Map's "
        //   "communicator and node must equal the supplied communicator resp. "
        //   "node.");

        // Process 0 will read in the matrix dimensions from the file,
        // and broadcast them to all ranks in the given communicator.
        // There are only 2 dimensions in the matrix, but we use the
        // third element of the Tuple to encode the banner's reported
        // data type: "real" == 0, "complex" == 1, and "integer" == 0
        // (same as "real").  We don't allow pattern matrices (i.e.,
        // graphs) since they only make sense for sparse data.
        Tuple<GO, 3> dims;
        dims[0] = 0;
        dims[1] = 0;
        dims[2] = 0;

        // Current line number in the input stream.  Only valid on
        // Proc 0.  Various calls will modify this depending on the
        // number of lines that are read from the input stream.
        size_t lineNumber = 1;

        // Capture errors and their messages on Proc 0.
        std::ostringstream exMsg;
        int localBannerReadSuccess = 1;
        int localDimsReadSuccess = 1;

        // Only Proc 0 gets to read matrix data from the input stream.
        if (myRank == 0) {
          if (! binary) {
            if (debug) {
              *err << myRank << ": readDenseImpl: Reading banner line (dense)" << endl;
            }

            // The "Banner" tells you whether the input stream
            // represents a dense matrix, the symmetry type of the
            // matrix, and the type of the data it contains.
            RCP<const Banner> pBanner;
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            } catch (std::exception& e) {
              exMsg << e.what ();
              localBannerReadSuccess = 0;
            }
            // Make sure the input stream is the right kind of data.
            if (localBannerReadSuccess) {
              if (pBanner->matrixType () != "array") {
                exMsg << "The Matrix Market file does not contain dense matrix "
                  "data.  Its banner (first) line says that its matrix type is \""
                      << pBanner->matrixType () << "\", rather that the required "
                  "\"array\".";
                localBannerReadSuccess = 0;
              } else if (pBanner->dataType() == "pattern") {
                exMsg << "The Matrix Market file's banner (first) "
                  "line claims that the matrix's data type is \"pattern\".  This does "
                  "not make sense for a dense matrix, yet the file reports the matrix "
                  "as dense.  The only valid data types for a dense matrix are "
                  "\"real\", \"complex\", and \"integer\".";
                localBannerReadSuccess = 0;
              } else {
                // Encode the data type reported by the Banner as the
                // third element of the dimensions Tuple.
                dims[2] = encodeDataType (pBanner->dataType ());
              }
            } // if we successfully read the banner line

            // At this point, we've successfully read the banner line.
            // Now read the dimensions line.
            if (localBannerReadSuccess) {
              if (debug) {
                *err << myRank << ": readDenseImpl: Reading dimensions line (dense)" << endl;
              }
              // Keep reading lines from the input stream until we find
              // a non-comment line, or until we run out of lines.  The
              // latter is an error, since every "array" format Matrix
              // Market file must have a dimensions line after the
              // banner (even if the matrix has zero rows or columns, or
              // zero entries).
              std::string line;
              bool commentLine = true;

              while (commentLine) {
                // Test whether it is even valid to read from the input
                // stream wrapping the line.
                if (in.eof () || in.fail ()) {
                  exMsg << "Unable to get array dimensions line (at all) from line "
                        << lineNumber << " of input stream.  The input stream "
                        << "claims that it is "
                        << (in.eof() ? "at end-of-file." : "in a failed state.");
                  localDimsReadSuccess = 0;
                } else {
                  // Try to get the next line from the input stream.
                  if (getline (in, line)) {
                    ++lineNumber; // We did actually read a line.
                  }
                  // Is the current line a comment line?  Ignore start
                  // and size; they are only useful for reading the
                  // actual matrix entries.  (We could use them here as
                  // an optimization, but we've chosen not to.)
                  size_t start = 0, size = 0;
                  commentLine = checkCommentLine (line, start, size, lineNumber, tolerant);
                } // whether we failed to read the line at all
              } // while the line we just read is a comment line

              //
              // Get <numRows> <numCols> from the line we just read.
              //
              std::istringstream istr (line);

              // Test whether it is even valid to read from the input
              // stream wrapping the line.
              if (istr.eof () || istr.fail ()) {
                exMsg << "Unable to read any data from line " << lineNumber
                      << " of input; the line should contain the matrix dimensions "
                      << "\"<numRows> <numCols>\".";
                localDimsReadSuccess = 0;
              } else { // It's valid to read from the line.
                GO theNumRows = 0;
                istr >> theNumRows; // Read in the number of rows.
                if (istr.fail ()) {
                  exMsg << "Failed to get number of rows from line "
                        << lineNumber << " of input; the line should contains the "
                        << "matrix dimensions \"<numRows> <numCols>\".";
                  localDimsReadSuccess = 0;
                } else { // We successfully read the number of rows
                  dims[0] = theNumRows; // Save the number of rows
                  if (istr.eof ()) { // Do we still have data to read?
                    exMsg << "No more data after number of rows on line "
                          << lineNumber << " of input; the line should contain the "
                          << "matrix dimensions \"<numRows> <numCols>\".";
                    localDimsReadSuccess = 0;
                  } else { // Still data left to read; read in number of columns.
                    GO theNumCols = 0;
                    istr >> theNumCols; // Read in the number of columns
                    if (istr.fail ()) {
                      exMsg << "Failed to get number of columns from line "
                            << lineNumber << " of input; the line should contain "
                            << "the matrix dimensions \"<numRows> <numCols>\".";
                      localDimsReadSuccess = 0;
                    } else { // We successfully read the number of columns
                      dims[1] = theNumCols; // Save the number of columns
                    } // if istr.fail ()
                  } // if istr.eof ()
                } // if we read the number of rows
              } // if the input stream wrapping the dims line was (in)valid
            } // if we successfully read the banner line
          } // not in binary format
          else { // in binary format
            global_size_t numRows, numCols;
            in.read(reinterpret_cast<char*>(&numRows), sizeof(numRows));
            in.read(reinterpret_cast<char*>(&numCols), sizeof(numCols));
            dims[0] = Teuchos::as<GO>(numRows);
            dims[1] = Teuchos::as<GO>(numCols);
            if ((typeid(ST) == typeid(double)) || Teuchos::ScalarTraits<ST>::isOrdinal) {
              dims[2] = 0;
            } else if (Teuchos::ScalarTraits<ST>::isComplex) {
              dims[2] = 1;
            } else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                         "Unrecognized Matrix Market data type. "
                                         "We should never get here.  "
                                         "Please report this bug to the Tpetra developers.");
            }
            localDimsReadSuccess = true;
            localDimsReadSuccess = true;
          } // in binary format
        } // if (myRank == 0)

        // Broadcast the matrix dimensions, the encoded data type, and
        // whether or not Proc 0 succeeded in reading the banner and
        // dimensions.
        Tuple<GO, 5> bannerDimsReadResult;
        if (myRank == 0) {
          bannerDimsReadResult[0] = dims[0]; // numRows
          bannerDimsReadResult[1] = dims[1]; // numCols
          bannerDimsReadResult[2] = dims[2]; // encoded data type
          bannerDimsReadResult[3] = localBannerReadSuccess;
          bannerDimsReadResult[4] = localDimsReadSuccess;
        }
        // Broadcast matrix dimensions and the encoded data type from
        // Proc 0 to all the MPI processes.
        broadcast (*comm, 0, bannerDimsReadResult);

        TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[3] == 0, std::runtime_error,
          "Failed to read banner line: " << exMsg.str ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[4] == 0, std::runtime_error,
          "Failed to read matrix dimensions line: " << exMsg.str ());
        if (myRank != 0) {
          dims[0] = bannerDimsReadResult[0];
          dims[1] = bannerDimsReadResult[1];
          dims[2] = bannerDimsReadResult[2];
        }

        // Tpetra objects want the matrix dimensions in these types.
        const global_size_t numRows = static_cast<global_size_t> (dims[0]);
        const size_t numCols = static_cast<size_t> (dims[1]);

        // Make a "Proc 0 owns everything" Map that we will use to
        // read in the multivector entries in the correct order on
        // Proc 0.  This must be a collective
        RCP<const map_type> proc0Map; // "Proc 0 owns everything" Map
        if (map.is_null ()) {
          // The user didn't supply a Map.  Make a contiguous
          // distributed Map for them, using the read-in multivector
          // dimensions.
          map = createUniformContigMapWithNode<LO, GO, NT> (numRows, comm);
          const size_t localNumRows = (myRank == 0) ? numRows : 0;
          // At this point, map exists and has a nonnull node.
          proc0Map = createContigMapWithNode<LO, GO, NT> (numRows, localNumRows,
                                                          comm);
        }
        else { // The user supplied a Map.
          proc0Map = Details::computeGatherMap<map_type> (map, err, debug);
        }

        // Make a multivector X owned entirely by Proc 0.
        RCP<MV> X = createMultiVector<ST, LO, GO, NT> (proc0Map, numCols);

        //
        // On Proc 0, read the Matrix Market data from the input
        // stream into the multivector X.
        //
        int localReadDataSuccess = 1;
        if (myRank == 0) {
          if (! binary) {
            try {
              if (debug) {
                *err << myRank << ": readDenseImpl: Reading matrix data (dense)"
                     << endl;
              }

              // Make sure that we can get a 1-D view of X.
              TEUCHOS_TEST_FOR_EXCEPTION(
                                         ! X->isConstantStride (), std::logic_error,
                                         "Can't get a 1-D view of the entries of the MultiVector X on "
                                         "Process 0, because the stride between the columns of X is not "
                                         "constant.  This shouldn't happen because we just created X and "
                                         "haven't filled it in yet.  Please report this bug to the Tpetra "
                                         "developers.");

              // Get a writeable 1-D view of the entries of X.  Rank 0
              // owns all of them.  The view will expire at the end of
              // scope, so (if necessary) it will be written back to X
              // at this time.
              ArrayRCP<ST> X_view = X->get1dViewNonConst ();
              TEUCHOS_TEST_FOR_EXCEPTION(
                                         as<global_size_t> (X_view.size ()) < numRows * numCols,
                                         std::logic_error,
                                         "The view of X has size " << X_view.size() << " which is not enough to "
                                         "accommodate the expected number of entries numRows*numCols = "
                                         << numRows << "*" << numCols << " = " << numRows*numCols << ".  "
                                         "Please report this bug to the Tpetra developers.");
              const size_t stride = X->getStride ();

              // The third element of the dimensions Tuple encodes the data
              // type reported by the Banner: "real" == 0, "complex" == 1,
              // "integer" == 0 (same as "real"), "pattern" == 2.  We do not
              // allow dense matrices to be pattern matrices, so dims[2] ==
              // 0 or 1.  We've already checked for this above.
              const bool isComplex = (dims[2] == 1);
              size_type count = 0, curRow = 0, curCol = 0;

              std::string line;
              while (getline (in, line)) {
                ++lineNumber;
                // Is the current line a comment line?  If it's not,
                // line.substr(start,size) contains the data.
                size_t start = 0, size = 0;
                const bool commentLine =
                  checkCommentLine (line, start, size, lineNumber, tolerant);
                if (! commentLine) {
                  // Make sure we have room in which to put the new matrix
                  // entry.  We check this only after checking for a
                  // comment line, because there may be one or more
                  // comment lines at the end of the file.  In tolerant
                  // mode, we simply ignore any extra data.
                  if (count >= X_view.size()) {
                    if (tolerant) {
                      break;
                    }
                    else {
                      TEUCHOS_TEST_FOR_EXCEPTION(
                                                 count >= X_view.size(),
                                                 std::runtime_error,
                                                 "The Matrix Market input stream has more data in it than "
                                                 "its metadata reported.  Current line number is "
                                                 << lineNumber << ".");
                    }
                  }

                  // mfh 19 Dec 2012: Ignore everything up to the initial
                  // colon.  writeDense() has the option to print out the
                  // global row index in front of each entry, followed by
                  // a colon and space.
                  {
                    const size_t pos = line.substr (start, size).find (':');
                    if (pos != std::string::npos) {
                      start = pos+1;
                    }
                  }
                  std::istringstream istr (line.substr (start, size));
                  // Does the line contain anything at all?  Can we
                  // safely read from the input stream wrapping the
                  // line?
                  if (istr.eof() || istr.fail()) {
                    // In tolerant mode, simply ignore the line.
                    if (tolerant) {
                      break;
                    }
                    // We repeat the full test here so the exception
                    // message is more informative.
                    TEUCHOS_TEST_FOR_EXCEPTION(
                                               ! tolerant && (istr.eof() || istr.fail()),
                                               std::runtime_error,
                                               "Line " << lineNumber << " of the Matrix Market file is "
                                               "empty, or we cannot read from it for some other reason.");
                  }
                  // Current matrix entry to read in.
                  ST val = STS::zero();
                  // Real and imaginary parts of the current matrix entry.
                  // The imaginary part is zero if the matrix is real-valued.
                  MT real = STM::zero(), imag = STM::zero();

                  // isComplex refers to the input stream's data, not to
                  // the scalar type S.  It's OK to read real-valued
                  // data into a matrix storing complex-valued data; in
                  // that case, all entries' imaginary parts are zero.
                  if (isComplex) {
                    // STS::real() and STS::imag() return a copy of
                    // their respective components, not a writeable
                    // reference.  Otherwise we could just assign to
                    // them using the istream extraction operator (>>).
                    // That's why we have separate magnitude type "real"
                    // and "imag" variables.

                    // Attempt to read the real part of the current entry.
                    istr >> real;
                    if (istr.fail()) {
                      TEUCHOS_TEST_FOR_EXCEPTION(
                                                 ! tolerant && istr.eof(), std::runtime_error,
                                                 "Failed to get the real part of a complex-valued matrix "
                                                 "entry from line " << lineNumber << " of the Matrix Market "
                                                 "file.");
                      // In tolerant mode, just skip bad lines.
                      if (tolerant) {
                        break;
                      }
                    } else if (istr.eof()) {
                      TEUCHOS_TEST_FOR_EXCEPTION(
                                                 ! tolerant && istr.eof(), std::runtime_error,
                                                 "Missing imaginary part of a complex-valued matrix entry "
                                                 "on line " << lineNumber << " of the Matrix Market file.");
                      // In tolerant mode, let any missing imaginary part be 0.
                    } else {
                      // Attempt to read the imaginary part of the current
                      // matrix entry.
                      istr >> imag;
                      TEUCHOS_TEST_FOR_EXCEPTION(
                                                 ! tolerant && istr.fail(), std::runtime_error,
                                                 "Failed to get the imaginary part of a complex-valued "
                                                 "matrix entry from line " << lineNumber << " of the "
                                                 "Matrix Market file.");
                      // In tolerant mode, let any missing or corrupted
                      // imaginary part be 0.
                    }
                  } else { // Matrix Market file contains real-valued data.
                    // Attempt to read the current matrix entry.
                    istr >> real;
                    TEUCHOS_TEST_FOR_EXCEPTION(
                                               ! tolerant && istr.fail(), std::runtime_error,
                                               "Failed to get a real-valued matrix entry from line "
                                               << lineNumber << " of the Matrix Market file.");
                    // In tolerant mode, simply ignore the line if
                    // we failed to read a matrix entry.
                    if (istr.fail() && tolerant) {
                      break;
                    }
                  }
                  // In tolerant mode, we simply let pass through whatever
                  // data we got.
                  TEUCHOS_TEST_FOR_EXCEPTION(
                                             ! tolerant && istr.fail(), std::runtime_error,
                                             "Failed to read matrix data from line " << lineNumber
                                             << " of the Matrix Market file.");

                  // Assign val = ST(real, imag).
                  Teuchos::MatrixMarket::details::assignScalar<ST> (val, real, imag);

                  curRow = count % numRows;
                  curCol = count / numRows;
                  X_view[curRow + curCol*stride] = val;
                  ++count;
                } // if not a comment line
              } // while there are still lines in the file, get the next one

              TEUCHOS_TEST_FOR_EXCEPTION(
                                         ! tolerant && static_cast<global_size_t> (count) < numRows * numCols,
                                         std::runtime_error,
                                         "The Matrix Market metadata reports that the dense matrix is "
                                         << numRows <<  " x " << numCols << ", and thus has "
                                         << numRows*numCols << " total entries, but we only found " << count
                                         << " entr" << (count == 1 ? "y" : "ies") << " in the file.");
            } catch (std::exception& e) {
              exMsg << e.what ();
              localReadDataSuccess = 0;
            }
          } // not binary file
          else {
            if (debug) {
                *err << myRank << ": readDenseImpl: Reading matrix data (dense)"
                     << endl;
              }

              // Make sure that we can get a 1-D view of X.
              TEUCHOS_TEST_FOR_EXCEPTION(
                                         ! X->isConstantStride (), std::logic_error,
                                         "Can't get a 1-D view of the entries of the MultiVector X on "
                                         "Process 0, because the stride between the columns of X is not "
                                         "constant.  This shouldn't happen because we just created X and "
                                         "haven't filled it in yet.  Please report this bug to the Tpetra "
                                         "developers.");

              // Get a writeable 1-D view of the entries of X.  Rank 0
              // owns all of them.  The view will expire at the end of
              // scope, so (if necessary) it will be written back to X
              // at this time.
              auto X_view = X->getLocalViewHost (Access::OverwriteAll);

	      TEUCHOS_TEST_FOR_EXCEPTION(
                                         as<global_size_t> (X_view.extent(0)) < numRows,
                                         std::logic_error,
                                         "The view of X has " << X_view.extent(0) << " rows which is not enough to "
                                         "accommodate the expected number of entries numRows = "
                                         << numRows << ".  "
                                         "Please report this bug to the Tpetra developers.");
	      TEUCHOS_TEST_FOR_EXCEPTION(
                                         as<global_size_t> (X_view.extent(1)) < numCols,
                                         std::logic_error,
                                         "The view of X has " << X_view.extent(1) << " colums which is not enough to "
                                         "accommodate the expected number of entries numRows = "
                                         << numCols << ".  "
                                         "Please report this bug to the Tpetra developers.");

              for (size_t curRow = 0; curRow < numRows; ++curRow) {
                for (size_t curCol = 0; curCol < numCols; ++curCol) {
                  if (Teuchos::ScalarTraits<ST>::isOrdinal){
                    global_size_t val;
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    X_view(curRow, curCol) = val;
                  } else {
                    double val;
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    X_view(curRow, curCol) = val;
                  }
		}
              }
          } // binary
        } // if (myRank == 0)

        if (debug) {
          *err << myRank << ": readDenseImpl: done reading data" << endl;
        }

        // Synchronize on whether Proc 0 successfully read the data.
        int globalReadDataSuccess = localReadDataSuccess;
        broadcast (*comm, 0, outArg (globalReadDataSuccess));
        TEUCHOS_TEST_FOR_EXCEPTION(
          globalReadDataSuccess == 0, std::runtime_error,
          "Failed to read the multivector's data: " << exMsg.str ());

        // If there's only one MPI process and the user didn't supply
        // a Map (i.e., pMap is null), we're done.  Set pMap to the
        // Map used to distribute X, and return X.
        if (comm->getSize () == 1 && map.is_null ()) {
          map = proc0Map;
          if (! err.is_null ()) {
            err->popTab ();
          }
          if (debug) {
            *err << myRank << ": readDenseImpl: done" << endl;
          }
          if (! err.is_null ()) {
            err->popTab ();
          }
          return X;
        }

        if (debug) {
          *err << myRank << ": readDenseImpl: Creating target MV" << endl;
        }

        // Make a multivector Y with the distributed map pMap.
        RCP<MV> Y = createMultiVector<ST, LO, GO, NT> (map, numCols);

        if (debug) {
          *err << myRank << ": readDenseImpl: Creating Export" << endl;
        }

        // Make an Export object that will export X to Y.  First
        // argument is the source map, second argument is the target
        // map.
        Export<LO, GO, NT> exporter (proc0Map, map, err);

        if (debug) {
          *err << myRank << ": readDenseImpl: Exporting" << endl;
        }
        // Export X into Y.
        Y->doExport (*X, exporter, INSERT);

        if (! err.is_null ()) {
          err->popTab ();
        }
        if (debug) {
          *err << myRank << ": readDenseImpl: done" << endl;
        }
        if (! err.is_null ()) {
          err->popTab ();
        }

        // Y is distributed over all process(es) in the communicator.
        return Y;
      }


      template<class VectorScalarType>
      static Teuchos::RCP<Tpetra::Vector<VectorScalarType,
                                     local_ordinal_type,
                                     global_ordinal_type,
                                     node_type> >
      readVectorImpl (std::istream& in,
                      const trcp_tcomm_t& comm,
                      Teuchos::RCP<const map_type>& map,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::MatrixMarket::checkCommentLine;
        using Teuchos::as;
        using Teuchos::broadcast;
        using Teuchos::outArg;
        using Teuchos::RCP;
        using Teuchos::Tuple;
        using std::endl;
        typedef VectorScalarType ST;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;
        typedef Teuchos::ScalarTraits<ST> STS;
        typedef typename STS::magnitudeType MT;
        typedef Teuchos::ScalarTraits<MT> STM;
        typedef Tpetra::Vector<ST, LO, GO, NT> MV;

        // Rank 0 is the only (MPI) process allowed to read from the
        // input stream.
        const int myRank = comm->getRank ();

        if (! err.is_null ()) {
          err->pushTab ();
        }
        if (debug) {
          *err << myRank << ": readVectorImpl" << endl;
        }
        if (! err.is_null ()) {
          err->pushTab ();
        }

        // mfh 17 Feb 2013: It's not strictly necessary that the Comm
        // instances be identical and that the Node instances be
        // identical.  The essential condition is more complicated to
        // test and isn't the same for all Node types.  Thus, we just
        // leave it up to the user.

        // // If map is nonnull, check the precondition that its
        // // communicator resp. node equal comm resp. node.  Checking
        // // now avoids doing a lot of file reading before we detect the
        // // violated precondition.
        // TEUCHOS_TEST_FOR_EXCEPTION(
        //   ! map.is_null() && (map->getComm() != comm || map->getNode () != node,
        //   std::invalid_argument, "If you supply a nonnull Map, the Map's "
        //   "communicator and node must equal the supplied communicator resp. "
        //   "node.");

        // Process 0 will read in the matrix dimensions from the file,
        // and broadcast them to all ranks in the given communicator.
        // There are only 2 dimensions in the matrix, but we use the
        // third element of the Tuple to encode the banner's reported
        // data type: "real" == 0, "complex" == 1, and "integer" == 0
        // (same as "real").  We don't allow pattern matrices (i.e.,
        // graphs) since they only make sense for sparse data.
        Tuple<GO, 3> dims;
        dims[0] = 0;
        dims[1] = 0;

        // Current line number in the input stream.  Only valid on
        // Proc 0.  Various calls will modify this depending on the
        // number of lines that are read from the input stream.
        size_t lineNumber = 1;

        // Capture errors and their messages on Proc 0.
        std::ostringstream exMsg;
        int localBannerReadSuccess = 1;
        int localDimsReadSuccess = 1;

        // Only Proc 0 gets to read matrix data from the input stream.
        if (myRank == 0) {
          if (debug) {
            *err << myRank << ": readVectorImpl: Reading banner line (dense)" << endl;
          }

          // The "Banner" tells you whether the input stream
          // represents a dense matrix, the symmetry type of the
          // matrix, and the type of the data it contains.
          RCP<const Banner> pBanner;
          try {
            pBanner = readBanner (in, lineNumber, tolerant, debug);
          } catch (std::exception& e) {
            exMsg << e.what ();
            localBannerReadSuccess = 0;
          }
          // Make sure the input stream is the right kind of data.
          if (localBannerReadSuccess) {
            if (pBanner->matrixType () != "array") {
              exMsg << "The Matrix Market file does not contain dense matrix "
                "data.  Its banner (first) line says that its matrix type is \""
                << pBanner->matrixType () << "\", rather that the required "
                "\"array\".";
              localBannerReadSuccess = 0;
            } else if (pBanner->dataType() == "pattern") {
              exMsg << "The Matrix Market file's banner (first) "
                "line claims that the matrix's data type is \"pattern\".  This does "
                "not make sense for a dense matrix, yet the file reports the matrix "
                "as dense.  The only valid data types for a dense matrix are "
                "\"real\", \"complex\", and \"integer\".";
              localBannerReadSuccess = 0;
            } else {
              // Encode the data type reported by the Banner as the
              // third element of the dimensions Tuple.
              dims[2] = encodeDataType (pBanner->dataType ());
            }
          } // if we successfully read the banner line

          // At this point, we've successfully read the banner line.
          // Now read the dimensions line.
          if (localBannerReadSuccess) {
            if (debug) {
              *err << myRank << ": readVectorImpl: Reading dimensions line (dense)" << endl;
            }
            // Keep reading lines from the input stream until we find
            // a non-comment line, or until we run out of lines.  The
            // latter is an error, since every "array" format Matrix
            // Market file must have a dimensions line after the
            // banner (even if the matrix has zero rows or columns, or
            // zero entries).
            std::string line;
            bool commentLine = true;

            while (commentLine) {
              // Test whether it is even valid to read from the input
              // stream wrapping the line.
              if (in.eof () || in.fail ()) {
                exMsg << "Unable to get array dimensions line (at all) from line "
                      << lineNumber << " of input stream.  The input stream "
                      << "claims that it is "
                      << (in.eof() ? "at end-of-file." : "in a failed state.");
                localDimsReadSuccess = 0;
              } else {
                // Try to get the next line from the input stream.
                if (getline (in, line)) {
                  ++lineNumber; // We did actually read a line.
                }
                // Is the current line a comment line?  Ignore start
                // and size; they are only useful for reading the
                // actual matrix entries.  (We could use them here as
                // an optimization, but we've chosen not to.)
                size_t start = 0, size = 0;
                commentLine = checkCommentLine (line, start, size, lineNumber, tolerant);
              } // whether we failed to read the line at all
            } // while the line we just read is a comment line

            //
            // Get <numRows> <numCols> from the line we just read.
            //
            std::istringstream istr (line);

            // Test whether it is even valid to read from the input
            // stream wrapping the line.
            if (istr.eof () || istr.fail ()) {
              exMsg << "Unable to read any data from line " << lineNumber
                    << " of input; the line should contain the matrix dimensions "
                    << "\"<numRows> <numCols>\".";
              localDimsReadSuccess = 0;
            } else { // It's valid to read from the line.
              GO theNumRows = 0;
              istr >> theNumRows; // Read in the number of rows.
              if (istr.fail ()) {
                exMsg << "Failed to get number of rows from line "
                      << lineNumber << " of input; the line should contains the "
                      << "matrix dimensions \"<numRows> <numCols>\".";
                localDimsReadSuccess = 0;
              } else { // We successfully read the number of rows
                dims[0] = theNumRows; // Save the number of rows
                if (istr.eof ()) { // Do we still have data to read?
                  exMsg << "No more data after number of rows on line "
                        << lineNumber << " of input; the line should contain the "
                        << "matrix dimensions \"<numRows> <numCols>\".";
                  localDimsReadSuccess = 0;
                } else { // Still data left to read; read in number of columns.
                  GO theNumCols = 0;
                  istr >> theNumCols; // Read in the number of columns
                  if (istr.fail ()) {
                    exMsg << "Failed to get number of columns from line "
                          << lineNumber << " of input; the line should contain "
                          << "the matrix dimensions \"<numRows> <numCols>\".";
                    localDimsReadSuccess = 0;
                  } else { // We successfully read the number of columns
                    dims[1] = theNumCols; // Save the number of columns
                  } // if istr.fail ()
                } // if istr.eof ()
              } // if we read the number of rows
            } // if the input stream wrapping the dims line was (in)valid
          } // if we successfully read the banner line
        } // if (myRank == 0)

        // Check if file has a Vector
        if (dims[1]!=1) {
          exMsg << "File does not contain a 1D Vector";
          localDimsReadSuccess = 0;
        }

        // Broadcast the matrix dimensions, the encoded data type, and
        // whether or not Proc 0 succeeded in reading the banner and
        // dimensions.
        Tuple<GO, 5> bannerDimsReadResult;
        if (myRank == 0) {
          bannerDimsReadResult[0] = dims[0]; // numRows
          bannerDimsReadResult[1] = dims[1]; // numCols
          bannerDimsReadResult[2] = dims[2]; // encoded data type
          bannerDimsReadResult[3] = localBannerReadSuccess;
          bannerDimsReadResult[4] = localDimsReadSuccess;
        }

        // Broadcast matrix dimensions and the encoded data type from
        // Proc 0 to all the MPI processes.
        broadcast (*comm, 0, bannerDimsReadResult);

        TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[3] == 0, std::runtime_error,
          "Failed to read banner line: " << exMsg.str ());
        TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[4] == 0, std::runtime_error,
          "Failed to read matrix dimensions line: " << exMsg.str ());
        if (myRank != 0) {
          dims[0] = bannerDimsReadResult[0];
          dims[1] = bannerDimsReadResult[1];
          dims[2] = bannerDimsReadResult[2];
        }

        // Tpetra objects want the matrix dimensions in these types.
        const global_size_t numRows = static_cast<global_size_t> (dims[0]);
        const size_t numCols = static_cast<size_t> (dims[1]);

        // Make a "Proc 0 owns everything" Map that we will use to
        // read in the multivector entries in the correct order on
        // Proc 0.  This must be a collective
        RCP<const map_type> proc0Map; // "Proc 0 owns everything" Map
        if (map.is_null ()) {
          // The user didn't supply a Map.  Make a contiguous
          // distributed Map for them, using the read-in multivector
          // dimensions.
          map = createUniformContigMapWithNode<LO, GO, NT> (numRows, comm);
          const size_t localNumRows = (myRank == 0) ? numRows : 0;
          // At this point, map exists and has a nonnull node.
          proc0Map = createContigMapWithNode<LO, GO, NT> (numRows, localNumRows,
                                                          comm);
        }
        else { // The user supplied a Map.
          proc0Map = Details::computeGatherMap<map_type> (map, err, debug);
        }

        // Make a multivector X owned entirely by Proc 0.
        RCP<MV> X = createVector<ST, LO, GO, NT> (proc0Map);

        //
        // On Proc 0, read the Matrix Market data from the input
        // stream into the multivector X.
        //
        int localReadDataSuccess = 1;
        if (myRank == 0) {
          try {
            if (debug) {
              *err << myRank << ": readVectorImpl: Reading matrix data (dense)"
                   << endl;
            }

            // Make sure that we can get a 1-D view of X.
            TEUCHOS_TEST_FOR_EXCEPTION(
              ! X->isConstantStride (), std::logic_error,
              "Can't get a 1-D view of the entries of the MultiVector X on "
              "Process 0, because the stride between the columns of X is not "
              "constant.  This shouldn't happen because we just created X and "
              "haven't filled it in yet.  Please report this bug to the Tpetra "
              "developers.");

            // Get a writeable 1-D view of the entries of X.  Rank 0
            // owns all of them.  The view will expire at the end of
            // scope, so (if necessary) it will be written back to X
            // at this time.
            Teuchos::ArrayRCP<ST> X_view = X->get1dViewNonConst ();
            TEUCHOS_TEST_FOR_EXCEPTION(
              as<global_size_t> (X_view.size ()) < numRows * numCols,
              std::logic_error,
              "The view of X has size " << X_view << " which is not enough to "
              "accommodate the expected number of entries numRows*numCols = "
              << numRows << "*" << numCols << " = " << numRows*numCols << ".  "
              "Please report this bug to the Tpetra developers.");
            const size_t stride = X->getStride ();

            // The third element of the dimensions Tuple encodes the data
            // type reported by the Banner: "real" == 0, "complex" == 1,
            // "integer" == 0 (same as "real"), "pattern" == 2.  We do not
            // allow dense matrices to be pattern matrices, so dims[2] ==
            // 0 or 1.  We've already checked for this above.
            const bool isComplex = (dims[2] == 1);
            size_type count = 0, curRow = 0, curCol = 0;

            std::string line;
            while (getline (in, line)) {
              ++lineNumber;
              // Is the current line a comment line?  If it's not,
              // line.substr(start,size) contains the data.
              size_t start = 0, size = 0;
              const bool commentLine =
                checkCommentLine (line, start, size, lineNumber, tolerant);
              if (! commentLine) {
                // Make sure we have room in which to put the new matrix
                // entry.  We check this only after checking for a
                // comment line, because there may be one or more
                // comment lines at the end of the file.  In tolerant
                // mode, we simply ignore any extra data.
                if (count >= X_view.size()) {
                  if (tolerant) {
                    break;
                  }
                  else {
                    TEUCHOS_TEST_FOR_EXCEPTION(
                       count >= X_view.size(),
                       std::runtime_error,
                       "The Matrix Market input stream has more data in it than "
                       "its metadata reported.  Current line number is "
                       << lineNumber << ".");
                  }
                }

                // mfh 19 Dec 2012: Ignore everything up to the initial
                // colon.  writeDense() has the option to print out the
                // global row index in front of each entry, followed by
                // a colon and space.
                {
                  const size_t pos = line.substr (start, size).find (':');
                  if (pos != std::string::npos) {
                    start = pos+1;
                  }
                }
                std::istringstream istr (line.substr (start, size));
                // Does the line contain anything at all?  Can we
                // safely read from the input stream wrapping the
                // line?
                if (istr.eof() || istr.fail()) {
                  // In tolerant mode, simply ignore the line.
                  if (tolerant) {
                    break;
                  }
                  // We repeat the full test here so the exception
                  // message is more informative.
                  TEUCHOS_TEST_FOR_EXCEPTION(
                    ! tolerant && (istr.eof() || istr.fail()),
                    std::runtime_error,
                    "Line " << lineNumber << " of the Matrix Market file is "
                    "empty, or we cannot read from it for some other reason.");
                }
                // Current matrix entry to read in.
                ST val = STS::zero();
                // Real and imaginary parts of the current matrix entry.
                // The imaginary part is zero if the matrix is real-valued.
                MT real = STM::zero(), imag = STM::zero();

                // isComplex refers to the input stream's data, not to
                // the scalar type S.  It's OK to read real-valued
                // data into a matrix storing complex-valued data; in
                // that case, all entries' imaginary parts are zero.
                if (isComplex) {
                  // STS::real() and STS::imag() return a copy of
                  // their respective components, not a writeable
                  // reference.  Otherwise we could just assign to
                  // them using the istream extraction operator (>>).
                  // That's why we have separate magnitude type "real"
                  // and "imag" variables.

                  // Attempt to read the real part of the current entry.
                  istr >> real;
                  if (istr.fail()) {
                    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.eof(), std::runtime_error,
                      "Failed to get the real part of a complex-valued matrix "
                      "entry from line " << lineNumber << " of the Matrix Market "
                      "file.");
                    // In tolerant mode, just skip bad lines.
                    if (tolerant) {
                      break;
                    }
                  } else if (istr.eof()) {
                    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.eof(), std::runtime_error,
                      "Missing imaginary part of a complex-valued matrix entry "
                      "on line " << lineNumber << " of the Matrix Market file.");
                    // In tolerant mode, let any missing imaginary part be 0.
                  } else {
                    // Attempt to read the imaginary part of the current
                    // matrix entry.
                    istr >> imag;
                    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.fail(), std::runtime_error,
                      "Failed to get the imaginary part of a complex-valued "
                      "matrix entry from line " << lineNumber << " of the "
                      "Matrix Market file.");
                    // In tolerant mode, let any missing or corrupted
                    // imaginary part be 0.
                  }
                } else { // Matrix Market file contains real-valued data.
                  // Attempt to read the current matrix entry.
                  istr >> real;
                  TEUCHOS_TEST_FOR_EXCEPTION(
                    ! tolerant && istr.fail(), std::runtime_error,
                    "Failed to get a real-valued matrix entry from line "
                    << lineNumber << " of the Matrix Market file.");
                  // In tolerant mode, simply ignore the line if
                  // we failed to read a matrix entry.
                  if (istr.fail() && tolerant) {
                    break;
                  }
                }
                // In tolerant mode, we simply let pass through whatever
                // data we got.
                TEUCHOS_TEST_FOR_EXCEPTION(
                  ! tolerant && istr.fail(), std::runtime_error,
                  "Failed to read matrix data from line " << lineNumber
                  << " of the Matrix Market file.");

                // Assign val = ST(real, imag).
                Teuchos::MatrixMarket::details::assignScalar<ST> (val, real, imag);

                curRow = count % numRows;
                curCol = count / numRows;
                X_view[curRow + curCol*stride] = val;
                ++count;
              } // if not a comment line
            } // while there are still lines in the file, get the next one

            TEUCHOS_TEST_FOR_EXCEPTION(
              ! tolerant && static_cast<global_size_t> (count) < numRows * numCols,
              std::runtime_error,
              "The Matrix Market metadata reports that the dense matrix is "
              << numRows <<  " x " << numCols << ", and thus has "
              << numRows*numCols << " total entries, but we only found " << count
              << " entr" << (count == 1 ? "y" : "ies") << " in the file.");
          } catch (std::exception& e) {
            exMsg << e.what ();
            localReadDataSuccess = 0;
          }
        } // if (myRank == 0)

        if (debug) {
          *err << myRank << ": readVectorImpl: done reading data" << endl;
        }

        // Synchronize on whether Proc 0 successfully read the data.
        int globalReadDataSuccess = localReadDataSuccess;
        broadcast (*comm, 0, outArg (globalReadDataSuccess));
        TEUCHOS_TEST_FOR_EXCEPTION(
          globalReadDataSuccess == 0, std::runtime_error,
          "Failed to read the multivector's data: " << exMsg.str ());

        // If there's only one MPI process and the user didn't supply
        // a Map (i.e., pMap is null), we're done.  Set pMap to the
        // Map used to distribute X, and return X.
        if (comm->getSize () == 1 && map.is_null ()) {
          map = proc0Map;
          if (! err.is_null ()) {
            err->popTab ();
          }
          if (debug) {
            *err << myRank << ": readVectorImpl: done" << endl;
          }
          if (! err.is_null ()) {
            err->popTab ();
          }
          return X;
        }

        if (debug) {
          *err << myRank << ": readVectorImpl: Creating target MV" << endl;
        }

        // Make a multivector Y with the distributed map pMap.
        RCP<MV> Y = createVector<ST, LO, GO, NT> (map);

        if (debug) {
          *err << myRank << ": readVectorImpl: Creating Export" << endl;
        }

        // Make an Export object that will export X to Y.  First
        // argument is the source map, second argument is the target
        // map.
        Export<LO, GO, NT> exporter (proc0Map, map, err);

        if (debug) {
          *err << myRank << ": readVectorImpl: Exporting" << endl;
        }
        // Export X into Y.
        Y->doExport (*X, exporter, INSERT);

        if (! err.is_null ()) {
          err->popTab ();
        }
        if (debug) {
          *err << myRank << ": readVectorImpl: done" << endl;
        }
        if (! err.is_null ()) {
          err->popTab ();
        }

        // Y is distributed over all process(es) in the communicator.
        return Y;
      }

    public:
      /// \brief Read Map (as a MultiVector) from the given input stream.
      ///
      /// Read the given input stream on MPI Process 0 (with respect
      /// to the given communicator).  The stream should contain
      /// Matrix Market "array" format dense matrix data with two
      /// columns, as generated by Writer::writeMap() or
      /// Writer::writeMapFile().  Distribute the data from Process 0
      /// to all processes.  Return the resulting Map.
      ///
      /// \param in [in/out] Input stream of Matrix Market data from
      ///   which to read.  This is only accessed on Process 0 of the
      ///   given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the Map will be distributed.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      /// \param debug [in] If true, read in binary mode.
      static Teuchos::RCP<const map_type>
      readMap (std::istream& in,
               const trcp_tcomm_t& comm,
               const bool tolerant=false,
               const bool debug=false,
               const bool binary=false)
      {
        Teuchos::RCP<Teuchos::FancyOStream> err =
          Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
        return readMap (in, comm, err, tolerant, debug, binary);
      }


      /// \brief Read Map (as a MultiVector) from the given input
      ///   stream, with optional debugging output stream.
      ///
      /// \warning We make no promises about backwards compatibility
      ///   for this method.  It may disappear or its interface may
      ///   change at any time.
      ///
      /// Read the given input stream on MPI Process 0 (with respect
      /// to the given communicator).  The stream should contain
      /// Matrix Market "array" format dense matrix data with two
      /// columns, as generated by Writer::writeMap() or
      /// Writer::writeMapFile().  Distribute the data from Process 0
      /// to all processes.  Return the resulting Map.
      ///
      /// \param in [in/out] Input stream of Matrix Market data from
      ///   which to read.  This is only accessed on Process 0 of the
      ///   given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the Map will be distributed.
      /// \param err [in] Optional output stream for debugging output.
      ///   This is only referenced if \c debug is true.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] If true, write copious debugging output to
      ///   \c err on all processes in \c comm.
      /// \param debug [in] If true, read in binary mode.
      static Teuchos::RCP<const map_type>
      readMap (std::istream& in,
               const trcp_tcomm_t& comm,
               const Teuchos::RCP<Teuchos::FancyOStream>& err,
               const bool tolerant=false,
               const bool debug=false,
               const bool binary=false)
      {
        using Teuchos::arcp;
        using Teuchos::Array;
        using Teuchos::ArrayRCP;
        using Teuchos::as;
        using Teuchos::broadcast;
        using Teuchos::Comm;
        using Teuchos::CommRequest;
        using Teuchos::inOutArg;
        using Teuchos::ireceive;
        using Teuchos::outArg;
        using Teuchos::RCP;
        using Teuchos::receive;
        using Teuchos::reduceAll;
        using Teuchos::REDUCE_MIN;
        using Teuchos::isend;
        using Teuchos::SerialComm;
        using Teuchos::toString;
        using Teuchos::wait;
        using std::endl;
        typedef Tpetra::global_size_t GST;
        typedef ptrdiff_t int_type; // Can hold int and GO
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;
        typedef Tpetra::MultiVector<GO, LO, GO, NT> MV;

        const int numProcs = comm->getSize ();
        const int myRank = comm->getRank ();

        if (err.is_null ()) {
          err->pushTab ();
        }
        if (debug) {
          std::ostringstream os;
          os << myRank << ": readMap: " << endl;
          *err << os.str ();
        }
        if (err.is_null ()) {
          err->pushTab ();
        }

        // Tag for receive-size / send-size messages.  writeMap used
        // tags 1337 and 1338; we count up from there.
        const int sizeTag = 1339;
        // Tag for receive-data / send-data messages.
        const int dataTag = 1340;

        // These are for sends on Process 0, and for receives on all
        // other processes.  sizeReq is for the {receive,send}-size
        // message, and dataReq is for the message containing the
        // actual GIDs to belong to the receiving process.
        RCP<CommRequest<int> > sizeReq;
        RCP<CommRequest<int> > dataReq;

        // Each process will have to receive the number of GIDs to
        // expect.  Thus, we can post the receives now, and cancel
        // them if something should go wrong in the meantime.
        ArrayRCP<int_type> numGidsToRecv (1);
        numGidsToRecv[0] = 0;
        if (myRank != 0) {
          sizeReq = ireceive<int, int_type> (numGidsToRecv, 0, sizeTag, *comm);
        }

        int readSuccess = 1;
        std::ostringstream exMsg;
        RCP<MV> data; // Will only be valid on Proc 0
        if (myRank == 0) {
          // If we want to reuse readDenseImpl, we have to make a
          // communicator that only contains Proc 0.  Otherwise,
          // readDenseImpl will redistribute the data to all
          // processes.  While we eventually want that, neither we nor
          // readDenseImpl know the correct Map to use at the moment.
          // That depends on the second column of the multivector.
          RCP<const Comm<int> > proc0Comm (new SerialComm<int> ());
          try {
            RCP<const map_type> dataMap;
            // This is currently the only place where we use the
            // 'tolerant' argument.  Later, if we want to be clever,
            // we could have tolerant mode allow PIDs out of order.
            data = readDenseImpl<GO> (in, proc0Comm, dataMap, err, tolerant, debug, binary);
            (void) dataMap; // Silence "unused" warnings
            if (data.is_null ()) {
              readSuccess = 0;
              exMsg << "readDenseImpl() returned null." << endl;
            }
          } catch (std::exception& e) {
            readSuccess = 0;
            exMsg << e.what () << endl;
          }
        }

        // Map from PID to all the GIDs for that PID.
        // Only populated on Process 0.
        std::map<int, Array<GO> > pid2gids;

        // The index base must be the global minimum GID.
        // We will compute this on Process 0 and broadcast,
        // so that all processes can set up the Map.
        int_type globalNumGIDs = 0;

        // The index base must be the global minimum GID.
        // We will compute this on Process 0 and broadcast,
        // so that all processes can set up the Map.
        GO indexBase = 0;

        // Process 0: If the above read of the MultiVector succeeded,
        // extract the GIDs and PIDs into pid2gids, and find the
        // global min GID.
        if (myRank == 0 && readSuccess == 1) {
          if (data->getNumVectors () == 2) { // Map format 1.0
            ArrayRCP<const GO> GIDs = data->getData (0);
            ArrayRCP<const GO> PIDs = data->getData (1); // convert to int
            globalNumGIDs = GIDs.size ();

            // Start computing the global min GID, while collecting
            // the GIDs for each PID.
            if (globalNumGIDs > 0) {
              const int pid = static_cast<int> (PIDs[0]);

              if (pid < 0 || pid >= numProcs) {
                readSuccess = 0;
                exMsg << "Tpetra::MatrixMarket::readMap: "
                      << "Encountered invalid PID " << pid << "." << endl;
              }
              else {
                const GO gid = GIDs[0];
                pid2gids[pid].push_back (gid);
                indexBase = gid; // the current min GID
              }
            }
            if (readSuccess == 1) {
              // Collect the rest of the GIDs for each PID, and compute
              // the global min GID.
              for (size_type k = 1; k < globalNumGIDs; ++k) {
                const int pid = static_cast<int> (PIDs[k]);
                if (pid < 0 || pid >= numProcs) {
                  readSuccess = 0;
                  exMsg << "Tpetra::MatrixMarket::readMap: "
                        << "Encountered invalid PID " << pid << "." << endl;
                }
                else {
                  const int_type gid = GIDs[k];
                  pid2gids[pid].push_back (gid);
                  if (gid < indexBase) {
                    indexBase = gid; // the current min GID
                  }
                }
              }
            }
          }
          else if (data->getNumVectors () == 1) { // Map format 2.0
            if (data->getGlobalLength () % 2 != static_cast<GST> (0)) {
              readSuccess = 0;
              exMsg << "Tpetra::MatrixMarket::readMap: Input data has the "
                "wrong format (for Map format 2.0).  The global number of rows "
                "in the MultiVector must be even (divisible by 2)." << endl;
            }
            else {
              ArrayRCP<const GO> theData = data->getData (0);
              globalNumGIDs = static_cast<GO> (data->getGlobalLength ()) /
                static_cast<GO> (2);

              // Start computing the global min GID, while
              // collecting the GIDs for each PID.
              if (globalNumGIDs > 0) {
                const int pid = static_cast<int> (theData[1]);
                if (pid < 0 || pid >= numProcs) {
                  readSuccess = 0;
                  exMsg << "Tpetra::MatrixMarket::readMap: "
                        << "Encountered invalid PID " << pid << "." << endl;
                }
                else {
                  const GO gid = theData[0];
                  pid2gids[pid].push_back (gid);
                  indexBase = gid; // the current min GID
                }
              }
              // Collect the rest of the GIDs for each PID, and
              // compute the global min GID.
              for (int_type k = 1; k < globalNumGIDs; ++k) {
                const int pid = static_cast<int> (theData[2*k + 1]);
                if (pid < 0 || pid >= numProcs) {
                  readSuccess = 0;
                  exMsg << "Tpetra::MatrixMarket::readMap: "
                        << "Encountered invalid PID " << pid << "." << endl;
                }
                else {
                  const GO gid = theData[2*k];
                  pid2gids[pid].push_back (gid);
                  if (gid < indexBase) {
                    indexBase = gid; // the current min GID
                  }
                }
              } // for each GID
            } // if the amount of data is correct
          }
          else {
            readSuccess = 0;
            exMsg << "Tpetra::MatrixMarket::readMap: Input data must have "
              "either 1 column (for the new Map format 2.0) or 2 columns (for "
              "the old Map format 1.0).";
          }
        } // myRank is zero

        // Broadcast the indexBase, the global number of GIDs, and the
        // current success status.  Use int_type for all of these.
        {
          int_type readResults[3];
          readResults[0] = static_cast<int_type> (indexBase);
          readResults[1] = static_cast<int_type> (globalNumGIDs);
          readResults[2] = static_cast<int_type> (readSuccess);
          broadcast<int, int_type> (*comm, 0, 3, readResults);

          indexBase = static_cast<GO> (readResults[0]);
          globalNumGIDs = static_cast<int_type> (readResults[1]);
          readSuccess = static_cast<int> (readResults[2]);
        }

        // Unwinding the stack will invoke sizeReq's destructor, which
        // will cancel the receive-size request on all processes that
        // posted it.
        TEUCHOS_TEST_FOR_EXCEPTION(
          readSuccess != 1, std::runtime_error,
          "Tpetra::MatrixMarket::readMap: Reading the Map failed with the "
          "following exception message: " << exMsg.str ());

        if (myRank == 0) {
          // Proc 0: Send each process' number of GIDs to that process.
          for (int p = 1; p < numProcs; ++p) {
            ArrayRCP<int_type> numGidsToSend (1);

            auto it = pid2gids.find (p);
            if (it == pid2gids.end ()) {
              numGidsToSend[0] = 0;
            } else {
              numGidsToSend[0] = it->second.size ();
            }
            sizeReq = isend<int, int_type> (numGidsToSend, p, sizeTag, *comm);
            wait<int> (*comm, outArg (sizeReq));
          }
        }
        else {
          // Wait on the receive-size message to finish.
          wait<int> (*comm, outArg (sizeReq));
        }

        // Allocate / get the array for my GIDs.
        // Only Process 0 will have its actual GIDs at this point.
        ArrayRCP<GO> myGids;
        int_type myNumGids = 0;
        if (myRank == 0) {
          GO* myGidsRaw = NULL;

          typename std::map<int, Array<GO> >::iterator it = pid2gids.find (0);
          if (it != pid2gids.end ()) {
            myGidsRaw = it->second.getRawPtr ();
            myNumGids = it->second.size ();
            // Nonowning ArrayRCP just views the Array.
            myGids = arcp<GO> (myGidsRaw, 0, myNumGids, false);
          }
        }
        else { // myRank != 0
          myNumGids = numGidsToRecv[0];
          myGids = arcp<GO> (myNumGids);
        }

        if (myRank != 0) {
          // Post receive for data, now that we know how much data we
          // will receive.  Only post receive if my process actually
          // has nonzero GIDs.
          if (myNumGids > 0) {
            dataReq = ireceive<int, GO> (myGids, 0, dataTag, *comm);
          }
        }

        for (int p = 1; p < numProcs; ++p) {
          if (myRank == 0) {
            ArrayRCP<GO> sendGids; // to send to Process p
            GO* sendGidsRaw = NULL;
            int_type numSendGids = 0;

            typename std::map<int, Array<GO> >::iterator it = pid2gids.find (p);
            if (it != pid2gids.end ()) {
              numSendGids = it->second.size ();
              sendGidsRaw = it->second.getRawPtr ();
              sendGids = arcp<GO> (sendGidsRaw, 0, numSendGids, false);
            }
            // Only send if that process actually has nonzero GIDs.
            if (numSendGids > 0) {
              dataReq = isend<int, GO> (sendGids, p, dataTag, *comm);
            }
            wait<int> (*comm, outArg (dataReq));
          }
          else if (myRank == p) {
            // Wait on my receive of GIDs to finish.
            wait<int> (*comm, outArg (dataReq));
          }
        } // for each process rank p in 1, 2, ..., numProcs-1

        if (debug) {
          std::ostringstream os;
          os << myRank << ": readMap: creating Map" << endl;
          *err << os.str ();
        }
        const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
        RCP<const map_type> newMap;

        // Create the Map; test whether the constructor threw.  This
        // avoids deadlock and makes error reporting more readable.

        int lclSuccess = 1;
        int gblSuccess = 0; // output argument
        std::ostringstream errStrm;
        try {
          newMap = rcp (new map_type (INVALID, myGids (), indexBase, comm));
        }
        catch (std::exception& e) {
          lclSuccess = 0;
          errStrm << "Process " << comm->getRank () << " threw an exception: "
                  << e.what () << std::endl;
        }
        catch (...) {
          lclSuccess = 0;
          errStrm << "Process " << comm->getRank () << " threw an exception "
            "not a subclass of std::exception" << std::endl;
        }
        Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN,
                                      lclSuccess, Teuchos::outArg (gblSuccess));
        if (gblSuccess != 1) {
          Tpetra::Details::gathervPrint (std::cerr, errStrm.str (), *comm);
        }
        TEUCHOS_TEST_FOR_EXCEPTION(gblSuccess != 1, std::runtime_error, "Map constructor failed!");

        if (err.is_null ()) {
          err->popTab ();
        }
        if (debug) {
          std::ostringstream os;
          os << myRank << ": readMap: done" << endl;
          *err << os.str ();
        }
        if (err.is_null ()) {
          err->popTab ();
        }
        return newMap;
      }


    private:

      /// \brief Encode the Matrix Market data type as an int.
      ///
      /// We assume that the Banner has already validated the data
      /// type string, so the string is one of four canonical values.
      /// We encode these values as follows: "real" as 0, "complex" as
      /// 1, "integer" as 0 (same as "real", since the way we
      /// implement reading integer values is the same as the way we
      /// implement reading real values), and "pattern" as 2.
      ///
      /// We use this encoding for communicating the data type.
      static int
      encodeDataType (const std::string& dataType)
      {
        if (dataType == "real" || dataType == "integer") {
          return 0;
        } else if (dataType == "complex") {
          return 1;
        } else if (dataType == "pattern") {
          return 2;
        } else {
          // We should never get here, since Banner validates the
          // reported data type and ensures it is one of the accepted
          // values.
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Unrecognized Matrix Market data type \"" << dataType
            << "\".  We should never get here.  "
            "Please report this bug to the Tpetra developers.");
        }
      }

    public:

      /// \brief Read a Tpetra::CrsMatrix from a file per rank setup
      ///
      //! Function to read a one-per-rank collection of MatrixMarket files
      //! and assemble it into a single big matrix.  The code will try to minimize
      //! the number of ranks hammering on the file system at once, but we don't
      //! make any guarantees.
      /// \param filename_prefix  [in] File for rank I is filename_prefix + to_string(I) + filename_suffix
      /// \param filename_sufffix [in] File for rank I is filename_prefix + to_string(I) + filename_suffix
      /// \param rowMap [in] The Map over which to distribute rows
      ///   of the sparse matrix.  This must be nonnull.
      /// \param colMap [in/out] If nonnull: the Map over which to
      ///   distribute columns of the sparse matrix.  If null and if
      ///   callFillComplete is true, we create this for you.
      /// \param domainMap [in] The sparse matrix's domain Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param rangeMap [in] The sparse matrix's range Map.  This
      ///   must be nonnull.  It may equal (pointer equality) the row
      ///   Map, if that would be appropriate for this matrix.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static Teuchos::RCP<sparse_matrix_type>
      readSparsePerRank (const std::string& filename_prefix,
                         const std::string& filename_suffix,
                         const Teuchos::RCP<const map_type>& rowMap,
                         Teuchos::RCP<const map_type>& colMap,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const bool callFillComplete=true,
                         const bool tolerant=false,
                         const int ranksToReadAtOnce=8,
                         const bool debug=false)
      {
        using ST = scalar_type;
        using LO = local_ordinal_type;
        using GO = global_ordinal_type;
        using STS = typename Teuchos::ScalarTraits<ST>;
        using Teuchos::RCP;
        using Teuchos::ArrayRCP;
        using Teuchos::arcp;
        using Teuchos::rcp;
        
        // Sanity Checks
        // Fast checks for invalid input.  We can't check other
        // attributes of the Maps until we've read in the matrix
        // dimensions.
        TEUCHOS_TEST_FOR_EXCEPTION(
          rowMap.is_null (), std::invalid_argument,
          "Row Map must be nonnull.");
        Teuchos::RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();
        TEUCHOS_TEST_FOR_EXCEPTION
          (comm.is_null (), std::invalid_argument,
           "The input row map's communicator (Teuchos::Comm object) is null.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          rangeMap.is_null (), std::invalid_argument,
          "Range Map must be nonnull.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          domainMap.is_null (), std::invalid_argument,
          "Domain Map must be nonnull.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          domainMap->getComm().getRawPtr() != comm.getRawPtr(),
          std::invalid_argument,
          "The specified domain Map's communicator (domainMap->getComm())"
          "differs from row Map's communicator");
        TEUCHOS_TEST_FOR_EXCEPTION(
          rangeMap->getComm().getRawPtr() != comm.getRawPtr(),
          std::invalid_argument,
          "The specified range Map's communicator (rangeMap->getComm())"
          "differs from row Map's communicator");
        
        // Setup
        const int myRank  = comm->getRank();
        const int numProc = comm->getSize();
        std::string filename = filename_prefix + std::to_string(myRank) + filename_suffix;  

        // Bounds check the writing limits
        int rank_limit = std::min(std::max(ranksToReadAtOnce,1),numProc);
        
        // Data structures for constructor
        ArrayRCP<size_t> numEntriesPerRow;             
        ArrayRCP<size_t> rowPtr;
        ArrayRCP<global_ordinal_type> colInd;
        ArrayRCP<scalar_type> values;
        std::ostringstream errMsg;
        
        ///////////////////////////////////////////////
        // Start the reading of the banners to get
        // local row / nnz counts and then read the 
        // data.  We'll pack everything into big ol'
        // rowptr/colind/values ArrayRCPs
        bool success = true;
        int bannerIsCorrect = 1, readSuccess = 1;
        LO numRows, numCols, numNonzeros;
        for(int base_rank = 0; base_rank < numProc; base_rank += rank_limit) {
          int stop = std::min(base_rank+rank_limit,numProc);

          // Is my rank in this batch?
          if(base_rank <= myRank  && myRank < stop) {
            // My turn to read
            std::ifstream in(filename);
            using Teuchos::MatrixMarket::Banner;
            size_t lineNumber = 1;
            RCP<const Banner> pBanner;
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }
            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse matrix.
              // We validate on Proc 0, since it reads the Banner.
              
              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                  "\"coordinate\"-format sparse matrix in order to create a "
                  "Tpetra::CrsMatrix object from it, but the file's matrix "
                  "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                  "format sparse matrix in order to create a Tpetra::CrsMatrix "
                  "object from it, but the file's matrix type is \"array\" "
                  "instead.  That probably means the file contains dense matrix "
                  "data.";
              }
            }
              
            // Unpacked coordinate matrix dimensions
            using Teuchos::MatrixMarket::readCoordinateDimensions;
            success = readCoordinateDimensions (in, numRows, numCols,
                                                numNonzeros, lineNumber,
                                                tolerant);

            // Sanity checking of headers
            TEUCHOS_TEST_FOR_EXCEPTION(numRows != (LO)rowMap->getLocalNumElements(), std::invalid_argument,
                                       "# rows in file does not match rowmap.");
            TEUCHOS_TEST_FOR_EXCEPTION(!colMap.is_null() && numCols != (LO)colMap->getLocalNumElements(), std::invalid_argument,
                                       "# rows in file does not match colmap.");
            
            
            // Read the data
            typedef Teuchos::MatrixMarket::Raw::Adder<scalar_type,global_ordinal_type> raw_adder_type;
            bool tolerant_required = true;
            Teuchos::RCP<raw_adder_type> pRaw =
              Teuchos::rcp (new raw_adder_type (numRows,numCols,numNonzeros,tolerant_required,debug));
            RCP<adder_type> pAdder =  Teuchos::rcp (new adder_type (pRaw, pBanner->symmType ()));
            
            if (debug) {
              std::cerr << "-- Reading matrix data" << std::endl;
            }

            try {
              // Reader for "coordinate" format sparse matrix data.
              typedef Teuchos::MatrixMarket::CoordDataReader<adder_type,
                                                             global_ordinal_type, scalar_type, STS::isComplex> reader_type;
              reader_type reader (pAdder);
              
              // Read the sparse matrix entries.
              std::pair<bool, std::vector<size_t> > results = reader.read (in, lineNumber, tolerant_required, debug);


              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }

            ///////////////////////////////////////
            // Create the CSR Arrays
            typedef Teuchos::MatrixMarket::Raw::Element<scalar_type,global_ordinal_type> element_type;
            
            // Additively merge duplicate matrix entries.
            pAdder->getAdder()->merge ();
            
            // Get a temporary const view of the merged matrix entries.
            const std::vector<element_type>& entries =  pAdder->getAdder()->getEntries();
            
            // Number of matrix entries (after merging).
            const size_t numEntries = (size_t)entries.size();
            
            if (debug) {
              std::cerr << "----- Proc "<<myRank<<": Matrix has numRows=" << numRows
                        << " rows and numEntries=" << numEntries
                   << " entries." << std::endl;
            }
            

            // Make space for the CSR matrix data.  Converting to
            // CSR is easier if we fill numEntriesPerRow with zeros
            // at first.
            numEntriesPerRow = arcp<size_t> (numRows);
            std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 0);
            rowPtr = arcp<size_t> (numRows+1);
            std::fill (rowPtr.begin(), rowPtr.end(), 0);
            colInd = arcp<global_ordinal_type> (numEntries);
            values = arcp<scalar_type> (numEntries);

            // Convert from array-of-structs coordinate format to CSR
            // (compressed sparse row) format.
            global_ordinal_type l_prvRow = 0;
            size_t curPos = 0;
            LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
            rowPtr[0] = 0;
            LO indexBase = rowMap->getIndexBase();
            for (curPos = 0; curPos < numEntries; ++curPos) {
              const element_type& curEntry = entries[curPos];
              const global_ordinal_type curRow = curEntry.rowIndex() + indexBase;
              LO l_curRow = rowMap->getLocalElement(curRow);


              TEUCHOS_TEST_FOR_EXCEPTION(l_curRow == INVALID,std::logic_error,
                                         "Current global row "<< curRow << " is invalid.");
                  
              TEUCHOS_TEST_FOR_EXCEPTION(l_curRow < l_prvRow, std::logic_error,      
                                         "Row indices are out of order, even though they are supposed "
                                         "to be sorted.  curRow = " << l_curRow << ", prvRow = "
                                         << l_prvRow << ", at curPos = " << curPos << ".  Please report "
                                         "this bug to the Tpetra developers.");
              if (l_curRow > l_prvRow) {
                for (LO r = l_prvRow+1; r <= l_curRow; ++r) {
                  rowPtr[r] = curPos;
                }
                l_prvRow = l_curRow;
              }
              numEntriesPerRow[l_curRow]++;
              colInd[curPos] = curEntry.colIndex() + indexBase;
              values[curPos] = curEntry.value();

            }
            // rowPtr has one more entry than numEntriesPerRow.  The
            // last entry of rowPtr is the number of entries in
            // colInd and values.
            rowPtr[numRows] = numEntries;

          }// end base_rank <= myRank < stop

          // Barrier between batches to keep the filesystem happy
          comm->barrier();

        }//end outer rank loop


        // Call the matrix constructor and fill.  This isn't particularly efficient
        RCP<sparse_matrix_type> A;
        if(colMap.is_null()) {
          A=rcp(new sparse_matrix_type(rowMap,numEntriesPerRow()));
          for(size_t i=0; i<rowMap->getLocalNumElements(); i++) {
            GO g_row = rowMap->getGlobalElement(i);
            size_t start = rowPtr[i];
            size_t size  = rowPtr[i+1] - rowPtr[i];
            if(size>0)  {
              A->insertGlobalValues(g_row,size,&values[start],&colInd[start]);
            }
          }
        }
        else {
          throw std::runtime_error("Reading with a column map is not yet implemented");
        }       
        RCP<const map_type> myDomainMap = domainMap.is_null() ? rowMap : domainMap;
        RCP<const map_type> myRangeMap  = rangeMap.is_null() ? rowMap : rangeMap;

        A->fillComplete(myDomainMap,myRangeMap);

        if(!readSuccess)
          success = false;
        TEUCHOS_TEST_FOR_EXCEPTION(success == false, std::runtime_error,
                                   "Read failed.");        

        return A;
      }// end readSparsePerRank


    }; // class Reader

    /// \class Writer
    /// \brief Matrix Market file writer for CrsMatrix and MultiVector.
    /// \tparam SparseMatrixType A specialization of Tpetra::CrsMatrix.
    /// \author Mark Hoemmen
    ///
    /// The Matrix Market (see their <a
    /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
    /// details) defines a human-readable ASCII text file format for
    /// interchange of sparse and dense matrices.  This class defines
    /// methods for writing sparse and dense matrices to a Matrix
    /// Market file or input stream.
    ///
    /// All methods of this class assume that the file is only
    /// openable resp. the input stream is only writeable, on MPI
    /// Process 0 (with respect to the MPI communicator over which the
    /// given CrsMatrix or MultiVector is to be distributed).
    ///
    /// We define the MultiVector type accepted by writeDense() and
    /// writeDenseFile() using the scalar_type, local_ordinal_type,
    /// global_ordinal_type, and node_type typedefs in
    /// <tt>SparseMatrixType</tt>.  This ensures that the
    /// Tpetra::MultiVector objects returned by those methods have a
    /// type compatible with the Tpetra::CrsMatrix sparse matrices
    /// accepted by writeSparse() and writeSparseFile().  We do this
    /// because the typical use case of Matrix Market files in
    /// Trilinos is to test sparse matrix methods, which usually
    /// involves reading a sparse matrix A and perhaps also a dense
    /// right-hand side b.
    template<class SparseMatrixType>
    class Writer {
    public:
      //! Template parameter of this class; specialization of CrsMatrix.
      typedef SparseMatrixType sparse_matrix_type;
      typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;

      //! Type of the entries of the sparse matrix.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      //! Type of the local indices of the sparse matrix.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// \brief Type of indices as read from the Matrix Market file.
      ///
      /// Indices of the sparse matrix are stored as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      //! The Kokkos Node type; fourth template parameter of Tpetra::CrsMatrix.
      typedef typename SparseMatrixType::node_type node_type;

      //! Specialization of Tpetra::MultiVector that matches SparseMatrixType.
      typedef MultiVector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> multivector_type;
      //! Specialization of Tpetra::Map that matches SparseMatrixType.
      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      //! Specialization of Tpetra::CrsGraph that matches SparseMatrixType.
      typedef CrsGraph<local_ordinal_type, global_ordinal_type, node_type> crs_graph_type;

      typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>            operator_type;
      typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>         mv_type;

      //! Type of the MPI communicator.
      using trcp_tcomm_t = Teuchos::RCP<const Teuchos::Comm<int>>;

      /// \brief Print the sparse matrix in Matrix Market format, with
      ///   comments.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to the given
      /// file, using the Matrix Market "coordinate" format.  MPI Proc
      /// 0 is the only MPI process that opens or writes to the file.
      /// Include the matrix name and description in the comments
      /// section of the file (after the initial banner line, but
      /// before the matrix metadata and data).
      ///
      /// \param filename [in] Name of the file to which to write the
      ///   given sparse matrix.  The matrix is distributed, but only
      ///   Proc 0 opens the file and writes to it.
      ///
      /// \param matrix [in] The sparse matrix to write to the file.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output file.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output file.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      static void
      writeSparseFile (const std::string& filename,
                       const sparse_matrix_type& matrix,
                       const std::string& matrixName,
                       const std::string& matrixDescription,
                       const bool debug=false)
      {
        trcp_tcomm_t comm = matrix.getComm ();
        TEUCHOS_TEST_FOR_EXCEPTION
          (comm.is_null (), std::invalid_argument,
          "The input matrix's communicator (Teuchos::Comm object) is null.");
        const int myRank = comm->getRank ();

        auto out = Writer::openOutFileOnRankZero(comm, filename, myRank, true);

        writeSparse (out, matrix, matrixName, matrixDescription, debug);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeSparse() throws an
        // exception.
      }

      //! Only for backwards compatibility; prefer the overload above.
      static void
      writeSparseFile (const std::string& filename,
                       const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                       const std::string& matrixName,
                       const std::string& matrixDescription,
                       const bool debug=false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (pMatrix.is_null (), std::invalid_argument,
           "The input matrix is null.");
        writeSparseFile (filename, *pMatrix, matrixName,
                         matrixDescription, debug);
      }

      /// \brief Print the sparse matrix in Matrix Market format.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to the given
      /// file, using the Matrix Market "coordinate" format.  MPI Proc
      /// 0 is the only MPI process that opens or writes to the file.
      ///
      /// \param filename [in] Name of the file to which to write the
      ///   given sparse matrix.  The matrix is distributed, but only
      ///   Proc 0 opens the file and writes to it.
      ///
      /// \param matrix [in] The sparse matrix to write to the file.
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      static void
      writeSparseFile (const std::string& filename,
                       const sparse_matrix_type& matrix,
                       const bool debug=false)
      {
        writeSparseFile (filename, matrix, "", "", debug);
      }

      //! Only for backwards compatibility; prefer the overload above.
      static void
      writeSparseFile (const std::string& filename,
                       const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                       const bool debug=false)
      {
        writeSparseFile (filename, *pMatrix, "", "", debug);
      }

      /// \brief Print the sparse matrix in Matrix Market format, with
      ///   comments.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to an output
      /// stream, using the Matrix Market "coordinate" format.  MPI
      /// Proc 0 is the only MPI process that writes to the output
      /// stream.
      ///
      /// \param out [out] Name of the output stream to which to write
      ///   the given sparse matrix.  The matrix is distributed, but
      ///   only Proc 0 writes to the output stream.
      ///
      /// \param matrix [in] The sparse matrix to write to the given
      ///   output stream.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Process 0.  This will cause out-of-memory errors
      ///   if the matrix is too big to fit on one process.  This will
      ///   be fixed in the future.
      static void
      writeSparse (std::ostream& out,
                   const sparse_matrix_type& matrix,
                   const std::string& matrixName,
                   const std::string& matrixDescription,
                   const bool debug=false)
      {
        using Teuchos::ArrayView;
        using Teuchos::Comm;
        using Teuchos::FancyOStream;
        using Teuchos::getFancyOStream;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcpFromRef;
        using std::cerr;
        using std::endl;
        using ST = scalar_type;
        using LO = local_ordinal_type;
        using GO = global_ordinal_type;
        using STS = typename Teuchos::ScalarTraits<ST>;

        // Make the output stream write floating-point numbers in
        // scientific notation.  It will politely put the output
        // stream back to its state on input, when this scope
        // terminates.
        Teuchos::SetScientific<ST> sci (out);

        // Get the matrix's communicator.
        trcp_tcomm_t comm = matrix.getComm ();
        TEUCHOS_TEST_FOR_EXCEPTION(
          comm.is_null (), std::invalid_argument,
          "The input matrix's communicator (Teuchos::Comm object) is null.");
        const int myRank = comm->getRank ();

        // Optionally, make a stream for debugging output.
        RCP<FancyOStream> err =
          debug ? getFancyOStream (rcpFromRef (std::cerr)) : null;
        if (debug) {
          std::ostringstream os;
          os << myRank << ": writeSparse" << endl;
          *err << os.str ();
          comm->barrier ();
          os << "-- " << myRank << ": past barrier" << endl;
          *err << os.str ();
        }

        // Whether to print debugging output to stderr.
        const bool debugPrint = debug && myRank == 0;

        RCP<const map_type> rowMap = matrix.getRowMap ();
        RCP<const map_type> colMap = matrix.getColMap ();
        RCP<const map_type> domainMap = matrix.getDomainMap ();
        RCP<const map_type> rangeMap = matrix.getRangeMap ();

        const global_size_t numRows = rangeMap->getGlobalNumElements ();
        const global_size_t numCols = domainMap->getGlobalNumElements ();

        if (debug && myRank == 0) {
          std::ostringstream os;
          os << "-- Input sparse matrix is:"
             << "---- " << numRows << " x " << numCols << endl
             << "---- "
             << (matrix.isGloballyIndexed() ? "Globally" : "Locally")
             << " indexed." << endl
             << "---- Its row map has " << rowMap->getGlobalNumElements ()
             << " elements." << endl
             << "---- Its col map has " << colMap->getGlobalNumElements ()
             << " elements." << endl;
          *err << os.str ();
        }
        // Make the "gather" row map, where Proc 0 owns all rows and
        // the other procs own no rows.
        const size_t localNumRows = (myRank == 0) ? numRows : 0;
        if (debug) {
          std::ostringstream os;
          os << "-- " << myRank << ": making gatherRowMap" << endl;
          *err << os.str ();
        }
        RCP<const map_type> gatherRowMap =
          Details::computeGatherMap (rowMap, err, debug);

        // Since the matrix may in general be non-square, we need to
        // make a column map as well.  In this case, the column map
        // contains all the columns of the original matrix, because we
        // are gathering the whole matrix onto Proc 0.  We call
        // computeGatherMap to preserve the original order of column
        // indices over all the processes.
        const size_t localNumCols = (myRank == 0) ? numCols : 0;
        RCP<const map_type> gatherColMap =
          Details::computeGatherMap (colMap, err, debug);

        // Current map is the source map, gather map is the target map.
        typedef Import<LO, GO, node_type> import_type;
        import_type importer (rowMap, gatherRowMap);

        // Create a new CrsMatrix to hold the result of the import.
        // The constructor needs a column map as well as a row map,
        // for the case that the matrix is not square.
        RCP<sparse_matrix_type> newMatrix =
          rcp (new sparse_matrix_type (gatherRowMap, gatherColMap,
                                       static_cast<size_t> (0)));
        // Import the sparse matrix onto Proc 0.
        newMatrix->doImport (matrix, importer, INSERT);

        // fillComplete() needs the domain and range maps for the case
        // that the matrix is not square.
        {
          RCP<const map_type> gatherDomainMap =
            rcp (new map_type (numCols, localNumCols,
                               domainMap->getIndexBase (),
                               comm));
          RCP<const map_type> gatherRangeMap =
            rcp (new map_type (numRows, localNumRows,
                               rangeMap->getIndexBase (),
                               comm));
          newMatrix->fillComplete (gatherDomainMap, gatherRangeMap);
        }

        if (debugPrint) {
          cerr << "-- Output sparse matrix is:"
               << "---- " << newMatrix->getRangeMap ()->getGlobalNumElements ()
               << " x "
               << newMatrix->getDomainMap ()->getGlobalNumElements ()
               << " with "
               << newMatrix->getGlobalNumEntries () << " entries;" << endl
               << "---- "
               << (newMatrix->isGloballyIndexed () ? "Globally" : "Locally")
               << " indexed." << endl
               << "---- Its row map has "
               << newMatrix->getRowMap ()->getGlobalNumElements ()
               << " elements, with index base "
               << newMatrix->getRowMap ()->getIndexBase () << "." << endl
               << "---- Its col map has "
               << newMatrix->getColMap ()->getGlobalNumElements ()
               << " elements, with index base "
               << newMatrix->getColMap ()->getIndexBase () << "." << endl
               << "---- Element count of output matrix's column Map may differ "
               << "from that of the input matrix's column Map, if some columns "
               << "of the matrix contain no entries." << endl;
        }

        //
        // Print the metadata and the matrix entries on Rank 0.
        //
        if (myRank == 0) {
          // Print the Matrix Market banner line.  CrsMatrix stores
          // data nonsymmetrically ("general").  This implies that
          // readSparse() on a symmetrically stored input file,
          // followed by writeSparse() on the resulting sparse matrix,
          // will result in an output file with a different banner
          // line than the original input file.
          out << "%%MatrixMarket matrix coordinate "
              << (STS::isComplex ? "complex" : "real")
              << " general" << endl;

          // Print comments (the matrix name and / or description).
          if (matrixName != "") {
            printAsComment (out, matrixName);
          }
          if (matrixDescription != "") {
            printAsComment (out, matrixDescription);
          }

          // Print the Matrix Market header (# rows, # columns, #
          // nonzeros).  Use the range resp. domain map for the number
          // of rows resp. columns, since Tpetra::CrsMatrix uses the
          // column map for the number of columns.  That only
          // corresponds to the "linear-algebraic" number of columns
          // when the column map is uniquely owned (a.k.a. one-to-one),
          // which only happens if the matrix is (block) diagonal.
          out << newMatrix->getRangeMap ()->getGlobalNumElements () << " "
              << newMatrix->getDomainMap ()->getGlobalNumElements () << " "
              << newMatrix->getGlobalNumEntries () << endl;

          // The Matrix Market format expects one-based row and column
          // indices.  We'll convert the indices on output from
          // whatever index base they use to one-based indices.
          const GO rowIndexBase = gatherRowMap->getIndexBase ();
          const GO colIndexBase = newMatrix->getColMap()->getIndexBase ();
          //
          // Print the entries of the matrix.
          //
          // newMatrix can never be globally indexed, since we called
          // fillComplete() on it.  We include code for both cases
          // (globally or locally indexed) just in case that ever
          // changes.
          if (newMatrix->isGloballyIndexed()) {
            // We know that the "gather" row Map is contiguous, so we
            // don't need to get the list of GIDs.
            const GO minAllGlobalIndex = gatherRowMap->getMinAllGlobalIndex ();
            const GO maxAllGlobalIndex = gatherRowMap->getMaxAllGlobalIndex ();
            for (GO globalRowIndex = minAllGlobalIndex;
                 globalRowIndex <= maxAllGlobalIndex; // inclusive range
                 ++globalRowIndex) {
              typename sparse_matrix_type::global_inds_host_view_type ind;
              typename sparse_matrix_type::values_host_view_type val;
              newMatrix->getGlobalRowView (globalRowIndex, ind, val);
              for (size_t ii = 0; ii < ind.extent(0); ii++) {
                const GO globalColIndex = ind(ii);
                // Convert row and column indices to 1-based.
                // This works because the global index type is signed.
                out << (globalRowIndex + 1 - rowIndexBase) << " "
                    << (globalColIndex + 1 - colIndexBase) << " ";
                if (STS::isComplex) {
                  out << STS::real (val(ii)) << " " << STS::imag (val(ii));
                } else {
                  out << val(ii);
                }
                out << endl;
              } // For each entry in the current row
            } // For each row of the "gather" matrix
          }
          else { // newMatrix is locally indexed
            using OTG = Teuchos::OrdinalTraits<GO>;
            for (LO localRowIndex = gatherRowMap->getMinLocalIndex();
                 localRowIndex <= gatherRowMap->getMaxLocalIndex();
                 ++localRowIndex) {
              // Convert from local to global row index.
              const GO globalRowIndex =
                gatherRowMap->getGlobalElement (localRowIndex);
              TEUCHOS_TEST_FOR_EXCEPTION(
                globalRowIndex == OTG::invalid(), std::logic_error,
                "Failed to convert the supposed local row index "
                << localRowIndex << " into a global row index.  "
                "Please report this bug to the Tpetra developers.");
              typename sparse_matrix_type::local_inds_host_view_type ind;
              typename sparse_matrix_type::values_host_view_type val;
              newMatrix->getLocalRowView (localRowIndex, ind, val);
              for (size_t ii = 0; ii < ind.extent(0); ii++) {
                // Convert the column index from local to global.
                const GO globalColIndex =
                  newMatrix->getColMap()->getGlobalElement (ind(ii));
                TEUCHOS_TEST_FOR_EXCEPTION(
                  globalColIndex == OTG::invalid(), std::logic_error,
                  "On local row " << localRowIndex << " of the sparse matrix: "
                  "Failed to convert the supposed local column index "
                  << ind(ii) << " into a global column index.  Please report "
                  "this bug to the Tpetra developers.");
                // Convert row and column indices to 1-based.
                // This works because the global index type is signed.
                out << (globalRowIndex + 1 - rowIndexBase) << " "
                    << (globalColIndex + 1 - colIndexBase) << " ";
                if (STS::isComplex) {
                  out << STS::real (val(ii)) << " " << STS::imag (val(ii));
                } else {
                  out << val(ii);
                }
                out << endl;
              } // For each entry in the current row
            } // For each row of the "gather" matrix
          } // Whether the "gather" matrix is locally or globally indexed
        } // If my process' rank is 0
      }

      //! Only for backwards compatibility; prefer the overload above.
      static void
      writeSparse (std::ostream& out,
                   const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                   const std::string& matrixName,
                   const std::string& matrixDescription,
                   const bool debug=false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (pMatrix.is_null (), std::invalid_argument,
           "The input matrix is null.");
        writeSparse (out, *pMatrix, matrixName, matrixDescription, debug);
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given output stream.
      ///
      /// Write the given Tpetra::CrsGraph sparse graph to an output
      /// stream, using the Matrix Market "coordinate" format.
      ///
      /// \param out [out] Name of the output stream to which to write
      ///   the given sparse graph.  The graph is distributed, but
      ///   only Process 0 in the graph's communicator writes to the
      ///   output stream.  Thus, the output stream need only be valid
      ///   (writeable) on Process 0.
      ///
      /// \param graph [in] The sparse graph to write.
      ///
      /// \param graphName [in] Name of the graph, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param graphDescription [in] Graph description, to print in
      ///   the comments section of the output stream.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Process 0 of the graph's
      ///   communicator.  False (do NOT print) by default.
      ///
      /// \warning The current implementation gathers the whole graph
      ///   onto Process 0 in the graph's communicator.  This will
      ///   cause out-of-memory errors if the graph is too big to fit
      ///   on one process.  This will be fixed in the future.
      static void
      writeSparseGraph (std::ostream& out,
                        const crs_graph_type& graph,
                        const std::string& graphName,
                        const std::string& graphDescription,
                        const bool debug=false)
      {
        using Teuchos::ArrayView;
        using Teuchos::Comm;
        using Teuchos::FancyOStream;
        using Teuchos::getFancyOStream;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcpFromRef;
        using std::cerr;
        using std::endl;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;

        // Get the graph's communicator.  Processes on which the
        // graph's Map or communicator is null don't participate in
        // this operation.  This function shouldn't even be called on
        // those processes.
        auto rowMap = graph.getRowMap ();
        if (rowMap.is_null ()) {
          return;
        }
        auto comm = rowMap->getComm ();
        if (comm.is_null ()) {
          return;
        }
        const int myRank = comm->getRank ();

        // Optionally, make a stream for debugging output.
        RCP<FancyOStream> err =
          debug ? getFancyOStream (rcpFromRef (std::cerr)) : null;
        if (debug) {
          std::ostringstream os;
          os << myRank << ": writeSparseGraph" << endl;
          *err << os.str ();
          comm->barrier ();
          os << "-- " << myRank << ": past barrier" << endl;
          *err << os.str ();
        }

        // Whether to print debugging output to stderr.
        const bool debugPrint = debug && myRank == 0;

        // We've already gotten the rowMap above.
        auto colMap = graph.getColMap ();
        auto domainMap = graph.getDomainMap ();
        auto rangeMap = graph.getRangeMap ();

        const global_size_t numRows = rangeMap->getGlobalNumElements ();
        const global_size_t numCols = domainMap->getGlobalNumElements ();

        if (debug && myRank == 0) {
          std::ostringstream os;
          os << "-- Input sparse graph is:"
             << "---- " << numRows << " x " << numCols << " with "
             << graph.getGlobalNumEntries () << " entries;" << endl
             << "---- "
             << (graph.isGloballyIndexed () ? "Globally" : "Locally")
             << " indexed." << endl
             << "---- Its row Map has " << rowMap->getGlobalNumElements ()
             << " elements." << endl
             << "---- Its col Map has " << colMap->getGlobalNumElements ()
             << " elements." << endl;
          *err << os.str ();
        }
        // Make the "gather" row map, where Proc 0 owns all rows and
        // the other procs own no rows.
        const size_t localNumRows = (myRank == 0) ? numRows : 0;
        if (debug) {
          std::ostringstream os;
          os << "-- " << myRank << ": making gatherRowMap" << endl;
          *err << os.str ();
        }
        auto gatherRowMap = Details::computeGatherMap (rowMap, err, debug);

        // Since the graph may in general be non-square, we need to
        // make a column map as well.  In this case, the column map
        // contains all the columns of the original graph, because we
        // are gathering the whole graph onto Proc 0.  We call
        // computeGatherMap to preserve the original order of column
        // indices over all the processes.
        const size_t localNumCols = (myRank == 0) ? numCols : 0;
        auto gatherColMap = Details::computeGatherMap (colMap, err, debug);

        // Current map is the source map, gather map is the target map.
        Import<LO, GO, node_type> importer (rowMap, gatherRowMap);

        // Create a new CrsGraph to hold the result of the import.
        // The constructor needs a column map as well as a row map,
        // for the case that the graph is not square.
        crs_graph_type newGraph (gatherRowMap, gatherColMap,
                                 static_cast<size_t> (0));
        // Import the sparse graph onto Proc 0.
        newGraph.doImport (graph, importer, INSERT);

        // fillComplete() needs the domain and range maps for the case
        // that the graph is not square.
        {
          RCP<const map_type> gatherDomainMap =
            rcp (new map_type (numCols, localNumCols,
                               domainMap->getIndexBase (),
                               comm));
          RCP<const map_type> gatherRangeMap =
            rcp (new map_type (numRows, localNumRows,
                               rangeMap->getIndexBase (),
                               comm));
          newGraph.fillComplete (gatherDomainMap, gatherRangeMap);
        }

        if (debugPrint) {
          cerr << "-- Output sparse graph is:"
               << "---- " << newGraph.getRangeMap ()->getGlobalNumElements ()
               << " x "
               << newGraph.getDomainMap ()->getGlobalNumElements ()
               << " with "
               << newGraph.getGlobalNumEntries () << " entries;" << endl
               << "---- "
               << (newGraph.isGloballyIndexed () ? "Globally" : "Locally")
               << " indexed." << endl
               << "---- Its row map has "
               << newGraph.getRowMap ()->getGlobalNumElements ()
               << " elements, with index base "
               << newGraph.getRowMap ()->getIndexBase () << "." << endl
               << "---- Its col map has "
               << newGraph.getColMap ()->getGlobalNumElements ()
               << " elements, with index base "
               << newGraph.getColMap ()->getIndexBase () << "." << endl
               << "---- Element count of output graph's column Map may differ "
               << "from that of the input matrix's column Map, if some columns "
               << "of the matrix contain no entries." << endl;
        }

        //
        // Print the metadata and the graph entries on Process 0 of
        // the graph's communicator.
        //
        if (myRank == 0) {
          // Print the Matrix Market banner line.  CrsGraph stores
          // data nonsymmetrically ("general").  This implies that
          // readSparseGraph() on a symmetrically stored input file,
          // followed by writeSparseGraph() on the resulting sparse
          // graph, will result in an output file with a different
          // banner line than the original input file.
          out << "%%MatrixMarket matrix coordinate pattern general" << endl;

          // Print comments (the graph name and / or description).
          if (graphName != "") {
            printAsComment (out, graphName);
          }
          if (graphDescription != "") {
            printAsComment (out, graphDescription);
          }

          // Print the Matrix Market header (# rows, # columns, #
          // stored entries).  Use the range resp. domain map for the
          // number of rows resp. columns, since Tpetra::CrsGraph uses
          // the column map for the number of columns.  That only
          // corresponds to the "linear-algebraic" number of columns
          // when the column map is uniquely owned
          // (a.k.a. one-to-one), which only happens if the graph is
          // block diagonal (one block per process).
          out << newGraph.getRangeMap ()->getGlobalNumElements () << " "
              << newGraph.getDomainMap ()->getGlobalNumElements () << " "
              << newGraph.getGlobalNumEntries () << endl;

          // The Matrix Market format expects one-based row and column
          // indices.  We'll convert the indices on output from
          // whatever index base they use to one-based indices.
          const GO rowIndexBase = gatherRowMap->getIndexBase ();
          const GO colIndexBase = newGraph.getColMap()->getIndexBase ();
          //
          // Print the entries of the graph.
          //
          // newGraph can never be globally indexed, since we called
          // fillComplete() on it.  We include code for both cases
          // (globally or locally indexed) just in case that ever
          // changes.
          if (newGraph.isGloballyIndexed ()) {
            // We know that the "gather" row Map is contiguous, so we
            // don't need to get the list of GIDs.
            const GO minAllGlobalIndex = gatherRowMap->getMinAllGlobalIndex ();
            const GO maxAllGlobalIndex = gatherRowMap->getMaxAllGlobalIndex ();
            for (GO globalRowIndex = minAllGlobalIndex;
                 globalRowIndex <= maxAllGlobalIndex; // inclusive range
                 ++globalRowIndex) {
              typename crs_graph_type::global_inds_host_view_type ind;
              newGraph.getGlobalRowView (globalRowIndex, ind);
              for (size_t ii = 0; ii < ind.extent(0); ii++) {
                const GO globalColIndex = ind(ii);
                // Convert row and column indices to 1-based.
                // This works because the global index type is signed.
                out << (globalRowIndex + 1 - rowIndexBase) << " "
                    << (globalColIndex + 1 - colIndexBase) << " ";
                out << endl;
              } // For each entry in the current row
            } // For each row of the "gather" graph
          }
          else { // newGraph is locally indexed
            typedef Teuchos::OrdinalTraits<GO> OTG;
            for (LO localRowIndex = gatherRowMap->getMinLocalIndex ();
                 localRowIndex <= gatherRowMap->getMaxLocalIndex ();
                 ++localRowIndex) {
              // Convert from local to global row index.
              const GO globalRowIndex =
                gatherRowMap->getGlobalElement (localRowIndex);
              TEUCHOS_TEST_FOR_EXCEPTION
                (globalRowIndex == OTG::invalid (), std::logic_error, "Failed "
                 "to convert the supposed local row index " << localRowIndex <<
                 " into a global row index.  Please report this bug to the "
                 "Tpetra developers.");
              typename crs_graph_type::local_inds_host_view_type ind;
              newGraph.getLocalRowView (localRowIndex, ind);
              for (size_t ii = 0; ii < ind.extent(0); ii++) {
                // Convert the column index from local to global.
                const GO globalColIndex =
                  newGraph.getColMap ()->getGlobalElement (ind(ii));
                TEUCHOS_TEST_FOR_EXCEPTION(
                  globalColIndex == OTG::invalid(), std::logic_error,
                  "On local row " << localRowIndex << " of the sparse graph: "
                  "Failed to convert the supposed local column index "
                  << ind(ii) << " into a global column index.  Please report "
                  "this bug to the Tpetra developers.");
                // Convert row and column indices to 1-based.
                // This works because the global index type is signed.
                out << (globalRowIndex + 1 - rowIndexBase) << " "
                    << (globalColIndex + 1 - colIndexBase) << " ";
                out << endl;
              } // For each entry in the current row
            } // For each row of the "gather" graph
          } // Whether the "gather" graph is locally or globally indexed
        } // If my process' rank is 0
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given output stream, with no comments.
      ///
      /// See the above five-argument version of this function for
      /// full documentation.
      static void
      writeSparseGraph (std::ostream& out,
                        const crs_graph_type& graph,
                        const bool debug=false)
      {
        writeSparseGraph (out, graph, "", "", debug);
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given file (by filename).
      ///
      /// Write the given Tpetra::CrsGraph sparse graph to the given
      /// file, using the Matrix Market "coordinate" format.  Process
      /// 0 in the graph's communicator is the only MPI process that
      /// opens or writes to the file.  Include the graph name and
      /// description in the comments section of the file (after the
      /// initial banner line, but before the graph metadata and
      /// data).
      ///
      /// \param filename [in] Name of the file to which to write the
      ///   given sparse graph.  The graph is distributed, but only
      ///   Process 0 in the graph's communicator opens the file and
      ///   writes to it.
      ///
      /// \param graph [in] The sparse graph to write.
      ///
      /// \param graphName [in] Name of the graph, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param graphDescription [in] Graph description, to print in
      ///   the comments section of the output stream.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Process 0 of the graph's
      ///   communicator.  False (do NOT print) by default.
      ///
      /// \warning The current implementation gathers the whole graph
      ///   onto Process 0 in the graph's communicator.  This will
      ///   cause out-of-memory errors if the graph is too big to fit
      ///   on one process.  This will be fixed in the future.
      static void
      writeSparseGraphFile (const std::string& filename,
                            const crs_graph_type& graph,
                            const std::string& graphName,
                            const std::string& graphDescription,
                            const bool debug=false)
      {
        auto comm = graph.getComm ();
        if (comm.is_null ()) {
          // Processes on which the communicator is null shouldn't
          // even call this function.  The convention is that
          // processes on which the object's communicator is null do
          // not participate in collective operations involving the
          // object.
          return;
        }
        const int myRank = comm->getRank ();

        auto out = Writer::openOutFileOnRankZero(comm, filename, myRank, true);

        writeSparseGraph (out, graph, graphName, graphDescription, debug);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeSparseGraph() throws
        // an exception.
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given file (by filename), with no comments.
      ///
      /// See the above five-argument overload for full documentation.
      static void
      writeSparseGraphFile (const std::string& filename,
                            const crs_graph_type& graph,
                            const bool debug=false)
      {
        writeSparseGraphFile (filename, graph, "", "", debug);
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given file (by filename), taking the graph by Teuchos::RCP.
      ///
      /// This is just a convenience for users who don't want to
      /// remember to dereference the Teuchos::RCP.  For
      /// documentation, see the above overload of this function that
      /// takes the graph by const reference, rather than by
      /// Teuchos::RCP.
      static void
      writeSparseGraphFile (const std::string& filename,
                            const Teuchos::RCP<const crs_graph_type>& pGraph,
                            const std::string& graphName,
                            const std::string& graphDescription,
                            const bool debug=false)
      {
        writeSparseGraphFile (filename, *pGraph, graphName, graphDescription, debug);
      }

      /// \brief Print the sparse graph in Matrix Market format to the
      ///   given file (by filename), with no comments, taking the
      ///   graph by Teuchos::RCP.
      ///
      /// This is just a convenience for users who don't want to
      /// remember to dereference the Teuchos::RCP.  For
      /// documentation, see the above overload of this function that
      /// takes the graph by const reference, rather than by
      /// Teuchos::RCP.
      static void
      writeSparseGraphFile (const std::string& filename,
                            const Teuchos::RCP<const crs_graph_type>& pGraph,
                            const bool debug=false)
      {
        writeSparseGraphFile (filename, *pGraph, "", "", debug);
      }

      /// \brief Print the sparse matrix in Matrix Market format.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to an output
      /// stream, using the Matrix Market "coordinate" format.  MPI
      /// Proc 0 is the only MPI process that writes to the output
      /// stream.
      ///
      /// \param out [out] Name of the output stream to which to write
      ///   the given sparse matrix.  The matrix is distributed, but
      ///   only Proc 0 writes to the output stream.
      ///
      /// \param matrix [in] The sparse matrix to write to the given
      ///   output stream.
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeSparse (std::ostream& out,
                   const sparse_matrix_type& matrix,
                   const bool debug=false)
      {
        writeSparse (out, matrix, "", "", debug);
      }

      //! Only for backwards compatibility; prefer the overload above.
      static void
      writeSparse (std::ostream& out,
                   const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                   const bool debug=false)
      {
        writeSparse (out, *pMatrix, "", "", debug);
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   matrix name and description.
      ///
      /// Write the given Tpetra::MultiVector matrix to the given
      /// file, using the Matrix Market "array" format for dense
      /// matrices.  MPI Process 0 is the only MPI process that opens
      /// or writes to the file.
      ///
      /// This is the preferred overload of writeDenseFile.  It is
      /// used to implement all other overloads of writeDenseFile.
      ///
      /// \param filename [in] Name of the output file to create (on
      ///   MPI Proc 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output file.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output file.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output file.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \param err [out] If nonnull, print any error messages to it.
      ///
      /// \param dbg [out] If nonnull, print copious debugging output to it.
      static void
      writeDenseFile (const std::string& filename,
                      const multivector_type& X,
                      const std::string& matrixName,
                      const std::string& matrixDescription,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                      const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        trcp_tcomm_t comm = Writer::getComm(X.getMap());
        const int myRank = Writer::getRank(comm);

        auto out = Writer::openOutFileOnRankZero(comm, filename, myRank, true);

        writeDense (out, X, matrixName, matrixDescription, err, dbg);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeDense() throws an
        // exception.
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   matrix name and description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDenseFile().
      static void
      writeDenseFile (const std::string& filename,
                      const Teuchos::RCP<const multivector_type>& X,
                      const std::string& matrixName,
                      const std::string& matrixDescription,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                      const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          X.is_null (), std::invalid_argument, "Tpetra::MatrixMarket::"
          "writeDenseFile: The input MultiVector X is null.");
        writeDenseFile (filename, *X, matrixName, matrixDescription, err, dbg);
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   no matrix name or description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDenseFile().
      static void
      writeDenseFile (const std::string& filename,
                      const multivector_type& X,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                      const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        writeDenseFile (filename, X, "", "", err, dbg);
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   no matrix name or description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDenseFile().
      static void
      writeDenseFile (const std::string& filename,
                      const Teuchos::RCP<const multivector_type>& X,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                      const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          X.is_null (), std::invalid_argument, "Tpetra::MatrixMarket::"
          "writeDenseFile: The input MultiVector X is null.");
        writeDenseFile (filename, *X, err, dbg);
      }


      /// \brief Print the multivector in Matrix Market format, with
      ///   matrix name and description.
      ///
      /// Write the given Tpetra::MultiVector matrix to an output
      /// stream, using the Matrix Market "array" format for dense
      /// matrices.  MPI Process 0 in the given MultiVector's
      /// communicator is the only MPI process that may write to the
      /// given output stream.
      ///
      /// This is the preferred overload of writeDense().  It is used
      /// to implement all other overloads of writeDense(), and is
      /// also used to implement all overloads of writeDenseFile().
      ///
      /// \param out [out] The output stream to which to write (on MPI
      ///   Process 0 only).
      ///
      /// \param X [in] The Tpetra::MultiVector to write to the given
      ///   output file \c out.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \param err [out] If nonnull, print any error messages to it.
      ///
      /// \param dbg [out] If nonnull, print copious debugging output to it.
      static void
      writeDense (std::ostream& out,
                  const multivector_type& X,
                  const std::string& matrixName,
                  const std::string& matrixDescription,
                  const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                  const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        using Teuchos::Comm;
        using Teuchos::outArg;
        using Teuchos::REDUCE_MAX;
        using Teuchos::reduceAll;
        using std::endl;

        trcp_tcomm_t comm = Writer::getComm(X.getMap());
        const int myRank = Writer::getRank(comm);

        // If the caller provides a nonnull debug output stream, we
        // print debugging output to it.  This is a local thing; we
        // don't have to check across processes.
        const bool debug = ! dbg.is_null ();
        if (debug) {
          dbg->pushTab ();
          std::ostringstream os;
          os << myRank << ": writeDense" << endl;
          *dbg << os.str ();
          dbg->pushTab ();
        }
        // Print the Matrix Market header.
        writeDenseHeader (out, X, matrixName, matrixDescription, err, dbg);

        // Print each column one at a time.  This is a (perhaps)
        // temporary fix for Bug 6288.
        const size_t numVecs = X.getNumVectors ();
        for (size_t j = 0; j < numVecs; ++j) {
          writeDenseColumn (out, * (X.getVector (j)), err, dbg);
        }

        if (debug) {
          dbg->popTab ();
          std::ostringstream os;
          os << myRank << ": writeDense: Done" << endl;
          *dbg << os.str ();
          dbg->popTab ();
        }
      }

    private:
      /**
       * @brief Open a file only on rank zero, possibly throwing if the stream is invalid.
       *
       * @note On processes that are not the rank zero process, the stream is left uninitialized.
       */
      static std::ofstream openOutFileOnRankZero(
        const trcp_tcomm_t& comm,
        const std::string& filename, const int rank, const bool safe = true,
        const std::ios_base::openmode mode = std::ios_base::out
      ){
        // Placeholder for the output stream.
        std::ofstream out;

        // State that will make all ranks throw if the root rank wasn't able to open the stream (using @c int for broadcasting).
        int all_should_stop = 0;

        // Try to open the file and update the state.
        if(rank == 0) {
          out.open(filename, mode);
          all_should_stop = !out && safe;
        }

        // Broadcast the stream state and throw from all ranks if needed.
        if(comm) Teuchos::broadcast(*comm, 0, &all_should_stop);

        TEUCHOS_TEST_FOR_EXCEPTION(
          all_should_stop,
          std::runtime_error,
          "Could not open output file '" + filename + "' on root rank 0."
        );

        return out;
      }

      /// \brief Print the MultiVector's Matrix Market header.
      ///
      /// Write the given Tpetra::MultiVector's Matrix Market header
      /// to the given output stream.  The header includes metadata
      /// like the storage format and dimensions, but not the actual
      /// entries.  MPI Process 0 is the only MPI process that may
      /// write to the output stream.
      ///
      /// \param out [out] The output stream to which to write (on MPI
      ///   Process 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output file.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \param err [out] If nonnull, print any error messages to it.
      ///
      /// \param dbg [out] If nonnull, print copious debugging output to it.
      static void
      writeDenseHeader (std::ostream& out,
                        const multivector_type& X,
                        const std::string& matrixName,
                        const std::string& matrixDescription,
                        const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                        const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        using Teuchos::Comm;
        using Teuchos::outArg;
        using Teuchos::RCP;
        using Teuchos::REDUCE_MAX;
        using Teuchos::reduceAll;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;
        const char prefix[] = "Tpetra::MatrixMarket::writeDenseHeader: ";

        trcp_tcomm_t comm = Writer::getComm(X.getMap());
        const int myRank = Writer::getRank(comm);
        int lclErr = 0; // whether this MPI process has seen an error
        int gblErr = 0; // whether we know if some MPI process has seen an error

        // If the caller provides a nonnull debug output stream, we
        // print debugging output to it.  This is a local thing; we
        // don't have to check across processes.
        const bool debug = ! dbg.is_null ();

        if (debug) {
          dbg->pushTab ();
          std::ostringstream os;
          os << myRank << ": writeDenseHeader" << endl;
          *dbg << os.str ();
          dbg->pushTab ();
        }

        //
        // Process 0: Write the MatrixMarket header.
        //
        if (myRank == 0) {
          try {
            // Print the Matrix Market header.  MultiVector stores data
            // nonsymmetrically, hence "general" in the banner line.
            // Print first to a temporary string output stream, and then
            // write it to the main output stream, so that at least the
            // header output has transactional semantics.  We can't
            // guarantee transactional semantics for the whole output,
            // since that would not be memory scalable.  (This could be
            // done in the file system by using a temporary file; we
            // don't do this, but users could.)
            std::ostringstream hdr;
            {
              std::string dataType;
              if (STS::isComplex) {
                dataType = "complex";
              } else if (STS::isOrdinal) {
                dataType = "integer";
              } else {
                dataType = "real";
              }
              hdr << "%%MatrixMarket matrix array " << dataType << " general"
                  << endl;
            }

            // Print comments (the matrix name and / or description).
            if (matrixName != "") {
              printAsComment (hdr, matrixName);
            }
            if (matrixDescription != "") {
              printAsComment (hdr, matrixDescription);
            }
            // Print the Matrix Market dimensions header for dense matrices.
            hdr << X.getGlobalLength () << " " << X.getNumVectors () << endl;

            // Write the MatrixMarket header to the output stream.
            out << hdr.str ();
          } catch (std::exception& e) {
            if (! err.is_null ()) {
              *err << prefix << "While writing the Matrix Market header, "
                "Process 0 threw an exception: " << e.what () << endl;
            }
            lclErr = 1;
          }
        } // if I am Process 0

        // Establish global agreement on the error state.  It wouldn't
        // be good for other processes to keep going, if Process 0
        // finds out that it can't write to the given output stream.
        reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
        TEUCHOS_TEST_FOR_EXCEPTION(
          gblErr == 1, std::runtime_error, prefix << "Some error occurred "
          "which prevented this method from completing.");

        if (debug) {
          dbg->popTab ();
          *dbg << myRank << ": writeDenseHeader: Done" << endl;
          dbg->popTab ();
        }
      }

      /// \brief Print a single column of the given MultiVector in
      ///   Matrix Market format.
      ///
      /// Write the given Tpetra::MultiVector matrix (which must have
      /// only one column!) to an output stream, using the Matrix
      /// Market "array" format for dense matrices.  MPI Process 0 is
      /// the only MPI process that writes to the output stream.
      ///
      /// \param out [out] The output stream to which to write (on MPI
      ///   Process 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output file.
      ///
      /// \param err [out] If nonnull, print any error messages to it.
      ///
      /// \param dbg [out] If nonnull, print copious debugging output to it.
      static void
      writeDenseColumn (std::ostream& out,
                        const multivector_type& X,
                        const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                        const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        using Teuchos::arcp;
        using Teuchos::Array;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::Comm;
        using Teuchos::CommRequest;
        using Teuchos::ireceive;
        using Teuchos::isend;
        using Teuchos::outArg;
        using Teuchos::REDUCE_MAX;
        using Teuchos::reduceAll;
        using Teuchos::RCP;
        using Teuchos::TypeNameTraits;
        using Teuchos::wait;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        const Comm<int>& comm = * (X.getMap ()->getComm ());
        const int myRank = comm.getRank ();
        const int numProcs = comm.getSize ();
        int lclErr = 0; // whether this MPI process has seen an error
        int gblErr = 0; // whether we know if some MPI process has seen an error

        // If the caller provides a nonnull debug output stream, we
        // print debugging output to it.  This is a local thing; we
        // don't have to check across processes.
        const bool debug = ! dbg.is_null ();

        if (debug) {
          dbg->pushTab ();
          std::ostringstream os;
          os << myRank << ": writeDenseColumn" << endl;
          *dbg << os.str ();
          dbg->pushTab ();
        }

        // Make the output stream write floating-point numbers in
        // scientific notation.  It will politely put the output
        // stream back to its state on input, when this scope
        // terminates.
        Teuchos::SetScientific<scalar_type> sci (out);

        const size_t myNumRows = X.getLocalLength ();
        const size_t numCols = X.getNumVectors ();
        // Use a different tag for the "size" messages than for the
        // "data" messages, in order to help us debug any mix-ups.
        const int sizeTag = 1337;
        const int dataTag = 1338;

        // Process 0 pipelines nonblocking receives with file output.
        //
        // Constraints:
        //   - Process 0 can't post a receive for another process'
        //     actual data, until it posts and waits on the receive
        //     from that process with the amount of data to receive.
        //     (We could just post receives with a max data size, but
        //     I feel uncomfortable about that.)
        //   - The C++ standard library doesn't allow nonblocking
        //     output to an std::ostream.  (Thus, we have to start a
        //     receive or send before starting the write, and hope
        //     that MPI completes it in the background.)
        //
        // Process 0: Post receive-size receives from Processes 1 and 2.
        // Process 1: Post send-size send to Process 0.
        // Process 2: Post send-size send to Process 0.
        //
        // All processes: Pack my entries.
        //
        // Process 1:
        //   - Post send-data send to Process 0.
        //   - Wait on my send-size send to Process 0.
        //
        // Process 0:
        //   - Print MatrixMarket header.
        //   - Print my entries.
        //   - Wait on receive-size receive from Process 1.
        //   - Post receive-data receive from Process 1.
        //
        // For each process p = 1, 2, ... numProcs-1:
        //   If I am Process 0:
        //     - Post receive-size receive from Process p + 2
        //     - Wait on receive-size receive from Process p + 1
        //     - Post receive-data receive from Process p + 1
        //     - Wait on receive-data receive from Process p
        //     - Write data from Process p.
        //   Else if I am Process p:
        //     - Wait on my send-data send.
        //   Else if I am Process p+1:
        //     - Post send-data send to Process 0.
        //     - Wait on my send-size send.
        //   Else if I am Process p+2:
        //     - Post send-size send to Process 0.
        //
        // Pipelining has three goals here:
        //   1. Overlap communication (the receives) with file I/O
        //   2. Give Process 0 a chance to prepost some receives,
        //      before sends show up, by packing local data before
        //      posting sends
        //   3. Don't post _all_ receives or _all_ sends, because that
        //      wouldn't be memory scalable.  (Just because we can't
        //      see how much memory MPI consumes, doesn't mean that it
        //      doesn't consume any!)

        // These are used on every process.  sendReqSize[0] holds the
        // number of rows on this process, and sendReqBuf holds this
        // process' data.  Process 0 packs into sendReqBuf, but
        // doesn't send; it only uses that for printing.  All other
        // processes send both of these to Process 0.
        RCP<CommRequest<int> > sendReqSize, sendReqData;

        // These are used only on Process 0, for received data.  Keep
        // 3 of each, and treat the arrays as circular buffers.  When
        // receiving from Process p, the corresponding array index
        // here is p % 3.
        Array<ArrayRCP<size_t> > recvSizeBufs (3);
        Array<ArrayRCP<scalar_type> > recvDataBufs (3);
        Array<RCP<CommRequest<int> > > recvSizeReqs (3);
        Array<RCP<CommRequest<int> > > recvDataReqs (3);

        // Buffer for nonblocking send of the "send size."
        ArrayRCP<size_t> sendDataSize (1);
        sendDataSize[0] = myNumRows;

        if (myRank == 0) {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Post receive-size receives from "
              "Procs 1 and 2: tag = " << sizeTag << endl;
            *dbg << os.str ();
          }
          // Process 0: Post receive-size receives from Processes 1 and 2.
          recvSizeBufs[0].resize (1);
          // Set these three to an invalid value as a flag.  If we
          // don't get these messages, then the invalid value will
          // remain, so we can test for it.
          (recvSizeBufs[0])[0] = Teuchos::OrdinalTraits<size_t>::invalid ();
          recvSizeBufs[1].resize (1);
          (recvSizeBufs[1])[0] = Teuchos::OrdinalTraits<size_t>::invalid ();
          recvSizeBufs[2].resize (1);
          (recvSizeBufs[2])[0] = Teuchos::OrdinalTraits<size_t>::invalid ();
          if (numProcs > 1) {
            recvSizeReqs[1] =
              ireceive<int, size_t> (recvSizeBufs[1], 1, sizeTag, comm);
          }
          if (numProcs > 2) {
            recvSizeReqs[2] =
              ireceive<int, size_t> (recvSizeBufs[2], 2, sizeTag, comm);
          }
        }
        else if (myRank == 1 || myRank == 2) {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Post send-size send: size = "
               << sendDataSize[0] << ", tag = " << sizeTag << endl;
            *dbg << os.str ();
          }
          // Prime the pipeline by having Processes 1 and 2 start
          // their send-size sends.  We don't want _all_ the processes
          // to start their send-size sends, because that wouldn't be
          // memory scalable.
          sendReqSize = isend<int, size_t> (sendDataSize, 0, sizeTag, comm);
        }
        else {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Not posting my send-size send yet" << endl;
            *dbg << os.str ();
          }
        }

        //
        // Pack my entries, in column-major order.
        //
        if (debug) {
          std::ostringstream os;
          os << myRank << ": Pack my entries" << endl;
          *dbg << os.str ();
        }
        ArrayRCP<scalar_type> sendDataBuf;
        try {
          sendDataBuf = arcp<scalar_type> (myNumRows * numCols);
          X.get1dCopy (sendDataBuf (), myNumRows);
        }
        catch (std::exception& e) {
          lclErr = 1;
          if (! err.is_null ()) {
            std::ostringstream os;
            os << "Process " << myRank << ": Attempt to pack my MultiVector "
              "entries threw an exception: " << e.what () << endl;
            *err << os.str ();
          }
        }
        if (debug) {
          std::ostringstream os;
          os << myRank << ": Done packing my entries" << endl;
          *dbg << os.str ();
        }

        //
        // Process 1: post send-data send to Process 0.
        //
        if (myRank == 1) {
          if (debug) {
            *dbg << myRank << ": Post send-data send: tag = " << dataTag
                 << endl;
          }
          sendReqData = isend<int, scalar_type> (sendDataBuf, 0, dataTag, comm);
        }

        //
        // Process 0: Write my entries.
        //
        if (myRank == 0) {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Write my entries" << endl;
            *dbg << os.str ();
          }

          // Write Process 0's data to the output stream.
          // Matrix Market prints dense matrices in column-major order.
          const size_t printNumRows = myNumRows;
          ArrayView<const scalar_type> printData = sendDataBuf ();
          const size_t printStride = printNumRows;
          if (static_cast<size_t> (printData.size ()) < printStride * numCols) {
            lclErr = 1;
            if (! err.is_null ()) {
              std::ostringstream os;
              os << "Process " << myRank << ": My MultiVector data's size "
                 << printData.size () << " does not match my local dimensions "
                 << printStride << " x " << numCols << "." << endl;
              *err << os.str ();
            }
          }
          else {
            // Matrix Market dense format wants one number per line.
            // It wants each complex number as two real numbers (real
            // resp. imaginary parts) with a space between.
            for (size_t col = 0; col < numCols; ++col) {
              for (size_t row = 0; row < printNumRows; ++row) {
                if (STS::isComplex) {
                  out << STS::real (printData[row + col * printStride]) << " "
                      << STS::imag (printData[row + col * printStride]) << endl;
                } else {
                  out << printData[row + col * printStride] << endl;
                }
              }
            }
          }
        }

        if (myRank == 0) {
          // Wait on receive-size receive from Process 1.
          const int recvRank = 1;
          const int circBufInd = recvRank % 3;
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Wait on receive-size receive from Process "
               << recvRank << endl;
            *dbg << os.str ();
          }
          if (numProcs > 1) {
            wait<int> (comm, outArg (recvSizeReqs[circBufInd]));

            // We received the number of rows of data.  (The data
            // come in two columns.)
            size_t recvNumRows = (recvSizeBufs[circBufInd])[0];
            if (recvNumRows == Teuchos::OrdinalTraits<size_t>::invalid ()) {
              lclErr = 1;
              if (! err.is_null ()) {
                std::ostringstream os;
                os << myRank << ": Result of receive-size receive from Process "
                   << recvRank << " is Teuchos::OrdinalTraits<size_t>::invalid() "
                   << "= " << Teuchos::OrdinalTraits<size_t>::invalid () << ".  "
                  "This should never happen, and suggests that the receive never "
                  "got posted.  Please report this bug to the Tpetra developers."
                   << endl;
                *err << os.str ();
              }

              // If we're going to continue after error, set the
              // number of rows to receive to a reasonable size.  This
              // may cause MPI_ERR_TRUNCATE if the sending process is
              // sending more than 0 rows, but that's better than MPI
              // overflowing due to the huge positive value that is
              // Teuchos::OrdinalTraits<size_t>::invalid().
              recvNumRows = 0;
            }

            // Post receive-data receive from Process 1.
            recvDataBufs[circBufInd].resize (recvNumRows * numCols);
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post receive-data receive from Process "
                 << recvRank << ": tag = " << dataTag << ", buffer size = "
                 << recvDataBufs[circBufInd].size () << endl;
              *dbg << os.str ();
            }
            if (! recvSizeReqs[circBufInd].is_null ()) {
              lclErr = 1;
              if (! err.is_null ()) {
                std::ostringstream os;
                os << myRank << ": recvSizeReqs[" << circBufInd << "] is not "
                  "null, before posting the receive-data receive from Process "
                   << recvRank << ".  This should never happen.  Please report "
                  "this bug to the Tpetra developers." << endl;
                *err << os.str ();
              }
            }
            recvDataReqs[circBufInd] =
              ireceive<int, scalar_type> (recvDataBufs[circBufInd],
                                          recvRank, dataTag, comm);
          } // numProcs > 1
        }
        else if (myRank == 1) {
          // Wait on my send-size send.
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Wait on my send-size send" << endl;
            *dbg << os.str ();
          }
          wait<int> (comm, outArg (sendReqSize));
        }

        //
        // Pipeline loop
        //
        for (int p = 1; p < numProcs; ++p) {
          if (myRank == 0) {
            if (p + 2 < numProcs) {
              // Post receive-size receive from Process p + 2.
              const int recvRank = p + 2;
              const int circBufInd = recvRank % 3;
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Post receive-size receive from Process "
                   << recvRank << ": tag = " << sizeTag << endl;
                *dbg << os.str ();
              }
              if (! recvSizeReqs[circBufInd].is_null ()) {
                lclErr = 1;
                if (! err.is_null ()) {
                  std::ostringstream os;
                  os << myRank << ": recvSizeReqs[" << circBufInd << "] is not "
                     << "null, for the receive-size receive from Process "
                     << recvRank << "!  This may mean that this process never "
                     << "finished waiting for the receive from Process "
                     << (recvRank - 3) << "." << endl;
                  *err << os.str ();
                }
              }
              recvSizeReqs[circBufInd] =
                ireceive<int, size_t> (recvSizeBufs[circBufInd],
                                       recvRank, sizeTag, comm);
            }

            if (p + 1 < numProcs) {
              const int recvRank = p + 1;
              const int circBufInd = recvRank % 3;

              // Wait on receive-size receive from Process p + 1.
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Wait on receive-size receive from Process "
                   << recvRank << endl;
                *dbg << os.str ();
              }
              wait<int> (comm, outArg (recvSizeReqs[circBufInd]));

              // We received the number of rows of data.  (The data
              // come in two columns.)
              size_t recvNumRows = (recvSizeBufs[circBufInd])[0];
              if (recvNumRows == Teuchos::OrdinalTraits<size_t>::invalid ()) {
                lclErr = 1;
                if (! err.is_null ()) {
                  std::ostringstream os;
                  os << myRank << ": Result of receive-size receive from Process "
                     << recvRank << " is Teuchos::OrdinalTraits<size_t>::invalid() "
                     << "= " << Teuchos::OrdinalTraits<size_t>::invalid () << ".  "
                    "This should never happen, and suggests that the receive never "
                    "got posted.  Please report this bug to the Tpetra developers."
                     << endl;
                  *err << os.str ();
                }
                // If we're going to continue after error, set the
                // number of rows to receive to a reasonable size.
                // This may cause MPI_ERR_TRUNCATE if the sending
                // process sends more than 0 rows, but that's better
                // than MPI overflowing due to the huge positive value
                // Teuchos::OrdinalTraits<size_t>::invalid().
                recvNumRows = 0;
              }

              // Post receive-data receive from Process p + 1.
              recvDataBufs[circBufInd].resize (recvNumRows * numCols);
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Post receive-data receive from Process "
                   << recvRank << ": tag = " << dataTag << ", buffer size = "
                   << recvDataBufs[circBufInd].size () << endl;
                *dbg << os.str ();
              }
              if (! recvDataReqs[circBufInd].is_null ()) {
                lclErr = 1;
                if (! err.is_null ()) {
                  std::ostringstream os;
                  os << myRank << ": recvDataReqs[" << circBufInd << "] is not "
                     << "null, for the receive-data receive from Process "
                     << recvRank << "!  This may mean that this process never "
                     << "finished waiting for the receive from Process "
                     << (recvRank - 3) << "." << endl;
                  *err << os.str ();
                }
              }
              recvDataReqs[circBufInd] =
                ireceive<int, scalar_type> (recvDataBufs[circBufInd],
                                            recvRank, dataTag, comm);
            }

            // Wait on receive-data receive from Process p.
            const int recvRank = p;
            const int circBufInd = recvRank % 3;
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on receive-data receive from Process "
                 << recvRank << endl;
              *dbg << os.str ();
            }
            wait<int> (comm, outArg (recvDataReqs[circBufInd]));

            // Write Process p's data.  Number of rows lives in
            // recvSizeBufs[circBufInd], and the actual data live in
            // recvDataBufs[circBufInd].  Do this after posting receives,
            // in order to expose overlap of comm. with file I/O.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Write entries from Process " << recvRank
                 << endl;
              *dbg << os.str () << endl;
            }
            size_t printNumRows = (recvSizeBufs[circBufInd])[0];
            if (printNumRows == Teuchos::OrdinalTraits<size_t>::invalid ()) {
              lclErr = 1;
              if (! err.is_null ()) {
                std::ostringstream os;
                os << myRank << ": Result of receive-size receive from Process "
                   << recvRank << " was Teuchos::OrdinalTraits<size_t>::"
                  "invalid() = " << Teuchos::OrdinalTraits<size_t>::invalid ()
                   << ".  This should never happen, and suggests that its "
                  "receive-size receive was never posted.  "
                  "Please report this bug to the Tpetra developers." << endl;
                *err << os.str ();
              }
              // If we're going to continue after error, set the
              // number of rows to print to a reasonable size.
              printNumRows = 0;
            }
            if (printNumRows > 0 && recvDataBufs[circBufInd].is_null ()) {
              lclErr = 1;
              if (! err.is_null ()) {
                std::ostringstream os;
                os << myRank << ": Result of receive-size receive from Proc "
                   << recvRank << " was " << printNumRows << " > 0, but "
                  "recvDataBufs[" << circBufInd << "] is null.  This should "
                  "never happen.  Please report this bug to the Tpetra "
                  "developers." << endl;
                *err << os.str ();
              }
              // If we're going to continue after error, set the
              // number of rows to print to a reasonable size.
              printNumRows = 0;
            }

            // Write the received data to the output stream.
            // Matrix Market prints dense matrices in column-major order.
            ArrayView<const scalar_type> printData = (recvDataBufs[circBufInd]) ();
            const size_t printStride = printNumRows;
            // Matrix Market dense format wants one number per line.
            // It wants each complex number as two real numbers (real
            // resp. imaginary parts) with a space between.
            for (size_t col = 0; col < numCols; ++col) {
              for (size_t row = 0; row < printNumRows; ++row) {
                if (STS::isComplex) {
                  out << STS::real (printData[row + col * printStride]) << " "
                      << STS::imag (printData[row + col * printStride]) << endl;
                } else {
                  out << printData[row + col * printStride] << endl;
                }
              }
            }
          }
          else if (myRank == p) { // Process p
            // Wait on my send-data send.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on my send-data send" << endl;
              *dbg << os.str ();
            }
            wait<int> (comm, outArg (sendReqData));
          }
          else if (myRank == p + 1) { // Process p + 1
            // Post send-data send to Process 0.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post send-data send: tag = " << dataTag
                 << endl;
              *dbg << os.str ();
            }
            sendReqData = isend<int, scalar_type> (sendDataBuf, 0, dataTag, comm);
            // Wait on my send-size send.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on my send-size send" << endl;
              *dbg << os.str ();
            }
            wait<int> (comm, outArg (sendReqSize));
          }
          else if (myRank == p + 2) { // Process p + 2
            // Post send-size send to Process 0.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post send-size send: size = "
                 << sendDataSize[0] << ", tag = " << sizeTag << endl;
              *dbg << os.str ();
            }
            sendReqSize = isend<int, size_t> (sendDataSize, 0, sizeTag, comm);
          }
        }

        // Establish global agreement on the error state.
        reduceAll<int, int> (comm, REDUCE_MAX, lclErr, outArg (gblErr));
        TEUCHOS_TEST_FOR_EXCEPTION(
          gblErr == 1, std::runtime_error, "Tpetra::MatrixMarket::writeDense "
          "experienced some kind of error and was unable to complete.");

        if (debug) {
          dbg->popTab ();
          *dbg << myRank << ": writeDenseColumn: Done" << endl;
          dbg->popTab ();
        }
      }

    public:

      /// \brief Print the multivector in Matrix Market format, with
      ///   matrix name and or description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDense().
      static void
      writeDense (std::ostream& out,
                  const Teuchos::RCP<const multivector_type>& X,
                  const std::string& matrixName,
                  const std::string& matrixDescription,
                  const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                  const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          X.is_null (), std::invalid_argument, "Tpetra::MatrixMarket::"
          "writeDense: The input MultiVector X is null.");
        writeDense (out, *X, matrixName, matrixDescription, err, dbg);
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   no matrix name or description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDense().
      static void
      writeDense (std::ostream& out,
                  const multivector_type& X,
                  const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                  const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        writeDense (out, X, "", "", err, dbg);
      }

      /// \brief Print the multivector in Matrix Market format, with
      ///   no matrix name or description.
      ///
      /// See the documentation of the above six-argument version of
      /// writeDense().
      static void
      writeDense (std::ostream& out,
                  const Teuchos::RCP<const multivector_type>& X,
                  const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                  const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          X.is_null (), std::invalid_argument, "Tpetra::MatrixMarket::"
          "writeDense: The input MultiVector X is null.");
        writeDense (out, *X, "", "", err, dbg);
      }

      /// \brief Print the Map to the given output stream.
      ///
      /// \param out [out] Output stream to which to print.  This only
      ///   needs to be accessible on Process 0 in the Map's
      ///   communicator; no other process will do anything with it.
      ///
      /// \param map [in] The Map to print.
      ///
      /// \param debug [in] Whether to print copious debugging output
      ///   to stderr on <i>all</i> processes in the Map's
      ///   communicator.
      ///
      /// We print the Map in Matrix Market format as a dense
      /// nonsymmetric integer matrix with two columns.  The first
      /// column holds global indices (GIDs), and the second column
      /// holds process ranks (PIDs).  In any row of the matrix, the
      /// first entry is a GID, and the second is a PID that owns the
      /// GID.  Multiple PIDs may own the same GID, and the order of
      /// rows with respect to a given PID is significant.
      static void
      writeMap (std::ostream& out, const map_type& map, const bool debug=false)
      {
        Teuchos::RCP<Teuchos::FancyOStream> err =
          Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
        writeMap (out, map, err, debug);
      }

      /// \brief Print the Map to the given output stream \c out.
      ///
      /// This version of writeMap() comes with an extra debug output
      /// stream \c err, that is only used if \c debug is true.
      ///
      /// \warning We make no promises of backwards compatibility with
      ///   this method.  It may go away or its interface may change
      ///   at any time.
      static void
      writeMap (std::ostream& out,
                const map_type& map,
                const Teuchos::RCP<Teuchos::FancyOStream>& err,
                const bool debug=false)
      {
        using Teuchos::Array;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::Comm;
        using Teuchos::CommRequest;
        using Teuchos::ireceive;
        using Teuchos::isend;
        using Teuchos::RCP;
        using Teuchos::TypeNameTraits;
        using Teuchos::wait;
        using std::endl;
        typedef global_ordinal_type GO;
        typedef int pid_type;

        // Treat the Map as a 1-column "multivector."  This differs
        // from the previous two-column format, in which column 0 held
        // the GIDs, and column 1 held the corresponding PIDs.  It
        // differs because printing that format requires knowing the
        // entire first column -- that is, all the GIDs -- in advance.
        // Sending messages from each process one at a time saves
        // memory, but it means that Process 0 doesn't ever have all
        // the GIDs at once.
        //
        // We pack the entries as ptrdiff_t, since this should be the
        // biggest signed built-in integer type that can hold any GO
        // or pid_type (= int) quantity without overflow.  Test this
        // assumption at run time.
        typedef ptrdiff_t int_type;
        TEUCHOS_TEST_FOR_EXCEPTION(
          sizeof (GO) > sizeof (int_type), std::logic_error,
          "The global ordinal type GO=" << TypeNameTraits<GO>::name ()
          << " is too big for ptrdiff_t.  sizeof(GO) = " << sizeof (GO)
          << " > sizeof(ptrdiff_t) = " << sizeof (ptrdiff_t) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(
          sizeof (pid_type) > sizeof (int_type), std::logic_error,
          "The (MPI) process rank type pid_type=" <<
          TypeNameTraits<pid_type>::name () << " is too big for ptrdiff_t.  "
          "sizeof(pid_type) = " << sizeof (pid_type) << " > sizeof(ptrdiff_t)"
          " = " << sizeof (ptrdiff_t) << ".");

        const Comm<int>& comm = * (map.getComm ());
        const int myRank = comm.getRank ();
        const int numProcs = comm.getSize ();

        if (! err.is_null ()) {
          err->pushTab ();
        }
        if (debug) {
          std::ostringstream os;
          os << myRank << ": writeMap" << endl;
          *err << os.str ();
        }
        if (! err.is_null ()) {
          err->pushTab ();
        }

        const size_t myNumRows = map.getLocalNumElements ();
        // Use a different tag for the "size" messages than for the
        // "data" messages, in order to help us debug any mix-ups.
        const int sizeTag = 1337;
        const int dataTag = 1338;

        // Process 0 pipelines nonblocking receives with file output.
        //
        // Constraints:
        //   - Process 0 can't post a receive for another process'
        //     actual data, until it posts and waits on the receive
        //     from that process with the amount of data to receive.
        //     (We could just post receives with a max data size, but
        //     I feel uncomfortable about that.)
        //   - The C++ standard library doesn't allow nonblocking
        //     output to an std::ostream.
        //
        // Process 0: Post receive-size receives from Processes 1 and 2.
        // Process 1: Post send-size send to Process 0.
        // Process 2: Post send-size send to Process 0.
        //
        // All processes: Pack my GIDs and PIDs.
        //
        // Process 1:
        //   - Post send-data send to Process 0.
        //   - Wait on my send-size send to Process 0.
        //
        // Process 0:
        //   - Print MatrixMarket header.
        //   - Print my GIDs and PIDs.
        //   - Wait on receive-size receive from Process 1.
        //   - Post receive-data receive from Process 1.
        //
        // For each process p = 1, 2, ... numProcs-1:
        //   If I am Process 0:
        //     - Post receive-size receive from Process p + 2
        //     - Wait on receive-size receive from Process p + 1
        //     - Post receive-data receive from Process p + 1
        //     - Wait on receive-data receive from Process p
        //     - Write data from Process p.
        //   Else if I am Process p:
        //     - Wait on my send-data send.
        //   Else if I am Process p+1:
        //     - Post send-data send to Process 0.
        //     - Wait on my send-size send.
        //   Else if I am Process p+2:
        //     - Post send-size send to Process 0.
        //
        // Pipelining has three goals here:
        //   1. Overlap communication (the receives) with file I/O
        //   2. Give Process 0 a chance to prepost some receives,
        //      before sends show up, by packing local data before
        //      posting sends
        //   3. Don't post _all_ receives or _all_ sends, because that
        //      wouldn't be memory scalable.  (Just because we can't
        //      see how much memory MPI consumes, doesn't mean that it
        //      doesn't consume any!)

        // These are used on every process.  sendReqSize[0] holds the
        // number of rows on this process, and sendReqBuf holds this
        // process' data.  Process 0 packs into sendReqBuf, but
        // doesn't send; it only uses that for printing.  All other
        // processes send both of these to Process 0.
        RCP<CommRequest<int> > sendReqSize, sendReqData;

        // These are used only on Process 0, for received data.  Keep
        // 3 of each, and treat the arrays as circular buffers.  When
        // receiving from Process p, the corresponding array index
        // here is p % 3.
        Array<ArrayRCP<int_type> > recvSizeBufs (3);
        Array<ArrayRCP<int_type> > recvDataBufs (3);
        Array<RCP<CommRequest<int> > > recvSizeReqs (3);
        Array<RCP<CommRequest<int> > > recvDataReqs (3);

        // Buffer for nonblocking send of the "send size."
        ArrayRCP<int_type> sendDataSize (1);
        sendDataSize[0] = myNumRows;

        if (myRank == 0) {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Post receive-size receives from "
              "Procs 1 and 2: tag = " << sizeTag << endl;
            *err << os.str ();
          }
          // Process 0: Post receive-size receives from Processes 1 and 2.
          recvSizeBufs[0].resize (1);
          (recvSizeBufs[0])[0] = -1; // error flag
          recvSizeBufs[1].resize (1);
          (recvSizeBufs[1])[0] = -1; // error flag
          recvSizeBufs[2].resize (1);
          (recvSizeBufs[2])[0] = -1; // error flag
          if (numProcs > 1) {
            recvSizeReqs[1] =
              ireceive<int, int_type> (recvSizeBufs[1], 1, sizeTag, comm);
          }
          if (numProcs > 2) {
            recvSizeReqs[2] =
              ireceive<int, int_type> (recvSizeBufs[2], 2, sizeTag, comm);
          }
        }
        else if (myRank == 1 || myRank == 2) {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Post send-size send: size = "
               << sendDataSize[0] << ", tag = " << sizeTag << endl;
            *err << os.str ();
          }
          // Prime the pipeline by having Processes 1 and 2 start
          // their send-size sends.  We don't want _all_ the processes
          // to start their send-size sends, because that wouldn't be
          // memory scalable.
          sendReqSize = isend<int, int_type> (sendDataSize, 0, sizeTag, comm);
        }
        else {
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Not posting my send-size send yet" << endl;
            *err << os.str ();
          }
        }

        //
        // Pack my GIDs and PIDs.  Each (GID,PID) pair gets packed
        // consecutively, for better locality.
        //

        if (debug) {
          std::ostringstream os;
          os << myRank << ": Pack my GIDs and PIDs" << endl;
          *err << os.str ();
        }

        ArrayRCP<int_type> sendDataBuf (myNumRows * 2);

        if (map.isContiguous ()) {
          const int_type myMinGblIdx =
            static_cast<int_type> (map.getMinGlobalIndex ());
          for (size_t k = 0; k < myNumRows; ++k) {
            const int_type gid = myMinGblIdx + static_cast<int_type> (k);
            const int_type pid = static_cast<int_type> (myRank);
            sendDataBuf[2*k] = gid;
            sendDataBuf[2*k+1] = pid;
          }
        }
        else {
          ArrayView<const GO> myGblInds = map.getLocalElementList ();
          for (size_t k = 0; k < myNumRows; ++k) {
            const int_type gid = static_cast<int_type> (myGblInds[k]);
            const int_type pid = static_cast<int_type> (myRank);
            sendDataBuf[2*k] = gid;
            sendDataBuf[2*k+1] = pid;
          }
        }

        if (debug) {
          std::ostringstream os;
          os << myRank << ": Done packing my GIDs and PIDs" << endl;
          *err << os.str ();
        }

        if (myRank == 1) {
          // Process 1: post send-data send to Process 0.
          if (debug) {
            *err << myRank << ": Post send-data send: tag = " << dataTag
                 << endl;
          }
          sendReqData = isend<int, int_type> (sendDataBuf, 0, dataTag, comm);
        }

        if (myRank == 0) {
          if (debug) {
            *err << myRank << ": Write MatrixMarket header" << endl;
          }

          // Process 0: Write the MatrixMarket header.
          // Description section explains each column.
          std::ostringstream hdr;

          // Print the Matrix Market header.  MultiVector stores data
          // nonsymmetrically, hence "general" in the banner line.
          hdr << "%%MatrixMarket matrix array integer general" << endl
              << "% Format: Version 2.0" << endl
              << "%" << endl
              << "% This file encodes a Tpetra::Map." << endl
              << "% It is stored as a dense vector, with twice as many " << endl
              << "% entries as the global number of GIDs (global indices)." << endl
              << "% (GID, PID) pairs are stored contiguously, where the PID " << endl
              << "% is the rank of the process owning that GID." << endl
              << (2 * map.getGlobalNumElements ()) << " " << 1 << endl;
          out << hdr.str ();

          if (debug) {
            std::ostringstream os;
            os << myRank << ": Write my GIDs and PIDs" << endl;
            *err << os.str ();
          }

          // Write Process 0's data to the output stream.
          // Matrix Market prints dense matrices in column-major order.
          const int_type printNumRows = myNumRows;
          ArrayView<const int_type> printData = sendDataBuf ();
          for (int_type k = 0; k < printNumRows; ++k) {
            const int_type gid = printData[2*k];
            const int_type pid = printData[2*k+1];
            out << gid << endl << pid << endl;
          }
        }

        if (myRank == 0) {
          // Wait on receive-size receive from Process 1.
          const int recvRank = 1;
          const int circBufInd = recvRank % 3;
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Wait on receive-size receive from Process "
               << recvRank << endl;
            *err << os.str ();
          }
          if (numProcs > 1) {
            wait<int> (comm, outArg (recvSizeReqs[circBufInd]));

            // We received the number of rows of data.  (The data
            // come in two columns.)
            const int_type recvNumRows = (recvSizeBufs[circBufInd])[0];
            if (debug && recvNumRows == -1) {
              std::ostringstream os;
              os << myRank << ": Result of receive-size receive from Process "
                 << recvRank << " is -1.  This should never happen, and "
                "suggests that the receive never got posted.  Please report "
                "this bug to the Tpetra developers." << endl;
              *err << os.str ();
            }

            // Post receive-data receive from Process 1.
            recvDataBufs[circBufInd].resize (recvNumRows * 2);
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post receive-data receive from Process "
                 << recvRank << ": tag = " << dataTag << ", buffer size = "
                 << recvDataBufs[circBufInd].size () << endl;
              *err << os.str ();
            }
            if (! recvSizeReqs[circBufInd].is_null ()) {
              std::ostringstream os;
              os << myRank << ": recvSizeReqs[" << circBufInd << "] is not "
                "null, before posting the receive-data receive from Process "
                 << recvRank << ".  This should never happen.  Please report "
                "this bug to the Tpetra developers." << endl;
              *err << os.str ();
            }
            recvDataReqs[circBufInd] =
              ireceive<int, int_type> (recvDataBufs[circBufInd],
                                       recvRank, dataTag, comm);
          } // numProcs > 1
        }
        else if (myRank == 1) {
          // Wait on my send-size send.
          if (debug) {
            std::ostringstream os;
            os << myRank << ": Wait on my send-size send" << endl;
            *err << os.str ();
          }
          wait<int> (comm, outArg (sendReqSize));
        }

        //
        // Pipeline loop
        //
        for (int p = 1; p < numProcs; ++p) {
          if (myRank == 0) {
            if (p + 2 < numProcs) {
              // Post receive-size receive from Process p + 2.
              const int recvRank = p + 2;
              const int circBufInd = recvRank % 3;
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Post receive-size receive from Process "
                   << recvRank << ": tag = " << sizeTag << endl;
                *err << os.str ();
              }
              if (! recvSizeReqs[circBufInd].is_null ()) {
                std::ostringstream os;
                os << myRank << ": recvSizeReqs[" << circBufInd << "] is not "
                   << "null, for the receive-size receive from Process "
                   << recvRank << "!  This may mean that this process never "
                   << "finished waiting for the receive from Process "
                   << (recvRank - 3) << "." << endl;
                *err << os.str ();
              }
              recvSizeReqs[circBufInd] =
                ireceive<int, int_type> (recvSizeBufs[circBufInd],
                                         recvRank, sizeTag, comm);
            }

            if (p + 1 < numProcs) {
              const int recvRank = p + 1;
              const int circBufInd = recvRank % 3;

              // Wait on receive-size receive from Process p + 1.
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Wait on receive-size receive from Process "
                   << recvRank << endl;
                *err << os.str ();
              }
              wait<int> (comm, outArg (recvSizeReqs[circBufInd]));

              // We received the number of rows of data.  (The data
              // come in two columns.)
              const int_type recvNumRows = (recvSizeBufs[circBufInd])[0];
              if (debug && recvNumRows == -1) {
                std::ostringstream os;
                os << myRank << ": Result of receive-size receive from Process "
                   << recvRank << " is -1.  This should never happen, and "
                  "suggests that the receive never got posted.  Please report "
                  "this bug to the Tpetra developers." << endl;
                *err << os.str ();
              }

              // Post receive-data receive from Process p + 1.
              recvDataBufs[circBufInd].resize (recvNumRows * 2);
              if (debug) {
                std::ostringstream os;
                os << myRank << ": Post receive-data receive from Process "
                   << recvRank << ": tag = " << dataTag << ", buffer size = "
                   << recvDataBufs[circBufInd].size () << endl;
                *err << os.str ();
              }
              if (! recvDataReqs[circBufInd].is_null ()) {
                std::ostringstream os;
                os << myRank << ": recvDataReqs[" << circBufInd << "] is not "
                   << "null, for the receive-data receive from Process "
                   << recvRank << "!  This may mean that this process never "
                   << "finished waiting for the receive from Process "
                   << (recvRank - 3) << "." << endl;
                *err << os.str ();
              }
              recvDataReqs[circBufInd] =
                ireceive<int, int_type> (recvDataBufs[circBufInd],
                                         recvRank, dataTag, comm);
            }

            // Wait on receive-data receive from Process p.
            const int recvRank = p;
            const int circBufInd = recvRank % 3;
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on receive-data receive from Process "
                 << recvRank << endl;
              *err << os.str ();
            }
            wait<int> (comm, outArg (recvDataReqs[circBufInd]));

            // Write Process p's data.  Number of rows lives in
            // recvSizeBufs[circBufInd], and the actual data live in
            // recvDataBufs[circBufInd].  Do this after posting receives,
            // in order to expose overlap of comm. with file I/O.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Write GIDs and PIDs from Process "
                 << recvRank << endl;
              *err << os.str () << endl;
            }
            const int_type printNumRows = (recvSizeBufs[circBufInd])[0];
            if (debug && printNumRows == -1) {
              std::ostringstream os;
              os << myRank << ": Result of receive-size receive from Process "
                 << recvRank << " was -1.  This should never happen, and "
                "suggests that its receive-size receive was never posted.  "
                "Please report this bug to the Tpetra developers." << endl;
              *err << os.str ();
            }
            if (debug && printNumRows > 0 && recvDataBufs[circBufInd].is_null ()) {
              std::ostringstream os;
              os << myRank << ": Result of receive-size receive from Proc "
                 << recvRank << " was " << printNumRows << " > 0, but "
                "recvDataBufs[" << circBufInd << "] is null.  This should "
                "never happen.  Please report this bug to the Tpetra "
                "developers." << endl;
              *err << os.str ();
            }
            ArrayView<const int_type> printData = (recvDataBufs[circBufInd]) ();
            for (int_type k = 0; k < printNumRows; ++k) {
              const int_type gid = printData[2*k];
              const int_type pid = printData[2*k+1];
              out << gid << endl << pid << endl;
            }
          }
          else if (myRank == p) { // Process p
            // Wait on my send-data send.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on my send-data send" << endl;
              *err << os.str ();
            }
            wait<int> (comm, outArg (sendReqData));
          }
          else if (myRank == p + 1) { // Process p + 1
            // Post send-data send to Process 0.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post send-data send: tag = " << dataTag
                 << endl;
              *err << os.str ();
            }
            sendReqData = isend<int, int_type> (sendDataBuf, 0, dataTag, comm);
            // Wait on my send-size send.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Wait on my send-size send" << endl;
              *err << os.str ();
            }
            wait<int> (comm, outArg (sendReqSize));
          }
          else if (myRank == p + 2) { // Process p + 2
            // Post send-size send to Process 0.
            if (debug) {
              std::ostringstream os;
              os << myRank << ": Post send-size send: size = "
                 << sendDataSize[0] << ", tag = " << sizeTag << endl;
              *err << os.str ();
            }
            sendReqSize = isend<int, int_type> (sendDataSize, 0, sizeTag, comm);
          }
        }

        if (! err.is_null ()) {
          err->popTab ();
        }
        if (debug) {
          *err << myRank << ": writeMap: Done" << endl;
        }
        if (! err.is_null ()) {
          err->popTab ();
        }
      }

      //! Write the Map to the given file.
      static void
      writeMapFile (const std::string& filename,
                    const map_type& map)
      {
        const int myRank = map.getComm ()->getRank ();

        auto out = Writer::openOutFileOnRankZero(map.getComm(), filename, myRank, true);

        writeMap (out, map);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeDense() throws an
        // exception.
      }

    private:
      /// \brief Print the given possibly multiline string as a comment.
      ///
      /// If the string is empty, don't print anything (not even an
      /// empty line).  Otherwise, print each line of the string (one
      /// or more lines) as a comment in the comments section of a
      /// Matrix Market file (the part of the file after the initial
      /// banner line, but before the matrix's size metadata and
      /// data).
      ///
      /// \param out [out] The output string to which to print.  This
      ///   function is <i>not</i> a collective operation; whichever
      ///   MPI process calls it will print to the given output
      ///   stream.
      ///
      /// \param str [in] The string to print.  It consists of zero or
      ///   more lines.  If empty, nothing is printed, not even an
      ///   empty line.
      ///
      /// \note Printing comments is tricky only because the string
      ///   might contain newlines.  We have to ensure that all the
      ///   lines start with a comment character.  If they already do,
      ///   we print each line as is; otherwise, we append a comment
      ///   character and a space to each line.
      static void
      printAsComment (std::ostream& out, const std::string& str)
      {
        using std::endl;
        std::istringstream inpstream (str);
        std::string line;

        while (getline (inpstream, line)) {
          if (! line.empty()) {
            // Note that getline() doesn't store '\n', so we have to
            // append the endline ourselves.
            if (line[0] == '%') { // Line starts with a comment character.
              out << line << endl;
            }
            else { // Line doesn't start with a comment character.
              out << "%% " << line << endl;
            }
          }
        }
      }

    public:

      /// \brief Write a Tpetra::Operator to a file.
      ///
      /// This method works by applying the Operator to columns of the
      /// identity matrix.  As a result, it effectively turns the
      /// Operator into a dense matrix.  However, it writes the
      /// Operator in sparse matrix format.  As such, you may read it
      /// back in again using Reader::readSparseFile.
      ///
      /// Probing calls apply() on the input Operator, using a
      /// MultiVector with a small, fixed number of columns.  If you
      /// want to change the number of columns used, you must invoke
      /// the overload of this method that takes an input
      /// Teuchos::ParameterList (see below).
      ///
      /// \param fileName [in] The name of the file to which to write.
      ///   Only Process 0 in the input Operator's communicator will
      ///   write to the file.
      /// \param A [in] The input Tpetra::Operator to write.
      static void
      writeOperator(const std::string& fileName, operator_type const &A) {
        Teuchos::ParameterList pl;
        writeOperator(fileName, A, pl);
      }

      /// \brief Write a Tpetra::Operator to an output stream.
      ///
      /// This method works by applying the Operator to columns of the
      /// identity matrix.  As a result, it effectively turns the
      /// Operator into a dense matrix.  However, it writes the
      /// Operator in sparse matrix format.  As such, you may read it
      /// back in again using Reader::readSparseFile.
      ///
      /// Probing calls apply() on the input Operator, using a
      /// MultiVector with a small, fixed number of columns.  If you
      /// want to change the number of columns used, you must invoke
      /// the overload of this method that takes an input
      /// Teuchos::ParameterList (see below).
      ///
      /// \param out [in] Output stream to which to write.  Only
      ///   Process 0 in the input Operator's communicator will write
      ///   to the output stream.  Other processes will not write to
      ///   it or call any methods on it.  Thus, the stream need only
      ///   be valid on Process 0.
      /// \param A [in] The input Tpetra::Operator to write.
      static void
      writeOperator (std::ostream& out, const operator_type& A) {
        Teuchos::ParameterList pl;
        writeOperator (out, A, pl);
      }

      /// \brief Write a Tpetra::Operator to a file, with options.
      ///
      /// This method works by applying the Operator to columns of the
      /// identity matrix.  As a result, it effectively turns the
      /// Operator into a dense matrix.  However, it writes the
      /// Operator in sparse matrix format.  As such, you may read it
      /// back in again using Reader::readSparseFile.
      ///
      /// Probing calls apply() on the input Operator, using a
      /// MultiVector with a small, fixed number of columns.  You may
      /// set this number of columns in the input ParameterList.
      ///
      /// \param fileName [in] The name of the file to which to write.
      ///   Only Process 0 in the Operator's communicator will write
      ///   to the file.
      /// \param A [in] The input Tpetra::Operator to write.
      /// \param params [in] List of options.  An empty list means
      ///   "use default values of options."
      ///
      /// If you always want the default options, use the overload of
      /// this method above that takes two arguments (the filename and
      /// the Operator).  This three-argument overload lets the user
      /// pass in options.  The currently supported options are:
      ///
      /// <ul>
      /// <li> "probing size" (integer [10]): number of columns to use
      ///       in the probing MultiVector </li>
      /// <li> "precision" (integer [C++ default]): precision to use
      ///      when writing floating-point values </li>
      /// <li> "print MatrixMarket header" (boolean [true]): whether
      ///      to print the MatrixMarket header </li>
      /// <li> "zero-based indexing" (boolean [false]): print matrix
      ///      using zero-based indexing.  The Matrix Market format
      ///      uses one-based indexing, so setting this option to true
      ///      violates the Matrix Market format standard. </li>
      /// </ul>
      static void
      writeOperator (const std::string& fileName,
                     const operator_type& A,
                     const Teuchos::ParameterList& params)
      {
        std::ofstream out;
        std::string tmpFile = "__TMP__" + fileName;
        const int myRank = A.getDomainMap()->getComm()->getRank();
        bool precisionChanged=false;
        int  oldPrecision;
        // The number of nonzero entries in a Tpetra::Operator is
        // unknown until probing is completed.  In order to write a
        // MatrixMarket header, we write the matrix to a temporary
        // file.
        //
        // FIXME (mfh 23 May 2015) IT WASN'T MY IDEA TO WRITE TO A
        // TEMPORARY FILE.
        if (myRank==0) {
          if (std::ifstream(tmpFile))
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                                       "writeOperator: temporary file " << tmpFile << " already exists");
          out.open(tmpFile.c_str());
          if (params.isParameter("precision")) {
            oldPrecision = out.precision(params.get<int>("precision"));
            precisionChanged=true;
          }
        }

        const std::string header = writeOperatorImpl(out, A, params);

        if (myRank==0) {
          if (precisionChanged)
            out.precision(oldPrecision);
          out.close();
          out.open(fileName.c_str(), std::ios::binary);
          bool printMatrixMarketHeader = true;
          if (params.isParameter("print MatrixMarket header"))
            printMatrixMarketHeader = params.get<bool>("print MatrixMarket header");
          if (printMatrixMarketHeader && myRank == 0) {
            // Write header to final file.
            out << header;
          }
          // Append matrix from temporary to final file.
          std::ifstream src(tmpFile, std::ios_base::binary);
          out << src.rdbuf();
          src.close();
          // Delete the temporary file.
          remove(tmpFile.c_str());
        }
      }

      /// \brief Write a Tpetra::Operator to an output stream, with options.
      ///
      /// This method works by applying the Operator to columns of the
      /// identity matrix.  As a result, it effectively turns the
      /// Operator into a dense matrix.  However, it writes the
      /// Operator in sparse matrix format.  As such, you may read it
      /// back in again using Reader::readSparseFile.
      ///
      /// Probing calls apply() on the input Operator, using a
      /// MultiVector with a small, fixed number of columns.  You may
      /// set this number of columns in the input ParameterList.
      ///
      /// \param out [in] Output stream to which to write.  Only
      ///   Process 0 in the input Operator's communicator will write
      ///   to the output stream.  Other processes will not write to
      ///   it or call any methods on it.  Thus, the stream need only
      ///   be valid on Process 0.
      /// \param A [in] The input Tpetra::Operator to write.
      /// \param params [in] List of options.  An empty list means
      ///   "use default values of options."
      ///
      /// If you always want the default options, use the overload of
      /// this method above that takes two arguments (the filename and
      /// the Operator).  This three-argument overload lets the user
      /// pass in options.  The currently supported options are:
      ///
      /// <ul>
      /// <li> "probing size" (integer [10]): number of columns to use
      ///       in the probing MultiVector </li>
      /// <li> "precision" (integer [C++ default]): precision to use
      ///      when writing floating-point values </li>
      /// <li> "print MatrixMarket header" (boolean [true]): whether
      ///      to print the MatrixMarket header </li>
      /// <li> "zero-based indexing" (boolean [false]): print matrix
      ///      using zero-based indexing.  The Matrix Market format
      ///      uses one-based indexing, so setting this option to true
      ///      violates the Matrix Market format standard. </li>
      /// </ul>
      static void
      writeOperator (std::ostream& out,
                     const operator_type& A,
                     const Teuchos::ParameterList& params)
      {
        const int myRank = A.getDomainMap ()->getComm ()->getRank ();

        // The number of nonzero entries in a Tpetra::Operator is
        // unknown until probing is completed.  In order to write a
        // MatrixMarket header, we write the matrix to a temporary
        // output stream.
        //
        // NOTE (mfh 23 May 2015): Writing to a temporary output
        // stream may double the memory usage, depending on whether
        // 'out' is a file stream or an in-memory output stream (e.g.,
        // std::ostringstream).  It might be wise to use a temporary
        // file instead.  However, please look carefully at POSIX
        // functions for safe creation of temporary files.  Don't just
        // prepend "__TMP__" to the filename and hope for the best.
        // Furthermore, it should be valid to call the std::ostream
        // overload of this method even when Process 0 does not have
        // access to a file system.
        std::ostringstream tmpOut;
        if (myRank == 0) {
          if (params.isParameter ("precision") && params.isType<int> ("precision")) {
            (void) tmpOut.precision (params.get<int> ("precision"));
          }
        }

        const std::string header = writeOperatorImpl (tmpOut, A, params);

        if (myRank == 0) {
          bool printMatrixMarketHeader = true;
          if (params.isParameter ("print MatrixMarket header") &&
              params.isType<bool> ("print MatrixMarket header")) {
            printMatrixMarketHeader = params.get<bool> ("print MatrixMarket header");
          }
          if (printMatrixMarketHeader && myRank == 0) {
            out << header; // write header to final output stream
          }
          // Append matrix from temporary output stream to final output stream.
          //
          // NOTE (mfh 23 May 2015) This might use a lot of memory.
          // However, we should not use temporary files in this
          // method.  Since it does not access the file system (unlike
          // the overload that takes a file name), it should not
          // require the file system at all.
          //
          // If memory usage becomes a problem, one thing we could do
          // is write the entries of the Operator one column (or a few
          // columns) at a time.  The Matrix Market sparse format does
          // not impose an order on its entries, so it would be OK to
          // write them in that order.
          out << tmpOut.str ();
        }
      }

    private:

      /// \brief Implementation detail of writeOperator overloads.
      ///
      /// Use column probing to discover the entries of the Operator.
      /// Write them in Matrix Market sparse matrix format, on Process
      /// 0 only, to the output stream \c os.  Do NOT write the Matrix
      /// Market header line, but DO return it (unless the input
      /// options specify otherwise).
      static std::string
      writeOperatorImpl (std::ostream& os,
                         const operator_type& A,
                         const Teuchos::ParameterList& params)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::ArrayRCP;
        using Teuchos::Array;

        typedef local_ordinal_type                 LO;
        typedef global_ordinal_type                GO;
        typedef scalar_type                        Scalar;
        typedef Teuchos::OrdinalTraits<LO>         TLOT;
        typedef Teuchos::OrdinalTraits<GO>         TGOT;
        typedef Tpetra::Import<LO, GO, node_type>  import_type;
        typedef Tpetra::MultiVector<GO, LO, GO, node_type> mv_type_go;

        const map_type&                domainMap = *(A.getDomainMap());
        RCP<const map_type>            rangeMap = A.getRangeMap();
        trcp_tcomm_t                   comm = rangeMap->getComm();
        const int                      myRank = comm->getRank();
        const size_t                   numProcs = comm->getSize();

        size_t numMVs = 10;
        if (params.isParameter("probing size"))
          numMVs = params.get<int>("probing size");

        GO globalNnz = 0;
        GO minColGid = domainMap.getMinAllGlobalIndex();
        GO maxColGid = domainMap.getMaxAllGlobalIndex();
        // Rather than replicating the domainMap on all processors, we instead
        // iterate from the min GID to the max GID.  If the map is gappy,
        // there will be invalid GIDs, i.e., GIDs no one has.  This will require
        // unnecessary matvecs  against potentially zero vectors.
        GO numGlobElts = maxColGid - minColGid + TGOT::one();
        GO numChunks = numGlobElts / numMVs;
        GO rem = numGlobElts % numMVs;
        GO indexBase = rangeMap->getIndexBase();

        int offsetToUseInPrinting = 1 - indexBase; // default is 1-based indexing
        if (params.isParameter("zero-based indexing")) {
          if (params.get<bool>("zero-based indexing") == true)
            offsetToUseInPrinting = -indexBase; // If 0-based, use as-is. If 1-based, subtract 1.
        }

        // Create map that replicates the range map on pid 0 and is empty for all other pids
        size_t numLocalRangeEntries = rangeMap->getLocalNumElements();

        // Create contiguous source map
        RCP<const map_type> allGidsMap = rcp(new map_type(TGOT::invalid(), numLocalRangeEntries,
                                                              indexBase, comm));
        // Create vector based on above map.  Populate it with GIDs corresponding to this pid's GIDs in rangeMap.
        mv_type_go allGids(allGidsMap,1);
        Teuchos::ArrayRCP<GO> allGidsData = allGids.getDataNonConst(0);

        for (size_t i=0; i<numLocalRangeEntries; i++)
          allGidsData[i] = rangeMap->getGlobalElement(i);
        allGidsData = Teuchos::null;

        // Create target map that is nontrivial only on pid 0
        GO numTargetMapEntries=TGOT::zero();
        Teuchos::Array<GO> importGidList;
        if (myRank==0) {
          numTargetMapEntries = rangeMap->getGlobalNumElements();
          importGidList.reserve(numTargetMapEntries);
          for (GO j=0; j<numTargetMapEntries; ++j) importGidList.push_back(j + indexBase);
        } else {
          importGidList.reserve(numTargetMapEntries);
        }
        RCP<map_type> importGidMap = rcp(new map_type(TGOT::invalid(), importGidList(), indexBase, comm));

        // Import all rangeMap GIDs to pid 0
        import_type gidImporter(allGidsMap, importGidMap);
        mv_type_go importedGids(importGidMap, 1);
        importedGids.doImport(allGids, gidImporter, INSERT);

        // The following import map will be non-trivial only on pid 0.
        ArrayRCP<const GO> importedGidsData = importedGids.getData(0);
        RCP<const map_type> importMap = rcp(new map_type(TGOT::invalid(), importedGidsData(), indexBase, comm) );

        // Importer from original range map to pid 0
        import_type importer(rangeMap, importMap);
        // Target vector on pid 0
        RCP<mv_type> colsOnPid0 = rcp(new mv_type(importMap,numMVs));

        RCP<mv_type> ei = rcp(new mv_type(A.getDomainMap(),numMVs));    //probing vector
        RCP<mv_type> colsA = rcp(new mv_type(A.getRangeMap(),numMVs));  //columns of A revealed by probing

        Array<GO> globalColsArray, localColsArray;
        globalColsArray.reserve(numMVs);
        localColsArray.reserve(numMVs);

        ArrayRCP<ArrayRCP<Scalar> > eiData(numMVs);
        for (size_t i=0; i<numMVs; ++i)
          eiData[i] = ei->getDataNonConst(i);

        // //////////////////////////////////////
        // Discover A by chunks
        // //////////////////////////////////////
        for (GO k=0; k<numChunks; ++k) {
          for (size_t j=0; j<numMVs; ++j ) {
            //GO curGlobalCol = maxColGid - numMVs + j + TGOT::one();
            GO curGlobalCol = minColGid + k*numMVs + j;
            globalColsArray.push_back(curGlobalCol);
            //TODO  extract the g2l map outside of this loop loop
            LO curLocalCol = domainMap.getLocalElement(curGlobalCol);
            if (curLocalCol != TLOT::invalid()) {
              eiData[j][curLocalCol] = TGOT::one();
              localColsArray.push_back(curLocalCol);
            }
          }

          // drop host views before apply
          for (size_t i=0; i<numMVs; ++i)
            eiData[i] = Teuchos::null;
          // probe
          A.apply(*ei,*colsA);

          colsOnPid0->doImport(*colsA,importer,INSERT);

          if (myRank==0)
            globalNnz += writeColumns(os,*colsOnPid0, numMVs, importedGidsData(),
                                      globalColsArray, offsetToUseInPrinting);

          // reconstruct dropped eiData
          for (size_t i=0; i<numMVs; ++i)
            eiData[i] = ei->getDataNonConst(i);
          for (size_t j=0; j<numMVs; ++j ) {
            GO curGlobalCol = minColGid + k*numMVs + j;
            //TODO  extract the g2l map outside of this loop loop
            LO curLocalCol = domainMap.getLocalElement(curGlobalCol);
            if (curLocalCol != TLOT::invalid()) {
              eiData[j][curLocalCol] = TGOT::one();
            }
          }

          //zero out the ei's
          for (size_t j=0; j<numMVs; ++j ) {
            for (int i=0; i<localColsArray.size(); ++i)
              eiData[j][localColsArray[i]] = TGOT::zero();
          }
          globalColsArray.clear();
          localColsArray.clear();

        }

        // //////////////////////////////////////
        // Handle leftover part of A
        // //////////////////////////////////////
        if (rem > 0) {
          for (int j=0; j<rem; ++j ) {
            GO curGlobalCol = maxColGid - rem + j + TGOT::one();
            globalColsArray.push_back(curGlobalCol);
            //TODO  extract the g2l map outside of this loop loop
            LO curLocalCol = domainMap.getLocalElement(curGlobalCol);
            if (curLocalCol != TLOT::invalid()) {
              eiData[j][curLocalCol] = TGOT::one();
              localColsArray.push_back(curLocalCol);
            }
          }

          // drop host views before apply
          for (size_t i=0; i<numMVs; ++i)
            eiData[i] = Teuchos::null;
          // probe
          A.apply(*ei,*colsA);

          colsOnPid0->doImport(*colsA,importer,INSERT);
          if (myRank==0)
            globalNnz += writeColumns(os,*colsOnPid0, rem, importedGidsData(),
                                      globalColsArray, offsetToUseInPrinting);

          // reconstruct dropped eiData
          for (size_t i=0; i<numMVs; ++i)
            eiData[i] = ei->getDataNonConst(i);
          for (int j=0; j<rem; ++j ) {
            GO curGlobalCol = maxColGid - rem + j + TGOT::one();
            //TODO  extract the g2l map outside of this loop loop
            LO curLocalCol = domainMap.getLocalElement(curGlobalCol);
            if (curLocalCol != TLOT::invalid()) {
              eiData[j][curLocalCol] = TGOT::one();
            }
          }

          //zero out the ei's
          for (int j=0; j<rem; ++j ) {
            for (int i=0; i<localColsArray.size(); ++i)
              eiData[j][localColsArray[i]] = TGOT::zero();
          }
          globalColsArray.clear();
          localColsArray.clear();

        }

        // Return the Matrix Market header.  It includes the header
        // line (that starts with "%%"), some comments, and the triple
        // of matrix dimensions and number of nonzero entries.  We
        // don't actually print this here, because we don't know the
        // number of nonzero entries until after probing.
        std::ostringstream oss;
        if (myRank == 0) {
          oss << "%%MatrixMarket matrix coordinate ";
          if (Teuchos::ScalarTraits<typename operator_type::scalar_type>::isComplex) {
            oss << "complex";
          } else {
            oss << "real";
          }
          oss << " general" << std::endl;
          oss << "% Tpetra::Operator" << std::endl;
          std::time_t now = std::time(NULL);
          oss << "% time stamp: " << ctime(&now);
          oss << "% written from " << numProcs << " processes" << std::endl;
          size_t numRows = rangeMap->getGlobalNumElements();
          size_t numCols = domainMap.getGlobalNumElements();
          oss << numRows << " " << numCols << " " << globalNnz << std::endl;
        }

        return oss.str ();
      }

      static global_ordinal_type
      writeColumns(std::ostream& os, mv_type const &colsA, size_t const &numCols,
                   Teuchos::ArrayView<const global_ordinal_type> const &rowGids,
                   Teuchos::Array<global_ordinal_type> const &colsArray,
                   global_ordinal_type const & indexBase) {

      typedef global_ordinal_type           GO;
      typedef scalar_type                   Scalar;
      typedef Teuchos::ScalarTraits<Scalar> STS;

        GO nnz=0;
        const Scalar zero = STS::zero();
        const size_t numRows = colsA.getGlobalLength();
        for (size_t j=0; j<numCols; ++j) {
          Teuchos::ArrayRCP<const Scalar> const curCol = colsA.getData(j);
          const GO J = colsArray[j];
          for (size_t i=0; i<numRows; ++i) {
            const Scalar val = curCol[i];
            if (val!=zero) {
              os << rowGids[i]+indexBase << " " << J+indexBase << " " << val << std::endl;
              ++nnz;
            }
          }
        }

        return nnz;

      }

    public:


      /// \brief Write a Tpetra::CrsMatrix to a file per rank.
      ///
      //! Function to write a one-per-rank collection of MatrixMarket files
      //! and assemble it into a single big matrix.  The code will try to minimize
      //! the number of ranks hammering on the file system at once, but we don't
      //! make any guarantees.
      
      static
      void 
      writeSparsePerRank (const std::string& filename_prefix,
                          const std::string& filename_suffix,
                          const sparse_matrix_type& matrix,
                          const std::string& matrixName,
                          const std::string& matrixDescription,
                          const int ranksToWriteAtOnce=8,
                          const bool debug=false) {
        
        using ST = scalar_type;
        //using LO = local_ordinal_type;
        using GO = global_ordinal_type;
        using STS = typename Teuchos::ScalarTraits<ST>;
        using Teuchos::RCP;
        
        // Sanity Checks
        trcp_tcomm_t comm = matrix.getComm ();
        TEUCHOS_TEST_FOR_EXCEPTION
          (comm.is_null (), std::invalid_argument,
           "The input matrix's communicator (Teuchos::Comm object) is null.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (matrix.isGloballyIndexed() || !matrix.isFillComplete(), std::invalid_argument,
           "The input matrix must not be GloballyIndexed and must be fillComplete.");
        
        // Setup
        const int myRank  = comm->getRank();
        const int numProc = comm->getSize();
        std::string filename = filename_prefix + std::to_string(myRank) + filename_suffix;  
        RCP<const map_type> rowMap = matrix.getRowMap();
        RCP<const map_type> colMap = matrix.getColMap();  
        size_t local_nnz      = matrix.getLocalNumEntries();
        size_t local_num_rows = rowMap->getLocalNumElements();
        size_t local_num_cols = colMap->getLocalNumElements();
        const GO rowIndexBase = rowMap->getIndexBase();
        const GO colIndexBase = colMap->getIndexBase();
        
        // Bounds check the writing limits
        int rank_limit = std::min(std::max(ranksToWriteAtOnce,1),numProc);
        
        // Start the writing
        for(int base_rank = 0; base_rank < numProc; base_rank += rank_limit) {
          int stop = std::min(base_rank+rank_limit,numProc);

          if(base_rank <= myRank  && myRank < stop) {          
            // My turn to write
            std::ofstream out(filename);
            
            // MatrixMarket Header
            out << "%%MatrixMarket matrix coordinate "
                << (STS::isComplex ? "complex" : "real")
                << " general" << std::endl;
            
            // Print comments (the matrix name and / or description).
            if (matrixName != "") {
              printAsComment (out, matrixName);
            }
            if (matrixDescription != "") {
              printAsComment (out, matrixDescription);
            }
            
            // Print the Matrix Market header (# local rows, # local columns, #
            // local enonzeros).  This will *not* be read correctly by a generic matrix
            // market reader since we'll be writing out GIDs here and local row/col counts
            out << local_num_rows << " " << local_num_cols << " " << local_nnz <<std::endl;
            
            {
              // Make the output stream write floating-point numbers in
              // scientific notation.  It will politely put the output
              // stream back to its state on input, when this scope
              // terminates.
              Teuchos::SetScientific<ST> sci (out);
              
              for(size_t l_row = 0; l_row < local_num_rows; l_row++) { 
                GO g_row = rowMap->getGlobalElement(l_row);            
                
                typename sparse_matrix_type::local_inds_host_view_type indices;
                typename sparse_matrix_type::values_host_view_type values;
                matrix.getLocalRowView(l_row, indices, values);
                for (size_t ii = 0; ii < indices.extent(0); ii++) {
                  const GO g_col = colMap->getGlobalElement(indices(ii));
                  // Convert row and column indices to 1-based.
                  // This works because the global index type is signed.
                  out << (g_row + 1 - rowIndexBase) << " "
                      << (g_col + 1 - colIndexBase) << " ";
                  if (STS::isComplex) {
                    out << STS::real(values(ii)) << " " << STS::imag(values(ii));
                  } else {
                    out << values(ii);
                  }
                  out << std::endl;
                } // For each entry in the current row
              } // For each row of the matrix
            }// end Teuchos::SetScientfic scoping
            
            out.close();
          }// end if base_rank <= myRank < stop
          
          // Barrier after each writing "batch" to make sure we're not hammering the file system
          // too aggressively
          comm->barrier();
          
        }// end outer loop
         
      }// end writeSparsePerRank

      //! Return @p obj MPI communicator or @ref Teuchos::null.
      template <typename T>
      static inline trcp_tcomm_t getComm(const Teuchos::RCP<T>& obj)
      {
        return obj.is_null() ? Teuchos::null : obj->getComm();
      }

      //! Return MPI rank or 0.
      static inline int getRank(const trcp_tcomm_t& comm)
      {
        return comm.is_null() ? 0 : comm->getRank();
      }

    }; // class Writer

  } // namespace MatrixMarket
} // namespace Tpetra

#endif // __MatrixMarket_Tpetra_hpp
