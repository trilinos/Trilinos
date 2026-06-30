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

/// \class MatrixMarketReader
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
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MatrixMarketReader {
 public:
  //! This class' template parameter; a specialization of CrsMatrix.
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> SparseMatrixType;
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
                   node_type>
      sparse_graph_type;

  //! The MultiVector specialization associated with SparseMatrixType.
  typedef MultiVector<scalar_type,
                      local_ordinal_type,
                      global_ordinal_type,
                      node_type>
      multivector_type;

  //! The Vector specialization associated with SparseMatrixType.
  typedef Vector<scalar_type,
                 local_ordinal_type,
                 global_ordinal_type,
                 node_type>
      vector_type;

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
  makeRangeMap(const trcp_tcomm_t& pComm,
               const global_ordinal_type numRows);

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
  makeRowMap(const Teuchos::RCP<const map_type>& pRowMap,
             const trcp_tcomm_t& pComm,
             const global_ordinal_type numRows);

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
  makeDomainMap(const Teuchos::RCP<const map_type>& pRangeMap,
                const global_ordinal_type numRows,
                const global_ordinal_type numCols);

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
  distribute(Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
             Teuchos::ArrayRCP<size_t>& myRowPtr,
             Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
             Teuchos::ArrayRCP<scalar_type>& myValues,
             const Teuchos::RCP<const map_type>& pRowMap,
             Teuchos::ArrayRCP<size_t>& numEntriesPerRow,
             Teuchos::ArrayRCP<size_t>& rowPtr,
             Teuchos::ArrayRCP<global_ordinal_type>& colInd,
             Teuchos::ArrayRCP<scalar_type>& values,
             const bool debug = false);

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
  makeMatrix(Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
             Teuchos::ArrayRCP<size_t>& myRowPtr,
             Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
             Teuchos::ArrayRCP<scalar_type>& myValues,
             const Teuchos::RCP<const map_type>& pRowMap,
             const Teuchos::RCP<const map_type>& pRangeMap,
             const Teuchos::RCP<const map_type>& pDomainMap,
             const bool callFillComplete = true);

  /// \brief Variant of makeMatrix() that takes parameters for
  ///   CrsMatrix's constructor and for fillComplete().
  ///
  /// Each process inserts its data into the sparse matrix, and
  /// then all processes call fillComplete().
  static Teuchos::RCP<sparse_matrix_type>
  makeMatrix(Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
             Teuchos::ArrayRCP<size_t>& myRowPtr,
             Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
             Teuchos::ArrayRCP<scalar_type>& myValues,
             const Teuchos::RCP<const map_type>& pRowMap,
             const Teuchos::RCP<const map_type>& pRangeMap,
             const Teuchos::RCP<const map_type>& pDomainMap,
             const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
             const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams);

  /// \brief Variant of makeMatrix() that takes an optional column Map.
  ///
  /// This method computes \c colMap only if it is null on input,
  /// and if \c callFillComplete is true.
  static Teuchos::RCP<sparse_matrix_type>
  makeMatrix(Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
             Teuchos::ArrayRCP<size_t>& myRowPtr,
             Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
             Teuchos::ArrayRCP<scalar_type>& myValues,
             const Teuchos::RCP<const map_type>& rowMap,
             Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const bool callFillComplete = true);

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
  readBanner(std::istream& in,
             size_t& lineNumber,
             const bool tolerant    = false,
             const bool /* debug */ = false,
             const bool isGraph     = false);

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
  readCoordDims(std::istream& in,
                size_t& lineNumber,
                const Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                const trcp_tcomm_t& pComm,
                const bool tolerant    = false,
                const bool /* debug */ = false);

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
  typedef Teuchos::MatrixMarket::SymmetrizingAdder<Teuchos::MatrixMarket::Raw::Adder<scalar_type, global_ordinal_type>> adder_type;

  typedef Teuchos::MatrixMarket::SymmetrizingGraphAdder<Teuchos::MatrixMarket::Raw::GraphAdder<global_ordinal_type>> graph_adder_type;

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
  makeAdder(const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
            Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
            const Teuchos::Tuple<global_ordinal_type, 3>& dims,
            const bool tolerant = false,
            const bool debug    = false);

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
  makeGraphAdder(const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
                 Teuchos::RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                 const Teuchos::Tuple<global_ordinal_type, 3>& dims,
                 const bool tolerant = false,
                 const bool debug    = false);

  /// \brief Read sparse graph from the given Matrix Market input stream.
  static Teuchos::RCP<sparse_graph_type>
  readSparseGraphHelper(std::istream& in,
                        const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
                        const Teuchos::RCP<const map_type>& rowMap,
                        Teuchos::RCP<const map_type>& colMap,
                        const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                        const bool tolerant,
                        const bool debug);

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
  /// \param comm [in] Communicator containing all processor(s)
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
  readSparseGraphFile(const std::string& filename,
                      const trcp_tcomm_t& comm,
                      const bool callFillComplete = true,
                      const bool tolerant         = false,
                      const bool debug            = false);

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
  readSparseGraphFile(const std::string& filename,
                      const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
                      const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                      const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                      const bool tolerant = false,
                      const bool debug    = false);

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
  readSparseGraphFile(const std::string& filename,
                      const Teuchos::RCP<const map_type>& rowMap,
                      Teuchos::RCP<const map_type>& colMap,
                      const Teuchos::RCP<const map_type>& domainMap,
                      const Teuchos::RCP<const map_type>& rangeMap,
                      const bool callFillComplete = true,
                      const bool tolerant         = false,
                      const bool debug            = false);

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
  readSparseGraph(std::istream& in,
                  const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
                  const bool callFillComplete = true,
                  const bool tolerant         = false,
                  const bool debug            = false);

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
  readSparseGraph(std::istream& in,
                  const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
                  const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                  const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                  const bool tolerant = false,
                  const bool debug    = false);

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
  readSparseGraph(std::istream& in,
                  const Teuchos::RCP<const map_type>& rowMap,
                  Teuchos::RCP<const map_type>& colMap,
                  const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const bool callFillComplete = true,
                  const bool tolerant         = false,
                  const bool debug            = false);

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
  /// \param comm [in] Communicator containing all processor(s)
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
  readSparseFile(const std::string& filename,
                 const trcp_tcomm_t& comm,
                 const bool callFillComplete = true,
                 const bool tolerant         = false,
                 const bool debug            = false);

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
  /// \param comm [in] Communicator containing all process(es)
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
  readSparseFile(const std::string& filename,
                 const trcp_tcomm_t& comm,
                 const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                 const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
                 const bool tolerant = false,
                 const bool debug    = false);

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
  readSparseFile(const std::string& filename,
                 const Teuchos::RCP<const map_type>& rowMap,
                 Teuchos::RCP<const map_type>& colMap,
                 const Teuchos::RCP<const map_type>& domainMap,
                 const Teuchos::RCP<const map_type>& rangeMap,
                 const bool callFillComplete = true,
                 const bool tolerant         = false,
                 const bool debug            = false);

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
  readSparse(std::istream& in,
             const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
             const bool callFillComplete = true,
             const bool tolerant         = false,
             const bool debug            = false);

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
  readSparse(std::istream& in,
             const Teuchos::RCP<const Teuchos::Comm<int>>& pComm,
             const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
             const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams,
             const bool tolerant = false,
             const bool debug    = false);

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
  readSparse(std::istream& in,
             const Teuchos::RCP<const map_type>& rowMap,
             Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const bool callFillComplete = true,
             const bool tolerant         = false,
             const bool debug            = false);

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
  /// \param binary [in] If true, read in binary mode.
  static Teuchos::RCP<multivector_type>
  readDenseFile(const std::string& filename,
                const trcp_tcomm_t& comm,
                Teuchos::RCP<const map_type>& map,
                const bool tolerant = false,
                const bool debug    = false,
                const bool binary   = false);

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
      std::ios_base::openmode mode = std::ios_base::in);

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
  readVectorFile(const std::string& filename,
                 const trcp_tcomm_t& comm,
                 Teuchos::RCP<const map_type>& map,
                 const bool tolerant = false,
                 const bool debug    = false);

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
  /// \param binary [in] If true, read in binary mode.

  static Teuchos::RCP<multivector_type>
  readDense(std::istream& in,
            const trcp_tcomm_t& comm,
            Teuchos::RCP<const map_type>& map,
            const bool tolerant = false,
            const bool debug    = false,
            const bool binary   = false);

  //! Read Vector from the given Matrix Market input stream.
  static Teuchos::RCP<vector_type>
  readVector(std::istream& in,
             const trcp_tcomm_t& comm,
             Teuchos::RCP<const map_type>& map,
             const bool tolerant = false,
             const bool debug    = false);

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
  /// \param binary [in] If true, read in binary mode.
  static Teuchos::RCP<const map_type>
  readMapFile(const std::string& filename,
              const trcp_tcomm_t& comm,
              const bool tolerant = false,
              const bool debug    = false,
              const bool binary   = false);

 private:
  template <class MultiVectorScalarType>
  static Teuchos::RCP<Tpetra::MultiVector<MultiVectorScalarType,
                                          local_ordinal_type,
                                          global_ordinal_type,
                                          node_type>>
  readDenseImpl(std::istream& in,
                const trcp_tcomm_t& comm,
                Teuchos::RCP<const map_type>& map,
                const Teuchos::RCP<Teuchos::FancyOStream>& err,
                const bool tolerant = false,
                const bool debug    = false,
                const bool binary   = false);

  template <class VectorScalarType>
  static Teuchos::RCP<Tpetra::Vector<VectorScalarType,
                                     local_ordinal_type,
                                     global_ordinal_type,
                                     node_type>>
  readVectorImpl(std::istream& in,
                 const trcp_tcomm_t& comm,
                 Teuchos::RCP<const map_type>& map,
                 const Teuchos::RCP<Teuchos::FancyOStream>& err,
                 const bool tolerant = false,
                 const bool debug    = false);

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
  /// \param binary [in] If true, read in binary mode.
  static Teuchos::RCP<const map_type>
  readMap(std::istream& in,
          const trcp_tcomm_t& comm,
          const bool tolerant = false,
          const bool debug    = false,
          const bool binary   = false);

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
  /// \param binary [in] If true, read in binary mode.
  static Teuchos::RCP<const map_type>
  readMap(std::istream& in,
          const trcp_tcomm_t& comm,
          const Teuchos::RCP<Teuchos::FancyOStream>& err,
          const bool tolerant = false,
          const bool debug    = false,
          const bool binary   = false);

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
  encodeDataType(const std::string& dataType);

 public:
  /// \brief Read a Tpetra::CrsMatrix from a file per rank setup
  ///
  //! Function to read a one-per-rank collection of MatrixMarket files
  //! and assemble it into a single big matrix.  The code will try to minimize
  //! the number of ranks hammering on the file system at once, but we don't
  //! make any guarantees.
  /// \param filename_prefix  [in] File for rank I is filename_prefix + to_string(I) + filename_suffix
  /// \param filename_suffix [in] File for rank I is filename_prefix + to_string(I) + filename_suffix
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
  /// \param ranksToReadAtOnce [in] Number of ranks to read at once.
  /// \param debug [in] Whether to produce copious status output
  ///   useful for Tpetra developers, but probably not useful for
  ///   anyone else.
  static Teuchos::RCP<sparse_matrix_type>
  readSparsePerRank(const std::string& filename_prefix,
                    const std::string& filename_suffix,
                    const Teuchos::RCP<const map_type>& rowMap,
                    Teuchos::RCP<const map_type>& colMap,
                    const Teuchos::RCP<const map_type>& domainMap,
                    const Teuchos::RCP<const map_type>& rangeMap,
                    const bool callFillComplete = true,
                    const bool tolerant         = false,
                    const int ranksToReadAtOnce = 8,
                    const bool debug            = false);  // end readSparsePerRank

};  // class MatrixMarketReader

/// \class MatrixMarketWriter
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
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MatrixMarketWriter {
 public:
  //! Template parameter of this class; specialization of CrsMatrix.
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> SparseMatrixType;
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
                      node_type>
      multivector_type;
  //! Specialization of Tpetra::Map that matches SparseMatrixType.
  typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  //! Specialization of Tpetra::CrsGraph that matches SparseMatrixType.
  typedef CrsGraph<local_ordinal_type, global_ordinal_type, node_type> crs_graph_type;

  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> operator_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;

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
  writeSparseFile(const std::string& filename,
                  const sparse_matrix_type& matrix,
                  const std::string& matrixName,
                  const std::string& matrixDescription,
                  const bool debug = false);

  //! Only for backwards compatibility; prefer the overload above.
  static void
  writeSparseFile(const std::string& filename,
                  const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                  const std::string& matrixName,
                  const std::string& matrixDescription,
                  const bool debug = false);

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
  writeSparseFile(const std::string& filename,
                  const sparse_matrix_type& matrix,
                  const bool debug = false);

  //! Only for backwards compatibility; prefer the overload above.
  static void
  writeSparseFile(const std::string& filename,
                  const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
                  const bool debug = false);

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
  writeSparse(std::ostream& out,
              const sparse_matrix_type& matrix,
              const std::string& matrixName,
              const std::string& matrixDescription,
              const bool debug = false);

  //! Only for backwards compatibility; prefer the overload above.
  static void
  writeSparse(std::ostream& out,
              const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
              const std::string& matrixName,
              const std::string& matrixDescription,
              const bool debug = false);

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
  writeSparseGraph(std::ostream& out,
                   const crs_graph_type& graph,
                   const std::string& graphName,
                   const std::string& graphDescription,
                   const bool debug = false);

  /// \brief Print the sparse graph in Matrix Market format to the
  ///   given output stream, with no comments.
  ///
  /// See the above five-argument version of this function for
  /// full documentation.
  static void
  writeSparseGraph(std::ostream& out,
                   const crs_graph_type& graph,
                   const bool debug = false);

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
  writeSparseGraphFile(const std::string& filename,
                       const crs_graph_type& graph,
                       const std::string& graphName,
                       const std::string& graphDescription,
                       const bool debug = false);

  /// \brief Print the sparse graph in Matrix Market format to the
  ///   given file (by filename), with no comments.
  ///
  /// See the above five-argument overload for full documentation.
  static void
  writeSparseGraphFile(const std::string& filename,
                       const crs_graph_type& graph,
                       const bool debug = false);

  /// \brief Print the sparse graph in Matrix Market format to the
  ///   given file (by filename), taking the graph by Teuchos::RCP.
  ///
  /// This is just a convenience for users who don't want to
  /// remember to dereference the Teuchos::RCP.  For
  /// documentation, see the above overload of this function that
  /// takes the graph by const reference, rather than by
  /// Teuchos::RCP.
  static void
  writeSparseGraphFile(const std::string& filename,
                       const Teuchos::RCP<const crs_graph_type>& pGraph,
                       const std::string& graphName,
                       const std::string& graphDescription,
                       const bool debug = false);

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
  writeSparseGraphFile(const std::string& filename,
                       const Teuchos::RCP<const crs_graph_type>& pGraph,
                       const bool debug = false);

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
  writeSparse(std::ostream& out,
              const sparse_matrix_type& matrix,
              const bool debug = false);

  //! Only for backwards compatibility; prefer the overload above.
  static void
  writeSparse(std::ostream& out,
              const Teuchos::RCP<const sparse_matrix_type>& pMatrix,
              const bool debug = false);

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
  writeDenseFile(const std::string& filename,
                 const multivector_type& X,
                 const std::string& matrixName,
                 const std::string& matrixDescription,
                 const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                 const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

  /// \brief Print the multivector in Matrix Market format, with
  ///   matrix name and description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDenseFile().
  static void
  writeDenseFile(const std::string& filename,
                 const Teuchos::RCP<const multivector_type>& X,
                 const std::string& matrixName,
                 const std::string& matrixDescription,
                 const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                 const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

  /// \brief Print the multivector in Matrix Market format, with
  ///   no matrix name or description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDenseFile().
  static void
  writeDenseFile(const std::string& filename,
                 const multivector_type& X,
                 const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                 const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

  /// \brief Print the multivector in Matrix Market format, with
  ///   no matrix name or description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDenseFile().
  static void
  writeDenseFile(const std::string& filename,
                 const Teuchos::RCP<const multivector_type>& X,
                 const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                 const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

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
  writeDense(std::ostream& out,
             const multivector_type& X,
             const std::string& matrixName,
             const std::string& matrixDescription,
             const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
             const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

 private:
  /**
   * @brief Open a file only on rank zero, possibly throwing if the stream is invalid.
   *
   * @note On processes that are not the rank zero process, the stream is left uninitialized.
   */
  static std::ofstream openOutFileOnRankZero(
      const trcp_tcomm_t& comm,
      const std::string& filename, const int rank, const bool safe = true,
      const std::ios_base::openmode mode = std::ios_base::out);

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
  writeDenseHeader(std::ostream& out,
                   const multivector_type& X,
                   const std::string& matrixName,
                   const std::string& matrixDescription,
                   const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                   const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

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
  writeDenseColumn(std::ostream& out,
                   const multivector_type& X,
                   const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
                   const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

 public:
  /// \brief Print the multivector in Matrix Market format, with
  ///   matrix name and or description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDense().
  static void
  writeDense(std::ostream& out,
             const Teuchos::RCP<const multivector_type>& X,
             const std::string& matrixName,
             const std::string& matrixDescription,
             const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
             const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

  /// \brief Print the multivector in Matrix Market format, with
  ///   no matrix name or description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDense().
  static void
  writeDense(std::ostream& out,
             const multivector_type& X,
             const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
             const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

  /// \brief Print the multivector in Matrix Market format, with
  ///   no matrix name or description.
  ///
  /// See the documentation of the above six-argument version of
  /// writeDense().
  static void
  writeDense(std::ostream& out,
             const Teuchos::RCP<const multivector_type>& X,
             const Teuchos::RCP<Teuchos::FancyOStream>& err = Teuchos::null,
             const Teuchos::RCP<Teuchos::FancyOStream>& dbg = Teuchos::null);

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
  writeMap(std::ostream& out, const map_type& map, const bool debug = false);

  /// \brief Print the Map to the given output stream \c out.
  ///
  /// This version of writeMap() comes with an extra debug output
  /// stream \c err, that is only used if \c debug is true.
  ///
  /// \warning We make no promises of backwards compatibility with
  ///   this method.  It may go away or its interface may change
  ///   at any time.
  static void
  writeMap(std::ostream& out,
           const map_type& map,
           const Teuchos::RCP<Teuchos::FancyOStream>& err,
           const bool debug = false);

  //! Write the Map to the given file.
  static void
  writeMapFile(const std::string& filename,
               const map_type& map);

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
  printAsComment(std::ostream& out, const std::string& str);

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
  writeOperator(const std::string& fileName, operator_type const& A);

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
  writeOperator(std::ostream& out, const operator_type& A);

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
  writeOperator(const std::string& fileName,
                const operator_type& A,
                const Teuchos::ParameterList& params);

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
  writeOperator(std::ostream& out,
                const operator_type& A,
                const Teuchos::ParameterList& params);

 private:
  /// \brief Implementation detail of writeOperator overloads.
  ///
  /// Use column probing to discover the entries of the Operator.
  /// Write them in Matrix Market sparse matrix format, on Process
  /// 0 only, to the output stream \c os.  Do NOT write the Matrix
  /// Market header line, but DO return it (unless the input
  /// options specify otherwise).
  static std::string
  writeOperatorImpl(std::ostream& os,
                    const operator_type& A,
                    const Teuchos::ParameterList& params);

  static global_ordinal_type
  writeColumns(std::ostream& os, mv_type const& colsA, size_t const& numCols,
               Teuchos::ArrayView<const global_ordinal_type> const& rowGids,
               Teuchos::Array<global_ordinal_type> const& colsArray,
               global_ordinal_type const& indexBase);

 public:
  /// \brief Write a Tpetra::CrsMatrix to a file per rank.
  ///
  //! Function to write a one-per-rank collection of MatrixMarket files
  //! and assemble it into a single big matrix.  The code will try to minimize
  //! the number of ranks hammering on the file system at once, but we don't
  //! make any guarantees.

  static void
  writeSparsePerRank(const std::string& filename_prefix,
                     const std::string& filename_suffix,
                     const sparse_matrix_type& matrix,
                     const std::string& matrixName,
                     const std::string& matrixDescription,
                     const int ranksToWriteAtOnce = 8,
                     const bool debug             = false);  // end writeSparsePerRank

  //! Return @p obj MPI communicator or @ref Teuchos::null.
  template <typename T>
  static trcp_tcomm_t getComm(const Teuchos::RCP<T>& obj);

  //! Return MPI rank or 0.
  static int getRank(const trcp_tcomm_t& comm);

};  // class MatrixMarketWriter

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

template <class SparseMatrixType>
using Reader = MatrixMarketReader<typename SparseMatrixType::scalar_type,
                                  typename SparseMatrixType::local_ordinal_type,
                                  typename SparseMatrixType::global_ordinal_type,
                                  typename SparseMatrixType::node_type>;

template <class SparseMatrixType>
using Writer = MatrixMarketWriter<typename SparseMatrixType::scalar_type,
                                  typename SparseMatrixType::local_ordinal_type,
                                  typename SparseMatrixType::global_ordinal_type,
                                  typename SparseMatrixType::node_type>;

}  // namespace MatrixMarket
}  // namespace Tpetra

#endif  // __MatrixMarket_Tpetra_hpp
