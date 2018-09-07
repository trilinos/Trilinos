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

#ifndef TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DECL_HPP
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DECL_HPP

/// \file Tpetra_Experimental_BlockCrsMatrix_decl.hpp
/// \brief Declaration of Tpetra::Experimental::BlockCrsMatrix

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Experimental_BlockMultiVector.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"

namespace Tpetra {
namespace Experimental {
namespace Classes {

/// \class BlockCrsMatrix
/// \brief Sparse matrix whose entries are small dense square blocks,
///   all of the same dimensions.
///
/// \tparam Scalar The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike in
///   Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the documentation of
///   the first template parameter of Map for requirements.
/// \tparam GO The type of global indices.  See the documentation of
///   the second template parameter of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of the
///   third template parameter of Map for requirements.
///
/// Please read the documentation of BlockMultiVector first.
///
/// This class implements a sparse matrix whose entries are small
/// dense square blocks, all of the same dimensions.  The intended
/// application is to store the discretization of a partial
/// differential equation with multiple degrees of freedom per mesh
/// point, where all mesh points have the same number of degrees of
/// freedom.  This class stores values associated with the degrees of
/// freedom of a single mesh point contiguously, in a getBlockSize()
/// by getBlockSize() block, in row-major format.  The matrix's graph
/// represents the mesh points, with one entry per mesh point.  This
/// saves storage over using a CrsMatrix, which requires one graph
/// entry per matrix entry.
///
/// This class requires a fill-complete Tpetra::CrsGraph for
/// construction.  Thus, it has a row Map and a column Map already.
/// As a result, BlockCrsMatrix only needs to provide access using
/// local indices.  Access using local indices is faster anyway, since
/// conversion from global to local indices requires a hash table
/// lookup per index.  Users are responsible for converting from
/// global to local indices if necessary.  Please be aware that the
/// row Map and column Map may differ, so you may not use local row
/// and column indices interchangeably.
///
/// We reserve the right to change the block layout in the future.
/// Best practice is to use this class' little_block_type and
/// const_little_block_type typedefs to access blocks, since their
/// types offer access to entries in a layout-independent way.
/// These two typedefs are both Kokkos::View specializations.
///
/// Here is an example of how to fill into this object using
/// raw-pointer views.
///
/// \code
/// int err = 0;
/// // At least one entry, so &offsets[0] always makes sense.
/// Teuchos::Array<ptrdiff_t> offsets (1);
/// for (LO localRowInd = 0; localRowInd < localNumRows; ++localRowInd) {
///   // Get a view of the current row.
///   // You may modify the values, but not the column indices.
///   const LO* localColInds;
///   Scalar* vals;
///   LO numEntries;
///   err = A.getLocalRowView (localRowInd, localColInds, vals, numEntries);
///   if (err != 0) {
///     break;
///   }
///
///   // Modify the entries in the current row.
///   for (LO k = 0; k < numEntries; ++k) {
///     Scalar* const curBlock = vals[blockSize * blockSize * k];
///     // Blocks are stored in row-major format.
///     for (LO j = 0; j < blockSize; ++j) {
///       for (LO i = 0; i < blockSize; ++i) {
///         const Scalar curVal = &curBlock[i + j * blockSize];
///         // Some function f of the current value and mesh point
///         curBlock[i + j * blockSize] = f (curVal, localColInds[k], ...);
///       }
///     }
///   }
/// }
/// \endcode
///
template<class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LO = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GO = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node = ::Tpetra::Details::DefaultTypes::node_type>
class BlockCrsMatrix :
  virtual public ::Tpetra::RowMatrix<Scalar, LO, GO, Node>,
  virtual public ::Tpetra::DistObject<char, LO, GO, Node>
{
private:
  using dist_object_type = ::Tpetra::DistObject<char, LO, GO, Node>;
  using BMV = BlockMultiVector<Scalar, LO, GO, Node>;
  using STS = Teuchos::ScalarTraits<Scalar>;

protected:
  //! Implementation detail; tells
  typedef char packet_type;

public:
  //! \name Public typedefs
  //@{

  //! The type of entries in the matrix (that is, of each entry in each block).
  using scalar_type = Scalar;

  /// \brief The implementation type of entries in the matrix.
  ///
  /// Letting scalar_type and impl_scalar_type differ helps this class
  /// work correctly for Scalar types like std::complex<T>, which lack
  /// the necessary CUDA device macros and volatile overloads to work
  /// correctly with Kokkos.
  using impl_scalar_type = typename BMV::impl_scalar_type;

  //! The type of local indices.
  typedef LO local_ordinal_type;
  //! The type of global indices.
  typedef GO global_ordinal_type;
  /// \brief The Node type.
  ///
  /// Prefer device_type, execution_space, and memory_space (see
  /// below), which relate directly to Kokkos.
  typedef Node node_type;

  //! The Kokkos::Device specialization that this class uses.
  typedef typename Node::device_type device_type;
  //! The Kokkos execution space that this class uses.
  typedef typename device_type::execution_space execution_space;
  //! The Kokkos memory space that this class uses.
  typedef typename device_type::memory_space memory_space;

  //! The implementation of Map that this class uses.
  typedef ::Tpetra::Map<LO, GO, node_type> map_type;
  //! The implementation of MultiVector that this class uses.
  typedef ::Tpetra::MultiVector<Scalar, LO, GO, node_type> mv_type;
  //! The implementation of CrsGraph that this class uses.
  typedef ::Tpetra::CrsGraph<LO, GO, node_type> crs_graph_type;

  //! The type used to access nonconst matrix blocks.
  typedef Kokkos::View<impl_scalar_type**,
                       Kokkos::LayoutRight,
                       device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
          little_block_type;
  //! The type used to access const matrix blocks.
  typedef Kokkos::View<const impl_scalar_type**,
                       Kokkos::LayoutRight,
                       device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
          const_little_block_type;
  //! The type used to access nonconst vector blocks.
  typedef Kokkos::View<impl_scalar_type*,
                       Kokkos::LayoutRight,
                       device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
          little_vec_type;
  //! The type used to access const vector blocks.
  typedef Kokkos::View<const impl_scalar_type*,
                       Kokkos::LayoutRight,
                       device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
          const_little_vec_type;

  //@}
  //! \name Constructors and destructor
  //@{

  //! Default constructor: Makes an empty block matrix.
  BlockCrsMatrix ();

  /// \brief Constructor that takes a graph and a block size.
  ///
  /// The graph represents the mesh.  This constructor computes the
  /// point Maps corresponding to the given graph's domain and range
  /// Maps.  If you already have those point Maps, it is better to
  /// call the four-argument constructor.
  ///
  /// \param graph [in] A fill-complete graph.
  /// \param blockSize [in] Number of degrees of freedom per mesh point.
  BlockCrsMatrix (const crs_graph_type& graph, const LO blockSize);

  /// \brief Constructor that takes a graph, domain and range point
  ///   Maps, and a block size.
  ///
  /// The graph represents the mesh.  This constructor uses the given
  /// domain and range point Maps, rather than computing them.  The
  /// given point Maps must be the same as the above two-argument
  /// constructor would have computed.
  BlockCrsMatrix (const crs_graph_type& graph,
                  const map_type& domainPointMap,
                  const map_type& rangePointMap,
                  const LO blockSize);

  //! Destructor (declared virtual for memory safety).
  virtual ~BlockCrsMatrix () {}

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  //! Get the (point) domain Map of this matrix.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! Get the (point) range Map of this matrix.
  Teuchos::RCP<const map_type> getRangeMap () const;

  //! get the (mesh) map for the rows of this block matrix.
  Teuchos::RCP<const map_type> getRowMap () const;

  //! get the (mesh) map for the columns of this block matrix.
  Teuchos::RCP<const map_type> getColMap () const;

  //! get the global number of block rows
  global_size_t getGlobalNumRows() const;

  //! get the local number of block rows
  size_t getNodeNumRows() const;

  size_t getNodeMaxNumRowEntries() const;

  /// \brief For this matrix A, compute <tt>Y := beta * Y + alpha * Op(A) * X</tt>.
  ///
  /// Op(A) is A if mode is Teuchos::NO_TRANS, the transpose of A if
  /// mode is Teuchos::TRANS, and the conjugate transpose of A if mode
  /// is Teuchos::CONJ_TRANS.
  ///
  /// If alpha is zero, ignore X's entries on input; if beta is zero,
  /// ignore Y's entries on input.  This follows the BLAS convention,
  /// and only matters if X resp. Y have Inf or NaN entries.
  void
  apply (const mv_type& X,
         mv_type& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const;

  /// \brief Whether it is valid to apply the transpose or conjugate
  ///   transpose of this matrix.
  bool hasTransposeApply () const {
    // FIXME (mfh 04 May 2014) Transpose and conjugate transpose modes
    // are not implemented yet.  Fill in applyBlockTrans() to fix this.
    return false;
  }

  //! Set all matrix entries equal to \c alpha.
  void setAllToScalar (const Scalar& alpha);

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{

  //! One-line description of this object.
  std::string description () const;

  /// \brief Print a description of this object to the given output stream.
  ///
  /// \param out [out] Output stream to which to print.  Valid values
  ///   include Teuchos::VERB_DEFAULT, Teuchos::VERB_NONE,
  ///   Teuchos::VERB_LOW, Teuchos::VERB_MEDIUM, Teuchos::VERB_HIGH,
  ///   and Teuchos::VERB_EXTREME.
  ///
  /// \param verbLevel [in] Verbosity level at which to print.
  ///
  /// \warning If verbLevel is Teuchos::VERB_EXTREME, this method has
  ///   collective semantics over the matrix's communicator.
  ///
  /// The following pseudocode shows how to wrap your std::ostream
  /// object in a Teuchos::FancyOStream, and pass it into this method:
  /// \code
  /// Tpetra::Experimental::BlockCrsMatrix<...> A (...);
  /// // ...
  /// std::ostream& yourObject = ...;
  /// Teuchos::RCP<Teuchos::FancyOStream> wrappedStream =
  ///   Teuchos::getFancyOStream (Teuchos::rcpFromRef (yourObject));
  /// const Teuchos::EVerbosityLevel verbLevel = ...;
  /// A.describe (*wrappedStream, verbLevel);
  /// \endcode
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const;

  //@}
  //! \name Block operations
  //@{

  //! The number of degrees of freedom per mesh point.
  LO getBlockSize () const { return blockSize_; }

  //! Get the (mesh) graph.
  virtual Teuchos::RCP<const ::Tpetra::RowGraph<LO,GO,Node> > getGraph () const;

  const crs_graph_type & getCrsGraph () const { return graph_; }

  /// \brief Version of apply() that takes BlockMultiVector input and output.
  ///
  /// This method is deliberately not marked const, because it may do
  /// lazy initialization of temporary internal block multivectors.
  void
  applyBlock (const BlockMultiVector<Scalar, LO, GO, Node>& X,
              BlockMultiVector<Scalar, LO, GO, Node>& Y,
              Teuchos::ETransp mode = Teuchos::NO_TRANS,
              const Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
              const Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ());

  /// \brief Version of gaussSeidel(), with fewer requirements on X.
  ///
  /// Not Implemented
  void
  gaussSeidelCopy (MultiVector<Scalar,LO,GO,Node> &X,
                   const ::Tpetra::MultiVector<Scalar,LO,GO,Node> &B,
                   const ::Tpetra::MultiVector<Scalar,LO,GO,Node> &D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps,
                   const bool zeroInitialGuess) const;

  /// \brief Version of reorderedGaussSeidel(), with fewer requirements on X.
  ///
  /// Not Implemented
  void
  reorderedGaussSeidelCopy (::Tpetra::MultiVector<Scalar,LO,GO,Node>& X,
                            const ::Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
                            const ::Tpetra::MultiVector<Scalar,LO,GO,Node>& D,
                            const Teuchos::ArrayView<LO>& rowIndices,
                            const Scalar& dampingFactor,
                            const ESweepDirection direction,
                            const int numSweeps,
                            const bool zeroInitialGuess) const;

  /// \brief Local Gauss-Seidel solve, given a factorized diagonal
  ///
  /// \param Residual [in] The "residual" (right-hand side) block
  ///   (multi)vector
  /// \param Solution [in/out] On input: the initial guess / current
  ///   approximate solution.  On output: the new approximate
  ///   solution.
  /// \param D_inv [in] Block diagonal, the explicit inverse of this
  ///   matrix's block diagonal (possibly modified for algorithmic
  ///   reasons).
  /// \param omega [in] (S)SOR relaxation coefficient
  /// \param direction [in] Forward, Backward, or Symmetric.
  ///
  /// One may access block i in \c D_inv using the
  /// following code:
  /// \code
  /// auto D_ii = Kokkos::subview(D_inv, i, Kokkos::ALL(), Kokkos::ALL());
  /// \endcode
  /// The resulting block is b x b, where <tt>b = this->getBlockSize()</tt>.
  void
  localGaussSeidel (const BlockMultiVector<Scalar, LO, GO, Node>& Residual,
                          BlockMultiVector<Scalar, LO, GO, Node>& Solution,
                    const Kokkos::View<impl_scalar_type***, device_type,
                          Kokkos::MemoryUnmanaged>& D_inv,
                    const Scalar& omega,
                    const ESweepDirection direction) const;

  /// \brief Replace values at the given (mesh, i.e., block) column
  ///   indices, in the given (mesh, i.e., block) row.
  ///
  /// \param localRowInd [in] Local mesh (i.e., block) index of the
  ///   row in which to replace.
  ///
  /// \param colInds [in] Local mesh (i.e., block) column ind{ex,ices}
  ///   at which to replace values.  \c colInds[k] is the local column
  ///   index whose new values start at
  ///   <tt>vals[getBlockSize() * getBlockSize() * k]</tt>,
  ///   and \c colInds has length at least \c numColInds.  This method
  ///   will only access the first \c numColInds entries of \c colInds.
  ///
  /// \param vals [in] The new values to use at the given column
  ///   indices.  Values for each block are stored contiguously, in
  ///   row major layout, with no padding between rows or between
  ///   blocks.  Thus, if <tt>b = getBlockSize()</tt>, then
  ///   <tt>vals[k*b*b] .. vals[(k+1)*b*b-1]</tt> are the values to
  ///   use for block \c colInds[k].
  ///
  /// \param numColInds [in] The number of entries of \c colInds.
  ///
  /// \return The number of valid entries of \c colInds.  \c colInds[k]
  ///   is valid if and only if it is a valid local mesh (i.e., block)
  ///   column index.  This method succeeded if and only if the return
  ///   value equals the input argument \c numColInds.
  LO
  replaceLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const;

  /// \brief Sum into values at the given (mesh, i.e., block) column
  ///   indices, in the given (mesh, i.e., block) row.
  ///
  /// \param localRowInd [in] Local mesh (i.e., block) index of the
  ///   row in which to sum.
  ///
  /// \param colInds [in] Local mesh (i.e., block) column ind{ex,ices}
  ///   at which to sum.  \c colInds[k] is the local column index whose
  ///   new values start at
  ///   <tt>vals[getBlockSize() * getBlockSize() * k]</tt>,
  ///   and \c colInds has length at least \c numColInds.  This method
  ///   will only access the first \c numColInds entries of \c colInds.
  ///
  /// \param vals [in] The new values to sum in at the given column
  ///   indices.  Values for each block are stored contiguously, in
  ///   row major layout, with no padding between rows or between
  ///   blocks.  Thus, if <tt>b = getBlockSize()</tt>, then
  ///   <tt>vals[k*b*b] .. vals[(k+1)*b*b-1]</tt> are the values to
  ///   use for block \c colInds[k].
  ///
  /// \param numColInds [in] The number of entries of \c colInds.
  ///
  /// \return The number of valid entries of \c colInds.  \c colInds[k]
  ///   is valid if and only if it is a valid local mesh (i.e., block)
  ///   column index.  This method succeeded if and only if the return
  ///   value equals the input argument \c numColInds.
  LO
  sumIntoLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const;

  /// \brief Get a view of the (mesh, i.e., block) row, using local
  ///   (mesh, i.e., block) indices.
  ///
  /// This matrix has a graph, and we assume that the graph is fill
  /// complete on input to the matrix's constructor.  Thus, the matrix
  /// has a column Map, and it stores column indices as local indices.
  /// This means you can view the column indices as local indices
  /// directly.  However, you may <i>not</i> view them as global
  /// indices directly, since the column indices are not stored as
  /// global indices in the graph.
  ///
  /// \param localRowInd [in] Local (mesh, i.e., block) row index.
  ///
  /// \param colInds [out] If \c localRowInd is valid on the calling
  ///   process, then on output, this is a pointer to the local (mesh,
  ///   i.e., block) column indices in the given (mesh, i.e., block)
  ///   row.  If localRowInd is <i>not</i> valid, then this is
  ///   undefined.  (Please check the return value of this method.)
  ///
  /// \param vals [out] If \c localRowInd is valid on the calling
  ///   process, then on output, this is a pointer to the row's
  ///   values.  If localRowInd is <i>not</i> valid, then this is
  ///   undefined.  (Please check the return value of this method.)
  ///
  /// \param numInds [in] The number of (mesh, i.e., block) indices in
  ///   \c colInds on output.
  ///
  /// \return 0 if \c localRowInd is valid, else
  ///   <tt>Teuchos::OrdinalTraits<LO>::invalid()</tt>.
  LO
  getLocalRowView (const LO localRowInd,
                   const LO*& colInds,
                   Scalar*& vals,
                   LO& numInds) const;

  /// \brief Not implemented.
  void
  getLocalRowView (LO LocalRow,
                   Teuchos::ArrayView<const LO> &indices,
                   Teuchos::ArrayView<const Scalar> &values) const;

  /// \brief Not implemented.
  void
  getLocalRowCopy (LO LocalRow,
                   const Teuchos::ArrayView<LO> &Indices,
                   const Teuchos::ArrayView<Scalar> &Values,
                   size_t &NumEntries) const;

  little_block_type
  getLocalBlock (const LO localRowInd, const LO localColInd) const;

  /// \brief Get relative offsets corresponding to the given rows,
  ///   given by local row index.
  ///
  /// The point of this method is to precompute the results of
  /// searching for the offsets corresponding to the given column
  /// indices.  You may then reuse these search results in
  /// replaceLocalValuesByOffsets or sumIntoLocalValuesByOffsets.
  ///
  /// Offsets are block offsets; they are for column indices,
  /// <i>not</i> for values.
  ///
  /// \param localRowInd [in] Local index of the row.
  /// \param offsets [out] On output: relative offsets corresponding
  ///   to the given column indices.  Must have at least numColInds
  ///   entries.
  /// \param colInds [in] The local column indices for which to
  ///   compute offsets.  Must have at least numColInds entries.
  ///   This method will only read the first numColsInds entries.
  /// \param numColInds [in] Number of entries in colInds to read.
  ///
  /// \return The number of valid column indices in colInds.  This
  ///   method succeeded if and only if the return value equals the
  ///   input argument numColInds.
  LO
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[],
                      const LO colInds[],
                      const LO numColInds) const;

  /// \brief Like replaceLocalValues, but avoids computing row offsets.
  ///
  /// \return The number of valid column indices in colInds.  This
  ///   method succeeded if and only if the return value equals the
  ///   input argument numColInds.
  LO
  replaceLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const;

  /// \brief Like sumIntoLocalValues, but avoids computing row offsets.
  ///
  /// \return The number of valid column indices in colInds.  This
  ///   method succeeded if and only if the return value equals the
  ///   input argument numColInds.
  LO
  sumIntoLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const;

  /// \brief Return the number of entries in the given row on the
  ///   calling process.
  ///
  /// If the given local row index is invalid, this method (sensibly)
  /// returns zero, since the calling process trivially does not own
  /// any entries in that row.
  size_t getNumEntriesInLocalRow (const LO localRowInd) const;

  /// \brief Whether this object had an error on the calling process.
  ///
  /// Import and Export operations using this object as the target of
  /// the Import or Export may incur local errors, if some process
  /// encounters an LID in its list which is not a valid mesh row
  /// local index on that process.  In that case, we don't want to
  /// throw an exception, because not all processes may throw an
  /// exception; this can result in deadlock or put Tpetra in an
  /// incorrect state, due to lack of consistency across processes.
  /// Instead, we set a local error flag and ignore the incorrect
  /// data.  When unpacking, we do the same with invalid column
  /// indices.  If you want to check whether some process experienced
  /// an error, you must do a reduction or all-reduce over this flag.
  /// Every time you initiate a new Import or Export with this object
  /// as the target, we clear this flag.  (Note to developers: we
  /// clear it at the beginning of checkSizes().)
  bool localError () const {
    return *localError_;
  }

  /// \brief The current stream of error messages.
  ///
  /// This is only nonempty on the calling process if localError()
  /// returns true.  In that case, it stores a stream of
  /// human-readable, endline-separated error messages encountered
  /// during an Import or Export cycle.  Every time you initiate a new
  /// Import or Export with this object as the target, we clear this
  /// stream.  (Note to developers: we clear it at the beginning of
  /// checkSizes().)
  ///
  /// If you want to print this, you are responsible for ensuring that
  /// it is valid for the calling MPI process to print to whatever
  /// output stream you use.  On some MPI implementations, you may
  /// need to send the string to Process 0 for printing.
  std::string errorMessages () const {
    return (*errs_).is_null () ? std::string ("") : (*errs_)->str ();
  }

  /// \brief Get offsets of the diagonal entries in the matrix.
  ///
  /// \warning This method is only for expert users.
  /// \warning We make no promises about backwards compatibility
  ///   for this method.  It may disappear or change at any time.
  /// \warning This method must be called collectively.  We reserve
  ///   the right to do extra checking in a debug build that will
  ///   require collectives.
  ///
  /// \pre The matrix must be locally indexed (which means that it
  ///   has a column Map).
  /// \pre All diagonal entries of the matrix's graph must be
  ///   populated on this process.  Results are undefined otherwise.
  /// \post <tt>offsets.extent(0) == getNodeNumRows()</tt>
  ///
  /// This method creates an array of offsets of the local diagonal
  /// entries in the matrix.  This array is suitable for use in the
  /// two-argument version of getLocalDiagCopy().  However, its
  /// contents are not defined in any other context.  For example,
  /// you should not rely on \c offsets(i) being the index of the
  /// diagonal entry in the views returned by getLocalRowView().
  /// This may be the case, but it need not be.  (For example, we
  /// may choose to optimize the lookups down to the optimized
  /// storage level, in which case the offsets will be computed with
  /// respect to the underlying storage format, rather than with
  /// respect to the views.)
  ///
  /// If the matrix has a const ("static") graph, and if that graph
  /// is fill complete, then the offsets array remains valid through
  /// calls to fillComplete() and resumeFill().  "Invalidates" means
  /// that you must call this method again to recompute the offsets.
  void
  getLocalDiagOffsets (const Kokkos::View<size_t*, device_type,
                         Kokkos::MemoryUnmanaged>& offsets) const;

  /// \brief DEPRECATED overload of this method that writes offsets to
  ///   a Teuchos::ArrayRCP instead of a Kokkos::View.
  ///
  /// Please use the version of this method directly above, that
  /// writes offsets a Kokkos::View instead of to a Teuchos::ArrayRCP.
  void TPETRA_DEPRECATED
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const;

  /// \brief Variant of getLocalDiagCopy() that uses precomputed
  ///   offsets and puts diagonal blocks in a 3-D Kokkos::View.
  ///
  /// \param diag [out] On input: Must be preallocated, with
  ///   dimensions at least (number of diagonal blocks on the calling
  ///   process) x getBlockSize() x getBlockSize(). On output: the
  ///   diagonal blocks.  Leftmost index is "which block," then the
  ///   row index within a block, then the column index within a
  ///   block.
  ///
  /// This method uses the offsets of the diagonal entries, as
  /// precomputed by getLocalDiagOffsets(), to speed up copying the
  /// diagonal of the matrix.
  void
  getLocalDiagCopy (const Kokkos::View<impl_scalar_type***, device_type,
                                       Kokkos::MemoryUnmanaged>& diag,
                    const Kokkos::View<const size_t*, device_type,
                                       Kokkos::MemoryUnmanaged>& offsets) const;

  /// \brief Variant of getLocalDiagCopy() that uses precomputed
  ///   offsets and puts diagonal blocks in a 3-D Kokkos::View.
  ///
  /// \param diag [out] On input: Must be preallocated, with
  ///   dimensions at least (number of diagonal blocks on the calling
  ///   process) x getBlockSize() x getBlockSize(). On output: the
  ///   diagonal blocks.  Leftmost index is "which block," then the
  ///   row index within a block, then the column index within a
  ///   block.
  ///
  /// This method uses the offsets of the diagonal entries, as
  /// precomputed by getLocalDiagOffsets(), to speed up copying the
  /// diagonal of the matrix.
  void
  getLocalDiagCopy (const Kokkos::View<impl_scalar_type***, device_type,
                                       Kokkos::MemoryUnmanaged>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const;

  /// \brief Variant of getLocalDiagCopy() that uses precomputed offsets.
  ///
  /// \warning This overload of the method is DEPRECATED.  Call the
  ///   overload above that returns the diagonal blocks as a 3-D
  ///   Kokkos::View.
  ///
  /// This method uses the offsets of the diagonal entries, as
  /// precomputed by getLocalDiagOffsets(), to speed up copying the
  /// diagonal of the matrix.
  ///
  /// If the matrix has a const ("static") graph, and if that graph
  /// is fill complete, then the offsets array remains valid through
  /// calls to fillComplete() and resumeFill().
  void TPETRA_DEPRECATED
  getLocalDiagCopy (BlockCrsMatrix<Scalar,LO,GO,Node>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const;

protected:
  //! Like sumIntoLocalValues, but for the ABSMAX combine mode.
  LO
  absMaxLocalValues (const LO localRowInd,
                     const LO colInds[],
                     const Scalar vals[],
                     const LO numColInds) const;

  //! Like sumIntoLocalValuesByOffsets, but for the ABSMAX combine mode.
  LO
  absMaxLocalValuesByOffsets (const LO localRowInd,
                              const ptrdiff_t offsets[],
                              const Scalar vals[],
                              const LO numOffsets) const;

  /// \brief \name Implementation of Tpetra::DistObject.
  ///
  /// The methods here implement Tpetra::DistObject.  They let
  /// BlockMultiVector participate in Import and Export operations.
  /// Users don't have to worry about these methods.
  //@{

  virtual bool checkSizes (const ::Tpetra::SrcDistObject& source);

  virtual void
  copyAndPermute (const ::Tpetra::SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LO>& permuteToLIDs,
                  const Teuchos::ArrayView<const LO>& permuteFromLIDs);

  virtual void
  packAndPrepare (const ::Tpetra::SrcDistObject& source,
                  const Teuchos::ArrayView<const LO>& exportLIDs,
                  Teuchos::Array<packet_type>& exports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  ::Tpetra::Distributor& distor);

  virtual void
  unpackAndCombine (const Teuchos::ArrayView<const LO> &importLIDs,
                    const Teuchos::ArrayView<const packet_type> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    ::Tpetra::Distributor& distor,
                    ::Tpetra::CombineMode CM);
  //@}

private:
  //! The graph that describes the structure of this matrix.
  crs_graph_type graph_;
  Teuchos::RCP<crs_graph_type> graphRCP_;
  /// \brief The graph's row Map; the mesh row Map of this matrix.
  ///
  /// We keep this separately, not as an RCP, so that methods like
  /// replaceLocalValues and sumIntoLocalValues are thread safe.  (We
  /// could do this just by keeping the number of local indices in the
  /// mesh Map, but this more general approach will let us make
  /// replaceGlobalValues and sumIntoGlobalValues thread safe as
  /// well.)
  map_type rowMeshMap_;
  /// \brief The point Map version of the graph's domain Map.
  ///
  /// NOTE (mfh 16 May 2014) Since this is created at construction
  /// time, we don't have to worry about the lazy initialization
  /// issue, that we <i>do</i> have to worry about for X_colMap_ and
  /// Y_rowMap_.
  map_type domainPointMap_;
  /// \brief The point Map version of the graph's range Map.
  ///
  /// NOTE (mfh 16 May 2014) Since this is created at construction
  /// time, we don't have to worry about the lazy initialization
  /// issue, that we <i>do</i> have to worry about for X_colMap_ and
  /// Y_rowMap_.
  map_type rangePointMap_;
  //! The number of degrees of freedom per mesh point.
  LO blockSize_;

  /// \brief Host version of the graph's array of row offsets.
  ///
  /// The device version of this is already stored in the graph.  We
  /// need the host version here, because this class' interface needs
  /// to access it on host.  We don't want to assume UVM if we can
  /// avoid it.  Furthermore, we don't use Kokkos::DualView, because
  /// the graph is never supposed to change once a BlockCrsMatrix gets
  /// it.  The whole point of Kokkos::DualView is to help users
  /// synchronize changes between host and device versions of the
  /// data, but the data aren't changing in this case.  Furthermore,
  /// Kokkos::DualView has extra Views in it for the "modified" flags,
  /// and we don't want the (modest) overhead of creating and storing
  /// those.
  typename crs_graph_type::local_graph_type::row_map_type::HostMirror ptrHost_;

  /// \brief Host version of the graph's array of column indices.
  ///
  /// The device version of this is already stored in the graph.  We
  /// need the host version here, because this class' interface needs
  /// to access it on host.  See notes on ptrHost_ above.
  typename crs_graph_type::local_graph_type::entries_type::HostMirror indHost_;

  /// \brief The array of values in the matrix.
  ///
  /// Each blockSize_ x blockSize_ block of values is stored
  /// contiguously, in row major format, with no padding either inside
  /// a block or between blocks.
  typename Kokkos::DualView<impl_scalar_type*, device_type> val_;

  /// \brief Column Map block multivector (only initialized if needed).
  ///
  /// mfh 16 May 2014: This is a pointer to a pointer to BMV.  Ditto
  /// for Y_rowMap_ below.  This lets us do lazy initialization
  /// correctly with view semantics of BlockCrsMatrix.  All views of
  /// this BlockCrsMatrix have the same outer pointer.  That way, we
  /// can set the inner pointer in one view, and all other views will
  /// see it.  (Otherwise, other existing views of the BlockCrsMatrix
  /// would not get the benefit of lazy initialization.)
  ///
  /// The outer pointer is always nonull: It is always true that
  ///
  /// <tt>! X_colMap_.is_null()</tt>.
  ///
  /// However, the inner pointer starts out null and is lazily
  /// initialized in applyBlock().
  ///
  /// It's necessary to use a shared pointer (either Teuchos::RCP or
  /// std::shared_ptr) here, at least for the outer pointer type,
  /// because different views of the same block matrix will use the
  /// same BMV object.
  Teuchos::RCP<Teuchos::RCP<BMV> > X_colMap_;
  /// \brief Row Map block multivector (only initialized if needed).
  ///
  /// See the documentation of X_colMap_ above.
  Teuchos::RCP<Teuchos::RCP<BMV> > Y_rowMap_;

  /// \brief Use MV's importer instead of BMV's.
  ///
  /// Use MV's importer, which implements the "new" doTransfer interface,
  /// instead of BMV's importer, which implements the old. This probably entails
  /// inefficiency relative to a correct new-interface implementation in BMV;
  /// this workaround is a first step to getting a reasonably fast
  /// applyBlockNoTrans on the GPU.
  Teuchos::RCP<Teuchos::RCP<typename crs_graph_type::import_type> > pointImporter_;

  /// \brief Offset between blocks in the matrix.
  LO offsetPerBlock_;

  /// \brief Whether this object on the calling process is in an error state.
  ///
  /// See the documentation of localError() for details.
  ///
  /// The outer pointer is always nonnull.  Using a pointer rather
  /// than a \c bool value here ensures that all views of this object
  /// have access to the error state, because all views have the same
  /// (nonnull at construction) pointer.
  ///
  /// FIXME (mfh 16 Jun 2014) Use a Kokkos::View<bool, DeviceType>
  /// here, so that read access to localError_ will be thread safe.
  Teuchos::RCP<bool> localError_;

  /// \brief Stream of error messages.
  ///
  /// The outer pointer is always nonnull, but the inner pointer is
  /// only nonnull if localError_ is true.  Using a pointer to a
  /// pointer ensures that all views of this object have access to the
  /// error stream, because all views have the same (nonnull at
  /// construction) outer pointer.
  Teuchos::RCP<Teuchos::RCP<std::ostringstream> > errs_;

  //! Mark that a local error occurred, and get a stream for reporting it.
  std::ostream& markLocalErrorAndGetStream ();

  // //! Clear the local error state and stream.
  // void clearLocalErrorStateAndStream ();

public:
  //! \name Implementation of "dual view semantics"
  //@{

  //! Mark the matrix's values as modified in the given memory space.
  template<class MemorySpace>
  void modify ()
  {
    // It's legit to use a memory space, execution space, or
    // Kokkos::Device specialization as the template parameter of
    // Kokkos::DualView::modify.  That method just extracts the
    // memory_space typedef of its template parameter anyway.
    // However, insisting on a memory space avoids unnecessary
    // instantiations.
    val_.template modify<typename MemorySpace::memory_space> ();
  }

  //! Whether the matrix's values need sync'ing to the given memory space.
  template<class MemorySpace>
  bool need_sync () const
  {
    // It's legit to use a memory space, execution space, or
    // Kokkos::Device specialization as the template parameter of
    // Kokkos::DualView::need_sync.  That method just extracts the
    // memory_space typedef of its template parameter anyway.
    // However, insisting on a memory space avoids unnecessary
    // instantiations.
    return val_.template need_sync<typename MemorySpace::memory_space> ();
  }

  //! Sync the matrix's values <i>to</i> the given memory space.
  template<class MemorySpace>
  void sync ()
  {
    // It's legit to use a memory space, execution space, or
    // Kokkos::Device specialization as the template parameter of
    // Kokkos::DualView::sync.  That method just extracts the
    // memory_space typedef of its template parameter anyway.
    // However, insisting on a memory space avoids unnecessary
    // instantiations.
    val_.template sync<typename MemorySpace::memory_space> ();
  }

  /// \brief Get the host or device View of the matrix's values (\c val_).
  ///
  /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES OF
  ///   BACKWARDS COMPATIBILITY.  YOU HAVE BEEN WARNED.  CAVEAT
  ///   LECTOR!!!
  ///
  /// \tparam MemorySpace Memory space for which to get the
  ///   Kokkos::View.
  ///
  /// This is <i>not</i> const, because we reserve the right to do
  /// lazy allocation on device / "fast" memory.  Host / "slow" memory
  /// allocations generally are not lazy; that way, the host fill
  /// interface always works in a thread-parallel context without
  /// needing to synchronize on the allocation.
  template<class MemorySpace>
  auto getValues () -> decltype (val_.template view<typename MemorySpace::memory_space> ())
  {
    return val_.template view<typename MemorySpace::memory_space> ();
  }

  //@}

private:

  /// \brief Global sparse matrix-vector multiply for the transpose or
  ///   conjugate transpose cases.
  ///
  /// This method computes Y := beta*Y + alpha*Op(A)*X, where A is
  /// *this (the block matrix), Op(A) signifies either the transpose
  /// or the conjugate transpose of A, and X and Y are block
  /// multivectors.  The special cases alpha = 0 resp. beta = 0 have
  /// their usual BLAS meaning; this only matters if (A or X) resp. Y
  /// contain Inf or NaN values.
  void
  applyBlockTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                   BlockMultiVector<Scalar, LO, GO, Node>& Y,
                   const Teuchos::ETransp mode,
                   const Scalar alpha,
                   const Scalar beta);

  /// \brief Global sparse matrix-vector multiply for the non-transpose case.
  ///
  /// This method computes Y := beta*Y + alpha*A*X, where A is *this
  /// (the block matrix), and X and Y are block multivectors.  The
  /// special cases alpha = 0 resp. beta = 0 have their usual BLAS
  /// meaning; this only matters if (A or X) resp. Y contain Inf or
  /// NaN values.
  void
  applyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                     BlockMultiVector<Scalar, LO, GO, Node>& Y,
                     const Scalar alpha,
                     const Scalar beta);

  /// \brief Local sparse matrix-vector multiply for the non-transpose case.
  ///
  /// This method computes Y := beta*Y + alpha*A*X, where A is *this
  /// (the block matrix), and X and Y are block multivectors.  The
  /// special cases alpha = 0 resp. beta = 0 have their usual BLAS
  /// meaning; this only matters if (A or X) resp. Y contain Inf or
  /// NaN values.
  void
  localApplyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                          BlockMultiVector<Scalar, LO, GO, Node>& Y,
                          const Scalar alpha,
                          const Scalar beta);

  /// \brief Get the relative block offset of the given block.
  ///
  /// \param localRowIndex [in] Local index of the entry's row.
  /// \param colIndexToFind [in] Local index of the entry's column.
  /// \param hint [in] Relative offset hint (MUST be nonnegative).
  ///
  /// <i>Relative</i> offsets are relative to the current row, while
  /// <i>absolute</i> offsets are just direct indices into an array.
  /// For example, if <tt>k_abs</tt> is an absolute offset into the
  /// array of column indices <tt>indHost_</tt>, then one can use
  /// <tt>k_abs</tt> directly as <tt>indHost_[k_abs]</tt>.  If
  /// <tt>k_rel</tt> is a relative offset into <tt>indHost_</tt>, then
  /// one must know the current local row index in order to use
  /// <tt>k_rel</tt>.  For example:
  /// \code
  /// size_t k_abs = ptrHost_[curLocalRow] + k_rel; // absolute offset
  /// LO colInd = indHost_[k_abs];
  /// \endcode
  ///
  /// This method returns a relative block offset.  A <i>block</i>
  /// offset means a graph or mesh offset, vs. the <i>point</i> offset
  /// into the array of values \c val_.  A block offset is suitable
  /// for use in \c indHost_, but not in \c val_.  One must multiply a
  /// block offset by offsetPerBlock() in order to get the
  /// <i>point</i> offset into \c val_.
  ///
  /// The given "hint" is a relative block offset.  It can help avoid
  /// searches, for the common case of accessing several consecutive
  /// entries in the same row.
  ///
  /// Relative offsets may have type \c LO, since a row may not have
  /// more entries than the number of columns in the graph.  Absolute
  /// offsets <i>must</i> have type \c size_t or \c ptrdiff_t, since
  /// the total number of graph or matrix entries on a process may
  /// exceed the total number of rows and columns.
  ///
  /// \return Teuchos::OrdinalTraits<size_t>::invalid() if there is no
  ///   block at the given index pair; otherwise, the "absolute"
  ///   block offset of that block.
  LO
  findRelOffsetOfColumnIndex (const LO localRowIndex,
                              const LO colIndexToFind,
                              const LO hint = 0) const;

  /// \brief Number of entries consumed by each block in the matrix,
  ///   including padding; the stride between blocks.
  LO offsetPerBlock () const;

  const_little_block_type
  getConstLocalBlockFromInput (const impl_scalar_type* val, const size_t pointOffset) const;

  little_block_type
  getNonConstLocalBlockFromInput (impl_scalar_type* val, const size_t pointOffset) const;

  const_little_block_type
  getConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const;

  little_block_type
  getNonConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const;

  /// \c Block at the given local mesh row and relative (mesh) offset.
  ///
  /// Use this for 2-argument getLocalDiagCopy that writes to Kokkos::View.
  const_little_block_type
  getConstLocalBlockFromRelOffset (const LO lclMeshRow,
                                   const size_t relMeshOffset) const;

public:
  //! The communicator over which this matrix is distributed.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The Kokkos Node instance.
  virtual Teuchos::RCP<Node> getNode() const;

  //! The global number of columns of this matrix.
  virtual global_size_t getGlobalNumCols() const;

  virtual size_t getNodeNumCols() const;

  virtual GO getIndexBase() const;

  //! The global number of stored (structurally nonzero) entries.
  virtual global_size_t getGlobalNumEntries() const;

  //! The local number of stored (structurally nonzero) entries.
  virtual size_t getNodeNumEntries() const;

  /// \brief The current number of entries on the calling process in the specified global row.
  ///
  /// Note that if the row Map is overlapping, then the calling
  /// process might not necessarily store all the entries in the
  /// row.  Some other process might have the rest of the entries.
  ///
  /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
  ///   the specified global row does not belong to this graph, else
  ///   the number of entries.
  virtual size_t getNumEntriesInGlobalRow (GO globalRow) const;

  /// \brief Number of diagonal entries in the matrix's graph, over
  ///   all processes in the matrix's communicator.
  ///
  /// \pre <tt>! this->isFillActive()</tt>
  ///
  /// \warning This method is DEPRECATED.  DO NOT CALL IT.  It may
  ///   go away at any time.
  virtual global_size_t TPETRA_DEPRECATED getGlobalNumDiags() const;

  /// \brief Number of diagonal entries in the matrix's graph, on the
  ///   calling process.
  ///
  /// \pre <tt>! this->isFillActive()</tt>
  ///
  /// \warning This method is DEPRECATED.  DO NOT CALL IT.  It may
  ///   go away at any time.
  virtual size_t TPETRA_DEPRECATED getNodeNumDiags() const;

  //! The maximum number of entries across all rows/columns on all nodes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! Whether this matrix has a well-defined column map.
  virtual bool hasColMap() const;

  /// \brief Whether the matrix's graph is locally lower triangular.
  ///
  /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
  ///   and will DISAPPEAR VERY SOON per #2630.
  ///
  /// \pre <tt>! isFillActive()</tt>.
  ///   If fill is active, this method's behavior is undefined.
  ///
  /// \note This is entirely a local property.  That means this
  ///   method may return different results on different processes.
  virtual bool TPETRA_DEPRECATED isLowerTriangular () const;

  /// \brief Whether the matrix's graph is locally upper triangular.
  ///
  /// \warning DO NOT CALL THIS METHOD!  This method is DEPRECATED
  ///   and will DISAPPEAR VERY SOON per #2630.
  ///
  /// \pre <tt>! isFillActive()</tt>.
  ///   If fill is active, this method's behavior is undefined.
  ///
  /// \note This is entirely a local property.  That means this
  ///   method may return different results on different processes.
  virtual bool TPETRA_DEPRECATED isUpperTriangular () const;

  /// \brief Whether matrix indices are locally indexed.
  ///
  /// A RowMatrix may store column indices either as global indices
  /// (of type <tt>GO</tt>), or as local indices (of type
  /// <tt>LO</tt>).  In some cases (for example, if the
  /// column Map has not been computed), it is not possible to
  /// switch from global to local indices without extra work.
  /// Furthermore, some operations only work for one or the other
  /// case.
  virtual bool isLocallyIndexed() const;

  /// \brief Whether matrix indices are globally indexed.
  ///
  /// A RowMatrix may store column indices either as global indices
  /// (of type <tt>GO</tt>), or as local indices (of type
  /// <tt>LO</tt>).  In some cases (for example, if the
  /// column Map has not been computed), it is not possible to
  /// switch from global to local indices without extra work.
  /// Furthermore, some operations only work for one or the other
  /// case.
  virtual bool isGloballyIndexed() const;

  //! Whether fillComplete() has been called.
  virtual bool isFillComplete() const;

  //! Whether this object implements getLocalRowView() and getGlobalRowView().
  virtual bool supportsRowViews() const;

  //@}
  //! @name Extraction Methods
  //@{

  /// \brief Get a copy of the given global row's entries.
  ///
  /// This method only gets the entries in the given row that are
  /// stored on the calling process.  Note that if the matrix has an
  /// overlapping row Map, it is possible that the calling process
  /// does not store all the entries in that row.
  ///
  /// \param GlobalRow [in] Global index of the row.
  /// \param Indices [out] Global indices of the columns
  ///   corresponding to values.
  /// \param Values [out] Matrix values.
  /// \param NumEntries [out] Number of stored entries on the
  ///   calling process; length of Indices and Values.
  ///
  /// This method throws <tt>std::runtime_error</tt> if either
  /// Indices or Values is not large enough to hold the data
  /// associated with row GlobalRow. If GlobalRow does not belong to
  /// the calling process, then the method sets NumIndices to
  /// <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>, and does
  /// not modify Indices or Values.
  virtual void
  getGlobalRowCopy (GO GlobalRow,
                    const Teuchos::ArrayView<GO> &Indices,
                    const Teuchos::ArrayView<Scalar> &Values,
                    size_t& NumEntries) const;

  /// \brief Get a constant, nonpersisting, globally indexed view of
  ///   the given row of the matrix.
  ///
  /// The returned views of the column indices and values are not
  /// guaranteed to persist beyond the lifetime of <tt>this</tt>.
  /// Furthermore, some RowMatrix implementations allow changing the
  /// values, or the indices and values.  Any such changes
  /// invalidate the returned views.
  ///
  /// This method only gets the entries in the given row that are
  /// stored on the calling process.  Note that if the matrix has an
  /// overlapping row Map, it is possible that the calling process
  /// does not store all the entries in that row.
  ///
  /// \pre <tt>isGloballyIndexed () && supportsRowViews ()</tt>
  /// \post <tt>indices.size () == getNumEntriesInGlobalRow (GlobalRow)</tt>
  ///
  /// \param GlobalRow [in] Global index of the row.
  /// \param Indices [out] Global indices of the columns
  ///   corresponding to values.
  /// \param Values [out] Matrix values.
  ///
  /// If \c GlobalRow does not belong to this node, then \c indices
  /// is set to \c null.
  virtual void
  getGlobalRowView (GO GlobalRow,
                    Teuchos::ArrayView<const GO>& indices,
                    Teuchos::ArrayView<const Scalar>& values) const;

  /// \brief Get a copy of the diagonal entries, distributed by the row Map.
  ///
  /// On input, the Vector's Map must be the same as the row Map of the matrix.
  /// (That is, <tt>this->getRowMap ()->isSameAs (* (diag.getMap ())) == true</tt>.)
  ///
  /// On return, the entries of \c diag are filled with the diagonal
  /// entries of the matrix stored on this process.  Note that if
  /// the row Map is overlapping, multiple processes may own the
  /// same diagonal element.  You may combine these overlapping
  /// diagonal elements by doing an Export from the row Map Vector
  /// to a range Map Vector.
  virtual void getLocalDiagCopy (::Tpetra::Vector<Scalar,LO,GO,Node>& diag) const;

  //@}
  //! \name Mathematical methods
  //@{

  /**
   * \brief Scale the RowMatrix on the left with the given Vector x.
   *
   * On return, for all entries i,j in the matrix, \f$A(i,j) = x(i)*A(i,j)\f$.
   */
  virtual void leftScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& x);

  /**
   * \brief Scale the RowMatrix on the right with the given Vector x.
   *
   * On return, for all entries i,j in the matrix, \f$A(i,j) = x(j)*A(i,j)\f$.
   */
  virtual void rightScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& x);

  /// \brief The Frobenius norm of the matrix.
  ///
  /// This method computes and returns the Frobenius norm of the
  /// matrix.  The Frobenius norm \f$\|A\|_F\f$ for the matrix
  /// \f$A\f$ is defined as
  /// \f$\|A\|_F = \sqrt{ \sum_{i,j} |A(i,j)|^2 }\f$.
  /// It has the same value as the Euclidean norm of a vector made
  /// by stacking the columns of \f$A\f$.
  virtual typename ::Tpetra::RowMatrix<Scalar, LO, GO, Node>::mag_type
  getFrobeniusNorm () const;
  //@}
};

} // namespace Classes
} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DECL_HPP
