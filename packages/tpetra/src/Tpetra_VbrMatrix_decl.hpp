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

#ifndef TPETRA_VBRMATRIX_DECL_HPP
#define TPETRA_VBRMATRIX_DECL_HPP

/// \file Tpetra_VbrMatrix_decl.hpp
/// \brief Declarations for the class Tpetra::VbrMatrix.

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_VbrUtils.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

// Forward declarations, to avoid including files that we don't really
// need just for declaring VbrMatrix.  The #ifndef ... #endif prevents
// Doxygen from generating spurious documentation for the forward
// declarations.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class BlockMap;

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class BlockCrsGraph;

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class Vector;
} // namespace Tpetra

namespace Kokkos {
  template<class Scalar, class LocalOrdinal, class Node>
  class VbrMatrix;
} // namespace Kokkos

namespace Teuchos {
  template<class OrdinalType, class ScalarType>
  class SerialDenseMatrix;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

//! \brief VbrMatrix: Variable block row matrix.
/**
The VbrMatrix class has two significant 'states', distinguished by whether or not
storage has been optimized (packed) or not.

When the matrix is in the non-optimized-storage state, internal data
storage is in a non-contiguous data-structure that allows for
convenient insertion of data.

When the matrix is in the optimized-storage state, internal data is stored in
contiguous (packed) arrays. When in this state, existing entries may be updated
and replaced, but no new entries (indices and/or coefficients) may be inserted.
In other words, the sparsity pattern or structure of the matrix may not be
changed.

Use of the matrix as an Operator (performing matrix-vector multiplication) is
only allowed when it is in the optimized-storage state.

VbrMatrix has two constructors, one which leaves the matrix in the optimized-
storage state, and another which leaves the matrix in the non-optimized-storage
state.

When the VbrMatrix is constructed in the non-optimized-storage state, (and then
filled using methods such as setGlobalBlockEntry etc.), it can then be transformed
to the optimized-storage state by calling the method fillComplete().

Once in the optimized-storage state, the VbrMatrix can not be returned to the
non-optimized-storage state.
*/
template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Kokkos::DefaultNode::DefaultNodeType,
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::BlockSparseOps >
class VbrMatrix :
    public Tpetra::DistObject<char, LocalOrdinal, GlobalOrdinal, Node>,
    public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
public:
  typedef Scalar        scalar_type;
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;
  typedef LocalMatOps   mat_vec_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying the row-map and the max number of (block) non-zeros for all rows.
  /*! After this constructor completes, the VbrMatrix is in the non-packed,
    non-optimized-storage, isFillComplete()==false state.
    Block-entries (rectangular, dense submatrices) may be inserted using class
    methods such as setGlobalBlockEntry(...), declared below.
  */
  VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! Constructor specifying a pre-filled block-graph.
  /*! Constructing a VbrMatrix with a pre-filled graph means that the matrix will
      start out in the optimized-storage state, i.e., isFillComplete()==true.
      The graph provided to this constructor must be already filled.
      (If blkGraph->isFillComplete() == false, an exception is thrown.)

      Entries in the input BlockCrsGraph correspond to block-entries in the
      VbrMatrix. In other words, the VbrMatrix will have a block-row corresponding
      to each row in the graph, and a block-entry corresponding to each column-
      index in the graph.
  */
  VbrMatrix(const Teuchos::RCP<const BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& blkGraph);

  //! Destructor
  virtual ~VbrMatrix();

  //@}

  //! @name Advanced Mathematical operations
  //@{

  //! Multiply this matrix by a MultiVector.
  /*! \c X is required to be post-imported, i.e., described by the column map
      of the matrix. \c Y is required to be pre-exported, i.e., described by
      the row map of the matrix.
      See also the Operator::apply method which is implemented below.
  */
  template <class DomainScalar, class RangeScalar>
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

  //! Triangular Solve -- Matrix must be triangular.
  /*! Find X such that A*X = Y.
      \c X is required to be post-imported, i.e., described by the column map
      of the matrix. \c Y is required to be pre-exported, i.e., described by
      the row map of the matrix.

      Both \c X and \c Y are required to have constant stride.

      Note that if the diagonal block-entries are stored, they must be triangular.
      I.e., the matrix structure must be block-triangular, and any diagonal blocks
      must be "point"-triangular, meaning that coefficients on the "wrong" side of the
      point-diagonal must be zero.
  */
  template <class DomainScalar, class RangeScalar>
      void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;

  //@}

  //! @name Operator Methods
  //@{

  //! Returns the (point-entry) Map associated with the domain of this operator.
  /*! Note that this is a point-entry map, not a block-map.
  */
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

  //! Returns the (point-entry) Map associated with the range of this operator.
  /*! Note that this is a point-entry map, not a block-map.
  */
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{trans}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
   */
  void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
             MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
             Teuchos::ETransp trans = Teuchos::NO_TRANS,
             Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Triangular Solve -- Matrix must be triangular.
  /*! Find X such that A*X = Y.
      Both \c X and \c Y are required to have constant stride.
  */
  void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                    Teuchos::ETransp trans) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  /*!
    VbrMatrix does support transpose-apply. (This method returns true.)
  */
  bool hasTransposeApply() const;

  //@}

  //! @name Attribute Query Methods
  //@{

  //! Returns the block-row map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const;

  //! Returns the block-column map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const;

  //! Returns the block-domain map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const;

  //! Returns the block-range map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const;

  //! Returns the point-row map.
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointRowMap() const;

  //! Returns the point-column map.
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getPointColMap() const;

  //! Return true if fillComplete has been called, false otherwise.
  bool isFillComplete() const;
  //@}

  //! @name Insertion Methods
  //@{

  //! Set the specified scalar throughout the matrix.
  /*!
    This method may be called any time (before or after fillComplete()).
  */
  void putScalar(Scalar s);

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, it will be
    over-written (replaced) by the input block-entry.

    Note that if globalBlockRow is not owned by the local processor (as
    indicated by getBlockRowMap()) then the block-entry is held in
    temporary storage until globalAssemble() is called (which is called
    internally by fillComplete()) and then globalAssemble() performans
    the communication needed to move the data to the owning processor.

    This method may be called any time (before or after fillComplete()).
  */
  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<int,Scalar>& blockEntry);

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The coefficients of the specified block-entry will be
    over-written (replaced) by the input block-entry.
  */
  void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<int,Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    Note that if globalBlockRow is not owned by the local processor (as
    indicated by getBlockRowMap()) then the block-entry is held in
    temporary storage until globalAssemble() is called (which is called
    internally by fillComplete()) and then globalAssemble() performans
    the communication needed to move the data to the owning processor.

    This method may be called any time (before or after fillComplete()).
  */
  void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<int,Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The contents of the input block-entry will be added to the values that are
    already present in the matrix.
  */
  void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<int,Scalar>& blockEntry);

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, it will be
    over-written (replaced) by the input block-entry.

    Note that if globalBlockRow is not owned by the local processor (as
    indicated by getBlockRowMap()) then the block-entry is held in
    temporary storage until globalAssemble() is called (which is called
    internally by fillComplete()) and then globalAssemble() performans
    the communication needed to move the data to the owning processor.

    This method may be called any time (before or after fillComplete()).
  */
  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The coefficients for the specified block-entry will be
    over-written (replaced) by the input block-entry.
  */
  void setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    Note that if globalBlockRow is not owned by the local processor (as
    indicated by getBlockRowMap()) then the block-entry is held in
    temporary storage until globalAssemble() is called (which is called
    internally by fillComplete()) and then globalAssemble() performans
    the communication needed to move the data to the owning processor.

    This method may be called any time (before or after fillComplete()).
  */
  void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will throw an exception if fillComplete() has not yet been called,
    or if the specified block-entry doesn't already exist in the matrix.

    The contents of the input block-entry will be added to the values that are
    already present in the matrix.
  */
  void sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //@}

  //! @name Transformational Methods
  //@{

  //! Transition the matrix to the packed, optimized-storage state.
  /*!
    This method also sets the domain and range maps.
    This method internally calls globalAssemble().
  */
  void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap);

  //! Transition the matrix to the packed, optimized-storage state.
  /*!
    This method internally calls fillComplete(getBlockRowMap(),getBlockRowMap()).
  */
  void fillComplete();

  //! Communicate non-local contributions to the processors that own those contributions.
  void globalAssemble();
  //@}

  //! @name Extraction Methods
  //@{

  //! Returns a const read-only view of the block-entries in the specified row.
  /*! Can only be called if isFillComplete()==false
  */
  void getGlobalBlockRowView(GlobalOrdinal globalBlockRow,
                             LocalOrdinal& numPtRows,
                             Teuchos::ArrayView<const GlobalOrdinal>& blockCols,
                             Teuchos::Array<LocalOrdinal>& ptColsPerBlockCol,
                             Teuchos::Array<Teuchos::ArrayRCP<const Scalar> >& blockEntries) const;

  //! Returns a const read-only view of the block-entries in the specified row.
  /*! Can only be called if isFillComplete()==true
  */
  void getLocalBlockRowView(LocalOrdinal localBlockRow,
                            LocalOrdinal& numPtRows,
                            Teuchos::ArrayView<const LocalOrdinal>& blockCols,
                            Teuchos::Array<LocalOrdinal>& ptColsPerBlockCol,
                            Teuchos::ArrayRCP<const Scalar>& blockEntries) const;

  //! Returns a const read-only view of a block-entry.
  /*!
    The arguments numPtRows and numPtCols are set to the dimensions of the block-
    entry on output.
    The stride (LDA in Blas terminology) is equal to numPtRows.

    This method may be called any time (before or after fillComplete()), but will
    throw an exception if the specified block-entry doesn't already exist.
  */
  void getGlobalBlockEntryView(GlobalOrdinal globalBlockRow,
                               GlobalOrdinal globalBlockCol,
                               LocalOrdinal& numPtRows,
                               LocalOrdinal& numPtCols,
                               Teuchos::ArrayRCP<const Scalar>& blockEntry) const;

  //! Returns a non-const read-write view of a block-entry.
  /*! Creates the block-entry if it doesn't already exist, and if:
     - the arguments numPtRows and numPtCols are set on entry (and nonzero),
     - and if fillComplete() has not yet been called.

     Important Note: Be very careful managing the lifetime of this view.
       If fillComplete() has been called, and if you are running on a GPU,
       this view may be a copy of memory from the GPU, and your changes to the
       view won't be copied back to the GPU until your ArrayRCP is destroyed
       or set to Teuchos::null.
  */
  void getGlobalBlockEntryViewNonConst(GlobalOrdinal globalBlockRow,
                                       GlobalOrdinal globalBlockCol,
                                       LocalOrdinal& numPtRows,
                                       LocalOrdinal& numPtCols,
                                       Teuchos::ArrayRCP<Scalar>& blockEntry);

  //! Returns a const read-only view of a block-entry.
  /*!
    The arguments numPtRows and numPtCols are set to the dimensions of the block-
    entry on output.
    The stride (LDA in Blas terminology) is equal to numPtRows.
    Throws an exception if fillComplete() has not yet been called, or if the
    specified block-entry doesn't exist.

    This method may only be called after fillComplete() has been called, and will
    throw an exception if the specified block-entry doesn't already exist.
  */
  void getLocalBlockEntryView(LocalOrdinal localBlockRow,
                              LocalOrdinal localBlockCol,
                              LocalOrdinal& numPtRows,
                              LocalOrdinal& numPtCols,
                              Teuchos::ArrayRCP<const Scalar>& blockEntry) const;

  //! Returns a non-const read-write view of a block-entry.
  /*!
    The arguments numPtRows and numPtCols are set to the dimensions of the block-
    entry on output.
    The stride (LDA in Blas terminology) is equal to numPtRows.
    Throws an exception if fillComplete() has not yet been called, or if the
    specified block-entry doesn't exist.

     Important Note: Be very careful managing the lifetime of this view.
       If fillComplete() has been called, and if you are running on a GPU,
       this view may be a copy of memory from the GPU, and your changes to the
       view won't be copied back to the GPU until your ArrayRCP is destroyed
       or set to Teuchos::null.

    This method may only be called after fillComplete() has been called, and will
    throw an exception if the specified block-entry doesn't already exist.
  */
  void getLocalBlockEntryViewNonConst(LocalOrdinal localBlockRow,
                                      LocalOrdinal localBlockCol,
                                      LocalOrdinal& numPtRows,
                                      LocalOrdinal& numPtCols,
                                      Teuchos::ArrayRCP<Scalar>& blockEntry);

  //! Return a copy of the (point-entry) diagonal values.
  /*!
    Throws an exception if the input-vector's map is not the same as
    getBlockRowMap()->getPointMap().
  */
  void getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& diag) const;

  const Teuchos::RCP<const BlockCrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& getBlockCrsGraph() {return constBlkGraph_;}
  //@}

  //! @name Overridden from Teuchos::DistObject
  //@{

  bool checkSizes(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source);

  void copyAndPermute(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source, size_t numSameIDs, const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs, const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs);

  void packAndPrepare(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source, const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs, Teuchos::Array<char>& exports, const Teuchos::ArrayView<size_t>& numPacketsPerLID, size_t& constantNumPackets, Distributor& distor);

  void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal>& importLIDs, const Teuchos::ArrayView<const char>& imports, const Teuchos::ArrayView<size_t>& numPacketsPerLID, size_t constantNumPackets, Distributor& distor, CombineMode CM);

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{
  std::string description() const;

  /** \brief Print the object with some verbosity level to a FancyOStream object.
  */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
  //@}

 private:
  //private methods:

  Teuchos::RCP<Node> getNode() const;

  void updateImport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const;
  void updateExport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void createImporterExporter();
  void optimizeStorage();
  void fillLocalMatrix();
  void fillLocalMatVec();

  //private data members:

  //We hold two graph pointers, one const and the other non-const.
  //If a BlockCrsGraph is provided at construction, it is const and VbrMatrix
  //never changes it.
  //If a BlockCrsGraph is not provided at construction, VbrMatrix creates one
  //internally and fills it as the matrix is filled, up until fillComplete()
  //is called.
  //
  //blkGraph_ is either the internally created graph, or is null.
  //constBlkGraph_ is either a pointer to blkGraph_, or a pointer to the
  //graph provided at construction time.
  //
  Teuchos::RCP<BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node> > blkGraph_;
  Teuchos::RCP<const BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node> > constBlkGraph_;

  // We use a pointer so that the _decl header file only needs a
  // forward declaration, not an include.  We don't really need an RCP
  // here, since the local matrix object never gets shared outside the
  // class.  std::unique_ptr (C++11) would be more appropriate.
  Teuchos::RCP<Kokkos::VbrMatrix<Scalar,LocalOrdinal,Node> > lclMatrix_;

  //A variable-block-row matrix is represented by 6 arrays
  //in packed (contiguous storage) form. For a description of these
  //arrays, see the text at the bottom of this file.
  //(2 of those arrays, rptr and cptr, are represented by arrays in the
  //getBlockRowMap() and getBlockColMap() objects, and
  //another two of those arrays, bptr and bindx, are represented by arrays in
  //the BlockCrsGraph object.)
  //This is noted in the comments for rptr,cptr,bptr,bindx below.
  //
  //These arrays are handled as if they may point to memory that resides on
  //a separate device (e.g., a GPU). In other words, when the contents of these
  //arrays are manipulated, we use views or buffers obtained from the Node
  //object.
  Teuchos::ArrayRCP<Scalar> pbuf_values1D_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_indx_;

  LocalMatOps lclMatOps_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;
  Teuchos::RCP<Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node> > exporter_;
  mutable Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importedVec_;
  mutable Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > exportedVec_;

  typedef typename std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > RowGlobalCols;

  //We use an array-of-maps to represent the variable-block-row matrix in
  //un-packed '2D' form.
  //
  //This unpacked data is assumed to be resident in CPU (host) memory.
  //It doesn't make sense to copy this data back and forth to a separate
  //compute device (e.g., a GPU), since we don't support doing matrix-vector
  //products until after fillComplete is called, at which time contiguous
  //arrays are allocated on the device and matrix data is copied into them.
  Teuchos::RCP<Teuchos::Array<RowGlobalCols> > data_2D_;

  VbrUtils::VbrData<LocalOrdinal,GlobalOrdinal,Scalar> nonlocal_data_;

  bool is_fill_completed_;
  bool is_storage_optimized_;
};//class VbrMatrix

}//namespace Tpetra

//----------------------------------------------------------------------------
// Description of arrays representing the VBR format:
//
// (For more description, see this URL (valid as of 5/26/2010):
// http://docs.sun.com/source/817-0086-10/prog-sparse-support.html)
// ...and of course more can be found using google...
// The old Aztec manual was a great resource for this but I can't
// find a copy of that these days...
//
//
// Here is a brief description of the 6 arrays that are required to
// represent a VBR matrix in packed (contiguous-memory-storage) format:
//
// rptr: length num_block_rows + 1
//       rptr[i]: the pt-row corresponding to the i-th block-row
//       Note: rptr is getBlockRowMap()->getNodeFirstPointInBlocks().
//
// cptr: length num_distinct_block_cols + 1
//       cptr[j]: the pt-col corresponding to the j-th block-col
//       Note: cptr is getBlockColMap()->getNodeFirstPointInBlocks().
//
// bptr: length num_block_rows + 1
//       bptr[i]: location in bindx of the first nonzero block-entry
//                of the i-th block-row
//       Note: bptr is blkGraph_->getNodeRowOffsets();
//
// bindx: length num-nonzero-block-entries
//        bindx[j]: block-col-index of j-th block-entry
//        Note: bindx is blkGraph_->getNodePackedIndices();
//
// indx: length num-nonzero-block-entries + 1
//       indx[j]: location in vals of the beginning of the j-th
//       block-entry
//
// vals: length num-nonzero-scalar-entries
//
//
// Some example look-ups:
//
// int nbr = num_block_rows;
// int total_num_block_nonzeros = bptr[nbr];
// int total_num_scalar_nonzeros = indx[num_block_nonzeros];
//
// //get arrays for i-th block-row:
// int* bindx_i = &bindx[bptr[i]];
// double* vals_i = &val[indx[bptr[i]]];
// int num_block_nonzeros_in_row_i = bptr[i+1]-bptr[i];
//
//----------------------------------------------------------------------------

#endif //TPETRA_VBRMATRIX_DECL_HPP

