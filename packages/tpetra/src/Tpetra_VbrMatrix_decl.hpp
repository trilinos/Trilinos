//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_VBRMATRIX_DECL_HPP
#define TPETRA_VBRMATRIX_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_VbrMatrix.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerializationTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_CrsGraph.hpp"

/** \file Tpetra_VbrMatrix_decl.hpp

  Declarations for the class Tpetra::VbrMatrix.
*/
namespace Tpetra {

//! \brief VbrMatrix: Variable block row matrix.
/**
The VbrMatrix class has two significant 'states', distinguished by whether or not
storage has been optimized (packed) or not.
After construction and before fillComplete() has been called, the internal data
storage is in an un-packed, non-contiguous data-structure that allows for
convenient insertion of data.
After fillComplete() has been called, internal storage is in contiguous 'packed'
arrays and no new entries (indices and/or coefficients) may be inserted. Existing
entries may still be updated and replaced though. The structure or sparsity
pattern of the matrix is finalized when fillComplete() is called.
*/
template <class Scalar, 
          class LocalOrdinal  = int, 
          class GlobalOrdinal = LocalOrdinal, 
          class Node          = Kokkos::DefaultNode::DefaultNodeType, 
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::BlockSparseOps >
class VbrMatrix : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
 public:
  typedef Scalar        scalar_type;
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;
  typedef LocalMatOps   mat_vec_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying the row-map and the max number of (block) non-zeros for all rows.
  VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! Destructor
  virtual ~VbrMatrix();

  //@}

  //! @name Advanced Matrix-vector multiplication method

  //! Multiply this matrix by a MultiVector.
  /*! \c X is required to be post-imported, i.e., described by the column map
      of the matrix. \c Y is required to be pre-exported, i.e., described by
      the row map of the matrix.
      See also the Operator::apply method which is implemented below.
  */
  template <class DomainScalar, class RangeScalar>
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

  //@}

  //! @name Operator Methods
  //@{

    //! Returns the Map associated with the domain of this operator, which must be compatible with X.getMap().
    /*! Note that this is a point-entry map, not a block-map.
    */
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
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

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

  //@}

  //! @name Attribute Query Methods
  //@{

  //! Returns the block-row map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const;

  //! Returns the block-column map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const;

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

    This method may be called any time (before or after fillComplete()).
  */
  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    This method may be called any time (before or after fillComplete()).
  */
  void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry);

  //!Copy the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, it will be
    over-written (replaced) by the input block-entry.

    This method may be called any time (before or after fillComplete()).
  */
  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //!Add the contents of the input block-entry into the matrix.
  /*!
    This method will create the specified block-entry if it doesn't already exist,
    but only if fillComplete() has not yet been called.

    If the specified block-entry already exists in the matrix, the contents of the
    input block-entry will be added to the values that are already present.

    This method may be called any time (before or after fillComplete()).
  */
  void sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //@}

  //! @name Transformational Methods
  //@{

  void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap);

  void fillComplete();
  //@}

  //! @name Extraction Methods
  //@{

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

  //@}

 private:
  //private methods:

  Teuchos::RCP<Node> getNode() const;

  void updateImport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const;
  void updateExport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void optimizeStorage();
  void fillLocalMatrix();
  void fillLocalMatVec();

  //private data members:

  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRowMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkColMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkDomainMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRangeMap_;

  Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > blkGraph_;
  Kokkos::VbrMatrix<Scalar,LocalOrdinal,Node> lclMatrix_;

  //It takes 6 arrays to adequately represent a variable-block-row
  //matrix in packed (contiguous storage) form. For a description of these
  //arrays, see the text at the bottom of this file.
  //
  //These arrays are handled as if they may point to memory that resides on
  //a separate device (e.g., a GPU). In other words, when the contents of these
  //arrays are manipulated, we use views or buffers obtained from the Node object.
  Teuchos::ArrayRCP<Scalar> pbuf_values1D_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_rptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_cptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_bptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_bindx_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_indx_;

  LocalMatOps lclMatVec_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;
  Teuchos::RCP<Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node> > exporter_;
  mutable Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importedVec_;
  mutable Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > exportedVec_;

  typedef typename std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > MapGlobalArrayRCP;
  typedef typename std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > MapLocalArrayRCP;

  //We use 2 arrays (well, array-of-maps, array-of-array-of-arrays...) to
  //represent the variable-block-row matrix in un-packed '2D' form.
  //
  //Note that these arrays are assumed to be resident in CPU (host) memory.
  //It doesn't make sense to copy this kind of data back and forth to a separate
  //compute device (e.g., a GPU), since we don't support doing matrix-vector
  //products until after fillComplete is called, at which time contiguous
  //arrays are allocated on the device and matrix data is copied into them.
  Teuchos::RCP<Teuchos::Array<MapGlobalArrayRCP> > col_ind_2D_global_;
  Teuchos::RCP<Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > > > values2D_;

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
// rptr: length num_block_rows + 1
//       rptr[i]: the pt-row corresponding to the i-th block-row
//
// cptr: length num_distinct_block_cols + 1
//       cptr[j]: the pt-col corresponding to the j-th block-col
//
// bptr: length num_block_rows + 1
//       bptr[i]: location in bindx of the first nonzero block-entry
//                of the i-th block-row
//
// bindx: length num-nonzero-block-entries
//        bindx[j]: block-col-index of j-th block-entry
//
// indx: length num-nonzero-block-entries + 1
//       indx[j] location in vals of the beginning of the j-th
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

