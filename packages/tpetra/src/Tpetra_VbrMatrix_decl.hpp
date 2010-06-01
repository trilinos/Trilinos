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
#include <Kokkos_VbrMatrix.hpp>
#include <Kokkos_DefaultBlockSparseMultiply.hpp>
#include <Kokkos_DefaultSparseSolve.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerializationTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_CrsGraph.hpp"

namespace Tpetra {

//! \brief VbrMatrix: Variable block row matrix.
/**
This class is under construction, not yet usable.
*/
template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultBlockSparseMultiply<Scalar,LocalOrdinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,LocalOrdinal,Node> >
class VbrMatrix : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
 public:
  typedef Scalar        scalar_type;
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;
  typedef LocalMatVec   mat_vec_type;
  typedef LocalMatSolve mat_solve_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying the row-map and the max number of non-zeros for all rows.
  VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! Destructor
  virtual ~VbrMatrix();

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
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
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
  //@}

  //! @name Insertion/Removal Methods
  //@{

  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry);

  void setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry);

  //@}

  //! @name Transformational Methods
  //@{

  void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap);

  void fillComplete();
  //@}

  //! @name Extraction Methods
  //@{

  void getGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol,
                           LocalOrdinal& numPtRows, LocalOrdinal& numPtCols,
                           Teuchos::ArrayRCP<const Scalar>& blockEntry);

  //@}

 private:
  //private methods:

  Teuchos::ArrayRCP<Scalar> getGlobalBlockEntryView(GlobalOrdinal globalBlockRow,
                                                    GlobalOrdinal globalBlockCol,
                                                    size_t rowsPerBlock = 0,
                                                    size_t colsPerBlock = 0);
  Teuchos::RCP<Node> getNode();

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
  //matrix in packed (optimized) form. For a description of these
  //arrays, see the text at the bottom of this file.
  Teuchos::ArrayRCP<Scalar> pbuf_values1D_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_rptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_cptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_bptr_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_bindx_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_indx_;

  LocalMatVec lclMatVec_;

  typedef typename std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > MapGlobalArrayRCP;
  typedef typename std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > MapLocalArrayRCP;

  //We use 3 arrays (well, arrays-of-maps, arrays-of-arrays...) to
  //represent the variable-block-row matrix in un-packed '2D' form.
  Teuchos::Array<std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > > col_ind_2D_global_;
  Teuchos::Array<std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > > col_ind_2D_local_;
  Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > > pbuf_values2D_;

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

