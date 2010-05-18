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
#include <Kokkos_DefaultSparseMultiply.hpp>
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
template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultSparseMultiply<Scalar,LocalOrdinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,LocalOrdinal,Node> >
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

  //@}

  //! @name Transformational Methods
  //@{

  void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap, OptimizeOption opt = DoOptimizeStorage);

  void fillComplete(OptimizeOption opt = DoOptimizeStorage);

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

  //private data members:

  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRowMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkColMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkDomainMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRangeMap_;

  Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > blkGraph_;
  Kokkos::VbrMatrix<Scalar,Node> lclMatrix_;
  Teuchos::ArrayRCP<Scalar> pbuf_values1D_;

  LocalMatVec lclMatVec_;

  typedef typename std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > MapGlobalArrayRCP;
  typedef typename std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > MapLocalArrayRCP;

  Teuchos::Array<std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > > col_ind_2D_global_;
  Teuchos::Array<std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > > col_ind_2D_local_;
  Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > > pbuf_values2D_;
};//class VbrMatrix

}//namespace Tpetra

#endif //TPETRA_VBRMATRIX_DECL_HPP

