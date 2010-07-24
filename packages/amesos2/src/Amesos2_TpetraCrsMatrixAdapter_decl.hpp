// @HEADER
// ***********************************************************************
//
//                Amesos2: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

/**
  \file   Amesos2_TpetraCrsMatrixAdapter_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:13:28 CDT 2010

  \brief  Amesos2::MatrixAdapter specialization for the
          Tpetra::CrsMatrix class.
*/

#ifndef AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DECL_HPP
#define AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_MPIContainerComm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2_MatrixAdapter_decl.hpp"

namespace Amesos {


/**
 * \brief Amesos2 Matrix adapter for the Tpetra::CrsMatrix class.
 *
 * \tparam Scalar type for scalar values
 * \tparam LocalOrdinal the ordinal type for local index references
 * \tparam GlobalOrdinal the ordinal type for global index references
 * \tparam a Kokkos node type
 */
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
class MatrixAdapter<Tpetra::CrsMatrix<Scalar,
                                      LocalOrdinal,
                                      GlobalOrdinal,
                                      Node,
                                      LocalMatOps > >
{
public:

  // public type definitions
  typedef Scalar                           scalar_type;
  typedef LocalOrdinal                     local_ordinal_type;
  typedef GlobalOrdinal                    global_ordinal_type;
  typedef Node                             node_type;
  typedef typename Tpetra::global_size_t   global_size_type;
  typedef Tpetra::CrsMatrix<Scalar,
                            LocalOrdinal,
                            GlobalOrdinal,
                            Node,
                            LocalMatOps>   matrix_type;

  /// The name of this adapter class.
  static const char* name;


  MatrixAdapter();


  /// Copy constructor
  MatrixAdapter(const MatrixAdapter<matrix_type>& adapter);


  /**
   * \brief Initialize an adapter from a matrix RCP
   *
   * \param m An RCP pointing to the matrix which is to be wrapped.
   */
  MatrixAdapter(const Teuchos::RCP<matrix_type>& m);


  /// Checks whether this matrix is local to the calling node.
  inline bool isLocal() const
    {
      return( mat_->isLocallyIndexed() );
    }


  /// Returns the Teuchos::Comm object associated with this matrix.
  inline const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return( mat_->getComm() );
    }


  /// Get the number of matrix rows local to the calling node.
  inline size_t getLocalNumRows() const
    {
      return( mat_->getNodeNumRows() );
    }


  /// Get the number of matrix columns local to the calling node.
  inline size_t getLocalNumCols() const
    {
      return mat_->getNodeNumCols();
    }


  /// Get the number of global matrix rows
  inline global_size_type getGlobalNumRows() const
    {
      // return mat_->getGlobalNumRows();
      return mat_->getRowMap()->getMaxAllGlobalIndex() + 1;
    }


  /// Get the number of global matrix columns
  inline global_size_type getGlobalNumCols() const
    {
      // return mat_->getGlobalNumCols();
      return mat_->getColMap()->getMaxAllGlobalIndex() + 1;
    }


  /// Get the number of non-zero matrix entries for the calling node.
  inline size_t getLocalNNZ() const
    {
      return mat_->getNodeNumEntries();
    }


  /// Get the number of global non-zero matrix entries.
  inline global_size_type getGlobalNNZ() const
    {
      return mat_->getGlobalNumEntries();
    }


  /// Get the maximum number of non-zeros in any global row of this matrix.
  inline size_t getMaxNNZ() const
    {
      return mat_->getGlobalMaxNumRowEntries();
    }


  /// Get the row map for this matrix.
  inline
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getRowMap() const
    {
      return mat_->getRowMap();
    }


  /// Get the column map for this matrix.
  inline
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getColMap() const
    {
      return mat_->getColMap();
    }


  // TODO:  Two methods for getting the Domain and Range Maps are
  // used at least in the matrixShapeOK_impl() method for
  // Amesos2::Superlu.  If their only function is to eventually get
  // the size of the Domain and Range, we might want to provide
  // functions that return those numbers instead of returning the Maps
  // themselves, because the subsequent accessor methods for size
  // might be (probably are) inconsistent across the types.


  /**
   * \brief Gets a compressed-row storage summary of \c this
   *
   * Extracts a compressed-row storage format of the matrix and stores the
   * information in the user-supplied containers.
   *
   * \param [out] nzval will hold the values of the nonzero entries of \c this
   * \param [out] colind will hold the column indices of \c this for each row.
   * \param [out] rowptr is of size <tt>nrow + 1</tt> and <tt>rowptr[j]</tt>
   *              stores the location in \c nzval and \c colind which starts
   *              row \c j of \c this.  <tt>rowptr[nrow] = nnz</tt>, where \c
   *              nrow is the number of rows in this matrix.
   * \param [out] nnz is the number of nonzero entries in this matrix.
   * \param [in]  local If \c false, the processor with ID \c root will contain
   *              a representation of the global matrix.  If \c true, then each
   *              processor will end up with a CRS representation of the matrix
   *              rows that it owns.
   *
   * \exception std::length_error Thrown if \c nzval or \c colind is not
   * large enough to hold the global number of nonzero values.
   *
   * \exception std::length_error Thrown if \c rowptr is not at least
   * <tt>nrow + 1</tt> in size, the required size.
   *
   * \exception std::runtime_error Thrown if there is an error while extracting
   * row values from the underlying matrix.
   */
  void getCrs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<global_size_type> rowptr,
    size_t& nnz,
    bool local = false);


  /**
   * \brief Gets a compressed-column storage summary of \c this
   *
   * Extracts a compressed-column storage format of the matrix and stores the
   * information in the user-supplied containers.
   *
   * \param [out] nzval will hold the values of the nonzero entries of \c this
   * \param [out] rowind will hold the row indices of \c this for each column.
   * \param [out] colptr is of size <tt>ncol + 1</tt> and <tt>colptr[j]</tt>
   *              stores the location in \c nzval and \c rowind which starts
   *              column \c j of \c this.  <tt>colptr[ncol] = nnz</tt>, where \c
   *              ncol is the number of columns in this matrix.
   * \param [out] nnz is the number of nonzero entries in this matrix.
   * \param [in]  local If \c false, the processor with ID \c root will contain
   *              a representation of the global matrix.  If \c true, then each
   *              processor will end up with a CRS representation of the matrix
   *              rows that it owns.
   *
   * \exception std::length_error Thrown if \c nzval or \c rowind is not
   * large enough to hold the global number of nonzero values.
   *
   * \exception std::length_error Thrown if \c colptr is not at least
   * <tt>ncol + 1</tt> in size, the required size.
   *
   * \exception std::runtime_error Thrown if there is an error while extracting
   * row values from the underlying matrix.
   */
  void getCcs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<global_size_type> colptr,
    size_t& nnz,
    bool local = false);


  /**
   * \brief Get a compressed-row representation for all nodes.
   *
   * Like \c getCrs() but at the end, each node will have a copy of the CRS
   * representation of the matrix.
   *
   * \note the \c local parameter of \c getCrs() does not make sense in this
   * context, so it is left out.
   */
  void getCrsAll(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<global_size_type> rowptr,
    size_t& nnz);


  /**
   * \brief Get a compressed-column representation for all nodes.
   *
   * Like \c getCcs() but at the end, each node will have a copy of the CCS
   * representation of the matrix.
   *
   * \note the \c local parameter of \c getCcs() does not make sense in this
   * context, so it is left out.
   */
  void getCcsAll(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<global_size_type> colptr,
    size_t& nnz);


  /**
   * \brief Updates the underlying matrix assuming CRS input.
   *
   * Handles both serial and distributed matrices.
   *
   * \param nzval  The new nonzero values of the matrix
   * \param colind The new column indices
   * \param rowptr The new row start indices
   */
  void updateValuesCrs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<global_size_type> rowptr);


  /**
   * \brief Updates the underlying matrix assuming CCS input.
   *
   * Handles both serial and distributed matrices.
   *
   * \param nzval  The new nonzero values of the matrix
   * \param rowind The new row indices
   * \param colptr The new column start indices
   */
  void updateValuesCcs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<global_size_type> colptr);


  /// Get a short description of this adapter class
  std::string description() const;


  /// Print a description of this adapter to the Fancy Output Stream.
  void describe(
    Teuchos::FancyOStream& os,
    const Teuchos::EVerbosityLevel verbLevel) const;


private:

  /// The matrix this adapter wraps.
  Teuchos::RCP<matrix_type> mat_;

};                              // end class MatrixAdapter<Tpetra::CrsMatrix>


} // end namespace Amesos

#endif  // AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DECL_HPP
