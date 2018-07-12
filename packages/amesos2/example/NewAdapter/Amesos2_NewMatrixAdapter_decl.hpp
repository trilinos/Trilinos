// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_NewMatrixAdapter_decl.hpp
  \author John Doe <jd@sandia.gov>
  \date   today

  \brief  Amesos2::MatrixAdapter specialization for the
          NewMatrix class.
*/

#ifndef AMESOS2_NEWMATRIX_ADAPTER_DECL_HPP
#define AMESOS2_NEWMATRIX_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_MPIContainerComm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <NewMatrix.hpp>

#include "Amesos2_MatrixAdapter_decl.hpp"

namespace Amesos {


/**
 * \brief Amesos2 Matrix adapter for the NewMatrix class.
 *
 * TODO: template parameters can be added if needed. 
 */
template <>
class MatrixAdapter< NewMatrix >
{
public:

  /* TODO: Redefine the following types as needed */
  // public type definitions
  typedef double                   scalar_type;
  typedef int                      local_ordinal_type;
  typedef int                      global_ordinal_type;
  typedef size_t                   global_size_type;
  typedef Tpetra::Map<>::node_type node_type;
  typedef NewMatrix                matrix_type;

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
  inline bool isLocal() const;
  

  /// Returns the Teuchos::Comm object associated with this matrix.
  inline const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;


  /// Get the number of matrix rows local to the calling node.
  inline size_t getLocalNumRows() const;


  /// Get the number of natrix columns local to the calling node.
  inline size_t getLocalNumCols() const;


  /// Get the number of global matrix rows
  inline global_size_type getGlobalNumRows() const;


  /// Get the number of global matrix columns
  inline global_size_type getGlobalNumCols() const;


  /// Get the number of non-zero matrix entries for the calling node.
  inline size_t getLocalNNZ() const;


  /// Get the number of global non-zero matrix entries.
  inline global_size_type getGlobalNNZ() const;


  /// Get the maximum number of non-zeros in any global row of this matrix.
  inline size_t getMaxNNZ() const;


  /// Get the row map for this matrix.
  inline
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getRowMap() const;


  /// Get the column map for this matrix.
  inline
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getColMap() const;


  // TODO:  Two methods for getting the Domain and Range Maps are
  // used at least in the matrixShapeOK_impl() method for
  // Amesos2::Superlu.  If their only function is to eventually get
  // the size of the Domain and Range, we might want to provide
  // functions that return those numbers instead of returning the Maps
  // themselves, because the subsequent accessor methods for size
  // might be (probably are) inconsitent accross the types.


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
   * \param [in]  root Is the processor ID of the node that will end up with
   *              the CRS representation of the matrix.
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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> colind,
    const Teuchos::ArrayView<global_size_type> rowptr,
    size_t& nnz,
    bool local = false,
    int root = 0);


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
   * \param [in]  root Is the processor ID of the node that will end up with
   *              the CRS representation of the matrix.
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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> rowind,
    const Teuchos::ArrayView<global_size_type> colptr,
    size_t& nnz,
    bool local = false,
    int root = 0);


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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> colind,
    const Teuchos::ArrayView<global_size_type> rowptr,
    size_t& nnz,
    int root = 0);


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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> rowind,
    const Teuchos::ArrayView<global_size_type> colptr,
    size_t& nnz,
    int root = 0);


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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> colind,
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
    const Teuchos::ArrayView<scalar_type> nzval,
    const Teuchos::ArrayView<global_ordinal_type> rowind,
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

};                              // end class MatrixAdapter<NewMatrix>


} // end namespace Amesos

#endif  // AMESOS2_NEWMATRIX_ADAPTER_DECL_HPP
