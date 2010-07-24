/**
  \file   Amesos2_MatrixAdapter.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Feb  6 10:00:22 2010

  \brief  A templated adapter class for Matrix Objects.
          Specializations may be created for any matrix that needs to
          be adapted for use by Amesos2.
*/

#ifndef AMESOS2_MATRIX_ADAPTER_DECL_HPP
#define AMESOS2_MATRIX_ADAPTER_DECL_HPP

namespace Amesos {


/**
 * \brief A templated Matrix class adapter for Amesos2.
 *
 * Specializations of this templated class provide a unified interface
 * to Matrix types for Amesos2.  Any specializations are expected to
 * implement the following methods:
 *
 * <br>
 * <b>Implementation Requirements:</b>
 * <ul>
 * <li> Default constructor
 * \code MatrixAdapter<matrix_type>(); \endcode
 * </li>
 *
 * <li> Wrapper constructor.
 *
 * \code
 * MatrixAdapter<matrix_type>(const Teuchos::RCP<matrix_type>& mat);
 * \endcode
 * </li>
 *
 * <li> Copy constructor.
 *
 * \code
 * MatrixAdapter<matrix_type>(const MatrixAdapter<matrix_type>& adapter);
 * \endcode
 * </li>
 *
 * <li> Method to get locality of matrix, either globally or locally indexed.
 *
 * \code
 * bool isLocal() const;
 * \endcode
 * </li>
 *
 * <li> Method to get matrix communicator.
 *
 * \code
 * const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const;
 * \endcode
 * </li>
 *
 * <li> Methods to get the local and global number of rows and columns
 *
 * \code
 * size_t getLocalNumRows() const;
 *
 * size_t getLocalNumCols() const;
 *
 * global_size_type getGlobalNumRows() const;
 *
 * global_size_type getGlobalNumCols() const;
 * \endcode
 * </li>
 *
 * <li> Method to access number of nonzero entries.
 *
 * \code
 * size_t getLocalNNZ() const;
 *
 * global_size_type getGlobalNNZ() const;
 * \endcode
 * </li>
 *
 * <li> Method to get the maximum number of nonzeros in all rows.
 * \code
 * size_t getMaxNNZ() const;
 * \endcode
 * </li>
 *
 * <li> Map methods
 * \code
 * Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > getRowMap() const;
 *
 * Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > getColMap() const;
 * \endcode
 * </li>
 *
 * <li> Methods to get compressed-row and compressed-column representations of
 * the underlying matrix.
 *
 * \code
 * void getCrs(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> colind,
 *   const Teuchos::ArrayView<global_size_type> rowptr,
 *   size_t& nnz,
 *   bool local = false);
 *
 * void getCcs(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> rowind,
 *   const Teuchos::ArrayView<global_size_type> colptr,
 *   size_t& nnz,
 *   bool local = false);
 *
 * void getCrsAll(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> colind,
 *   const Teuchos::ArrayView<global_size_type> rowptr,
 *   size_t& nnz);
 *
 * void getCcsAll(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> rowind,
 *   const Teuchos::ArrayView<global_size_type> colptr,
 *   size_t& nnz)
 * \endcode
 * </li>
 *
 * <li> Methods to update the underlying matrix given CRS and CCS
 * representations.
 *
 * \code
 * void updateValuesCrs(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> colind,
 *   const Teuchos::ArrayView<global_size_type> rowptr);
 *
 * void updateValuesCcs(
 *   const Teuchos::ArrayView<scalar_type> nzval,
 *   const Teuchos::ArrayView<GO> rowind,
 *   const Teuchos::ArrayView<global_size_type> colptr);
 * \endcode
 *
 * <li> Get a description of this adapter
 * \code
 * std::string description() const;
 * \endcode
 * </li>
 *
 * <li>Print the matrix to the \c os output stream
 * \code
 * void describe(
 *   Teuchos::FancyOStream& os,
 *   const Teuchos::EVerbosityLevel verbLevel) const;
 * \endcode
 * </li>
 */
template <class MatrixType>
struct MatrixAdapter
{};

} // end namespace Amesos

#endif // AMESOS2_MATRIX_ADAPTER_DECL_HPP
