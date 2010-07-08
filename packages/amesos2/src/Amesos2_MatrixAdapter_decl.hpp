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

// TODO: Update required methods list

/** \brief A templated Matrix class adapter for Amesos2.
 *
 * Specializations of this templated class provide a unified interface
 * to Matrix types for Amesos2.  Any specializations are excpected to
 * implement the following methods:
 * 
 * <br><b>Implementation Requirements:</b>
 * <ul>
 * <li>Default constructor
 * \code MatrixAdapter<MatrixType>(); \endcode
 * </li>
 *
 * <li>Wrapper constructor
 * \code MatrixAdapter<MatrixType>(MatrixType* const mat); \endcode
 * </li>
 *
 * <li>Copy constructor
 * \code MatrixAdapter<MatrixType>(MatrixAdapter<MatrixType>& const); \endcode
 * </li>
 *
 * <li>Matrix scaling operation.  Scales each element by \c alpha
 * \code void scale(const scalar_type alpha); \endcode
 * </li>
 *
 * <li>Matrix add operation.  Replace each element in \c this with \f$alpha*this + beta*B\f$
 * \code void add(const scalar_type beta, MatrixAdapter<MatrixType> B, const scalar_type alpha); \endcode
 * </li>
 *
 * <li>Method to get locality of matrix, either globally or locally indexed
 * \code bool isLocal() \endcode
 * </li>
 *
 * <li>Method to get matrix communicator
 * \code virtual const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const; \endcode
 * </li>
 *
 * <li>Methods to get the local and global number of rows and columns
 * \code local_ordinal_type  getLocalNumRows();
 * local_ordinal_type  getLocalNumCols();
 * global_ordinal_type getGlobalNumRows();
 * global_ordinal_type getGlobalNumCols(); \endcode
 * </li>
 *
 * <li>Method to access number of nonzero entries
 * \code local_ordinal_type  getLocalNNZ();
 * global_ordinal_type getGlobalNNZ(); \endcode
 * </li>
 *
 * <li>Method to get the maximum number of nonzeros in all rows.
 * \code global_ordinal_type getMaxNNZ(); \endcode
 * </li>
 *
 * <li>Map methods
 * \code Teuchos::RCP<Map> getRowMap();
 *  Teuchos::RCP<Map> getRangeMap();
 *	Teuchos::RCP<Map> getDomainMap(); \endcode
 * </li>
 *
 * <li>Row access methods
 * \code void getLocalRowCopy(local_ordinal_type localRow,
 *                            const Teuchos::ArrayView<local_ordinal_type> &indices,
 *                            const Teuchos::ArrayView<scalar_type> &values,
 *                            size_t &numEntries) const;
 * void getGlobalRowCopy(global_ordinal_type globalRow,
 *                       const Teuchos::ArrayView<global_ordinal_type> &indices,
 *                       const Teuchos::ArrayView<scalar_type> &values,
 *                       size_t &numEntries) const; \endcode
 * </li>
 *
 * <li>Method to convert Matrix to serial
 * \code matrix_type serial(); \endcode
 * </li>
 *
 * <li>Get a description of this adapter
 * \code std::string description(); \endcode
 * </li>
 * 
 * <li>Print the matrix to the \c os output stream
 * \code void print(Teuchos::FancyOStream& os) const; \endcode
 * </li>
 */
template <class MatrixType>
struct MatrixAdapter {};

//   /** \brief Default constuctor throws compiler error
//    *
//    * Without specialization, the Amesos2::MatrixAdapter simply throws a
//    * compiler error upon instantiation.
//    */
//   MatrixAdapter()
//     {}
// };


} // end namespace Amesos

#endif // AMESOS2_MATRIX_ADAPTER_DECL_HPP
