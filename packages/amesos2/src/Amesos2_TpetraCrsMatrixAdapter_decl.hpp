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
 * \tparam Scalar type for scalar values
 * \tparam LocalOrdinal the ordinal type for local index references
 * \tparam GlobalOrdinal the ordinal type for global index references
 * \tparam a Kokkos node type
 */
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node
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

  static const char* name;


  MatrixAdapter();


  /// Copy constructor
  MatrixAdapter(const MatrixAdapter<matrix_type>& adapter);

  
  MatrixAdapter(matrix_type* const m);


  MatrixAdapter(const Teuchos::RCP<matrix_type>& m);


  /** \brief Scales values of \c this by a factor of \c alpha
   *
   * \param [in] alpha scalar factor
   */
  MatrixAdapter<matrix_type>& scale(const Scalar alpha);


  /** \brief Updates the values of \c this to \f$this = alpha*this + beta*B\f$
   *
   * \param [in] beta scalar coefficient of \c B
   * \param [in] B additive MultiVector
   * \param [in] alpha scalar coefficient of \c this
   */
  /* TODO: May not need this function */
  MatrixAdapter<matrix_type>& update(const Scalar beta,
    const MatrixAdapter<matrix_type>& B,
    const Scalar alpha);
  

  bool isLocal();
  

  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;


  size_t getLocalNumRows() const;
    

  size_t getLocalNumCols() const;
  

  global_size_type getGlobalNumRows() const;
  

  global_size_type getGlobalNumCols() const;
  

  size_t getLocalNNZ() const;
  

  global_size_type getGlobalNNZ() const;
  

  size_t getMaxNNZ() const;
  

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const;

  // TODO: The following two methods for Domain and Range Maps are
  // used at least in the matrixShapeOK_impl() method for
  // Amesos2::Superlu.  If their only function is to eventually get
  // the size of the Domain and Range, we might want to provide
  // functions that return those numbers instead of returning the Maps
  // themselves, because the subsequent accessor methods for size
  // might be (probably are) inconsitent accross the types.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >getColMap() const;

  // Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal> > getDomainMap() const;

  /** \brief Gets a compressed-row storage summary of \c this
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
   * \exception std::runtime_error Thrown if \c nzval or \c colind is not
   * large enough to hold the global number of nonzero values.
   *
   * \exception std::runtime_error Thrown if \c rowptr is not at least
   * <tt>nrow + 1</tt> in size, the required size.
   */
  void getCrs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<global_size_type> rowptr,
    size_t& nnz,
    bool local = false,
    int root = 0);
  

  /** \brief Gets a compressed-column storage summary of \c this
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
   * \exception std::runtime_error Thrown if \c nzval or \c rowind is not
   * large enough to hold the global number of nonzero values.
   *
   * \exception std::runtime_error Thrown if \c colptr is not at least
   * <tt>ncol + 1</tt> in size, the required size.
   */
  void getCcs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<global_size_type> colptr,
    size_t& nnz,
    bool local = false,
    int root = 0);


  /**
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
    size_t& nnz,
    int root = 0);


  /** \brief Updates the underlying matrix assuming CRS input.
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
  

  /** \brief Updates the underlying matrix assuming CCS input.
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
  Teuchos::RCP<matrix_type> mat_;

};                              // end class MatrixAdapter<Tpetra::CrsMatrix>


} // end namespace Amesos

#endif  // AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DECL_HPP
