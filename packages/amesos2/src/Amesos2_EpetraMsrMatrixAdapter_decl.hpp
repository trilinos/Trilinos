/**
  \file   Amesos2_EpetraMsrMatrixAdapter_decl.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul 22 10:11:36 CDT 2010
  
  \brief  Amesos2 matrix adapter for the Epetra_MsrMatrix class.
*/

#ifndef AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DECL_HPP
#define AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DECL_HPP

#include <Epetra_MsrMatrix.h>

#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_EpetraRowMatrixAdapter.hpp"

namespace Amesos {


/**
 * \brief Amesos2 matrix adapter for the Epetra_MsrMatrix class.
 *
 * This class inherits all methods and functionality from the
 * MatrixAdapter<Epetra_RowMatrix> class, but overrides the \c matrix_type
 * typedef and the \c name attribute.
 */
template<>
class MatrixAdapter< Epetra_MsrMatrix > : public MatrixAdapter< Epetra_RowMatrix >
{
public:
  // override matrix typedef
  typedef Epetra_MsrMatrix matrix_type;

  /// The name of this adapter class.
  static const char* name;


  /// Default constructor
  MatrixAdapter();


  /// Copy constructor
  MatrixAdapter(const MatrixAdapter<matrix_type>& adapter);


  /**
   * \brief Initialize an adapter from a matrix RCP
   *
   * \param m An RCP pointing to the matrix which is to be wrapped.
   */
  MatrixAdapter(const Teuchos::RCP<matrix_type>& m);

};


} // end namespace Amesos

#endif  // AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DECL_HPP
