/**
  \file   Amesos2_EpetraCrsMatrixAdapter_decl.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Tue Jul 20 11:36:54 2010
  
  \brief  Amesos2 matrix adapter for the Epetra_CrsMatrix class.
*/

#ifndef AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DECL_HPP
#define AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DECL_HPP

#include <Epetra_CrsMatrix.h>

#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_EpetraRowMatrixAdapter.hpp"

namespace Amesos {


/**
 * \brief Amesos2 matrix adapter for the Epetra_CrsMatrix class.
 *
 * This class inherits all methods and functionality from the
 * MatrixAdapter<Epetra_RowMatrix> class, but overrides the \c matrix_type
 * typedef and the \c name attribute.
 */
template<>
class MatrixAdapter< Epetra_CrsMatrix > : public MatrixAdapter< Epetra_RowMatrix >
{
public:
  // override matrix typedef
  typedef Epetra_CrsMatrix matrix_type;

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

#endif  // AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DECL_HPP
