/**
  \file   Amesos2_EpetraMsrMatrixAdapter_def.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul 22 10:13:09 CDT 2010
  
  \brief  Amesos2 MatrixAdapter<Epetra_MsrMatrix> definitions.
*/

#ifndef AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DEF_HPP

namespace Amesos {


MatrixAdapter<Epetra_MsrMatrix>::MatrixAdapter()
  : MatrixAdapter<Epetra_RowMatrix>()
{ }


MatrixAdapter<Epetra_MsrMatrix>::MatrixAdapter(
  const MatrixAdapter<Epetra_MsrMatrix>& adapter)
  : MatrixAdapter<Epetra_RowMatrix>(adapter)
{ }


MatrixAdapter<Epetra_MsrMatrix>::MatrixAdapter(
  const RCP<Epetra_MsrMatrix>& m)
  : MatrixAdapter<Epetra_RowMatrix>(m)
{ }


const char* MatrixAdapter<Epetra_MsrMatrix>::name
= "Amesos2 adapter for Epetra_MsrMatrix";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_MSRMATRIX_ADAPTER_DEF_HPP
