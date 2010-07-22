/**
  \file   Amesos2_EpetraCrsMatrixAdapter_def.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Tue Jul 20 11:43:44 2010
  
  \brief  Amesos2 MatrixAdapter<Epetra_CrsMatrix> definitions.
*/

#ifndef AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DEF_HPP

namespace Amesos {


MatrixAdapter<Epetra_CrsMatrix>::MatrixAdapter()
  : MatrixAdapter<Epetra_RowMatrix>()
{ }


MatrixAdapter<Epetra_CrsMatrix>::MatrixAdapter(
  const MatrixAdapter<Epetra_CrsMatrix>& adapter)
  : MatrixAdapter<Epetra_RowMatrix>(adapter)
{ }


MatrixAdapter<Epetra_CrsMatrix>::MatrixAdapter(
  const RCP<Epetra_CrsMatrix>& m)
  : MatrixAdapter<Epetra_RowMatrix>(m)
{ }


const char* MatrixAdapter<Epetra_CrsMatrix>::name
= "Amesos2 adapter for Epetra_CrsMatrix";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_CRSMATRIX_ADAPTER_DEF_HPP
