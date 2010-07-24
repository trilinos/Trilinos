/**
  \file   Amesos2_EpetraJadMatrixAdapter_def.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul 22 10:28:32 CDT 2010
  
  \brief  Amesos2 MatrixAdapter<Epetra_JadMatrix> definitions.
*/

#ifndef AMESOS2_EPETRA_JADMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_JADMATRIX_ADAPTER_DEF_HPP

namespace Amesos {


MatrixAdapter<Epetra_JadMatrix>::MatrixAdapter()
  : MatrixAdapter<Epetra_RowMatrix>()
{ }


MatrixAdapter<Epetra_JadMatrix>::MatrixAdapter(
  const MatrixAdapter<Epetra_JadMatrix>& adapter)
  : MatrixAdapter<Epetra_RowMatrix>(adapter)
{ }


MatrixAdapter<Epetra_JadMatrix>::MatrixAdapter(
  const RCP<Epetra_JadMatrix>& m)
  : MatrixAdapter<Epetra_RowMatrix>(m)
{ }


const char* MatrixAdapter<Epetra_JadMatrix>::name
= "Amesos2 adapter for Epetra_JadMatrix";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_JADMATRIX_ADAPTER_DEF_HPP
