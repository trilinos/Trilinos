/**
  \file   Amesos2_EpetraVbrMatrixAdapter_def.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul 22 10:21:04 CDT 2010
  
  \brief  Amesos2 MatrixAdapter<Epetra_VbrMatrix> definitions.
*/

#ifndef AMESOS2_EPETRA_VBRMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_VBRMATRIX_ADAPTER_DEF_HPP

namespace Amesos {


MatrixAdapter<Epetra_VbrMatrix>::MatrixAdapter()
  : MatrixAdapter<Epetra_RowMatrix>()
{ }


MatrixAdapter<Epetra_VbrMatrix>::MatrixAdapter(
  const MatrixAdapter<Epetra_VbrMatrix>& adapter)
  : MatrixAdapter<Epetra_RowMatrix>(adapter)
{ }


MatrixAdapter<Epetra_VbrMatrix>::MatrixAdapter(
  const RCP<Epetra_VbrMatrix>& m)
  : MatrixAdapter<Epetra_RowMatrix>(m)
{ }


const char* MatrixAdapter<Epetra_VbrMatrix>::name
= "Amesos2 adapter for Epetra_VbrMatrix";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_VBRMATRIX_ADAPTER_DEF_HPP
