/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DETAILS_TPETRA_ROWGRAPH_DEF_HPP
#define IFPACK2_DETIALS_TPETRA_ROWGRAPH_DEF_HPP
#include "Ifpack2_Details_Tpetra_RowGraph_decl.hpp"
#include <vector>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif
namespace Ifpack2 {
namespace Details {

//Tpetra_RowGraph (const Teuchos::RCP<const MatrixType& Matrix) :
//Tpetra_RowGraph (const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Matrix) :

//==========================================================================
template<class MatrixType>
Tpetra_RowGraph<MatrixType>::
Tpetra_RowGraph (const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Matrix) :
  A_(Matrix)
{ 
}

//=========================================================================
template<class MatrixType>
Tpetra_RowGraph<MatrixType>::~Tpetra_RowGraph() { }

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & Tpetra_RowGraph<MatrixType>::getComm() const
{
  return A_->getComm ();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP <typename MatrixType::node_type> Tpetra_RowGraph<MatrixType>::getNode() const
{
  return A_->getNode();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Tpetra_RowGraph<MatrixType>::getRowMap() const
{
  return A_->getRowMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Tpetra_RowGraph<MatrixType>::getColMap() const
{
  return A_->getRowMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Tpetra_RowGraph<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Tpetra_RowGraph<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
global_size_t Tpetra_RowGraph<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}

//==========================================================================
template<class MatrixType>
global_size_t Tpetra_RowGraph<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows();
}

//==========================================================================

template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNodeNumCols() const
{
  return A_->getNodeNumCols();
}

//==========================================================================
template<class MatrixType>
typename MatrixType::global_ordinal_type Tpetra_RowGraph<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

//==========================================================================
template<class MatrixType>
global_size_t Tpetra_RowGraph<MatrixType>::getGlobalNumEntries() const
{
  return A_->getGlobalNumEntries();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNodeNumEntries() const
{
  return A_->getGlobalNumEntries();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNumEntriesInGlobalRow(global_ordinal_type globalRow) const
{
  return A_->getNumEntriesInGlobalRow(globalRow);
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNumEntriesInLocalRow(local_ordinal_type localRow) const
{
  return A_->getNumEntriesInLocalRow(localRow);
}

//==========================================================================
template<class MatrixType>
global_size_t Tpetra_RowGraph<MatrixType>::getGlobalNumDiags() const
{
  return A_->getGlobalNumDiags();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return A_->getGlobalMaxNumRowEntries();
}

//==========================================================================
template<class MatrixType>
size_t Tpetra_RowGraph<MatrixType>::getNodeMaxNumRowEntries() const
{
  return A_->getNodeMaxNumRowEntries();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::hasColMap() const
{
  return A_->hasColMap();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::isLocallyIndexed() const
{
  return A_->isLocallyIndexed();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::isGloballyIndexed() const
{
  return A_->isGloballyIndexed();
}

//==========================================================================
template<class MatrixType>
bool Tpetra_RowGraph<MatrixType>::isFillComplete() const
{
  return A_->isFillComplete();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::Import<typename MatrixType::local_ordinal_type,
                                  typename MatrixType::global_ordinal_type,
                                  typename MatrixType::node_type> > Tpetra_RowGraph<MatrixType>::getImporter() const
{
  throw std::runtime_error("Ifpack2::Tpetra_RowGraph: does not support getImporter()");
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::Export<typename MatrixType::local_ordinal_type,
                                  typename MatrixType::global_ordinal_type,
                                  typename MatrixType::node_type> > Tpetra_RowGraph<MatrixType>::getExporter() const
{
  throw std::runtime_error("Ifpack2::Tpetra_RowGraph: does not support getExporter()");
}

//==========================================================================
template<class MatrixType>
void Tpetra_RowGraph<MatrixType>::getGlobalRowCopy(global_ordinal_type GlobalRow,
                                                  const Teuchos::ArrayView<global_ordinal_type> &Indices,
                                                   size_t &NumEntries) const
{
  Teuchos::Array<scalar_type> Values(Indices.size());
  A_->getGlobalRowCopy(GlobalRow,Indices,Values (),NumEntries);
}

//==========================================================================
template<class MatrixType>
void Tpetra_RowGraph<MatrixType>::getLocalRowCopy(local_ordinal_type LocalRow,
                                              const Teuchos::ArrayView<local_ordinal_type> &Indices,
                                              size_t &NumEntries) const
{
  Teuchos::Array<scalar_type> Values(Indices.size());
  A_->getLocalRowCopy(LocalRow,Indices,Values (),NumEntries);
}

}// namespace Details
}// namespace Ifpack2

#endif
