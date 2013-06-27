/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_TPETRA_ROWGRAPH_DEF_HPP
#define IFPACK2_TPETRA_ROWGRAPH_DEF_HPP
#include "Ifpack2_Tpetra_RowGraph_decl.hpp"
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


}// namespace Ifpack2

#endif
