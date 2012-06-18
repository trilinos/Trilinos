// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP

#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::RowMatrixTransposer(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& origMatrix)
  : origMatrix_(origMatrix), comm_(origMatrix.getComm()), indexBase_(origMatrix_.getIndexBase()) {}

template<class Scalar, 
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::~RowMatrixTransposer() {}

template<class Scalar, 
  class LocalOrdinal,
  class GlobalOrdinal, 
  class Node, 
  class SpMatOps>
RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::createTranspose (
  const OptimizeOption optimizeTranspose,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap)
{
  optimizeTranspose_ = optimizeTranspose;
  transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>(
    transposeRowMap == null? origMatrix_.getDomainMap() : transposeRowMap, 0));
  ArrayView<const LocalOrdinal> localIndicies;
  ArrayView<const Scalar> localValues;
  ArrayView<const GlobalOrdinal> myGlobalRows = origMatrix_.getRowMap()->getNodeElementList();
  size_t numEntriesInRow;
  Array<GlobalOrdinal> rowNum(1);
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > origColMap = origMatrix_.getColMap();
  for(
    size_t i = Teuchos::OrdinalTraits<size_t>::zero();
    i < origMatrix_.getNodeNumRows();
    ++i
  )
  {
    rowNum[0] = origMatrix_.getRowMap()->getGlobalElement(i);
    numEntriesInRow = origMatrix_.getNumEntriesInLocalRow(i);
    origMatrix_.getLocalRowView(i, localIndicies, localValues);
    for(size_t j=0; j<numEntriesInRow; ++j){
      transposeMatrix_->insertGlobalValues(origColMap->getGlobalElement(localIndicies[j]), rowNum(0,1), localValues(j,1));
    }
  }

  RCP<ParameterList> params = parameterList();
  if (optimizeTranspose_ == DoOptimizeStorage) params->set("Optimize Storage",true);
  else                                         params->set("Optimize Storage",false);
  transposeMatrix_->fillComplete(origMatrix_.getRangeMap(), origMatrix_.getDomainMap(), params);

  return transposeMatrix_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;


}

#endif
