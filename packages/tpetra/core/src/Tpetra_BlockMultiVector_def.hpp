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

#ifndef TPETRA_BLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_BLOCKMULTIVECTOR_DEF_HPP

#include <Tpetra_ConfigDefs.hpp>
#include "Tpetra_BlockMultiVector_decl.hpp"

namespace Tpetra {

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::BlockMultiVector(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap, size_t NumVectors, bool zeroOut)
 : MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blockMap->getPointMap(), NumVectors, zeroOut),
   blockMap_(blockMap)
{
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
LocalOrdinal BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalPointIndex(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset) const
{
  LocalOrdinal LBID = blockMap_->getLocalBlockID(globalBlockRow);

  TEUCHOS_TEST_FOR_EXCEPTION( LBID == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error, "Tpetra::BlockMultiVector::getLocalPointIndex ERROR: specified globalBlockRow not found in local block-map.");

  LocalOrdinal blkSize = blockMap_->getLocalBlockSize(LBID);

  TEUCHOS_TEST_FOR_EXCEPTION( blockOffset >= blkSize, std::runtime_error,
     "Tpetra::BlockMultiVector::getLocalPointIndex ERROR: specified blockOffset >= blockSize.");

  LocalOrdinal pointIndex = blockMap_->getFirstLocalPointInLocalBlock(LBID);

  return pointIndex + blockOffset;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value)
{
  LocalOrdinal pointIndex = getLocalPointIndex(globalBlockRow, blockOffset);

  replaceLocalValue(pointIndex, vectorIndex, value);
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue(LocalOrdinal localBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value)
{
  LocalOrdinal pointIndex = blockMap_->getFirstLocalPointInLocalBlock(localBlockRow);

  replaceLocalValue(pointIndex+blockOffset, vectorIndex, value);
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value)
{
  LocalOrdinal pointIndex = getLocalPointIndex(globalBlockRow, blockOffset);

  sumIntoLocalValue(pointIndex, vectorIndex, value);
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue(LocalOrdinal localBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value)
{
  LocalOrdinal pointIndex = blockMap_->getFirstLocalPointInLocalBlock(localBlockRow);

  sumIntoLocalValue(pointIndex+blockOffset, vectorIndex, value);
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_BLOCKMULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class BlockMultiVector< SCALAR , LO , GO , NODE >;

}//namespace Tpetra

#endif // ! TPETRA_BLOCKMULTIVECTOR_DEF_HPP
