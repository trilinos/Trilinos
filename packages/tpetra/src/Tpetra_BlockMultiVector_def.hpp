// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_BLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_BLOCKMULTIVECTOR_DEF_HPP

#include "Tpetra_ConfigDefs.hpp"

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

  TEST_FOR_EXCEPTION( LBID == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error, "Tpetra::BlockMultiVector::getLocalPointIndex ERROR: specified globalBlockRow not found in local block-map.");

  LocalOrdinal blkSize = blockMap_->getLocalBlockSize(LBID);

  TEST_FOR_EXCEPTION( blockOffset >= blkSize, std::runtime_error,
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

#endif

