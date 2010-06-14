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

#ifndef TPETRA_BLOCKMULTIVECTOR_DECL_HPP
#define TPETRA_BLOCKMULTIVECTOR_DECL_HPP

#include "Tpetra_BlockMap.hpp"
#include "Tpetra_MultiVector.hpp"

/** \file Tpetra_BlockMultiVector_decl.hpp

  Declarations for the class Tpetra::BlockMultiVector.
*/
namespace Tpetra {

/** \brief Block-entry specialization of Tpetra::MultiVector.

  This class inherits (is-a) Tpetra::MultiVector, adding block-entry
  functionality for referencing/accessing data.
*/
template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
class BlockMultiVector : public MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
 public:
  typedef Scalar        scalar_type;
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;

  //! @name Constructor/Destructor Methods
  //@{

  BlockMultiVector(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap, size_t NumVectors, bool zeroOut=true);

  //! Destructor
  ~BlockMultiVector(){}

  //@}

  //! @name Attribute Queries
  //@{

  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& getBlockMap() const
  { return blockMap_; }

  //@}

  //! @name Post-construction modification routines
  //@{ 

  //using-declaration avoids warnings about hiding inherited methods:
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue;

  //! Replace current value at the specified (globalBlockRow, blockOffset, vectorIndex) location with specified value.
  void replaceGlobalValue(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value);

  //using-declaration avoids warnings about hiding inherited methods:
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue;

  //! Replace current value at the specified (localBlockRow, blockOffset, vectorIndex) location with specified value.
  void replaceLocalValue(LocalOrdinal localBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value);

  //using-declaration avoids warnings about hiding inherited methods:
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue;

  //! Adds specified value to existing value at the specified (globalBlockRow, blockOffset, vectorIndex) location.
  void sumIntoGlobalValue(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value);

  //using-declaration avoids warnings about hiding inherited methods:
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue;

  //! Adds specified value to existing value at the specified (localBlockRow, blockOffset, vectorIndex) location with specified value.
  void sumIntoLocalValue(LocalOrdinal localBlockRow, LocalOrdinal blockOffset, size_t vectorIndex, const Scalar &value);

  //@}

 private:
  LocalOrdinal getLocalPointIndex(GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset) const;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blockMap_;
};//class BlockMultiVector
}//namespace Tpetra

#endif

