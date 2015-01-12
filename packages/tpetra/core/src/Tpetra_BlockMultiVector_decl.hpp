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

#ifndef TPETRA_BLOCKMULTIVECTOR_DECL_HPP
#define TPETRA_BLOCKMULTIVECTOR_DECL_HPP

/// \file Tpetra_BlockMultiVector_decl.hpp
/// \brief Declarations for the class Tpetra::BlockMultiVector.

#include <Tpetra_ConfigDefs.hpp>

#ifndef HAVE_TPETRA_CLASSIC_VBR
#  error "It is an error to include this file if VBR (variable-block-size) sparse matrix support is disabled in Tpetra.  If you would like to enable VBR support, please reconfigure Trilinos with the CMake option Tpetra_ENABLE_CLASSIC_VBR set to ON, and rebuild Trilinos."
#else

#include "Tpetra_BlockMap.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra {

/// \class BlockMultiVector
/// \brief Block-entry specialization of Tpetra::MultiVector.
///
/// This class inherits (is-a) Tpetra::MultiVector, adding block-entry
/// functionality for referencing/accessing data.
///
/// \warning This class is DEPRECATED.  There are known outstanding
///   bugs with the current implementations of variable-block-size
///   sparse matrices and related classes in Tpetra.
template <class Scalar = Details::DefaultTypes::scalar_type,
          class LocalOrdinal = Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = Details::DefaultTypes::global_ordinal_type,
          class Node = Details::DefaultTypes::node_type>
class TPETRA_DEPRECATED BlockMultiVector :
  public MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
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

  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > getBlockMap() const
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
  LocalOrdinal getLocalPointIndex (GlobalOrdinal globalBlockRow, LocalOrdinal blockOffset) const;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blockMap_;
};

} // namespace Tpetra

#endif // ! HAVE_TPETRA_CLASSIC_VBR
#endif // ! TPETRA_BLOCKMULTIVECTOR_DECL_HPP
