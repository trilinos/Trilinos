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

#ifndef TPETRA_EXPERIMENTAL_BLOCKVECTOR_DEF_HPP
#define TPETRA_EXPERIMENTAL_BLOCKVECTOR_DEF_HPP

#include "Tpetra_Experimental_BlockVector_decl.hpp"

namespace Tpetra {
namespace Experimental {

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const map_type& meshMap, const LO blockSize) :
    base_type (meshMap, blockSize, 1)
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const map_type& meshMap,
               const map_type& pointMap,
               const LO blockSize) :
    base_type (meshMap, pointMap, blockSize, 1)
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const mv_type& X_mv,
               const map_type& meshMap,
               const LO blockSize) :
    base_type (X_mv, meshMap, blockSize)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      X_mv.getNumVectors () != 1, std::invalid_argument,
      "Tpetra::Experimental::BlockVector: Input MultiVector has "
      << X_mv.getNumVectors () << " != 1 columns.");
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector () : base_type () {}

  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::vec_type
  BlockVector<Scalar, LO, GO, Node>::getVectorView () {
    Teuchos::RCP<vec_type> vPtr = this->mv_.getVectorNonConst (0);
    vPtr->setCopyOrView (Teuchos::View);
    return *vPtr; // shallow copy
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowIndex, const Scalar vals[]) const {
    return ((const base_type*) this)->replaceLocalValues (localRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  replaceGlobalValues (const GO globalRowIndex, const Scalar vals[]) const {
    return ((const base_type*) this)->replaceGlobalValues (globalRowIndex, 0, vals);
  }


  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  sumIntoLocalValues (const LO localRowIndex, const Scalar vals[]) const {
    return ((const base_type*) this)->sumIntoLocalValues (localRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  sumIntoGlobalValues (const GO globalRowIndex, const Scalar vals[]) const {
    return ((const base_type*) this)->sumIntoLocalValues (globalRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  getLocalRowView (const LO localRowIndex, Scalar*& vals) const {
    return ((const base_type*) this)->getLocalRowView (localRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  getGlobalRowView (const GO globalRowIndex, Scalar*& vals) const {
    return ((const base_type*) this)->getGlobalRowView (globalRowIndex, 0, vals);
  }

  /// \brief Get a view of the degrees of freedom at the given mesh point.
  ///
  /// \warning This method's interface may change or disappear at any
  ///   time.  Please do not rely on it in your code yet.
  ///
  /// The preferred way to refer to little_vec_type is to get it from
  /// BlockVector's typedef.  This is because different
  /// specializations of BlockVector reserve the right to use
  /// different types to implement little_vec_type.  This gives us a
  /// porting strategy to move from "classic" Tpetra to the Kokkos
  /// refactor version.
  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::little_vec_type
  BlockVector<Scalar, LO, GO, Node>::
  getLocalBlock (const LO localRowIndex) const
  {
    if (! this->isValidLocalMeshIndex (localRowIndex)) {
      return little_vec_type (NULL, 0, 0);
    } else {
      const LO strideX = this->getStrideX ();
      const size_t blockSize = this->getBlockSize ();
      const size_t offset = localRowIndex * blockSize * strideX;
      return little_vec_type (this->getRawPtr () + offset, blockSize, strideX);
    }
  }

} // namespace Experimental
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_EXPERIMENTAL_BLOCKVECTOR_INSTANT(S,LO,GO,NODE) \
  template class Experimental::BlockVector< S, LO, GO, NODE >;

#endif // TPETRA_EXPERIMENTAL_BLOCKVECTOR_DEF_HPP
