// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKVECTOR_DEF_HPP
#define TPETRA_BLOCKVECTOR_DEF_HPP

namespace Tpetra {

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector () :
    base_type ()
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& in,
               const Teuchos::DataAccess copyOrView) :
    base_type (in, copyOrView)
  {}

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
      "Tpetra::BlockVector: Input MultiVector has "
      << X_mv.getNumVectors () << " != 1 columns.");
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const vec_type& X_vec,
               const map_type& meshMap,
               const LO blockSize) :
    base_type (X_vec, meshMap, blockSize)
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& X,
               const map_type& newMeshMap,
               const map_type& newPointMap,
               const size_t offset) :
    base_type (X, newMeshMap, newPointMap, offset)
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockVector<Scalar, LO, GO, Node>::
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& X,
               const map_type& newMeshMap,
               const size_t offset) :
    base_type (X, newMeshMap, offset)
  {}

  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::vec_type
  BlockVector<Scalar, LO, GO, Node>::getVectorView () {
    Teuchos::RCP<vec_type> vPtr = this->mv_.getVectorNonConst (0);
    return *vPtr; // shallow copy
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowIndex, const Scalar vals[]) {
    return ((base_type*) this)->replaceLocalValues (localRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  replaceGlobalValues (const GO globalRowIndex, const Scalar vals[]) {
    return ((base_type*) this)->replaceGlobalValues (globalRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  sumIntoLocalValues (const LO localRowIndex, const Scalar vals[]) {
    return ((base_type*) this)->sumIntoLocalValues (localRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockVector<Scalar, LO, GO, Node>::
  sumIntoGlobalValues (const GO globalRowIndex, const Scalar vals[]) {
    return ((base_type*) this)->sumIntoLocalValues (globalRowIndex, 0, vals);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::const_little_host_vec_type
  BlockVector<Scalar, LO, GO, Node>::
  getLocalBlockHost (const LO localRowIndex, Access::ReadOnlyStruct) const
  {
    return ((const base_type*) this)->getLocalBlockHost(localRowIndex, 0, 
                                                        Access::ReadOnly);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::little_host_vec_type
  BlockVector<Scalar, LO, GO, Node>::
  getLocalBlockHost (const LO localRowIndex, Access::ReadWriteStruct)
  {
    return ((base_type*) this)->getLocalBlockHost(localRowIndex, 0, 
                                                  Access::ReadWrite);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockVector<Scalar, LO, GO, Node>::little_host_vec_type
  BlockVector<Scalar, LO, GO, Node>::
  getLocalBlockHost (const LO localRowIndex, Access::OverwriteAllStruct)
  {
    return ((base_type*) this)->getLocalBlockHost(localRowIndex, 0, 
                                                  Access::OverwriteAll);
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_BLOCKVECTOR_INSTANT(S,LO,GO,NODE) \
  template class BlockVector< S, LO, GO, NODE >; 

#endif // TPETRA_BLOCKVECTOR_DEF_HPP
