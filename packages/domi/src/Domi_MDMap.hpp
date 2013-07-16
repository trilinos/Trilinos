// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_MDMAP_HPP
#define DOMI_MDMAP_HPP

// Teuchos includes
#include "Teuchos_Tuple.hpp"

// Kokkos includes
#include "Kokkos_DefaultNode.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_Slice.hpp"
#include "Domi_MDComm.hpp"

namespace Domi
{

template< class LocalOrd,
          class GlobalOrd = LocalOrd,
          class Node = Kokkos::DefaultNode::DefaultNodeType >
class MDMap
{
public:

  MDMap(const MDCommRCP mdComm,
        const Teuchos::ArrayView< GlobalOrd > & dimensions,
        const Teuchos::ArrayView< int > & halos =
          Teuchos::ArrayView< int >(),
        const Teuchos::ArrayView< int > & ghosts =
          Teuchos::ArrayView< int >(),
        const Teuchos::ArrayView< int > & periodic =
          Teuchos::ArrayView< int >(),
        const EStorageOrder storageOrder = DEFAULT_ORDER,
        const Teuchos::RCP< Node > & node =
          Kokkos::DefaultNode::getDefaultNode());
  MDMap(const Teuchos::RCP< MDMap< LocalOrd, GlobalOrd, Node > > parent,
        const Teuchos::ArrayView< Slice > slices);
  ~MDMap();
  /** \name MDComm pass-through methods */
  //@{
  bool onSubcommunicator() const;
  TeuchosCommRCP getTeuchosComm() const;
  int getNumDims() const;
  int getAxisCommSize(int axis) const;
  int getAxisRank(int axis) const;
  int getLowerNeighbor(int axis) const;
  int getUpperNeighbor(int axis) const;
  //@}
  // GlobalOrd getGlobalDim(int axis, bool withGhosts=false) const;
  // LocalOrd getLocalDim(int axis, bool withHalos=false) const;
  // Slice getBounds(int axis, bool withHalos=false) const;
  // bool hasHalos() const;
  // int getLowerHalo(int axis) const;
  // int getUpperHalo(int axis) const;
  // int getHaloSize(int axis) const;
  // int getGhostSize(int axis) const;
  // bool getPeriodic(int axis) const;
  // EStorageOrder getStorageOrder() const;
  // GlobalOrd getGlobalStride(int axis) const;
  // LocalOrd getLocalStride(int axis) const;
  // // Axis map methods should go here ...
  // // Conversion to Epetra, Tpetra or Xpetra Maps should go here ...
  // Teuchos::Array< GlobalOrd >
  // getGlobalAxisIndex(GlobalOrd globalIndex) const;
  // Teuchos::Array< LocalOrd >
  // getLocalAxisIndex(LocalOrd LocalIndex) const;
  // GlobalOrd getGlobalElement(LocalOrd localIndex) const;
  // GlobalOrd
  // getGlobalElement(const Teuchos::ArrayView< GlobalOrd > globalAxisIndex) const;
  // LocalOrd getLocalElement(GlobalOrd globalIndex) const;
  // LocalOrd
  // getLocalElement(const Teuchos::ArrayView< LocalOrd > localAxisIndex) const;

private:

  typedef Teuchos::Tuple< int, 2 > halo_t;

  MDCommRCP                   _mdComm;
  Teuchos::Array< GlobalOrd > _globalDims;
  Teuchos::Array< Slice >     _localBounds;
  Teuchos::Array< int >       _haloSizes;
  Teuchos::Array< halo_t >    _halos;
  Teuchos::Array< int >       _ghostSizes;
  Teuchos::Array< int >       _periodic;
  Teuchos::Array< GlobalOrd > _globalStrides;
  Teuchos::Array< LocalOrd >  _localStrides;
  EStorageOrder               _storageOrder;
  Teuchos::RCP< Node >        _node;
};

////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const MDCommRCP mdComm,
      const Teuchos::ArrayView< GlobalOrd > & dimensions,
      const Teuchos::ArrayView< int > & halos,
      const Teuchos::ArrayView< int > & ghosts,
      const Teuchos::ArrayView< int > & periodic,
      const EStorageOrder storageOrder,
      const Teuchos::RCP< Node > & node) :
  _mdComm(mdComm),
  _globalDims(dimensions),
  _localBounds(),
  _ghostSizes(mdComm->getNumDims(), 0),
  _periodic(mdComm->getNumDims(), 0),
  _globalStrides(mdComm->getNumDims(), 1),
  _localStrides(mdComm->getNumDims(), 1),
  _haloSizes(mdComm->getNumDims(), 0),
  _halos(),
  _storageOrder(storageOrder),
  _node(node)
{
  // Temporarily store the number of dimensions
  int numDims = mdComm->getNumDims();

  // Check the global dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");

  // Copy the ghost sizes
  for (int axis = 0; axis < ghosts.size() && axis < numDims; ++axis)
  {
    _ghostSizes[axis] = ghosts[axis];
  }

  // Copy the halo sizes and set the halos
  for (int axis = 0; axis < halos.size() && axis < numDims; ++axis)
  {
    _haloSizes[axis] = halos[axis];
    int lower, upper;
    if (mdComm->getLowerNeighbor(axis) == -1)
      lower = _ghostSizes[axis];
    else
      lower = _haloSizes[axis];
    if (mdComm->getUpperNeighbor(axis) == -1)
      upper = _ghostSizes[axis];
    else
      upper = _haloSizes[axis];
    _halos.push_back(Teuchos::tuple(lower, upper));
  }

  // Copy the periodic flags
  for (int axis = 0; axis < periodic.size() && axis < numDims; ++axis)
  {
    _periodic[axis] = periodic[axis];
  }

  // Compute the global strides
  if (_storageOrder == FIRST_INDEX_FASTEST)
    for (int axis = 1; axis < numDims; ++axis)
    {
      _globalStrides[axis] = _globalStrides[axis-1] * _globalDims[axis-1];
    }
  else
    for (int axis = numDims - 2; axis >= 0; ++axis)
    {
      _globalStrides[axis] = _globalStrides[axis+1] * _globalDims[axis+1];
    }
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const Teuchos::RCP< MDMap< LocalOrd, GlobalOrd, Node > > parent,
      const Teuchos::ArrayView< Slice > slices)
{

}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::~MDMap()
{
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
bool
MDMap< LocalOrd, GlobalOrd, Node >::onSubcommunicator() const
{
  return _mdComm->onSubcommunicator();
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
TeuchosCommRCP
MDMap< LocalOrd, GlobalOrd, Node >::getTeuchosComm() const
{
  return _mdComm->getTeuchosComm();
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getNumDims() const
{
  return _mdComm->getNumDims();
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getAxisCommSize(int axis) const
{
  return _mdComm->getAxisSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getAxisRank(int axis) const
{
  return _mdComm->getAxisRank(axis);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getLowerNeighbor(int axis) const
{
  return _mdComm->getLowerNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperNeighbor(int axis) const
{
  return _mdComm->getUpperNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// GlobalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getGlobalDim(int axis,
//              bool withGhosts) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// LocalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getLocalDim(int axis,
//             bool withHalos) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// Slice
// MDMap< LocalOrd, GlobalOrd, Node >::
// getBounds(int axis,
//           bool withHalos) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// int
// MDMap< LocalOrd, GlobalOrd, Node >::getLowerHalo(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// int
// MDMap< LocalOrd, GlobalOrd, Node >::getUpperHalo(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// int
// MDMap< LocalOrd, GlobalOrd, Node >::getHaloSize(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// int
// MDMap< LocalOrd, GlobalOrd, Node >::getGhostSize(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// bool
// MDMap< LocalOrd, GlobalOrd, Node >::getPeriodic(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// EStorageOrder
// MDMap< LocalOrd, GlobalOrd, Node >::getStorageOrder() const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// GlobalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::getGlobalStride(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// LocalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::getLocalStride(int axis) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// Teuchos::Array< GlobalOrd >
// MDMap< LocalOrd, GlobalOrd, Node >::
// getGlobalAxisIndex(GlobalOrd globalIndex) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// Teuchos::Array< LocalOrd >
// MDMap< LocalOrd, GlobalOrd, Node >::
// getLocalAxisIndex(LocalOrd localIndex) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// GlobalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getGlobalElement(LocalOrd localIndex) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// GlobalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getGlobalElement(const Teuchos::ArrayView< GlobalOrd > globalAxisIndex) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// LocalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getLocalElement(GlobalOrd globalIndex) const
// {

// }

// ////////////////////////////////////////////////////////////////////////

// template< class LocalOrd, class GlobalOrd, class Node >
// LocalOrd
// MDMap< LocalOrd, GlobalOrd, Node >::
// getLocalElement(const Teuchos::ArrayView< LocalOrd > localAxisIndex) const
// {

// }

}

#endif
