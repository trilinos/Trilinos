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
  bool isPeriodic(int axis) const;
  int getAxisRank(int axis) const;
  int getLowerNeighbor(int axis) const;
  int getUpperNeighbor(int axis) const;
  //@}
  GlobalOrd getGlobalDim(int axis, bool withGhosts=false) const;
  LocalOrd getLocalDim(int axis, bool withHalos=false) const;
  Slice getGlobalAxisBounds(int axis, bool withGhosts=false) const;
  Slice getLocalAxisBounds(int axis, bool withHalos=false) const;
  bool hasHalos() const;
  int getLowerHalo(int axis) const;
  int getUpperHalo(int axis) const;
  int getHaloSize(int axis) const;
  int getGhostSize(int axis) const;
  EStorageOrder getStorageOrder() const;
  // GlobalOrd getGlobalStride(int axis) const;
  // LocalOrd getLocalStride(int axis) const;
  // // Axis map methods should go here ...
  // // Conversion to Epetra, Tpetra or Xpetra Maps should go here ...

  Teuchos::Array< GlobalOrd >
  getGlobalAxisIndex(GlobalOrd globalIndex) const;

  Teuchos::Array< LocalOrd >
  getLocalAxisIndex(LocalOrd LocalIndex) const;

  GlobalOrd getGlobalIndex(LocalOrd localIndex) const;

  GlobalOrd
  getGlobalIndex(const Teuchos::ArrayView< GlobalOrd > globalAxisIndex) const;

  LocalOrd getLocalIndex(GlobalOrd globalIndex) const;

  LocalOrd
  getLocalIndex(const Teuchos::ArrayView< LocalOrd > localAxisIndex) const;

private:

  typedef Teuchos::Tuple< int, 2 > halo_t;

  void computeAxisBounds();

  MDCommRCP                   _mdComm;
  Teuchos::Array< GlobalOrd > _globalDims;
  Teuchos::Array< LocalOrd >  _localDims;
  Teuchos::Array< Slice >     _globalAxisBounds;
  Teuchos::Array< Slice >     _localAxisBounds;
  Teuchos::Array< int >       _haloSizes;
  Teuchos::Array< halo_t >    _halos;
  Teuchos::Array< int >       _ghostSizes;
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
      const EStorageOrder storageOrder,
      const Teuchos::RCP< Node > & node) :
  _mdComm(mdComm),
  _globalDims(dimensions),
  _localDims(mdComm->getNumDims(), 0),
  _globalAxisBounds(),
  _localAxisBounds(),
  _ghostSizes(mdComm->getNumDims(), 0),
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
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < halos.size())
      _haloSizes[axis] = halos[axis];
    int lower, upper;
    if (getLowerNeighbor(axis) == -1)
      lower = _ghostSizes[axis];
    else
      lower = _haloSizes[axis];
    if (getUpperNeighbor(axis) == -1)
      upper = _ghostSizes[axis];
    else
      upper = _haloSizes[axis];
    _halos.push_back(Teuchos::tuple(lower, upper));
  }

  // Compute the axis bounds
  computeAxisBounds();

  // Compute the strides
  if (_storageOrder == FIRST_INDEX_FASTEST)
    for (int axis = 1; axis < numDims; ++axis)
    {
      _globalStrides[axis] = _globalStrides[axis-1] * (_globalDims[axis-1] +
                                                       2*_ghostSizes[axis-1]);
      _localStrides[axis]  = _localStrides[axis-1]  * _localDims[axis-1];
    }
  else
    for (int axis = numDims - 2; axis >= 0; ++axis)
    {
      _globalStrides[axis] = _globalStrides[axis+1] * (_globalDims[axis+1] +
                                                       2*_ghostSizes[axis+1]);
      _localStrides[axis]  = _localStrides[axis+1]  * _localDims[axis+1];
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
void
MDMap< LocalOrd, GlobalOrd, Node >::computeAxisBounds()
{
  // Initialization
  int numDims = getNumDims();

  // Decompose the multi-dimensional domain
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Get the communicator info for this axis
    int axisCommSize = getAxisCommSize(axis);
    int axisRank     = getAxisRank(axis);

    // First estimates assuming even division of global dimensions by
    // the number of processors along this axis, and ignoring ghosts
    // and halos.
    _localDims[axis]    = _globalDims[axis] / axisCommSize;
    GlobalOrd axisStart = axisRank * _localDims[axis];

    // Adjustments for non-zero remainder.  Compute the remainder
    // using the mod operator.  If the remainder is > 0, then add an
    // element to the appropriate number of processors with the
    // highest axis ranks.  Note that this is the opposite of the
    // standard Tpetra::Map constructor (which adds an elements to the
    // lowest processor ranks), and provides better balance for finite
    // differencing systems with staggered data location.
    GlobalOrd remainder = _globalDims[axis] % axisCommSize;
    if (axisCommSize - axisRank - 1 < remainder)
    {
      ++_localDims[axis];
      axisStart += (remainder - axisCommSize + axisRank);
    }

    // Global adjustment for ghost points
    axisStart += _ghostSizes[axis];

    // Compute and store the global axis bounds
    _globalAxisBounds.push_back(ConcreteSlice(axisStart,
                                              axisStart + _localDims[axis]));

    // Local adjustment for halos.  Note that _halos should already be
    // corrected to use ghost sizes for processors on the boundaries
    _localDims[axis] += _halos[axis][0] + _halos[axis][1];
    axisStart        -= _halos[axis][0];

    // Compute and store the axis bounds
    _localAxisBounds.push_back(ConcreteSlice(_localDims[axis]));
  }
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
bool
MDMap< LocalOrd, GlobalOrd, Node >::isPeriodic(int axis) const
{
  return _mdComm->isPeriodic(axis);
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

template< class LocalOrd, class GlobalOrd, class Node >
GlobalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalDim(int axis,
             bool withGhosts) const
{
  if (withGhosts)
    return _globalDims[axis] + 2 * _ghostSizes[axis];
  else
    return _globalDims[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
LocalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalDim(int axis,
            bool withHalos) const
{
  if (withHalos)
    return _localDims[axis];
  else
    return _localDims[axis] - _halos[axis][0] - _halos[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalAxisBounds(int axis,
                    bool withGhosts) const
{
  if (withGhosts)
  {
    GlobalOrd start = _globalAxisBounds[axis].start();
    GlobalOrd stop  = _globalAxisBounds[axis].stop();
    if (getAxisRank(axis) == 0)
      start -= _ghostSizes[axis];
    if (getAxisRank(axis) == getAxisCommSize(axis)-1)
      stop += _ghostSizes[axis];
    return ConcreteSlice(start,stop);
  }
  else
    return _globalAxisBounds[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalAxisBounds(int axis,
                   bool withHalos) const
{
  if (withHalos)
    return _localAxisBounds[axis];
  else
  {
    LocalOrd start = _localAxisBounds[axis].start() + _halos[axis][0];
    LocalOrd stop  = _localAxisBounds[axis].stop()  - _halos[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
bool
MDMap< LocalOrd, GlobalOrd, Node >::hasHalos() const
{
  bool result = false;
  for (int axis = 0; axis < getNumDims(); ++axis)
    if (_halos[axis][0] + _halos[axis][1]) result = true;
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getLowerHalo(int axis) const
{
  return _halos[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperHalo(int axis) const
{
  return _halos[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getHaloSize(int axis) const
{
  return _haloSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getGhostSize(int axis) const
{
  return _ghostSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
EStorageOrder
MDMap< LocalOrd, GlobalOrd, Node >::getStorageOrder() const
{
  return _storageOrder;
}

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

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::Array< GlobalOrd >
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalAxisIndex(GlobalOrd globalIndex) const
{
  int numDims = getNumDims();
  Teuchos::Array< GlobalOrd > result(numDims);
  GlobalOrd index = globalIndex;
  if (_storageOrder == LAST_INDEX_FASTEST)
  {
    for (int axis = 0; axis < numDims-1; ++axis)
    {
      result[axis] = index / _globalStrides[axis];
      index        = index % _globalStrides[axis];
    }
    result[numDims-1] = index;
  }
  else
  {
    for (int axis = numDims-1; axis > 0; --axis)
    {
      result[axis] = index / _globalStrides[axis];
      index        = index % _globalStrides[axis];
    }
    result[0] = index;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::Array< LocalOrd >
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalAxisIndex(LocalOrd localIndex) const
{
  int numDims = getNumDims();
  Teuchos::Array< LocalOrd > result(numDims);
  LocalOrd index = localIndex;
  if (_storageOrder == LAST_INDEX_FASTEST)
  {
    for (int axis = 0; axis < numDims-1; ++axis)
    {
      result[axis] = index / _localStrides[axis];
      index        = index % _localStrides[axis];
    }
    result[numDims-1] = index;
  }
  else
  {
    for (int axis = numDims-1; axis > 0; --axis)
    {
      result[axis] = index / _localStrides[axis];
      index        = index % _localStrides[axis];
    }
    result[0] = index;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
GlobalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalIndex(LocalOrd localIndex) const
{
  Teuchos::Array< LocalOrd > localAxisIndex = getLocalAxisIndex(localIndex);
  GlobalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    GlobalOrd globalAxisIndex =
      localAxisIndex[axis] + _globalAxisBounds[axis].start() - _halos[axis][0];
    result += globalAxisIndex * _globalStrides[axis];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
GlobalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalIndex(const Teuchos::ArrayView< GlobalOrd > globalAxisIndex) const
{
  GlobalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
    result += globalAxisIndex[axis] * _globalStrides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
LocalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalIndex(GlobalOrd globalIndex) const
{
  Teuchos::Array< GlobalOrd > globalAxisIndex =
    getGlobalAxisIndex(globalIndex);
  LocalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    LocalOrd localAxisIndex =
      globalAxisIndex[axis] - _globalAxisBounds[axis].start() + _halos[axis][0];
    TEUCHOS_TEST_FOR_EXCEPTION(
      (localAxisIndex < 0 || localAxisIndex >= _localDims[axis]),
      RangeError,
      "global index not on local processor")
    result += localAxisIndex * _localStrides[axis];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
LocalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalIndex(const Teuchos::ArrayView< LocalOrd > localAxisIndex) const
{
  LocalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
    result += localAxisIndex[axis] * _localStrides[axis];
  return result;
}

}

#endif
