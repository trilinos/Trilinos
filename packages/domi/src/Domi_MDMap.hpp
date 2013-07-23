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

/** \brief Multi-dimensional map
 *
 * The <tt>MDMap</tt> class is intended to perform the functions of
 * Epetra or Tpetra maps, except that the data is multi-dimensional.
 * Rather than taking a <tt>Teuchos::Comm</tt> object in its
 * constructor, an <tt>MDMap</tt> object takes a <tt>Domi::MDComm</tt>
 * object.  An <tt>MDComm</tt> object is ordered in a structured way,
 * with processors assigned to each axis, whose product is the total
 * number of processors of its underlying communicator.
 *
 * At a bare minimum, an <tt>MDMap</tt> constructor requires an
 * <tt>MDComm</tt> and an array of dimensions indicating the size and
 * shape of the data structure(s) the <tt>MDMap</tt> will describe.
 * Each dimension will be decomposed along each axis, according to the
 * number of processors along that axis in the <tt>MDComm</tt>, in as
 * even a fashion as is possible.
 *
 * Each axis may be flagged as periodic when constructing the
 * <tt>MDComm</tt>.  This attribute is transferred to the
 * <tt>MDMap</tt> at construction.
 *
 * There are two distinct exterior buffer concepts employed by
 * <tt>MDMaps</tt>, halos and ghosts.  Halos refer to local extensions
 * of the data structure used to receive communication from
 * neighboring processors.  Ghosts refer to points added to the
 * exterior of the global domain that are useful, for example, in
 * certain finite differencing algorithms.  They are treated
 * separately from the dimensions so that the decomposition will be
 * performed relative to the dimensions only, and not influenced by
 * the ghost points.
 *
 * From a global perspective, the only points that are visible are
 * ghost points.  Halos are internal and essentially invisible.  From
 * a local perspective, the exterior buffers are refered to as the
 * halos, although if a processor is on a domain boundary, the buffer
 * will actually be ghost points and not communication halos.
 *
 * The sizes of ghost and halo buffers can be different along each
 * axis, and the sizes of ghosts and halos can be different from each
 * other.
 *
 * Whe you want to loop over the elements of a data structure
 * described by an <tt>MDMap</tt>, you can obtain the local looping
 * bounds along a given axis with the <tt>getLocalAxisBounds(int
 * axis)</tt> method.  This returns a <tt>Slice</tt> object, whose
 * <tt>start()</tt> and <tt>stop()</tt> methods return the loop bounds
 * (<tt>stop()</tt> is non-inclusive).  This method also takes a
 * boolean flag indicating whether the bounds should include the halo
 * points, which defaults to false.
 *
 * There is a corresponding <tt>getGlobalAxisBounds(int axis)</tt>
 * method that returns global looping bounds.  While it is highly
 * unlikely that a user will ever actually loop over global bounds,
 * the method can be useful for obtaining size information.  This
 * method also has a boolean flag argument to control inclusion of
 * ghost points in these bounds.
 */
template< class LocalOrd,
          class GlobalOrd = LocalOrd,
          class Node = Kokkos::DefaultNode::DefaultNodeType >
class MDMap
{
public:

  /** \name Constructors and destructor */
  //@{

  /** \brief Main constructor
   */
  MDMap(const MDCommRCP mdComm,
        const Teuchos::ArrayView< GlobalOrd > & dimensions,
        const Teuchos::ArrayView< int > & halos =
          Teuchos::ArrayView< int >(),
        const Teuchos::ArrayView< int > & ghosts =
          Teuchos::ArrayView< int >(),
        const EStorageOrder storageOrder = DEFAULT_ORDER,
        const Teuchos::RCP< Node > & node =
          Kokkos::DefaultNode::getDefaultNode());

  /** \brief Parent/slice constructor
   */
  MDMap(const Teuchos::RCP< MDMap< LocalOrd, GlobalOrd, Node > > parent,
        const Teuchos::ArrayView< Slice > slices);

  /** \brief Destructor
   */
  ~MDMap();

  //@}

  /** \name MDComm pass-through methods */
  //@{

  /** \brief Query whether this processor is on the sub-communicator
   *
   * Sub-communicators are formed when a parent MDComm is sliced by
   * using the (parent,slices) constructor.  For a full communicator,
   * this method will always return true.
   */
  bool onSubcommunicator() const;

  /** \brief Get the Teuchos communicator
   *
   * Note that if the communicator is not a full communicator, i.e. a
   * sub-communicator, that the underlying Comm pointer may be NULL,
   * depending on this processor's rank.
   */
  TeuchosCommRCP getTeuchosComm() const;

  /** \brief Get the number of dimensions
   *
   * This method will return 0 if the communicator is a
   * sub-communicator and this processor does not belong to the
   * sub-communicator.
   */
  int getNumDims() const;

  /** \brief Get the size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  int getAxisCommSize(int axis) const;

  /** \brief Return the periodic flag for the given axis.
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  bool isPeriodic(int axis) const;

  /** \brief Get the rank of the lower neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  int getAxisRank(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   *
   * If the periodic flag for the given axis is set, the returned
   * lower neighbor will be the highest rank of the highest axis rank
   * processor along this axis.
   */
  int getLowerNeighbor(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * If the periodic flag for the given axis is set, the lower
   * neighbor will be the rank of the zero axis rank processor along
   * this axis.
   */
  int getUpperNeighbor(int axis) const;

  //@}

  /** \name Dimensions and looping bounds */
  //@{

  /** \brief Get the global dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withGhosts [in] specify whether the dimension should
   *        include ghost points or not
   */
  GlobalOrd getGlobalDim(int axis, bool withGhosts=false) const;

  /** \brief Get the local dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withHalos [in] specify whether the dimension should
   *        include halos or not
   */
  LocalOrd getLocalDim(int axis, bool withHalos=false) const;

  /** \brief Get the global loop bounds along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withGhosts [in] specify whether the dimension should
   *        include ghost points or not
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  Slice getGlobalAxisBounds(int axis, bool withGhosts=false) const;

  /** \brief Get the local dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withHalos [in] specify whether the dimension should
   *        include halos or not
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  Slice getLocalAxisBounds(int axis, bool withHalos=false) const;

  //@}

  /** \name Halos, ghost points, and storage order */
  //@{

  /** \brief Return true if there are any exterior buffers stored
   *         locally
   *
   * Note that it is not as simple as whether there were halo points
   * specified in the constructor.  Suppose non-zero halos were
   * specified, but the problem is on one processor and no ghost
   * points were specified: this method will return false.  Similarly,
   * if no halos were specified but ghost points were, then boundary
   * processors will have exterior buffers and this method will return
   * true.
   */
  bool hasHalos() const;

  /** \brief Get the size of the lower halo along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that on processors that contain lower domain boundaries, the
   * returned value will be the ghost size.
   */
  int getLowerHalo(int axis) const;

  /** \brief Get the size of the upper halo along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that on processors that contain upper domain boundaries, the
   * returned value will be the ghost size.
   */
  int getUpperHalo(int axis) const;

  /** \brief Get the halo size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the halos along the given axis at the
   * time of construction, regardless of the processor's position
   * relative to the domain boundary.
   */
  int getHaloSize(int axis) const;

  /** \brief Get the ghost size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getGhostSize(int axis) const;

  /** \brief Get the storage order
   */
  EStorageOrder getStorageOrder() const;

  //@}

  // GlobalOrd getGlobalStride(int axis) const;
  // LocalOrd getLocalStride(int axis) const;
  // // Axis map methods should go here ...
  // // Conversion to Epetra, Tpetra or Xpetra Maps should go here ...

  /** \name GLobal/local index conversion methods */
  //@{

  /** \brief Convert a global index to an array of global axis indexes
   *
   * \param globalIndex [in] a unique 1D global identifier
   */
  Teuchos::Array< GlobalOrd >
  getGlobalAxisIndex(GlobalOrd globalIndex) const;

  /** \brief Convert a local index to an array of local axis indexes
   *
   * \param localIndex [in] a unique 1D local identifier
   */
  Teuchos::Array< LocalOrd >
  getLocalAxisIndex(LocalOrd LocalIndex) const;

  /** \brief Convert a local index to a global index
   *
   * \param localIndex [in] a unique 1D local identifier
   */
  GlobalOrd getGlobalIndex(LocalOrd localIndex) const;

  /** \brief convert an array of global axis indexes to a global index
   *
   * \param globalAxisIndex [in] a multi-dimensional global axis index
   */
  GlobalOrd
  getGlobalIndex(const Teuchos::ArrayView< GlobalOrd > globalAxisIndex) const;

  /** \brief Convert a global index to a local index
   *
   * \param globalIndex [in] a unique 1D global identifier
   *
   * This method can throw a Domi::RangeError if the global index is
   * not on the current processor.
   */
  LocalOrd getLocalIndex(GlobalOrd globalIndex) const;

  /** \brief Convert an array of local axis indexes to a local index
   *
   * \param localAxisIndex [in] a multi-dimensional local axis index
   */
  LocalOrd
  getLocalIndex(const Teuchos::ArrayView< LocalOrd > localAxisIndex) const;

  //@}

private:

  // Typedef for the halo type, a tuple that stores lower and upper
  // halo sizes, which can be different due to processors on domain
  // boundaries.
  typedef Teuchos::Tuple< int, 2 > halo_t;

  // A private method for computing the axis bounds, after the
  // dimensions, halos and ghost points have been properly assigned.
  void computeAxisBounds();

  // The underlying multi-dimensional communicator.
  MDCommRCP _mdComm;

  // The size of the global dimensions along each axis.  This does NOT
  // include the values of the ghost sizes.
  Teuchos::Array< GlobalOrd > _globalDims;

  // The size of the local dimensions along each axis.  This DOES
  // include the values of the halos.
  Teuchos::Array< LocalOrd > _localDims;

  // The global loop bounds along each axis, stored as an array of
  // Slices.  These bounds do NOT include the ghost points.
  Teuchos::Array< Slice > _globalAxisBounds;

  // The local loop bounds along each axis, stored as an array of
  // Slices.  These bounds DO include the halos.
  Teuchos::Array< Slice > _localAxisBounds;

  // The halo sizes that were specified at construction, one value
  // along each axis.
  Teuchos::Array< int > _haloSizes;

  // The actual halos stored on this processor along each axis.  The
  // lower and upper halo can be different due to the processor
  // position on the boundary of a domain.
  Teuchos::Array< halo_t > _halos;

  // The size of the ghost points along each axis.
  Teuchos::Array< int > _ghostSizes;

  // The global stride between adjacent elements.  This quantity is
  // "virtual" as it does not describe actual memory layot, but it is
  // useful for index conversion algorithms.
  Teuchos::Array< GlobalOrd > _globalStrides;

  // The local stride between adjecent elements in memory.
  Teuchos::Array< LocalOrd > _localStrides;

  // The storage order
  EStorageOrder _storageOrder;

  // The Kokkos node type
  Teuchos::RCP< Node > _node;
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
  return _mdComm->getAxisCommSize(axis);
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
