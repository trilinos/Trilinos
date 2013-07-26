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
 * axis, and the sizes of ghosts and halos on any given axis can be
 * different from each other.
 *
 * When you want to loop over the elements of a data structure
 * described by an <tt>MDMap</tt>, you can obtain the local looping
 * bounds along a given axis with the <tt>getLocalBounds(int
 * axis)</tt> method.  This returns a <tt>Slice</tt> object, whose
 * <tt>start()</tt> and <tt>stop()</tt> methods return the loop bounds
 * (<tt>stop()</tt> is non-inclusive).  This method also takes a
 * boolean flag indicating whether the bounds should include the halo
 * points, which defaults to false.
 *
 * There are two methods for obtaining global loop bounds,
 * <tt>getGlobalBounds(int axis)</tt> and <tt>getGlobalRankBounds(int
 * axis)</tt>.  The first returns the global bounds for the full data
 * structure and the second returns the bounds on this processor.
 *
 * Note that all indexes start at zero, whether global or local,
 * whether unique 1D IDs or multidimensional indexes, whether ghosts
 * or halos are present.  This allows negative indexes to represent
 * reverse indexing from the end of a dimension.
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
        const Teuchos::ArrayView< Slice > & slices,
        const Teuchos::ArrayView< int > & ghosts =
        Teuchos::ArrayView< int >());

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

  /** \brief Get the bounds of the global problem
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withGhosts [in] specify whether the bounds should include
   *        ghost points or not
   */
  Slice getGlobalBounds(int axis, bool withGhosts=false) const;

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
  Slice getGlobalRankBounds(int axis, bool withGhosts=false) const;

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
  Slice getLocalBounds(int axis, bool withHalos=false) const;

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

  /** \brief Get the size of the lower ghost along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getLowerGhost(int axis) const;

  /** \brief Get the size of the upper ghost along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getUpperGhost(int axis) const;

  /** \brief Get the ghost size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the ghosts along the given axis at the
   * time of construction, regardless of whether a sub-map reduced
   * these values.
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

  // Typedef for the exteriior buffer type, a tuple that stores lower
  // and upper ghost or halo sizes, which can be different due to
  // processors on domain boundaries.
  typedef Teuchos::Tuple< int, 2 > ext_buffer;

  // A private method for computing the bounds and local dimensions,
  // after the global dimensions, halos and ghost points have been
  // properly assigned.
  void computeBounds();

  // The underlying multi-dimensional communicator.
  MDCommRCP _mdComm;

  // The size of the global dimensions along each axis.  This includes
  // the values of the ghost sizes.
  Teuchos::Array< GlobalOrd > _globalDims;

  // Store the start and stop indexes of the global problem along each
  // axis.  This includes the values of the ghost sizes.
  Teuchos::Array< Slice > _globalBounds;

  // The minumum 1D index of the global data structure, including
  // ghost points.  This will only be non-zero on a sub-map.
  GlobalOrd _globalMin;

  // The maximum 1D index of the global data structure, including
  // ghost points.
  GlobalOrd _globalMax;

  // The size of the local dimensions along each axis.  This includes
  // the values of the halos.
  Teuchos::Array< LocalOrd > _localDims;

  // The maximum 1D index of the local data structure, including halo
  // points.
  LocalOrd _localMax;

  // The minimum 1D index of the local data structure, including halo
  // points.
  LocalOrd _localMin;

  // A double array of Slices that stores the starting and stopping
  // indexes for each axis processor along each axis.  This data
  // structure allows any processor to know the map bounds on any
  // other processor.  These bounds do NOT include the ghost points.
  // The outer array will have a size equal to the number of
  // dimensions.  Each inner array will have a size equal to the
  // number of axis processors along that axis.
  Teuchos::Array< Teuchos::Array< Slice > > _globalRankBounds;

  // The local loop bounds along each axis, stored as an array of
  // Slices.  These bounds DO include the halos.
  Teuchos::Array< Slice > _localBounds;

  // The halo sizes that were specified at construction, one value
  // along each axis.
  Teuchos::Array< int > _haloSizes;

  // The actual halos stored on this processor along each axis.  The
  // lower and upper halos can be different due to the processor
  // position on the boundary of a domain.
  Teuchos::Array< ext_buffer > _halos;

  // The size of the ghost points along each axis.
  Teuchos::Array< int > _ghostSizes;

  // The actual ghost points stored along each axis.  For a full
  // communicator, the lower and upper ghosts are both the same as the
  // corresponding value in _ghostSizes.  However, the introduction of
  // sub-maps creates the possibility of different upper and lower
  // ghost values.
  Teuchos::Array< ext_buffer > _ghosts;

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
  _globalDims(mdComm->getNumDims()),
  _globalBounds(),
  _globalMin(0),
  _localDims(mdComm->getNumDims(), 0),
  _globalRankBounds(mdComm->getNumDims()),
  _localBounds(),
  _ghostSizes(mdComm->getNumDims(), 0),
  _ghosts(),
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

  // Copy the ghost sizes and compute the global dimensions and bounds
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < ghosts.size())
      _ghostSizes[axis] = ghosts[axis];
    _ghosts.push_back(Teuchos::tuple(_ghostSizes[axis], _ghostSizes[axis]));
    _globalDims[axis] = dimensions[axis] + 2*_ghostSizes[axis];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

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

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  computeBounds();
  _localMax = computeSize(_localDims());

  // Compute the strides
  _globalStrides = computeStrides(_globalDims, _storageOrder);
  _localStrides  = computeStrides(_localDims , _storageOrder);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const Teuchos::RCP< MDMap< LocalOrd, GlobalOrd, Node > > parent,
      const Teuchos::ArrayView< Slice > & slices,
      const Teuchos::ArrayView< int > & ghosts) :
  _globalDims(parent->getNumDims()),
  _globalBounds(),
  _globalMin(0),
  _globalMax(0),
  _localDims(parent->getNumDims(), 0),
  _localMin(0),
  _localMax(0),
  _globalRankBounds(parent->_globalRankBounds),
  _localBounds(),
  _ghostSizes(parent->getNumDims(), 0),
  _ghosts(),
  _haloSizes(parent->_haloSizes),
  _halos(),
  _globalStrides(parent->_globalStrides),
  _localStrides(parent->_localStrides),
  _storageOrder(parent->getStorageOrder()),
  _node(parent->_node)
{
  // Temporarily store the number of dimensions
  int numDims = parent->getNumDims();

  // Sanity check on dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    (slices.size() != numDims),
    InvalidArgument,
    "number of slices = " << slices.size() << " != parent MDMap number of "
    "dimensions = " << numDims);

  // Copy the ghost sizes and set initial values for _ghosts
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < ghosts.size())
      _ghostSizes[axis] = ghosts[axis];
    _ghosts[axis][0] = _ghostSizes[axis];
    _ghosts[axis][1] = _ghostSizes[axis];
  }

  // Convert the slices to concrete, add the ghost points, and store
  // in _globalBounds, altering _ghosts if necessary.  Compute
  // _globalDims, _globalMax, and _globalMin.
  for (int axis = 0; axis < numDims; ++axis)
  {
    Slice bounds =
      slices[axis].bounds(parent->getGlobalBounds(axis).stop());
    GlobalOrd start = bounds.start() - _ghostSizes[axis];
    if (start < 0)
    {
      _ghosts[axis][0] = bounds.start();
      start = 0;
    }
    GlobalOrd stop = bounds.stop() + _ghostSizes[axis];
    if (stop > parent->getGlobalDim(axis,true))
    {
      _ghosts[axis][1] = parent->getGlobalDim(axis,true) - bounds.stop();
      stop = parent->getGlobalDim(axis,true);
    }
    _globalBounds.push_back(ConcreteSlice(start,stop));
    _globalDims[axis] = stop - start;
    _globalMin +=  start   * _globalStrides[axis];
    _globalMax += (stop-1) * _globalStrides[axis];
  }

  // Build the array of slices for the MDComm sub-communicator
  // constructor
  Teuchos::Array< Slice > axisRankSlices;
  for (int axis = 0; axis < numDims; ++axis)
  {
    int start = -1;
    int stop  = -1;
    for (int axisRank = 0; axisRank < getAxisCommSize(); ++axisRank)
    {
      if ((_globalRankBounds[axis][axisRank].start() <=
           _globalBounds[axis].start()) &&
          (_globalBounds[axis].start() <
           _globalRankBounds[axis][axisRank].stop()))
        start = axisRank;
      if ((_globalRankBounds[axis][axisRank].start() <
           _globalBounds[axis].stop()) &&
          (_globalBounds[axis].stop() <=
           _globalRankBounds[axis][axisRank].stop()))
        stop = axisRank;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      (start == -1 || stop == -1),
      InvalidArgument,
      "error computing axis rank slices");
    axisRankSlices.push_back(ConcreteSlice(start,stop));
  }
  // Construct the MDComm sub-communicator
  _mdComm = Teuchos::rcp(new MDComm(parent->_mdComm, axisRankSlices));

  // We now have a sub-communicator, and should only construct this
  // MDMap if this processor is on it.  If this processor is off the
  // communicator, then we should clear many of the data members.
  if (_mdComm->onSubcommunicator())
  {
    for (int axis = 0; axis < numDims; ++axis)
    {
      int axisRank = getAxisRank(axis);
      LocalOrd start = parent->_localBounds[axis].start();
      if (_globalBounds[axis].start() >
          _globalRankBounds[axis][axisRank].start())
      {
        start = _globalBounds[axis].start() -
          _globalRankBounds[axis][axisRank].start();
      }
      else
      {
        if (_globalRankBounds[axis][axisRank].start() -
            _globalBounds[axis].start() < _halos[axis][0])
        {
          _halos[axis][0] = _globalRankBounds[axis][axisRank].start() -
            _globalBounds[axis].start();
        }
      }
      LocalOrd stop = parent->_localBounds[axis].stop();
      if (_globalBounds[axis].stop() <
          _globalRankBounds[axis][axisRank].stop())
      {
        stop = _globalBounds[axis].stop() -
          _globalRankBounds[axis][axisRank].start();
      }
      else
      {
        if (_globalBounds[axis].stop() -
            _globalRankBounds[axis][axisRank].stop() < _halos[axis][1])
        {
          _halos[axis][1] = _globalBounds[axis].stop() -
            _globalRankBounds[axis][axisRank].stop();
        }
      }
      _localBounds.push_back(ConcreteSlice(start,stop));
      _localDims[axis] = stop - start;
      _localMin +=  start   * _localStrides[axis];
      _localMax += (stop-1) * _localStrides[axis];
    }
  }
  else
  {
    _localDims.clear();
    _localMin = 0;
    _localMax = 0;
    _localBounds.clear();
    _haloSizes.clear();
    _halos.clear();
    _ghostSizes.clear();
    _ghosts.clear();
    _localStrides.clear();
  }
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withGhosts)
    return _globalDims[axis];
  else
    return _globalDims[axis] - _ghosts[axis][0] - _ghosts[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalBounds(int axis,
                bool withGhosts) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withGhosts)
    return _globalBounds[axis];
  else
  {
    GlobalOrd start = _globalBounds[axis].start() + _ghosts[axis][0];
    GlobalOrd stop  = _globalBounds[axis].stop()  - _ghosts[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
LocalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalDim(int axis,
            bool withHalos) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withHalos)
    return _localDims[axis];
  else
    return _localDims[axis] - _halos[axis][0] - _halos[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalRankBounds(int axis,
                    bool withGhosts) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  int axisRank = getAxisRank(axis);
  if (withGhosts)
  {
    GlobalOrd start = _globalRankBounds[axis][axisRank].start();
    GlobalOrd stop  = _globalRankBounds[axis][axisRank].stop();
    if (getAxisRank(axis) == 0)
      start -= _ghosts[axis][0];
    if (getAxisRank(axis) == getAxisCommSize(axis)-1)
      stop += _ghosts[axis][1];
    return ConcreteSlice(start,stop);
  }
  else
    return _globalRankBounds[axis][axisRank];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalBounds(int axis,
               bool withHalos) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withHalos)
    return _localBounds[axis];
  else
  {
    LocalOrd start = _localBounds[axis].start() + _halos[axis][0];
    LocalOrd stop  = _localBounds[axis].stop()  - _halos[axis][1];
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _halos[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperHalo(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _halos[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getHaloSize(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _haloSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getLowerGhost(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _ghosts[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperGhost(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _ghosts[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getGhostSize(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((globalIndex < 0) || (globalIndex >= _globalMax)),
    RangeError,
    "invalid global index = " << globalIndex << " (global size = " <<
    _globalMax << ")");
#endif
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((localIndex < 0) || (localIndex >= _localMax)),
    RangeError,
    "invalid local index = " << localIndex << " (local size = " <<
    _localMax << ")");
#endif
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((localIndex < 0) || (localIndex >= _localMax)),
    RangeError,
    "invalid local index = " << localIndex << " (local size = " <<
    _localMax << ")");
#endif
  Teuchos::Array< LocalOrd > localAxisIndex = getLocalAxisIndex(localIndex);
  GlobalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    GlobalOrd globalAxisIndex = localAxisIndex[axis] +
      _globalRankBounds[axis][getAxisRank(axis)].start() - _halos[axis][0];
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (globalAxisIndex.size() != getNumDims()),
    InvalidArgument,
    "globalAxisIndex has " << globalAxisIndex.size() << " entries; expecting "
    << getNumDims());
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((globalAxisIndex[axis] < 0) ||
       (globalAxisIndex[axis] >= _globalDims[axis])),
      RangeError,
      "invalid globalAxisIndex[" << axis << "] = " << globalAxisIndex[axis] <<
      " (global dimension = " << _globalDims[axis] << ")");
  }
#endif
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((globalIndex < 0) || (globalIndex >= _globalMax)),
    RangeError,
    "invalid global index = " << globalIndex << " (global size = " <<
    _globalMax << ")");
#endif
  Teuchos::Array< GlobalOrd > globalAxisIndex =
    getGlobalAxisIndex(globalIndex);
  LocalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    LocalOrd localAxisIndex = globalAxisIndex[axis] -
      _globalRankBounds[axis][getAxisRank(axis)].start() + _halos[axis][0];
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
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (localAxisIndex.size() != getNumDims()),
    InvalidArgument,
    "localAxisIndex has " << localAxisIndex.size() << " entries; expecting "
    << getNumDims());
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((localAxisIndex[axis] < 0) ||
       (localAxisIndex[axis] >= _localDims[axis])),
      RangeError,
      "invalid localAxisIndex[" << axis << "] = " << localAxisIndex[axis] <<
      " (local dimension = " << _localDims[axis] << ")");
  }
#endif
  LocalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
    result += localAxisIndex[axis] * _localStrides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
void
MDMap< LocalOrd, GlobalOrd, Node >::computeBounds()
{
  // Initialization
  int numDims = getNumDims();

  // Decompose the multi-dimensional domain
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Get the communicator info for this axis
    int axisCommSize = getAxisCommSize(axis);
    for (int axisRank = 0; axisRank < axisCommSize; ++axisRank)
    {
      // First estimates assuming even division of global dimensions
      // by the number of processors along this axis, and ignoring
      // ghosts and halos.
      LocalOrd  localDim  = (_globalDims[axis] - 2*_ghostSizes[axis]) /
                            axisCommSize;
      GlobalOrd axisStart = axisRank * localDim;

      // Adjustments for non-zero remainder.  Compute the remainder
      // using the mod operator.  If the remainder is > 0, then add an
      // element to the appropriate number of processors with the
      // highest axis ranks.  Note that this is the opposite of the
      // standard Tpetra::Map constructor (which adds an elements to
      // the lowest processor ranks), and provides better balance for
      // finite differencing systems with staggered data location.
      GlobalOrd remainder = (_globalDims[axis] - 2*_ghostSizes[axis]) %
                            axisCommSize;
      if (axisCommSize - axisRank - 1 < remainder)
      {
        ++localDim;
        axisStart += (remainder - axisCommSize + axisRank);
      }

      // Global adjustment for ghost points
      axisStart += _ghostSizes[axis];

      // Compute and store the global axis bounds
      _globalRankBounds[axis].push_back(
        ConcreteSlice(axisStart, axisStart + localDim));

      // Set _localDims[axis] and _localBounds[axis] only if
      // axisRank equals the axis rank of this processor
      if (axisRank == getAxisRank(axis))
      {
        // Local adjustment for halos.  Note that _halos should
        // already be corrected to use ghost sizes for processors on
        // the boundaries
        _localDims[axis] = localDim + _halos[axis][0] + _halos[axis][1];

        // Compute and store the axis bounds
        _localBounds.push_back(ConcreteSlice(_localDims[axis]));
      }
    }
  }
}

}

#endif
