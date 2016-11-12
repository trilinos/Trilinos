// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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

// System includes
#include <algorithm>
#include <limits>
#include <sstream>

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_DefaultNode.hpp"
#include "Domi_Utils.hpp"
#include "Domi_getValidParameters.hpp"
#include "Domi_Slice.hpp"
#include "Domi_MDComm.hpp"
#include "Domi_MDArray.hpp"

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Tuple.hpp"

// Kokkos includes
#include "Kokkos_Core.hpp"

#ifdef HAVE_EPETRA
#include "Epetra_Map.h"
#endif

#ifdef HAVE_TPETRA
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#endif

namespace Domi
{

/** \brief Multi-dimensional map
 *
 * The <tt>MDMap</tt> class is intended to perform the functions of
 * Epetra or Tpetra maps, except that the data is multi-dimensional.
 * Rather than taking a <tt>Teuchos::Comm</tt> object in its
 * constructor, an <tt>MDMap</tt> object takes a <tt>Domi::MDComm</tt>
 * object.  An <tt>MDComm</tt> object is ordered in a structured way,
 * with processors assigned to each axis, whose product of lengths is
 * the total number of processors of its underlying communicator.
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
 * There are two distinct buffer padding concepts employed by
 * <tt>MDMaps</tt>, communication padding and boundary padding.
 * Communication padding refers to local extensions of the data
 * structure used to receive communication from neighboring
 * processors.  Boundary padding refers to points added to the
 * exterior of the global domain that are useful, for example, in
 * certain finite differencing algorithms.  They are treated
 * separately from the dimensions so that the decomposition will be
 * performed relative to the dimensions only, and not influenced by
 * the size of the boundary padding.
 *
 * From a global perspective, the only padding points that are visible
 * are boundary padding points.  Communication padding is internal and
 * essentially invisible from this perspective.  From a local
 * perspective, a local buffer can have padding on all sides.  Some of
 * this padding could be boundary padding, some could be communication
 * padding, and some could be both.
 *
 * The sizes of boundary and communication padding can be different
 * along each axis, and the sizes of boundary and communication
 * padding on any given axis can be different from each other.
 *
 * When you want to loop over the elements of a data structure
 * described by an <tt>MDMap</tt>, you can obtain the local looping
 * bounds along a given axis with the <tt>getLocalBounds(int
 * axis)</tt> method.  This returns a <tt>Slice</tt> object, whose
 * <tt>start()</tt> and <tt>stop()</tt> methods return the loop bounds
 * (<tt>stop()</tt> is non-inclusive).  This method also takes a
 * boolean flag indicating whether the bounds should include the
 * communication padding, which defaults to false.
 *
 * There are two methods for obtaining global loop bounds,
 * <tt>getGlobalBounds(int axis)</tt> and <tt>getGlobalRankBounds(int
 * axis)</tt>.  The first returns the global bounds for the full data
 * structure and the second returns the bounds on this processor.
 *
 * Note that all indexes start at zero, whether global or local,
 * whether unique 1D IDs or multidimensional indexes, whether boundary
 * or communication padding is present.  This allows negative indexes
 * to represent reverse indexing from the end of a dimension.
 */
template< class Node = DefaultNode::DefaultNodeType >
class MDMap
{
public:

  /** \name Constructors and destructor */
  //@{

  /** \brief Main constructor
   *
   * \param mdComm [in] an RCP of an MDComm (multi-dimensional
   *        communicator), on which this MDMap will be built.
   *
   * \param dimensions [in] the dimensions of the map along each axis.
   *        This array must be the same number of dimensions as the
   *        MDComm.
   *
   * \param commPad [in] the number of indexes in the communication
   *        padding along each axis.  If this array is less than the
   *        number of dimensions, unspecified communication padding
   *        will be set to zero.
   *
   * \param bndryPad [in] the number of indexes in the boundary
   *        padding along each axis.  If this array is less than the
   *        number of dimensions, unspecified unspecified boundary
   *        padding sizes will be set to zero.
   *
   * \param replicatedBoundary [in] An array of ints which are simple
   *        flags denoting whether each axis contains replicated
   *        boundary points (RBPs). RBPs pertain only to periodic
   *        axes. As an example, consider a map that describes the
   *        decomposition of a longitude coordinate. If the global
   *        lower boundary of this domain is longitude = 0 degrees,
   *        and the global upper boundary of this domain is longitude
   *        = 360 degrees, then this axis has a replicated
   *        boundary. Note that communication of the boundary padding
   *        will skip these boundary points, and it is up to the user
   *        to ensure that their values are equal or compatible.
   *
   * \param layout [in] the storage order of the map
   */
  MDMap(const Teuchos::RCP< const MDComm > mdComm,
        const Teuchos::ArrayView< const dim_type > & dimensions,
        const Teuchos::ArrayView< const int > & commPad =
          Teuchos::ArrayView< const int >(),
        const Teuchos::ArrayView< const int > & bndryPad =
          Teuchos::ArrayView< const int >(),
        const Teuchos::ArrayView< const int > & replicatedBoundary =
          Teuchos::ArrayView< const int >(),
        const Layout layout = DEFAULT_ORDER);

  /** \brief Constructor with ParameterList
   *
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scroling="no" frameborder="no">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   *
   * This constructor uses the Teuchos::DefaultComm
   */
  MDMap(Teuchos::ParameterList & plist);

  /** \brief Constructor with Teuchos::Comm and ParameterList
   *
   * \param teuchosComm [in] The Teuchos Communicator.  Note that an
   *        MDComm will be constructed from the information in plist.
   *
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   */
  MDMap(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
        Teuchos::ParameterList & plist);

  /** \brief Constructor with MDComm and ParameterList
   *
   * \param mdComm [in] an RCP of an MDComm (multi-dimensional
   *        communicator), on which this MDMap will be built.
   *
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   */
  MDMap(const Teuchos::RCP< const MDComm > mdComm,
        Teuchos::ParameterList & plist);

  /** \brief Constructor with global bounds for this processor
   *
   * \param mdComm [in] an RCP of an MDComm (multi-dimensional
   *        communicator), on which the MDMap will be built
   *
   * \param myGlobalBounds [in] an array of Slices, one for each axis,
   *        that represent the global indexes of the bounds on this
   *        processor, excluding padding
   *
   * \param replicatedBoundary [in] An array of ints which are simple
   *        flags denoting whether each axis contains replicated
   *        boundary points (RBPs). RBPs pertain only to periodic
   *        axes. As an example, consider a map that describes the
   *        decomposition of a longitude coordinate. If the global
   *        lower boundary of this domain is longitude = 0 degrees,
   *        and the global upper boundary of this domain is longitude
   *        = 360 degrees, then this axis has a replicated
   *        boundary. Note that communication of the boundary padding
   *        will skip these boundary points, and it is up to the user
   *        to ensure that their values are equal or compatible.
   *
   * \param layout [in] the storage order of the map
   */
  MDMap(const Teuchos::RCP< const MDComm > mdComm,
        const Teuchos::ArrayView< Slice > & myGlobalBounds,
        const Teuchos::ArrayView< padding_type > & padding =
          Teuchos::ArrayView< padding_type >(),
        const Teuchos::ArrayView< const int > & replicatedBoundary =
          Teuchos::ArrayView< const int >(),
        const Layout layout = DEFAULT_ORDER);

  /** \brief Copy constructor
   *
   * \param source [in] the source MDMap to be copied
   */
  MDMap(const MDMap< Node > & source);

  /** \brief Parent/single global ordinal sub-map constructor
   *
   * \param parent [in] an MDMap, from which this sub-map will be
   *        derived.
   *
   * \param axis [in] the axis to which the slice argument applies
   *
   * \param index [in] a global ordinal that defines the sub-map.
   *
   * This constructor will return an MDMap that is one dimension less
   * than the dimensions of the parent MDMap (unless the parent MDMap
   * is already one dimension).
   */
  MDMap(const MDMap< Node > & parent,
        int axis,
        dim_type index);

  /** \brief Parent/single slice sub-map constructor
   *
   * \param parent [in] an MDMap, from which this sub-map will be
   *        derived.
   *
   * \param axis [in] the axis to which the slice argument applies
   *
   * \param slice [in] a Slice of global axis indexes that
   *        defines the sub-map.  This slice must not include
   *        indexes from the boundary padding along each axis.
   *
   * \param bndryPad [in] The size of the boundary padding along the
   *        altered axis of the new sub-map.  This may include indexes
   *        from the boundary padding of the parent MDMap, but it does
   *        not have to.
   *
   * This constructor will return an MDMap that is the same number of
   * dimensions as the parent MDMap, but with a smaller dimension
   * along the given axis (unless the given Slice represents the
   * entire axis).
   */
  MDMap(const MDMap< Node > & parent,
        int axis,
        const Slice & slice,
        int bndryPad = 0);

  /** \brief Parent/array of slices sub-map constructor
   *
   * \param parent [in] an MDMap, from which this sub-map will be
   *        derived.
   *
   * \param slices [in] an array of Slices of global axis indexes that
   *        defines the sub-map.  These slices must not include
   *        indexes from the boundary padding along each axis.
   *
   * \param bndryPad [in] The boundary padding of the new sub-map.
   *        These may include indexes from the boundary padding of the
   *        parent MDMap, but they do not have to.
   */
  MDMap(const MDMap< Node > & parent,
        const Teuchos::ArrayView< Slice > & slices,
        const Teuchos::ArrayView< int > & bndryPad =
          Teuchos::ArrayView< int >());

  /** \brief Destructor
   */
  ~MDMap();

  //@}

  /** \name MDMap standard operators */
  //@{

  /** \brief Assignment operator
   *
   * \param source [in] source MDComm
   */
  MDMap< Node > & operator=(const MDMap< Node > & source);

  //@}

  /** \name MDComm accessor and pass-through methods */
  //@{

  /** \brief Access the underlying MDComm object
   *
   * Return a reference counted pointer to the Domi::MDComm object
   * upon which this MDMap is built.
   */
  inline Teuchos::RCP< const MDComm > getMDComm() const;

  /** \brief Query whether this processor is on the sub-communicator
   *
   * Sub-communicators are formed when a parent MDComm is sliced by
   * using the (parent,ordinal) or (parent,slices) constructor.  For a
   * full communicator, this method will always return true.
   */
  inline bool onSubcommunicator() const;

  /** \brief Get the Teuchos communicator
   *
   * Note that if the communicator is not a full communicator, i.e. a
   * sub-communicator, that the underlying Comm pointer may be NULL,
   * depending on this processor's rank.
   */
  inline Teuchos::RCP< const Teuchos::Comm< int > > getTeuchosComm() const;

  /** \brief Get the number of dimensions
   *
   * This method will return 0 if the communicator is a
   * sub-communicator and this processor does not belong to the
   * sub-communicator.
   */
  inline int numDims() const;

  /** \brief Get the communicator sizes along each axis
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  inline Teuchos::Array< int > getCommDims() const;

  /** \brief Get the communicator size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  inline int getCommDim(int axis) const;

  /** \brief Return the periodic flag for the given axis.
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  inline bool isPeriodic(int axis) const;

  /** \brief Get the axis rank of this processor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
 inline  int getCommIndex(int axis) const;

  /** \brief Get the rank of the lower neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   *
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is zero, then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is the highest axis rank processor along
   * this axis, then the returned lower neighbor will be zero.
   */
  inline int getLowerNeighbor(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   *
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is the highest axis rank processor along
   * this axis, then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is zero, then the returned lower
   * neighbor will be the highest axis rank processor along this axis.
   */
  inline int getUpperNeighbor(int axis) const;

  //@}

  /** \name Dimensions and looping bounds */
  //@{

  /** \brief Get an array of the the global dimensions, including
   *         boundary padding
   */
  Teuchos::Array< dim_type > getGlobalDims() const;

  /** \brief Get the global dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the dimension should
   *        include boundary padding or not
   */
  dim_type getGlobalDim(int axis,
                        bool withBndryPad=false) const;

  /** \brief Get the bounds of the global problem
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the bounds should
   *        include boundary padding or not
   */
  Slice getGlobalBounds(int axis,
                        bool withBndryPad=false) const;

  /** \brief Get the global loop bounds along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the dimension should
   *        include boundary padding or not
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  Slice getGlobalRankBounds(int axis,
                            bool withBndryPad=false) const;

  /** \brief Get an array of the local dimensions, including padding
   */
  Teuchos::Array< dim_type > getLocalDims() const;

  /** \brief Get the local dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withPad [in] specify whether the dimension should include
   *        padding or not
   */
  dim_type getLocalDim(int axis,
                       bool withPad=false) const;

  /** \brief Get the local loop bounds along every axis
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.  For this method, padding is included in the bounds.
   */
  Teuchos::ArrayView< const Slice > getLocalBounds() const;

  /** \brief Get the local loop bounds along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withPad [in] specify whether the dimension should include
   *        padding or not
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  Slice getLocalBounds(int axis,
                       bool withPad=false) const;

  /** \brief Get the local interior loop bounds along the specified
   *         axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Local interior loop bounds are the same as local loop bounds
   * without padding, except that for non-periodic axes the global end
   * points of the given axis are excluded. For periodic axes, the
   * local interior loop bounds are exactly the same as local loop
   * bounds without padding.
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  Slice getLocalInteriorBounds(int axis) const;

  //@}

  /** \name Storage order, communication and boundary padding */
  //@{

  /** \brief Return true if there is any padding stored locally
   *
   * Note that it is not as simple as whether there were communication
   * padding specified in the constructor.  Suppose non-zero
   * communication padding was specified, but the problem is on one
   * processor and no boundary padding was specified: this method will
   * return false.  Similarly, if no communication padding was
   * specified but boundary padding was, then boundary processors will
   * have padding and this method will return true.
   */
  bool hasPadding() const;

  /** \brief Get the size of the lower padding along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that the returned padding can be either communication
   * padding or boundary padding as appropriate.
   */
  int getLowerPadSize(int axis) const;

  /** \brief Get the size of the upper padding along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that the returned padding can be either communication
   * padding or boundary padding as appropriate.
   */
  int getUpperPadSize(int axis) const;

  /** \brief Get the communication padding size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the communication padding along the
   * given axis at the time of construction, regardless of the
   * processor's position relative to the domain boundary.
   */
  int getCommPadSize(int axis) const;

  /** \brief Get the size of the lower boundary padding along the
   *         given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getLowerBndryPad(int axis) const;

  /** \brief Get the size of the upper boundary padding along the
   *         given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getUpperBndryPad(int axis) const;

  /** \brief Get the boundary padding sizes along each axis
   *
   * This returns an array of the boundary padding along each axis at
   * the time of construction, regardless of whether a subsequent
   * sub-map reduced these values.
   */
  Teuchos::Array< int > getBndryPadSizes() const;

  /** \brief Get the boundary padding size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the boundary padding along the given
   * axis at the time of construction, regardless of whether a sub-map
   * reduced these values.
   */
  int getBndryPadSize(int axis) const;

  /* \brief Return whether given local index is in the padding
   *
   * \param index [in] An array of indexes ((i) for 1D, (i,j) for 2D,
   *        (i,j,k) for 3D, etc)
   */
  bool isPad(const Teuchos::ArrayView< dim_type > & index) const;

  /* \brief Return whether given local index is in the communication
   *        padding
   *
   * \param index [in] An array of indexes ((i) for 1D, (i,j) for 2D,
   *        (i,j,k) for 3D, etc)
   */
  bool isCommPad(const Teuchos::ArrayView< dim_type > & index) const;

  /* \brief Return whether given local index is in the boundary
   *        padding
   *
   * \param index [in] An array of indexes ((i) for 1D, (i,j) for 2D,
   *        (i,j,k) for 3D, etc)
   */
  bool isBndryPad(const Teuchos::ArrayView< dim_type > & index) const;

  /** \brief Return whether the given axis supports a replicated boundary
   */
  bool isReplicatedBoundary(int axis) const;

  /** \brief Get the storage order
   */
  Layout getLayout() const;

  /** \brief Get the Kokkos node
   */
  Teuchos::RCP< Node > getNode() const;

  //@}

  /** \name Conversions to other Maps */
  //@{

  /** \brief Return an array of axis maps
   *
   * An axis map is a 1D map along a given axis.
   */
  Teuchos::ArrayView< Teuchos::RCP< const Domi::MDMap< Node > > >
  getAxisMaps() const;

  /** \brief Return an axis map for the given axis
   *
   * An axis map is a 1D map along a given axis.
   */
  Teuchos::RCP< const Domi::MDMap< Node > > getAxisMap(int axis) const;

  /** \brief Return an RCP to a new MDMap that is a simple
   *         augmentation of this MDMap
   *
   * \param leadingDim [in] If this value is greater than 0, then add
   *        a dimension to the new MDMap that comes before the
   *        existing dimensions and has this many degrees of freedom.
   *
   * \param trailingDim [in] If this value is greater than 0, then
   *        add a dimension to the new MDMap that comes after the
   *        existing dimensions and has this many degrees of freedom.
   *
   * Note that the new dimension(s) will not be distributed; i.e. the
   * commDim for the new leading dimension (if requested) will be 1
   * and the commDim for the trailing dimension (if requested) will be
   * 1.
   */
  Teuchos::RCP< const MDMap< Node > >
  getAugmentedMDMap(const dim_type leadingDim,
                    const dim_type trailingDim=0) const;

#ifdef HAVE_EPETRA

  /** \brief Return an RCP to an Epetra_Map that is equivalent to this
   *         MDMap
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  Teuchos::RCP< const Epetra_Map >
  getEpetraMap(bool withCommPad=true) const;

  /** \brief Return an RCP to an Epetra_Map that represents the
   *         decomposition of this MDMap along the given axis
   *
   * \param axis [in] the requested axis
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  Teuchos::RCP< const Epetra_Map >
  getEpetraAxisMap(int axis,
                   bool withCommPad=true) const;

#endif

#ifdef HAVE_TPETRA

  /** \brief Return an RCP to a Tpetra::Map that is equivalent to this
   *         MDMap, specifying the new LocalOrdinal, GlobalOrdinal and
   *         Node types
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  template< class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node2 = Node >
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node2 > >
  getTpetraMap(bool withCommPad=true) const;

  /** \brief Return an RCP to a Tpetra::Map that represents the
   *         decomposition of this MDMap along the given axis,
   *         specifying the new LocalOrdinal, GlobalOrdinal and Node
   *         types
   *
   * \param axis [in] the requested axis
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  template< class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node2 = Node >
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node2 > >
  getTpetraAxisMap(int axis,
                   bool withCommPad=true) const;

#endif

  //@}

  /** \name Global/local index and ID conversion methods */
  //@{

  /** \brief Convert a global ID to a global index
   *
   * \param globalID [in] a unique 1D global identifier
   */
  Teuchos::Array< dim_type >
  getGlobalIndex(size_type globalID) const;

  /** \brief Convert a local ID to a local index
   *
   * \param localID [in] a unique 1D local identifier
   */
  Teuchos::Array< dim_type >
  getLocalIndex(size_type localID) const;

  /** \brief Convert a local ID to a global ID
   *
   * \param localID [in] a unique 1D local identifier
   */
  size_type getGlobalID(size_type localID) const;

  /** \brief convert a global index to a global ID
   *
   * \param globalIndex [in] a global index
   */
  size_type
  getGlobalID(const Teuchos::ArrayView< dim_type > & globalIndex) const;

  /** \brief Convert a global ID to a local ID
   *
   * \param globalID [in] a unique 1D global identifier
   *
   * This method can throw a Domi::RangeError if the global ID is
   * not on the current processor.
   */
  size_type getLocalID(size_type globalID) const;

  /** \brief Convert a local index to a local ID
   *
   * \param localIndex [in] a local index
   */
  size_type
  getLocalID(const Teuchos::ArrayView< dim_type > & localIndex) const;

  //@}

  /** \name Boolean comparison methods */
  //@{

  /** \brief True if two MDMaps are "compatible"
   *
   * \param mdMap [in] MDMap to compare against
   *
   * Two MDMaps are considered "compatible" if all of the following
   * are true:
   * <ol>
   *   <li> The commDims of their underlying MDComms are
   *        identical.</li>
   *   <li> Their global dimensions are identical.</li>
   *   <li> Their local dimensions, not including padding, are
   *        identical.</li>
   * <ol>
   */
  bool isCompatible(const MDMap< Node > & mdMap) const;

  /** \brief True if two MDMaps are "identical"
   *
   * \param mdMap [in] MDMap to compare against
   *
   * \param verbose [in] set to one to print why MDMaps are not the
   *        same
   *
   * Two MDMaps are considered "identical" if all of the following
   * are true:
   * <ol>
   *   <li> They are compatible.</li>
   *   <li> Their underlying communicators are <i>congruent</i> (have
   *        the same number of processes, in the same order: this
   *        corresponds to the \c MPI_IDENT or \c MPI_CONGRUENT return
   *        values of MPI_Comm_compare).</li>
   *   <li> Both are distributed or not distributed.</li>
   *   <li> Their global bounds identical on each process.</li>
   *   <li> Their local dimensions, including padding, are
   *        identical.</li>
   * <ol>
   */
  bool isSameAs(const MDMap< Node > & mdMap,
                const int verbose = 0) const;

  /** \brief True if there are no stride gaps on any processor
   *
   * An MDMap constructed from a communicator and dimensions will
   * always be contiguous.  An MDMap that is a slice of a parent MDMap
   * will generally be non-contiguous, with some exceptions.  There
   * are cases where some local data is contiguous and some is not,
   * but this method returns True only if all processes' local data is
   * contiguous.
   */
  bool isContiguous() const;

  //@}

private:

  // A private method for computing the bounds and local dimensions,
  // after the global dimensions, communication and boundary padding
  // have been properly assigned.
  void computeBounds();

  // The underlying multi-dimensional communicator.
  Teuchos::RCP< const MDComm > _mdComm;

  // The size of the global dimensions along each axis.  This includes
  // the values of the boundary padding.
  Teuchos::Array< dim_type > _globalDims;

  // Store the start and stop indexes of the global problem along each
  // axis.  This includes the values of the boundary padding.
  Teuchos::Array< Slice > _globalBounds;

  // A double array of Slices that stores the starting and stopping
  // indexes for each axis processor along each axis.  This data
  // structure allows any processor to know the map bounds on any
  // other processor.  These bounds do NOT include the boundary or
  // communication padding.  The outer array will have a size equal to
  // the number of dimensions.  Each inner array will have a size
  // equal to the number of axis processors along that axis.
  Teuchos::Array< Teuchos::Array< Slice > > _globalRankBounds;

  // The global stride between adjacent IDs.  This quantity is
  // "virtual" as it does not describe actual memory layout, but it is
  // useful for index conversion algorithms.
  Teuchos::Array< size_type > _globalStrides;

  // The minumum global ID of the global data structure, including
  // boundary padding.  This will only be non-zero on a sub-map.
  size_type _globalMin;

  // The maximum global ID of the global data structure, including
  // boundary padding.
  size_type _globalMax;

  // The size of the local dimensions along each axis.  This includes
  // the values of the padding.
  Teuchos::Array< dim_type > _localDims;

  // The local loop bounds along each axis, stored as an array of
  // Slices.  These bounds DO include the padding.  By definition, the
  // start() attribute of these Slices will all be zero by definition.
  // This may seem wasteful (we could store an array of dim_type
  // rather than an array of Slice), but this storage is convenient
  // when it comes to the getLocalBounds() method.
  Teuchos::Array< Slice > _localBounds;

  // The local stride between adjacent elements in memory.
  Teuchos::Array< size_type > _localStrides;

  // The minimum local ID of the local data structure, including
  // padding.
  size_type _localMin;

  // The maximum local ID of the local data structure, including
  // padding.
  size_type _localMax;

  // The communication padding that was specified at construction, one
  // value along each axis.
  Teuchos::Array< int > _commPadSizes;

  // The actual padding stored on this processor along each axis.  The
  // padding can be either communication padding or boundary padding
  // based on the processor position on the boundary of a domain.
  Teuchos::Array< padding_type > _pad;

  // The size of the boundary padding along each axis.
  Teuchos::Array< int > _bndryPadSizes;

  // The actual boundary padding stored along each axis.  For a full
  // communicator, the lower and upper boundary padding are both the
  // same as the corresponding value in _bndryPadSizes.  However, the
  // introduction of sub-maps creates the possibility of different
  // upper and lower boundary padding values.  If an axis is periodic,
  // then these values will be set to the communication padding.
  Teuchos::Array< padding_type > _bndryPad;

  // An array of flags, one for each axis, that specifies whether the
  // endpoints of a periodic boundary are replicated (1 or true) or
  // unique (0 or false). These flags only have meaning for axes that
  // are periodic. The arrray type is int, beacause Array< bool > is
  // specialized and sometimes difficult to work with.
  Teuchos::Array< int > _replicatedBoundary;

  // The storage order
  Layout _layout;

  // An array of axis maps that correspond to the full MDMap.  An axis
  // map describes the map along a single axis. This member is mutable
  // because it is logically const but does not get constructed until
  // requested by the user.
  mutable Teuchos::Array< Teuchos::RCP< const MDMap< Node > > > _axisMaps;

  // Kokkos node
  Teuchos::RCP< Node > _node;

#ifdef HAVE_EPETRA
  // An RCP pointing to an Epetra_Map that is equivalent to this
  // MDMap, including communication padding.  It is mutable because we
  // do not construct it until it asked for by a get method that is
  // const.
  mutable Teuchos::RCP< const Epetra_Map > _epetraMap;

  // An RCP pointing to an Epetra_Map that is equivalent to this
  // MDMap, excluding communication padding.  The returned Epetra_Map
  // thus indicates what IDs are owned by this processor.  It is
  // mutable because we do not construct it until it asked for by a
  // get method that is const.
  mutable Teuchos::RCP< const Epetra_Map > _epetraOwnMap;

  // An ArrayRCP that stores Epetra_Maps that represent the
  // distribution of indexes along each axis, including communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array< Teuchos::RCP< const Epetra_Map > > _epetraAxisMaps;

  // An ArrayRCP that stores Epetra_Maps that represent the
  // distribution of indexes along each axis, excluding communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array< Teuchos::RCP< const Epetra_Map > > _epetraAxisOwnMaps;
#endif

};

////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const Teuchos::RCP< const MDComm > mdComm,
      const Teuchos::ArrayView< const dim_type > & dimensions,
      const Teuchos::ArrayView< const int > & commPad,
      const Teuchos::ArrayView< const int > & bndryPad,
      const Teuchos::ArrayView< const int > & replicatedBoundary,
      const Layout layout) :
  _mdComm(mdComm),
  _globalDims(mdComm->numDims()),
  _globalBounds(),
  _globalRankBounds(mdComm->numDims()),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(mdComm->numDims(), 0),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(),
  _commPadSizes(mdComm->numDims(), 0),
  _pad(),
  _bndryPadSizes(mdComm->numDims(), 0),
  _bndryPad(),
  _replicatedBoundary(createArrayOfInts(mdComm->numDims(),
                                        replicatedBoundary)),
  _layout(layout)
{
  // Temporarily store the number of dimensions
  int numDims = mdComm->numDims();

  // Check the global dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");

  // Copy the communication and boundary padding sizes, compute the
  // global dimensions and bounds, and the actual padding
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPad.size() ) _commPadSizes[ axis] = commPad[ axis];
    if (axis < bndryPad.size()) _bndryPadSizes[axis] = bndryPad[axis];
    if (_mdComm->isPeriodic(axis))
      _bndryPad.push_back(Teuchos::tuple(_commPadSizes[axis],
                                         _commPadSizes[axis]));
    else
      _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                         _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + _bndryPad[axis][0] +
                        _bndryPad[axis][1];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
    int lower, upper;
    if (getLowerNeighbor(axis) == -1)
      lower = _bndryPadSizes[axis];
    else
      lower = _commPadSizes[axis];
    if (getUpperNeighbor(axis) == -1)
      upper = _bndryPadSizes[axis];
    else
      upper = _commPadSizes[axis];
    _pad.push_back(Teuchos::tuple(lower, upper));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Compute _globalRankBounds, _localBounds, and _localDims.  Then
  // compute the local size
  computeBounds();
  _localMax = computeSize(_localDims());

  // Compute the global and local strides
  _globalStrides = computeStrides< size_type, dim_type >(_globalDims, _layout);
  _localStrides  = computeStrides< size_type, dim_type >(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(Teuchos::ParameterList & plist) :
  _mdComm(Teuchos::rcp(new MDComm(plist))),
  _globalDims(),
  _globalBounds(),
  _globalRankBounds(),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(0),
  _commPadSizes(),
  _pad(),
  _bndryPadSizes(),
  _bndryPad(),
  _replicatedBoundary(),
  _layout()
{
  // Note that the call to the MDComm constructor in the constructor
  // initialization list will validate the ParameterList, so we don't
  // have to do it again here.

  // Temporarily store the number of dimensions
  int numDims = _mdComm->numDims();

  // Check the global dimensions
  Teuchos::Array< dim_type > dimensions =
    plist.get("dimensions", Teuchos::Array< dim_type >());
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");

  // Initialize _bndryPadSizes, _commPadSizes and _globalDims from the
  // ParameterList
  int commPad  = plist.get("communication pad size", int(0));
  int bndryPad = plist.get("boundary pad size"     , int(0));
  Teuchos::Array< int > commPads  =
    plist.get("communication pad sizes", Teuchos::Array< int >());
  Teuchos::Array< int > bndryPads =
    plist.get("boundary pad sizes"     , Teuchos::Array< int >());
  _commPadSizes.resize( numDims);
  _bndryPadSizes.resize(numDims);
  _globalDims.resize(   numDims);

  // Copy the communication and boundary padding sizes, compute the
  // global dimensions and bounds, and the actual padding
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPads.size() ) _commPadSizes[ axis] = commPads[ axis];
    else                         _commPadSizes[ axis] = commPad;
    if (axis < bndryPads.size()) _bndryPadSizes[axis] = bndryPads[axis];
    else                         _bndryPadSizes[axis] = bndryPad;
    if (_mdComm->isPeriodic(axis))
      _bndryPad.push_back(Teuchos::tuple(_commPadSizes[axis],
                                         _commPadSizes[axis]));
    else
      _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                         _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + _bndryPad[axis][0] +
                        _bndryPad[axis][1];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
    int lower, upper;
    if (getLowerNeighbor(axis) == -1)
      lower = _bndryPadSizes[axis];
    else
      lower = _commPadSizes[axis];
    if (getUpperNeighbor(axis) == -1)
      upper = _bndryPadSizes[axis];
    else
      upper = _commPadSizes[axis];
    _pad.push_back(Teuchos::tuple(lower, upper));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  _globalRankBounds.resize(numDims);
  _localDims.resize(numDims);
  computeBounds();
  _localMax = computeSize(_localDims());

  // Set the replicated boundary flags along each axis
  Teuchos::Array< int > repBndry = plist.get("replicated boundary",
                                             Teuchos::Array< int >());
  _replicatedBoundary = createArrayOfInts(numDims, repBndry);

  // Set the layout
  std::string layout = plist.get("layout", "DEFAULT");
  std::transform(layout.begin(), layout.end(), layout.begin(), ::toupper);
  if (layout == "C ORDER")
    _layout = C_ORDER;
  else if (layout == "FORTRAN ORDER")
    _layout = FORTRAN_ORDER;
  else if (layout == "ROW MAJOR")
    _layout = ROW_MAJOR;
  else if (layout == "COLUMN MAJOR")
    _layout = COLUMN_MAJOR;
  else if (layout == "LAST INDEX FASTEST")
    _layout = LAST_INDEX_FASTEST;
  else if (layout == "FIRST INDEX FASTEST")
    _layout = FIRST_INDEX_FASTEST;
  else
    _layout = DEFAULT_ORDER;

  // Compute the strides
  _globalStrides = computeStrides< size_type, dim_type >(_globalDims, _layout);
  _localStrides  = computeStrides< size_type, dim_type >(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
      Teuchos::ParameterList & plist) :
  _mdComm(Teuchos::rcp(new MDComm(teuchosComm, plist))),
  _globalDims(),
  _globalBounds(),
  _globalRankBounds(),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(0),
  _commPadSizes(),
  _pad(),
  _bndryPadSizes(),
  _bndryPad(),
  _replicatedBoundary(),
  _layout()
{
  // Note that the call to the MDComm constructor in the constructor
  // initialization list will validate the ParameterList, so we don't
  // have to do it again here.

  // Temporarily store the number of dimensions
  int numDims = _mdComm->numDims();

  // Check the global dimensions
  Teuchos::Array< dim_type > dimensions =
    plist.get("dimensions", Teuchos::Array< dim_type >());
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");

  // Initialize _bndryPadSizes, _commPadSizes and _globalDims from the
  // ParameterList
  int commPad  = plist.get("communication pad size", int(0));
  int bndryPad = plist.get("boundary pad size"     , int(0));
  Teuchos::Array< int > commPads  =
    plist.get("communication pad sizes", Teuchos::Array< int >());
  Teuchos::Array< int > bndryPads =
    plist.get("boundary pad sizes"     , Teuchos::Array< int >());
  _commPadSizes.resize( numDims);
  _bndryPadSizes.resize(numDims);
  _globalDims.resize(   numDims);

  // Copy the communication and boundary padding sizes, compute the
  // global dimensions and bounds, and the actual padding
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPads.size() ) _commPadSizes[ axis] = commPads[ axis];
    else                         _commPadSizes[ axis] = commPad;
    if (axis < bndryPads.size()) _bndryPadSizes[axis] = bndryPads[axis];
    else                         _bndryPadSizes[axis] = bndryPad;
    if (_mdComm->isPeriodic(axis))
      _bndryPad.push_back(Teuchos::tuple(_commPadSizes[axis],
                                         _commPadSizes[axis]));
    else
      _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                         _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + _bndryPad[axis][0] +
                        _bndryPad[axis][1];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
    int lower, upper;
    if (getLowerNeighbor(axis) == -1)
      lower = _bndryPadSizes[axis];
    else
      lower = _commPadSizes[axis];
    if (getUpperNeighbor(axis) == -1)
      upper = _bndryPadSizes[axis];
    else
      upper = _commPadSizes[axis];
    _pad.push_back(Teuchos::tuple(lower, upper));
  }
  // std::cout << "MDMap constructor: _commPadSizes = " << _commPadSizes
  //           << ", bndryPadSizes = " << _bndryPadSizes << ", _bndryPad = "
  //           << _bndryPad << ", _pad = " << _pad << ", _globalDims = "
  //           << _globalDims << ", _globalBounds = " << _globalBounds
  //           << std::endl;

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  _globalRankBounds.resize(numDims);
  _localDims.resize(numDims);
  computeBounds();
  _localMax = computeSize(_localDims());

  // Set the replicated boundary flags along each axis
  Teuchos::Array< int > repBndry = plist.get("replicated boundary",
                                             Teuchos::Array< int >());
  _replicatedBoundary = createArrayOfInts(numDims, repBndry);

  // Set the layout
  std::string layout = plist.get("layout", "DEFAULT");
  std::transform(layout.begin(), layout.end(), layout.begin(), ::toupper);
  if (layout == "C ORDER")
    _layout = C_ORDER;
  else if (layout == "FORTRAN ORDER")
    _layout = FORTRAN_ORDER;
  else if (layout == "ROW MAJOR")
    _layout = ROW_MAJOR;
  else if (layout == "COLUMN MAJOR")
    _layout = COLUMN_MAJOR;
  else if (layout == "LAST INDEX FASTEST")
    _layout = LAST_INDEX_FASTEST;
  else if (layout == "FIRST INDEX FASTEST")
    _layout = FIRST_INDEX_FASTEST;
  else
    _layout = DEFAULT_ORDER;

  // Compute the strides
  _globalStrides = computeStrides< size_type, dim_type >(_globalDims, _layout);
  _localStrides  = computeStrides< size_type, dim_type >(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(Teuchos::RCP< const MDComm > mdComm,
      Teuchos::ParameterList & plist) :
  _mdComm(mdComm),
  _globalDims(mdComm->numDims()),
  _globalBounds(),
  _globalRankBounds(mdComm->numDims()),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(mdComm->numDims(), 0),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(),
  _commPadSizes(mdComm->numDims(), 0),
  _pad(),
  _bndryPadSizes(mdComm->numDims(), 0),
  _bndryPad(),
  _replicatedBoundary(),
  _layout()
{
  // Validate the ParameterList
  plist.validateParameters(*getValidParameters());

  // Temporarily store the number of dimensions
  int numDims = _mdComm->numDims();

  // Check the global dimensions
  Teuchos::Array< dim_type > dimensions =
    plist.get("dimensions", Teuchos::Array< dim_type >());
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Number of dimensions does not match MDComm number of dimensions");

  // Initialize _bndryPadSizes, _commPadSizes and _globalDims from the
  // ParameterList
  int commPad  = plist.get("communication pad size", int(0));
  int bndryPad = plist.get("boundary pad size"     , int(0));
  Teuchos::Array< int > commPads  =
    plist.get("communication pad sizes", Teuchos::Array< int >());
  Teuchos::Array< int > bndryPads =
    plist.get("boundary pad sizes"     , Teuchos::Array< int >());
  _commPadSizes.resize( numDims);
  _bndryPadSizes.resize(numDims);
  _globalDims.resize(   numDims);

  // Copy the communication and boundary padding sizes, compute the
  // global dimensions and bounds, and the actual padding
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPads.size() ) _commPadSizes[ axis] = commPads[ axis];
    else                         _commPadSizes[ axis] = commPad;
    if (axis < bndryPads.size()) _bndryPadSizes[axis] = bndryPads[axis];
    else                         _bndryPadSizes[axis] = bndryPad;
    if (_mdComm->isPeriodic(axis))
      _bndryPad.push_back(Teuchos::tuple(_commPadSizes[axis],
                                         _commPadSizes[axis]));
    else
      _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                         _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + _bndryPad[axis][0] +
                        _bndryPad[axis][1];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
    int lower, upper;
    if (getLowerNeighbor(axis) == -1)
      lower = _bndryPadSizes[axis];
    else
      lower = _commPadSizes[axis];
    if (getUpperNeighbor(axis) == -1)
      upper = _bndryPadSizes[axis];
    else
      upper = _commPadSizes[axis];
    _pad.push_back(Teuchos::tuple(lower, upper));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  _globalRankBounds.resize(numDims);
  _localDims.resize(numDims);
  computeBounds();
  _localMax = computeSize(_localDims());

  // Set the replicated boundary flags along each axis
  Teuchos::Array< int > repBndry = plist.get("replicated boundary",
                                             Teuchos::Array< int >());
  _replicatedBoundary = createArrayOfInts(numDims, repBndry);

  // Set the layout
  std::string layout = plist.get("layout", "DEFAULT");
  std::transform(layout.begin(), layout.end(), layout.begin(), ::toupper);
  if (layout == "C ORDER")
    _layout = C_ORDER;
  else if (layout == "FORTRAN ORDER")
    _layout = FORTRAN_ORDER;
  else if (layout == "ROW MAJOR")
    _layout = ROW_MAJOR;
  else if (layout == "COLUMN MAJOR")
    _layout = COLUMN_MAJOR;
  else if (layout == "LAST INDEX FASTEST")
    _layout = LAST_INDEX_FASTEST;
  else if (layout == "FIRST INDEX FASTEST")
    _layout = FIRST_INDEX_FASTEST;
  else
    _layout = DEFAULT_ORDER;

  // Compute the strides
  _globalStrides = computeStrides< size_type, dim_type >(_globalDims, _layout);
  _localStrides  = computeStrides< size_type, dim_type >(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const Teuchos::RCP< const MDComm > mdComm,
      const Teuchos::ArrayView< Slice > & myGlobalBounds,
      const Teuchos::ArrayView< padding_type > & padding,
      const Teuchos::ArrayView< const int > & replicatedBoundary,
      const Layout layout) :
  _mdComm(mdComm),
  _globalDims(mdComm->numDims()),
  _globalBounds(mdComm->numDims()),
  _globalRankBounds(mdComm->numDims()),
  _globalStrides(mdComm->numDims()),
  _globalMin(0),
  _globalMax(0),
  _localDims(mdComm->numDims(), 0),
  _localBounds(mdComm->numDims()),
  _localStrides(mdComm->numDims()),
  _localMin(0),
  _localMax(0),
  _commPadSizes(mdComm->numDims(), 0),
  _pad(mdComm->numDims(), Teuchos::tuple(0,0)),
  _bndryPadSizes(mdComm->numDims(), 0),
  _bndryPad(mdComm->numDims()),
  _replicatedBoundary(createArrayOfInts(mdComm->numDims(),
                                        replicatedBoundary)),
  _layout(layout)
{
  // Check that myGlobalBounds is the correct size
  int numDims = _mdComm->numDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    myGlobalBounds.size() < numDims,
    InvalidArgument,
    "MDMap: myGlobalBounds too small");

  // Copy the padding to _pad
  int maxAxis = std::min(numDims, (int)padding.size());
  for (int axis = 0; axis < maxAxis; ++axis)
    _pad[axis] = padding[axis];

  // All of the required info for the MDMap is contained in the
  // myGlobalBounds and padding arguments, but it is distributed.  We
  // will do a gather so that each processor has the global bounds
  // data from each other processor.

  // Resize _globalRankBounds
  for (int axis = 0; axis < numDims; ++axis)
    _globalRankBounds[axis].resize(_mdComm->getCommDim(axis));

  // Allocate and initialize the communication buffers
  int numProc = _mdComm->getTeuchosComm()->getSize();
  MDArray< dim_type > sendBuffer(Teuchos::tuple(5,numDims),
                                 FIRST_INDEX_FASTEST);
  MDArray< dim_type > recvBuffer(Teuchos::tuple(5,numDims,numProc),
                                 FIRST_INDEX_FASTEST);
  for (int axis = 0; axis < numDims; ++axis)
  {
    sendBuffer(0,axis) = mdComm->getCommIndex(axis);
    sendBuffer(1,axis) = myGlobalBounds[axis].start();
    sendBuffer(2,axis) = myGlobalBounds[axis].stop();
    sendBuffer(3,axis) = _pad[axis][0];
    sendBuffer(4,axis) = _pad[axis][1];
  }

  // Perform the gather all
  Teuchos::gatherAll(*(_mdComm->getTeuchosComm()),
                     (int)sendBuffer.size(),
                     sendBuffer.getRawPtr(),
                     (int)recvBuffer.size(),
                     recvBuffer.getRawPtr());

  // Extract _globalRankBounds and _bndryPad.  Because of the
  // structure, there will be duplicate Slices and padding, and we
  // will check to make sure they are the expected equivalent values.
  for (int axis = 0; axis < numDims; ++axis)
  {
    for (int commIndex = 0; commIndex < _mdComm->getCommDim(axis); ++commIndex)
    {
      Slice bounds;
      padding_type pad;
      for (int rank = 0; rank < numProc; ++rank)
      {
        if (recvBuffer(0,axis,rank) == commIndex)
        {
          dim_type start = recvBuffer(1,axis,rank);
          dim_type stop  = recvBuffer(2,axis,rank);
          int      loPad = recvBuffer(3,axis,rank);
          int      hiPad = recvBuffer(4,axis,rank);
          if (bounds.stop() == Slice::Default)
          {
            bounds = Slice(start, stop);
            pad[0] = loPad;
            pad[1] = hiPad;
          }
          else
          {
            TEUCHOS_TEST_FOR_EXCEPTION(
              (bounds.start() != start) || (bounds.stop() != stop),
              BoundsError,
              "Global rank bounds mismatch: bounds = " << bounds <<
              ", (start,stop) = (" << start << "," << stop << ")");
            TEUCHOS_TEST_FOR_EXCEPTION(
              (pad[0] != loPad) || (pad[1] != hiPad),
              BoundsError,
              "Padding value mismatch: pad = " << pad << ", (loPad,hiPad) = ("
              << loPad << "," << hiPad << ")");
          }
        }
      }

      // Extract the _bndryPad data
      if (commIndex == 0                         ) _bndryPad[axis][0] = pad[0];
      if (commIndex == mdComm->getCommDim(axis)-1) _bndryPad[axis][1] = pad[1];

      // Extract the verified _globalRankBounds
      _globalRankBounds[axis][commIndex] = bounds;
    }
  }

  // Check the sanity of _globalRankBounds
  for (int axis = 0; axis < numDims; ++axis)
  {
    for (int commIndex = 1; commIndex < _mdComm->getCommDim(axis); ++commIndex)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        _globalRankBounds[axis][commIndex-1].stop() !=
          _globalRankBounds[axis][commIndex  ].start(),
        MDMapNoncontiguousError,
        "Global rank bounds are not contiguous");
    }
  }

  // Set the global data
  for (int axis = 0; axis < numDims; ++axis)
  {
    int commSize = _mdComm->getCommDim(axis);
    dim_type start =
      _globalRankBounds[axis][0         ].start() - _pad[axis][0];
    dim_type stop  =
      _globalRankBounds[axis][commSize-1].stop()  + _pad[axis][1];
    _globalDims[axis]   = stop - start;
    _globalBounds[axis] = Slice(start, stop);
  }
  _globalStrides = computeStrides< size_type, dim_type >(_globalDims, _layout);

  // Set the global min and max
  for (int axis = 0; axis < numDims; ++axis)
  {
    _globalMin += _globalBounds[axis].start() * _globalStrides[axis];
    _globalMax += _globalBounds[axis].stop()  * _globalStrides[axis];
  }

  // Set the local data
  for (int axis = 0; axis < numDims; ++axis)
  {
    int commIndex = _mdComm->getCommIndex(axis);
    dim_type start =
      _globalRankBounds[axis][commIndex].start() - _pad[axis][0];
    dim_type stop  =
      _globalRankBounds[axis][commIndex].stop()  + _pad[axis][1];
    _localDims[axis]   = stop - start;
    _localBounds[axis] = Slice(stop - start);
  }
  _localStrides = computeStrides< size_type, dim_type >(_localDims, _layout);

  // Compute the local max
  for (int axis = 0; axis < numDims; ++axis)
    _localMax += (_localDims[axis] - 1) * _localStrides[axis];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const MDMap< Node > & source) :
  _mdComm(source._mdComm),
  _globalDims(source._globalDims),
  _globalBounds(source._globalBounds),
  _globalRankBounds(source._globalRankBounds),
  _globalStrides(source._globalStrides),
  _globalMin(source._globalMin),
  _globalMax(source._globalMax),
  _localDims(source._localDims),
  _localBounds(source._localBounds),
  _localStrides(source._localStrides),
  _localMin(source._localMin),
  _localMax(source._localMax),
  _commPadSizes(source._commPadSizes),
  _pad(source._pad),
  _bndryPadSizes(source._bndryPadSizes),
  _bndryPad(source._bndryPad),
  _replicatedBoundary(source._replicatedBoundary),
  _layout(source._layout)
{
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const MDMap< Node > & parent,
      int axis,
      dim_type index) :
  _mdComm(parent._mdComm),
  _globalDims(),
  _globalBounds(),
  _globalRankBounds(),
  _globalStrides(),
  _globalMin(),
  _globalMax(),
  _localDims(),
  _localBounds(),
  _localStrides(),
  _localMin(),
  _localMax(),
  _commPadSizes(),
  _pad(),
  _bndryPadSizes(),
  _bndryPad(),
  _replicatedBoundary(),
  _layout(parent._layout)
{
  int rank = parent.getTeuchosComm()->getRank();
  if (parent.onSubcommunicator())
  {
    int numDims = parent.numDims();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for communicator with " <<
        numDims << " dimensions");

    Slice globalBounds = parent.getGlobalBounds(axis,true);
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((index < globalBounds.start()) || (index >= globalBounds.stop())),
      RangeError,
      "index = " << index  << " is invalid for MDMap axis " <<
      axis << " with bounds " << globalBounds);

    // Determine the axis rank for the processor on which the index
    // lives, and construct the MDComm
    int thisAxisRank = -1;
    for (int axisRank = 0; axisRank < parent.getCommDim(axis); ++axisRank)
      if (index >= parent._globalRankBounds[axis][axisRank].start() &&
          index < parent._globalRankBounds[axis][axisRank].stop())
        thisAxisRank = axisRank;
    TEUCHOS_TEST_FOR_EXCEPTION(
      (thisAxisRank == -1),
      InvalidArgument,
      "error computing axis rank for sub-communicator");
    _mdComm = Teuchos::rcp(new MDComm(*(parent._mdComm), axis, thisAxisRank));
  }

  // There are now two ways for this processor to be off the
  // sub-communicator: (1) it came in that way, or (2) it is not on
  // the new sub-communicator.  Either way, this will be reflected in
  // the new _mdComm, so we check it now.
  if (_mdComm->onSubcommunicator())
  {
    int numDims = parent.numDims();
    if (numDims == 1)
    {
      _globalDims.push_back(1);
      _globalBounds.push_back(ConcreteSlice(index,index+1));
      Teuchos::Array< Slice > bounds(1);
      bounds[0] = ConcreteSlice(index, index+1);
      _globalRankBounds.push_back(bounds);
      _globalStrides.push_back(1);
      _globalMin = index * parent._globalStrides[axis];
      _globalMax = _globalMin;
      _localDims.push_back(1);
      _localBounds.push_back(ConcreteSlice(0,1));
      _localStrides.push_back(1);
      _localMin = parent._localMin +
        (index - parent._globalRankBounds[axis][0].start()) *
        parent._localStrides[axis];
      _localMax = _localMin + 1;
      _commPadSizes.push_back(0);
      _pad.push_back(Teuchos::tuple(0,0));
      _bndryPadSizes.push_back(0);
      _bndryPad.push_back(Teuchos::tuple(0,0));
      _replicatedBoundary.push_back(0);
    }
    else
    {
      _globalMin = parent._globalMin;
      _globalMax = parent._globalMax;
      _localMin  = parent._localMin;
      _localMax  = parent._localMax;
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        if (myAxis != axis)
        {
          _globalDims.push_back(parent._globalDims[myAxis]);
          _globalBounds.push_back(parent._globalBounds[myAxis]);
          _globalRankBounds.push_back(parent._globalRankBounds[myAxis]);
          _globalStrides.push_back(parent._globalStrides[myAxis]);
          _localDims.push_back(parent._localDims[myAxis]);
          _localBounds.push_back(parent._localBounds[myAxis]);
          _localStrides.push_back(parent._localStrides[myAxis]);
          _commPadSizes.push_back(parent._commPadSizes[myAxis]);
          _pad.push_back(parent._pad[myAxis]);
          _bndryPadSizes.push_back(parent._bndryPadSizes[myAxis]);
          _bndryPad.push_back(parent._bndryPad[myAxis]);
          _replicatedBoundary.push_back(parent._replicatedBoundary[myAxis]);
        }
        else
        {
          int axisRank = parent.getCommIndex(axis);
          _globalMin += index * parent._globalStrides[axis];
          _globalMax -= (parent._globalBounds[axis].stop() - index) *
            parent._globalStrides[axis];
          _localMin += (index-parent._globalRankBounds[axis][axisRank].start())
            * parent._localStrides[axis];
          _localMax -= (parent._globalRankBounds[axis][axisRank].stop()-index-1)
            * parent._localStrides[axis];
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const MDMap< Node > & parent,
      int axis,
      const Slice & slice,
      int bndryPad) :
  _mdComm(parent._mdComm),
  _globalDims(parent._globalDims),
  _globalBounds(parent._globalBounds),
  _globalRankBounds(parent._globalRankBounds),
  _globalStrides(parent._globalStrides),
  _globalMin(parent._globalMin),
  _globalMax(parent._globalMax),
  _localDims(parent._localDims),
  _localBounds(parent._localBounds),
  _localStrides(parent._localStrides),
  _localMin(parent._localMin),
  _localMax(parent._localMax),
  _commPadSizes(parent._commPadSizes),
  _pad(parent._pad),
  _bndryPadSizes(parent._bndryPadSizes),
  _bndryPad(parent._bndryPad),
  _replicatedBoundary(parent._replicatedBoundary),
  _layout(parent._layout)
{
  if (parent.onSubcommunicator())
  {
    // Temporarily store the number of dimensions
    int numDims = parent.numDims();

    // Sanity check
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for MDMap with " <<
        numDims << " dimensions");

    // Convert the slice to concrete and check
    Slice bounds =
      slice.bounds(parent.getGlobalBounds(axis,true).stop());
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((bounds.start() < parent.getGlobalBounds(axis).start()) ||
       (bounds.stop() > parent.getGlobalBounds(axis).stop())),
      RangeError,
      "Slice along axis " << axis << " is " << bounds << " but must be within "
      << parent.getGlobalBounds(axis));
    TEUCHOS_TEST_FOR_EXCEPTION(
      (bounds.stop() == bounds.start()),
      RangeError,
      "Slice along axis " << axis << " has length zero");

    // Copy the boundary padding sizes and set initial values for
    // _bndryPad
    _bndryPadSizes[axis] = bndryPad;
    _bndryPad[axis][0]   = bndryPad;
    _bndryPad[axis][1]   = bndryPad;

    // Adjust _globalRankBounds
    for (int axisRank = 0; axisRank < parent.getCommDim(axis);
         ++axisRank)
    {
      dim_type start = _globalRankBounds[axis][axisRank].start();
      dim_type stop  = _globalRankBounds[axis][axisRank].stop();
      if (start < bounds.start()) start = bounds.start();
      if (stop  > bounds.stop() ) stop  = bounds.stop();
      _globalRankBounds[axis][axisRank] = ConcreteSlice(start, stop);
    }

    // Alter _bndryPad if necessary
    dim_type start = bounds.start() - _bndryPadSizes[axis];
    if (start < 0)
    {
      _bndryPad[axis][0] = bounds.start();
      start = 0;
    }
    dim_type stop = bounds.stop() + _bndryPadSizes[axis];
    if (stop > parent.getGlobalBounds(axis,true).stop())
    {
      _bndryPad[axis][1] = parent.getGlobalBounds(axis,true).stop() -
        bounds.stop();
      stop = parent.getGlobalBounds(axis,true).stop();
    }

    // Compute _globalBounds, _globalDims, _globalMax, and _globalMin
    _globalBounds[axis] = ConcreteSlice(start,stop);
    _globalDims[axis]   = stop - start;
    _globalMin         += start * _globalStrides[axis];
    _globalMax         -= (parent.getGlobalDim(axis,true) - stop) *
                           _globalStrides[axis];

    // Alter _replicatedBoundary
    if ((parent.getGlobalBounds(axis,true).start() < _globalBounds[axis].start())
        || (parent.getGlobalBounds(axis,true).stop() > _globalBounds[axis].stop()))
      _replicatedBoundary[axis] = 0;

    // Build the slice for the MDComm sub-communicator constructor
    int pStart = -1;
    int pStop  = -1;
    for (int axisRank = 0; axisRank < parent.getCommDim(axis);
         ++axisRank)
    {
      if ((_globalRankBounds[axis][axisRank].start() - _bndryPad[axis][0]
           <= _globalBounds[axis].start()) &&
          (_globalBounds[axis].start() <
           _globalRankBounds[axis][axisRank].stop() + _bndryPad[axis][1]))
        if (pStart == -1) pStart = axisRank;
      if ((_globalRankBounds[axis][axisRank].start() - _bndryPad[axis][0]
           < _globalBounds[axis].stop()) &&
          (_globalBounds[axis].stop() <=
           _globalRankBounds[axis][axisRank].stop() + _bndryPad[axis][1]))
        pStop = axisRank+1;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      (pStart == -1 || pStop == -1),
      InvalidArgument,
      "error computing axis rank slice");
    Slice axisRankSlice = ConcreteSlice(pStart,pStop);

    // Construct the MDComm sub-communicator
    _mdComm = Teuchos::rcp(new MDComm(*(parent._mdComm), axis, axisRankSlice));

    // We now have a sub-communicator, and should only construct this
    // MDMap if this processor is on it.  If this processor is off the
    // communicator, then we clear many of the data members.
    if (_mdComm->onSubcommunicator())
    {
      // Fix _pad, if needed
      int parentAxisRank = parent.getCommIndex(axis);
      int myAxisRank     = _mdComm->getCommIndex(axis);
      if (myAxisRank == 0)
        _pad[axis][0] = _bndryPad[axis][0];
      if (myAxisRank == _mdComm->getCommDim(axis)-1)
        _pad[axis][1] = _bndryPad[axis][1];

      // Compute the local start and stop indexes.  Note that
      // _globalRankBounds has an axis-rank dimension, and that it
      // still uses the parent's commDim, not the new ones.  We will
      // fix this later.
      dim_type start = (_globalRankBounds[axis][parentAxisRank].start() -
                        _pad[axis][0]) -
                       (parent._globalRankBounds[axis][parentAxisRank].start() -
                        parent._pad[axis][0]);
      dim_type stop  = (_globalRankBounds[axis][parentAxisRank].stop() +
                        _pad[axis][1]) -
                       (parent._globalRankBounds[axis][parentAxisRank].start() -
                        parent._pad[axis][0]);

      // Compute the local bounds, dims, min, and max
      _localBounds[axis] = ConcreteSlice(stop - start);
      _localDims[axis]   = stop - start;
      _localMin         += start * _localStrides[axis];
      _localMax         -= (parent.getLocalBounds(axis,true).stop() - stop) *
                            _localStrides[axis];

      // The new sub-communicator may have fewer processors than the
      // parent communicator, so we need to fix _globalRankBounds
      Teuchos::Array< Slice > newRankBounds;
      for (int axisRank = 0; axisRank < parent.getCommDim(axis);
           ++axisRank)
        if ((axisRank >= axisRankSlice.start()) &&
            (axisRank <  axisRankSlice.stop() )   )
          newRankBounds.push_back(_globalRankBounds[axis][axisRank]);
      _globalRankBounds[axis] = newRankBounds;
    }
    else
    {
      _localDims.clear();
      _localMin = 0;
      _localMax = 0;
      _localBounds.clear();
      _commPadSizes.clear();
      _pad.clear();
      _bndryPadSizes.clear();
      _bndryPad.clear();
      _localStrides.clear();
    }
  }
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::
MDMap(const MDMap< Node > & parent,
      const Teuchos::ArrayView< Slice > & slices,
      const Teuchos::ArrayView< int > & bndryPad)
{
  // Temporarily store the number of dimensions
  int numDims = parent.numDims();

  // Sanity check on dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    (slices.size() != numDims),
    InvalidArgument,
    "number of slices = " << slices.size() << " != parent MDMap number of "
    "dimensions = " << numDims);

  // Apply the single-Slice constructor to each axis in succession
  MDMap< Node > tempMDMap1(parent);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int bndryPadding = (axis < bndryPad.size()) ? bndryPad[axis] : 0;
    MDMap< Node > tempMDMap2(tempMDMap1,
                             axis,
                             slices[axis],
                             bndryPadding);
    tempMDMap1 = tempMDMap2;
  }
  *this = tempMDMap1;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node >::~MDMap()
{
}

////////////////////////////////////////////////////////////////////////

template< class Node >
MDMap< Node > &
MDMap< Node >::operator=(const MDMap< Node > & source)
{
  _mdComm           = source._mdComm;
  _globalDims       = source._globalDims;
  _globalBounds     = source._globalBounds;
  _globalRankBounds = source._globalRankBounds;
  _globalStrides    = source._globalStrides;
  _globalMin        = source._globalMin;
  _globalMax        = source._globalMax;
  _localDims        = source._localDims;
  _localBounds      = source._localBounds;
  _localStrides     = source._localStrides;
  _localMin         = source._localMin;
  _localMax         = source._localMax;
  _commPadSizes     = source._commPadSizes;
  _pad              = source._pad;
  _bndryPadSizes    = source._bndryPadSizes;
  _bndryPad         = source._bndryPad;
  _layout           = source._layout;
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::RCP< const MDComm >
MDMap< Node >::getMDComm() const
{
  return _mdComm;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::onSubcommunicator() const
{
  return _mdComm->onSubcommunicator();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::RCP< const Teuchos::Comm< int > >
MDMap< Node >::getTeuchosComm() const
{
  return _mdComm->getTeuchosComm();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::numDims() const
{
  return _mdComm->numDims();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< int >
MDMap< Node >::getCommDims() const
{
  return _mdComm->getCommDims();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getCommDim(int axis) const
{
  return _mdComm->getCommDim(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::isPeriodic(int axis) const
{
  return _mdComm->isPeriodic(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getCommIndex(int axis) const
{
  return _mdComm->getCommIndex(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getLowerNeighbor(int axis) const
{
  return _mdComm->getLowerNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getUpperNeighbor(int axis) const
{
  return _mdComm->getUpperNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< dim_type >
MDMap< Node >::
getGlobalDims() const
{
  return _globalDims;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
dim_type
MDMap< Node >::
getGlobalDim(int axis,
             bool withBndryPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (withBndryPad)
    return _globalDims[axis];
  else
    return _globalDims[axis] - _bndryPad[axis][0] - _bndryPad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Slice
MDMap< Node >::
getGlobalBounds(int axis,
                bool withBndryPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (withBndryPad)
    return _globalBounds[axis];
  else
  {
    dim_type start = _globalBounds[axis].start() + _bndryPad[axis][0];
    dim_type stop  = _globalBounds[axis].stop()  - _bndryPad[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class Node >
dim_type
MDMap< Node >::
getLocalDim(int axis,
            bool withPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (withPad)
    return _localDims[axis];
  else
    return _localDims[axis] - _pad[axis][0] - _pad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< dim_type >
MDMap< Node >::
getLocalDims() const
{
  return _localDims;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Slice
MDMap< Node >::
getGlobalRankBounds(int axis,
                    bool withBndryPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  int axisRank = getCommIndex(axis);
  if (withBndryPad)
  {
    dim_type start = _globalRankBounds[axis][axisRank].start();
    dim_type stop  = _globalRankBounds[axis][axisRank].stop();
    if (getCommIndex(axis) == 0)
      start -= _bndryPad[axis][0];
    if (getCommIndex(axis) == getCommDim(axis)-1)
      stop += _bndryPad[axis][1];
    return ConcreteSlice(start,stop);
  }
  else
    return _globalRankBounds[axis][axisRank];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::ArrayView< const Slice >
MDMap< Node >::
getLocalBounds() const
{
  return _localBounds();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Slice
MDMap< Node >::
getLocalBounds(int axis,
               bool withPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (withPad)
    return _localBounds[axis];
  else
  {
    dim_type start = _localBounds[axis].start() + _pad[axis][0];
    dim_type stop  = _localBounds[axis].stop()  - _pad[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Slice
MDMap< Node >::
getLocalInteriorBounds(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  dim_type start = _localBounds[axis].start() + _pad[axis][0];
  dim_type stop  = _localBounds[axis].stop()  - _pad[axis][1];
  if (_mdComm->getLowerNeighbor(axis) == -1) ++start;
  if (_mdComm->getUpperNeighbor(axis) == -1) --stop;
  return ConcreteSlice(start, stop);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::hasPadding() const
{
  bool result = false;
  for (int axis = 0; axis < numDims(); ++axis)
    if (_pad[axis][0] + _pad[axis][1]) result = true;
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getLowerPadSize(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _pad[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getUpperPadSize(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _pad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getCommPadSize(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _commPadSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getLowerBndryPad(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _bndryPad[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getUpperBndryPad(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _bndryPad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< int >
MDMap< Node >::getBndryPadSizes() const
{
  return _bndryPadSizes;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
int
MDMap< Node >::getBndryPadSize(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _bndryPadSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::
isPad(const Teuchos::ArrayView< dim_type > & index) const
{
  bool result = false;
  for (int axis = 0; axis < numDims(); ++axis)
  {
    if (index[axis] < getLowerPadSize(axis))
      result = true;
    if (index[axis] >= getLocalDim(axis,true) - getUpperPadSize(axis))
      result = true;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::
isCommPad(const Teuchos::ArrayView< dim_type > & index) const
{
  bool result = false;
  for (int axis = 0; axis < numDims(); ++axis)
  {
    // Check the ranks of the lower and upper neighbor processors.  If
    // either of these values is non-negative, then we are on a
    // processor that contains communication padding
    if (getLowerNeighbor(axis) >= 0)
    {
      if (index[axis] < getLowerPadSize(axis))
        result = true;
    }
    if (getUpperNeighbor(axis) >= 0)
    {
      if (index[axis] >= getLocalDim(axis,true) - getUpperPadSize(axis))
        result = true;
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::
isBndryPad(const Teuchos::ArrayView< dim_type > & index) const
{
  bool result = false;
  for (int axis = 0; axis < numDims(); ++axis)
  {
    // Check the ranks of the lower and upper neighbor processors.  If
    // either of these values is -1, then we are on a processor that
    // contains a boundary
    if (getLowerNeighbor(axis) == -1)
    {
      if (index[axis] < getLowerPadSize(axis))
        result = true;
    }
    if (getUpperNeighbor(axis) == -1)
    {
      if (index[axis] >= getLocalDim(axis,true) - getUpperPadSize(axis))
        result = true;
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::isReplicatedBoundary(int axis) const
{
  return _mdComm->isPeriodic(axis) && bool(_replicatedBoundary[axis]);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Layout
MDMap< Node >::getLayout() const
{
  return _layout;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::RCP< Node >
MDMap< Node >::getNode() const
{
  return _node;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::ArrayView< Teuchos::RCP< const MDMap< Node > > >
MDMap< Node >::getAxisMaps() const
{
  if (_axisMaps.size() == 0) _axisMaps.resize(numDims());
  for (int axis = 0; axis < numDims(); ++axis)
    if (_axisMaps[axis].is_null()) getAxisMap(axis);
  return _axisMaps();
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::RCP< const MDMap< Node > >
MDMap< Node >::getAxisMap(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (_axisMaps.size() == 0) _axisMaps.resize(numDims());
  if (_axisMaps[axis].is_null())
  {
    Teuchos::RCP< const MDComm > axisComm = _mdComm->getAxisComm(axis);
    Domi::dim_type axisDim = _globalDims[axis] - 2*_bndryPadSizes[axis];
    _axisMaps[axis] =
      Teuchos::rcp(new MDMap(axisComm,
                             Teuchos::tuple(axisDim),
                             Teuchos::tuple(_commPadSizes[axis]),
                             Teuchos::tuple(_bndryPadSizes[axis]),
                             Teuchos::tuple(_replicatedBoundary[axis]),
                             _layout));
  }
  return _axisMaps[axis];
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::RCP< const MDMap< Node > >
MDMap< Node >::getAugmentedMDMap(const dim_type leadingDim,
                                 const dim_type trailingDim) const
{
  // Construct the new MDMap
  MDMap< Node > * newMdMap = new MDMap< Node >(*this);

  // Compute the number of new dimensions
  int newNumDims = 0;
  if (leadingDim  > 0) ++newNumDims;
  if (trailingDim > 0) ++newNumDims;

  // Trivial result
  if (newNumDims == 0) return Teuchos::rcp(newMdMap);

  // Compute the new MDComm
  int oldNumDims = numDims();
  Teuchos::Array< int > newCommDims(oldNumDims);
  Teuchos::Array< int > newPeriodic(oldNumDims);
  Teuchos::Array< int > newReplicatedBndry(oldNumDims);
  
  for (int axis = 0; axis < oldNumDims; ++axis)
  {
    newCommDims[axis]        = getCommDim(axis);
    newPeriodic[axis]        = int(isPeriodic(axis));
    newReplicatedBndry[axis] = int(isReplicatedBoundary(axis));
  }
  if (leadingDim > 0)
  {
    newCommDims.insert(newCommDims.begin(),1);
    newPeriodic.insert(newPeriodic.begin(),0);
    newReplicatedBndry.insert(newReplicatedBndry.begin(),0);
  }
  if (trailingDim > 0)
  {
    newCommDims.push_back(1);
    newPeriodic.push_back(0);
    newReplicatedBndry.push_back(0);
  }
  newMdMap->_mdComm = Teuchos::rcp(new MDComm(getTeuchosComm(),
                                              newCommDims,
                                              newPeriodic));
  newMdMap->_replicatedBoundary = newReplicatedBndry;

  // Adjust new MDMap arrays for a new leading dimension
  Slice slice = Slice(leadingDim);
  padding_type pad(Teuchos::tuple(0,0));
  if (leadingDim > 0)
  {
    newMdMap->_globalDims.insert(newMdMap->_globalDims.begin(), leadingDim);
    newMdMap->_globalBounds.insert(newMdMap->_globalBounds.begin(), slice);
    newMdMap->_globalRankBounds.insert(newMdMap->_globalRankBounds.begin(),
                                      Teuchos::Array< Slice >(1,slice));
    newMdMap->_localDims.insert(newMdMap->_localDims.begin(), leadingDim);
    newMdMap->_localBounds.insert(newMdMap->_localBounds.begin(), slice);
    newMdMap->_commPadSizes.insert(newMdMap->_commPadSizes.begin(),0);
    newMdMap->_pad.insert(newMdMap->_pad.begin(), pad);
    newMdMap->_bndryPadSizes.insert(newMdMap->_bndryPadSizes.begin(),0);
    newMdMap->_bndryPad.insert(newMdMap->_bndryPad.begin(), pad);
  }

  // Adjust new MDMap arrays for a new trailing dimension
  slice = Slice(trailingDim);
  if (trailingDim > 0)
  {
    newMdMap->_globalDims.push_back(trailingDim);
    newMdMap->_globalBounds.push_back(slice);
    newMdMap->_globalRankBounds.push_back(Teuchos::Array< Slice >(1,slice));
    newMdMap->_localDims.push_back(trailingDim);
    newMdMap->_localBounds.push_back(slice);
    newMdMap->_commPadSizes.push_back(0);
    newMdMap->_pad.push_back(pad);
    newMdMap->_bndryPadSizes.push_back(0);
    newMdMap->_bndryPad.push_back(pad);
  }

  // Compute the new stride related data
  newMdMap->_globalStrides =
    computeStrides< size_type, dim_type >(newMdMap->_globalDims,
                                          newMdMap->_layout);
  newMdMap->_localStrides =
    computeStrides< size_type, dim_type >(newMdMap->_localDims,
                                          newMdMap->_layout);
  newMdMap->_globalMin = 0;
  newMdMap->_globalMax = 0;
  newMdMap->_localMin  = 0;
  newMdMap->_localMax  = 0;
  for (int axis = 0; axis < oldNumDims + newNumDims; ++axis)
  {
    newMdMap->_globalMin += newMdMap->_globalBounds[axis].start() *
                            newMdMap->_globalStrides[axis];
    newMdMap->_globalMax += newMdMap->_globalBounds[axis].stop() *
                            newMdMap->_globalStrides[axis];
    newMdMap->_localMin  += newMdMap->_localBounds[axis].start() *
                            newMdMap->_localStrides[axis];
    newMdMap->_localMax  += newMdMap->_localBounds[axis].stop() *
                            newMdMap->_localStrides[axis];
  }

  // Return the result
  return Teuchos::rcp(newMdMap);
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA

template< class Node >
Teuchos::RCP< const Epetra_Map >
MDMap< Node >::getEpetraMap(bool withCommPad) const
{
  if (withCommPad)
  {
    if (_epetraMap.is_null())
    {
      // Check if the maximum global ID is larger than what an int can
      // hold (because Epetra uses int ordinals)
      TEUCHOS_TEST_FOR_EXCEPTION(
        computeSize(_globalDims) - 1 > std::numeric_limits< int >::max(),
        MapOrdinalError,
        "The maximum global ID of this MDMap is too large for an Epetra_Map");

      // Allocate the myElements MDArray and the index array
      int num_dims = numDims();
      Teuchos::Array<dim_type> localDims(num_dims);
      for (int axis = 0; axis < num_dims; ++axis)
        localDims[axis] = _localDims[axis];
      MDArray<int> myElements(localDims);
      Teuchos::Array<int> index(num_dims);

      // Iterate over the local MDArray and assign global IDs
      for (MDArray<int>::iterator it = myElements.begin();
           it != myElements.end(); ++it)
      {
        int globalID = 0;
        for (int axis = 0; axis < num_dims; ++axis)
        {
          int axisRank = getCommIndex(axis);
          int start    = _globalRankBounds[axis][axisRank].start() -
                         _pad[axis][0];
          globalID += (start + it.index(axis)) * _globalStrides[axis];
        }
        *it = globalID;
      }

      // Construct the Epetra_Map
      Teuchos::RCP< const Epetra_Comm > epetraComm = _mdComm->getEpetraComm();
      _epetraMap = Teuchos::rcp(new Epetra_Map(-1,
                                               myElements.size(),
                                               myElements.getRawPtr(),
                                               0,
                                               *epetraComm));
    }
    return _epetraMap;
  }
  else
  {
    if (_epetraOwnMap.is_null())
    {
      // Check if the maximum global ID is larger than what an int can
      // hold (because Epetra uses int ordinals)
      if (computeSize(_globalDims) - 1 > std::numeric_limits< int >::max())
        throw MapOrdinalError("The maximum global ID of this MDMap is too "
                              "large for an Epetra_Map");

      // Allocate the myElements MDArray and the index array
      int num_dims = numDims();
      Teuchos::Array<int> index(num_dims);
      Teuchos::Array<dim_type> myDims(num_dims);
      for (int axis = 0; axis < num_dims; ++axis)
      {
        myDims[axis] = _localDims[axis] - _pad[axis][0] - _pad[axis][1];
        int axisRank = getCommIndex(axis);
        if (axisRank == 0)
          myDims[axis] += _bndryPad[axis][0];
        if (axisRank == getCommDim(axis)-1)
          myDims[axis] += _bndryPad[axis][1];
      }
      MDArray<int> myElements(myDims());

      // Iterate over the local MDArray and assign global IDs
      for (MDArray<int>::iterator it = myElements.begin();
           it != myElements.end(); ++it)
      {
        int globalID = 0;
          for (int axis = 0; axis < num_dims; ++axis)
          {
            int axisRank = getCommIndex(axis);
            int start    = _globalRankBounds[axis][axisRank].start();
            if (axisRank == 0)
              start -= _bndryPad[axis][0];
            if (axisRank == getCommDim(axis)-1)
              start += _bndryPad[axis][1];
            globalID += (start + it.index(axis)) * _globalStrides[axis];
          }
      }

      // Construct the Epetra_Map
      Teuchos::RCP< const Epetra_Comm > epetraComm = _mdComm->getEpetraComm();
      _epetraOwnMap = Teuchos::rcp(new Epetra_Map(-1,
                                                  myElements.size(),
                                                  myElements.getRawPtr(),
                                                  0,
                                                  *epetraComm));
    }
    return _epetraOwnMap;
  }
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA

template< class Node >
Teuchos::RCP< const Epetra_Map >
MDMap< Node >::
getEpetraAxisMap(int axis,
                 bool withCommPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if ((withCommPad     && (_epetraAxisMaps.size()    == 0)) ||
      (not withCommPad && (_epetraAxisOwnMaps.size() == 0)))
  {
    int num_dims = numDims();
    Teuchos::RCP< const Epetra_Comm > epetraComm = _mdComm->getEpetraComm();
    for (int axis=0; axis < num_dims; ++axis)
    {
      Teuchos::Array<int> elements(getLocalDim(axis, withCommPad));
      int start = getGlobalRankBounds(axis,true).start();
      if (withCommPad && (getCommIndex(axis) != 0)) start -= _pad[axis][0];
      for (int i = 0; i < elements.size(); ++i)
        elements[i] = i + start;
      if (withCommPad)
      {
        _epetraAxisMaps.push_back(
          Teuchos::rcp(new Epetra_Map(-1,
                                      elements.size(),
                                      elements.getRawPtr(),
                                      0,
                                      *epetraComm)));
      }
      else
      {
        _epetraAxisOwnMaps.push_back(
          Teuchos::rcp(new Epetra_Map(-1,
                                      elements.size(),
                                      elements.getRawPtr(),
                                      0,
                                      *epetraComm)));
      }
    }
  }

  if (withCommPad)
    return _epetraAxisMaps[axis];
  else
    return _epetraAxisOwnMaps[axis];
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class Node >
template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node2 >
Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node2 > >
MDMap< Node >::getTpetraMap(bool withCommPad) const
{
  if (withCommPad)
  {
    // Allocate the elementsMDArray and the index array
    int num_dims = numDims();
    Teuchos::Array<dim_type> localDims(num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
      localDims[axis] = _localDims[axis];
    MDArray< GlobalOrdinal > elementMDArray(localDims);
    Teuchos::Array< LocalOrdinal > index(num_dims);

    // Iterate over the local MDArray and assign global IDs
    for (typename MDArray< GlobalOrdinal >::iterator it = elementMDArray.begin();
         it != elementMDArray.end(); ++it)
    {
      GlobalOrdinal globalID = 0;
      for (int axis = 0; axis < num_dims; ++axis)
      {
        int axisRank        = getCommIndex(axis);
        GlobalOrdinal start = _globalRankBounds[axis][axisRank].start() -
                              _pad[axis][0];
        globalID += (start + it.index(axis)) * _globalStrides[axis];
      }
      *it = globalID;
    }

    // Return the Tpetra::Map
    const Teuchos::Array< GlobalOrdinal > & myElements =
      elementMDArray.array();
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      _mdComm->getTeuchosComm();
    return
      Teuchos::rcp(new Tpetra::Map< LocalOrdinal,
                                    GlobalOrdinal,
                                    Node2 >(Teuchos::OrdinalTraits<
                                              Tpetra::global_size_t>::invalid(),
                                            myElements(),
                                            0,
                                            teuchosComm));
  }
  else
  {
    // Allocate the elementMDArray MDArray and the index array
    int num_dims = numDims();
    Teuchos::Array< LocalOrdinal > index(num_dims);
    Teuchos::Array< dim_type >     myDims(num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      myDims[axis] =
        _localDims[axis] - _pad[axis][0] - _pad[axis][1];
      int axisRank = getCommIndex(axis);
      if (axisRank == 0)
        myDims[axis] += _bndryPad[axis][0];
      if (axisRank == getCommDim(axis)-1)
        myDims[axis] += _bndryPad[axis][1];
    }
    MDArray< GlobalOrdinal > elementMDArray(myDims());

    // Iterate over the local MDArray and assign global IDs
    for (typename MDArray< GlobalOrdinal >::iterator it = elementMDArray.begin();
         it != elementMDArray.end(); ++it)
    {
      GlobalOrdinal globalID = 0;
      for (int axis = 0; axis < num_dims; ++axis)
      {
        int axisRank        = getCommIndex(axis);
        GlobalOrdinal start = _globalRankBounds[axis][axisRank].start();
        if (axisRank == 0)
          start -= _bndryPad[axis][0];
        if (axisRank == getCommDim(axis)-1)
          start += _bndryPad[axis][1];
        globalID += (start + it.index(axis)) * _globalStrides[axis];
      }
    }

    // Return the Tpetra::Map
    const Teuchos::Array< GlobalOrdinal> & myElements =
      elementMDArray.array();
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      _mdComm->getTeuchosComm();
    return
      Teuchos::rcp(new Tpetra::Map< LocalOrdinal,
                                    GlobalOrdinal,
                                    Node >(Teuchos::OrdinalTraits<
                                             Tpetra::global_size_t>::invalid(),
                                           myElements(),
                                           0,
                                           teuchosComm));
  }
}
#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class Node >
template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node2 >
Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node2 > >
MDMap< Node >::
getTpetraAxisMap(int axis,
                 bool withCommPad) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  int num_dims = numDims();
  Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
    _mdComm->getTeuchosComm();
  Teuchos::Array< GlobalOrdinal > elements(getLocalDim(axis,withCommPad));
  GlobalOrdinal start = getGlobalRankBounds(axis,true).start();
  if (withCommPad && (getCommIndex(axis) != 0)) start -= _pad[axis][0];
  for (LocalOrdinal i = 0; i < elements.size(); ++i)
    elements[i] = i + start;
  return Teuchos::rcp(new Tpetra::Map< LocalOrdinal,
                                       GlobalOrdinal,
                                       Node >(Teuchos::OrdinalTraits<
                                                   Tpetra::global_size_t>::invalid(),
                                                   elements,
                                                   0,
                                                   teuchosComm));
}
#endif

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< dim_type >
MDMap< Node >::
getGlobalIndex(size_type globalID) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((globalID < _globalMin) || (globalID >= _globalMax)),
    RangeError,
    "invalid global index = " << globalID << " (should be between " <<
    _globalMin << " and " << _globalMax << ")");
#endif
  int num_dims = numDims();
  Teuchos::Array< dim_type > result(num_dims);
  size_type index = globalID;
  if (_layout == LAST_INDEX_FASTEST)
  {
    for (int axis = 0; axis < num_dims-1; ++axis)
    {
      result[axis] = index / _globalStrides[axis];
      index        = index % _globalStrides[axis];
    }
    result[num_dims-1] = index;
  }
  else
  {
    for (int axis = num_dims-1; axis > 0; --axis)
    {
      result[axis] = index / _globalStrides[axis];
      index        = index % _globalStrides[axis];
    }
    result[0] = index;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
Teuchos::Array< dim_type >
MDMap< Node >::
getLocalIndex(size_type localID) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((localID < _localMin) || (localID >= _localMax)),
    RangeError,
    "invalid local index = " << localID << " (should be between " <<
    _localMin << " and " << _localMax << ")");
#endif
  int num_dims = numDims();
  Teuchos::Array< dim_type > result(num_dims);
  size_type index = localID;
  if (_layout == LAST_INDEX_FASTEST)
  {
    for (int axis = 0; axis < num_dims-1; ++axis)
    {
      result[axis] = index / _localStrides[axis];
      index        = index % _localStrides[axis];
    }
    result[num_dims-1] = index;
  }
  else
  {
    for (int axis = num_dims-1; axis > 0; --axis)
    {
      result[axis] = index / _localStrides[axis];
      index        = index % _localStrides[axis];
    }
    result[0] = index;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
size_type
MDMap< Node >::
getGlobalID(size_type localID) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((localID < 0) || (localID >= _localMax)),
    RangeError,
    "invalid local index = " << localID << " (local size = " <<
    _localMax << ")");
#endif
  Teuchos::Array< dim_type > localIndex = getLocalIndex(localID);
  size_type result = 0;
  for (int axis = 0; axis < numDims(); ++axis)
  {
    dim_type globalIndex = localIndex[axis] +
      _globalRankBounds[axis][getCommIndex(axis)].start() - _pad[axis][0];
    result += globalIndex * _globalStrides[axis];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
size_type
MDMap< Node >::
getGlobalID(const Teuchos::ArrayView< dim_type > & globalIndex) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (globalIndex.size() != numDims()),
    InvalidArgument,
    "globalIndex has " << globalIndex.size() << " entries; expecting "
    << numDims());
  for (int axis = 0; axis < numDims(); ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((globalIndex[axis] < 0) ||
       (globalIndex[axis] >= _globalDims[axis])),
      RangeError,
      "invalid globalIndex[" << axis << "] = " << globalIndex[axis] <<
      " (global dimension = " << _globalDims[axis] << ")");
  }
#endif
  size_type result = 0;
  for (int axis = 0; axis < numDims(); ++axis)
    result += globalIndex[axis] * _globalStrides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
size_type
MDMap< Node >::
getLocalID(size_type globalID) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((globalID < _globalMin) || (globalID >= _globalMax)),
    RangeError,
    "invalid global index = " << globalID << " (should be between " <<
    _globalMin << " and " << _globalMax << ")");
#endif
  Teuchos::Array< dim_type > globalIndex =
    getGlobalIndex(globalID);
  size_type result = 0;
  for (int axis = 0; axis < numDims(); ++axis)
  {
    dim_type localIndex = globalIndex[axis] -
      _globalRankBounds[axis][getCommIndex(axis)].start() + _pad[axis][0];
    TEUCHOS_TEST_FOR_EXCEPTION(
      (localIndex < 0 || localIndex >= _localDims[axis]),
      RangeError,
      "global index not on local processor")
    result += localIndex * _localStrides[axis];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
size_type
MDMap< Node >::
getLocalID(const Teuchos::ArrayView< dim_type > & localIndex) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (localIndex.size() != numDims()),
    InvalidArgument,
    "localIndex has " << localIndex.size() << " entries; expecting "
    << numDims());
  for (int axis = 0; axis < numDims(); ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((localIndex[axis] < 0) ||
       (localIndex[axis] >= _localDims[axis])),
      RangeError,
      "invalid localIndex[" << axis << "] = " << localIndex[axis] <<
      " (local dimension = " << _localDims[axis] << ")");
  }
#endif
  size_type result = 0;
  for (int axis = 0; axis < numDims(); ++axis)
    result += localIndex[axis] * _localStrides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::isCompatible(const MDMap< Node > & mdMap) const
{
  // Trivial comparison.  We assume that if the object pointers match
  // on this processor, then they match on all processors
  if (this == &mdMap) return true;

  // Check the number of dimensions.  Assume this check produces the
  // same result on all processors
  int num_dims = numDims();
  if (num_dims != mdMap.numDims()) return false;

  // Check the commDims.  Assume this check produces the same result
  // on all processes
  for (int axis = 0; axis < num_dims; ++axis)
    if (getCommDim(axis) != mdMap.getCommDim(axis)) return false;

  // Check the global dimensions.  Assume this check produces the same
  // result on all processes
  if (_globalDims != mdMap._globalDims) return false;

  // Check the local dimensions.  This needs to be checked locally on
  // each processor and then the results communicated to obtain global
  // result
  int localResult  = 1;
  int globalResult = 1;
  for (int axis = 0; axis < num_dims; ++axis)
    if (getLocalDim(axis,false) != mdMap.getLocalDim(axis,false))
      localResult = 0;
  Teuchos::reduceAll(*(getTeuchosComm()),
                     Teuchos::REDUCE_MIN,
                     1,
                     &localResult,
                     &globalResult);

  // Return the result
  return bool(globalResult);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::isSameAs(const MDMap< Node > & mdMap,
                        const int verbose) const
{
  // Trivial comparison.  We assume that if the object pointers match
  // on this processor, then they match on all processors
  if (this == &mdMap) return true;

  // Start by setting a local result to true.  We will perform a
  // number of tests, and if they fail, the local result will be set
  // to false.  At the end, we will perform a global reduction to
  // obtain the global result.
  int localResult = 1;
  Teuchos::RCP< const Teuchos::Comm< int > > comm = getTeuchosComm();
  int rank = comm->getRank();

  // Check if MDMaps are compatible.
  if (! isCompatible(mdMap))
  {
    localResult = 0;
    if (verbose)
      std::cout << rank << ": MDMaps are incompatible" << std::endl;
  }

  // Check that underlying communicators are the same size
  if (comm->getSize() != mdMap.getTeuchosComm()->getSize())
  {
    localResult = 0;
    if (verbose)
      std::cout << rank << ": this Comm size = " << comm->getSize() << " != "
                << mdMap.getTeuchosComm()->getSize() << std::endl;
  }

  // Check that underlying communicators have the same rank
  if (rank != mdMap.getTeuchosComm()->getRank())
  {
    localResult = 0;
    if (verbose)
      std::cout << rank << ": this Comm rank = " << rank << " != "
                << mdMap.getTeuchosComm()->getRank() << std::endl;
  }

  // Check the global bounds.
  if (_globalBounds != mdMap._globalBounds)
  {
    localResult = 0;
    if (verbose)
      std::cout << rank << ": global bounds " << _globalBounds << " != "
                << mdMap._globalBounds << std::endl;
  }

  // Check the local dimensions.
  if (_localDims != mdMap._localDims)
  {
    localResult = 0;
    if (verbose)
      std::cout << rank << ": local dimensions " << _localDims << " != "
                << mdMap._localDims << std::endl;
  }

  // Obtain the global result
  int globalResult = 1;
  Teuchos::reduceAll(*(getTeuchosComm()),
                     Teuchos::REDUCE_MIN,
                     1,
                     &localResult,
                     &globalResult);

  // Return the result
  return bool(globalResult);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
bool
MDMap< Node >::isContiguous() const
{
  // Compute the local strides if they were contiguous
  Teuchos::Array< size_type > contiguousStrides =
    computeStrides< size_type, dim_type >(_localDims, _layout);

  // Compute the local result: 0 = contiguous, 1 = non-contiguous
  int localResult = int(_localStrides != contiguousStrides);

  // Compute the global result
  int globalResult = 0;
  Teuchos::reduceAll(*(_mdComm->getTeuchosComm()),
                     Teuchos::REDUCE_SUM,
                     1,
                     &localResult,
                     &globalResult);
  return (globalResult == 0);
}

////////////////////////////////////////////////////////////////////////

template< class Node >
void
MDMap< Node >::computeBounds()
{
  // Initialization
  int num_dims = numDims();

  // Decompose the multi-dimensional domain
  for (int axis = 0; axis < num_dims; ++axis)
  {
    // Get the communicator info for this axis
    int commDim = getCommDim(axis);
    for (int axisRank = 0; axisRank < commDim; ++axisRank)
    {
      // First estimates assuming even division of global dimensions
      // by the number of processors along this axis, and ignoring
      // communication and boundary padding.
      dim_type  localDim  = (_globalDims[axis] - _bndryPad[axis][0] -
                             _bndryPad[axis][1]) / commDim;
      dim_type axisStart = axisRank * localDim;

      // Adjustments for non-zero remainder.  Compute the remainder
      // using the mod operator.  If the remainder is > 0, then add an
      // element to the appropriate number of processors with the
      // highest axis ranks.  Note that this is the opposite of the
      // standard Tpetra::Map constructor (which adds an elements to
      // the lowest processor ranks), and provides better balance for
      // finite differencing systems with staggered data location.
      dim_type remainder = (_globalDims[axis] - _bndryPad[axis][0] -
                            _bndryPad[axis][1]) % commDim;
      if (commDim - axisRank - 1 < remainder)
      {
        ++localDim;
        axisStart += (remainder - commDim + axisRank);
      }

      // Global adjustment for boundary padding
      axisStart += _bndryPad[axis][0];

      // Compute and store the global axis bounds
      _globalRankBounds[axis].push_back(
        ConcreteSlice(axisStart, axisStart + localDim));

      // Set _localDims[axis] and _localBounds[axis] only if
      // axisRank equals the axis rank of this processor
      if (axisRank == getCommIndex(axis))
      {
        // Local adjustment for padding.  Note that _pad should
        // already be corrected to be either the communication padding
        // or boundary padding as appropriate
        _localDims[axis] = localDim + _pad[axis][0] + _pad[axis][1];

        // Compute and store the axis bounds
        _localBounds.push_back(ConcreteSlice(_localDims[axis]));
      }
    }
  }
}

}

#endif
