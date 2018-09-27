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
typedef Tpetra::Details::DefaultTypes::local_ordinal_type  TpetraLOType;
typedef Tpetra::Details::DefaultTypes::global_ordinal_type TpetraGOType;
typedef Tpetra::Details::DefaultTypes::node_type           TpetraNodeType;
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
   * \param padding [in] an array of padding_type (a 2-tuple of
   *        integers) specifying the local padding along each
   *        axis. Since this is local, the MDMap constructor can
   *        determine from context whether each pad refers to
   *        communication padding or boundary padding.
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
  MDMap(const MDMap & source);

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
  MDMap(const MDMap & parent,
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
  MDMap(const MDMap & parent,
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
  MDMap(const MDMap & parent,
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
  MDMap & operator=(const MDMap & source);

  //@}

  /** \name MDComm accessor and pass-through methods */
  //@{

  /** \brief Get the Teuchos communicator
   *
   * Note that if the communicator is not a full communicator, i.e. a
   * sub-communicator, that the underlying Comm pointer may be NULL,
   * depending on this processor's rank.
   */
  inline Teuchos::RCP< const Teuchos::Comm< int > > getTeuchosComm() const;

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
  inline Teuchos::Array< dim_type > getGlobalDims() const;

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

  //@}

  /** \name Conversions to other Maps */
  //@{

  /** \brief Return an array of axis maps
   *
   * An axis map is a 1D map along a given axis.
   */
  Teuchos::ArrayView< Teuchos::RCP< const Domi::MDMap > >
  getAxisMaps() const;

  /** \brief Return an axis map for the given axis
   *
   * An axis map is a 1D map along a given axis.
   */
  Teuchos::RCP< const Domi::MDMap > getAxisMap(int axis) const;

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
  Teuchos::RCP< const MDMap >
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
   *         MDMap, specifying the new LocalOrdinal and GlobalOrdinal
   *         types
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  template< class LocalOrdinal  = TpetraLOType,
            class GlobalOrdinal = TpetraGOType,
            class Node          = TpetraNodeType>
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >
  getTpetraMap(bool withCommPad=true) const;

  /** \brief Return an RCP to a Tpetra::Map that represents the
   *         decomposition of this MDMap along the given axis,
   *         specifying the new LocalOrdinal and GlobalOrdinal
   *         types
   *
   * \param axis [in] the requested axis
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  template< class LocalOrdinal  = TpetraLOType,
            class GlobalOrdinal = TpetraGOType,
            class Node          = TpetraNodeType >
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >
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
  bool isCompatible(const MDMap & mdMap) const;

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
  bool isSameAs(const MDMap & mdMap,
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
  mutable Teuchos::Array< Teuchos::RCP< const MDMap > > _axisMaps;

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

////////////////////////////
// Inline implementations //
////////////////////////////

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Teuchos::Comm< int > >
MDMap::getTeuchosComm() const
{
  return _mdComm->getTeuchosComm();
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const MDComm >
MDMap::getMDComm() const
{
  return _mdComm;
}

////////////////////////////////////////////////////////////////////////

bool
MDMap::onSubcommunicator() const
{
  return _mdComm->onSubcommunicator();
}

////////////////////////////////////////////////////////////////////////

int
MDMap::numDims() const
{
  return _mdComm->numDims();
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
MDMap::getCommDims() const
{
  return _mdComm->getCommDims();
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getCommDim(int axis) const
{
  return _mdComm->getCommDim(axis);
}

////////////////////////////////////////////////////////////////////////

bool
MDMap::isPeriodic(int axis) const
{
  return _mdComm->isPeriodic(axis);
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getCommIndex(int axis) const
{
  return _mdComm->getCommIndex(axis);
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getLowerNeighbor(int axis) const
{
  return _mdComm->getLowerNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getUpperNeighbor(int axis) const
{
  return _mdComm->getUpperNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< dim_type >
MDMap::
getGlobalDims() const
{
  return _globalDims;
}

////////////////////////////////////////////////////////////////////////

//////////////////////////////////////
// Templated method implementations //
//////////////////////////////////////

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node >
Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >
MDMap::getTpetraMap(bool withCommPad) const
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
                                    Node >(Teuchos::OrdinalTraits< Tpetra::global_size_t >::invalid(),
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
                                    Node >(Teuchos::OrdinalTraits< Tpetra::global_size_t>::invalid(),
                                           myElements(),
                                           0,
                                           teuchosComm));
  }
}
#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node >
Teuchos::RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >
MDMap::
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
                                       Node >(Teuchos::OrdinalTraits< Tpetra::global_size_t>::invalid(),
                                              elements,
                                              0,
                                              teuchosComm));
}
#endif

////////////////////////////////////////////////////////////////////////

}

#endif
