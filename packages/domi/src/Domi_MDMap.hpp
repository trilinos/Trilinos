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

// System includes
#include <limits>

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_getValidParameters.hpp"
#include "Domi_Slice.hpp"
#include "Domi_MDComm.hpp"
#include "Domi_MDArray.hpp"

// Teuchos includes
#include "Teuchos_Tuple.hpp"

// Kokkos includes
#include "Kokkos_DefaultNode.hpp"

#ifdef HAVE_EPETRA
#include "Epetra_Map.h"
#endif

#ifdef HAVE_TPETRA
#include "Tpetra_Map.hpp"
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
 * From a global perspective, the only points that are visible are
 * boundary padding.  Communication padding is internal and
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
template< class LocalOrd,
          class GlobalOrd = LocalOrd,
          class Node = Kokkos::DefaultNode::DefaultNodeType >
class MDMap
{
public:

  /** \brief Adopt the MDArray <tt>sizetype</tt> type
   */
  typedef typename MDArray< LocalOrd >::size_type size_type;

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
   * \param layout [in] the storage order of the map
   *
   * \param node [in] the Kokkos node of the map
   */
  MDMap(const MDCommRCP mdComm,
        const Teuchos::ArrayView< GlobalOrd > & dimensions,
        const Teuchos::ArrayView< int > & commPad =
          Teuchos::ArrayView< int >(),
        const Teuchos::ArrayView< int > & bndryPad =
          Teuchos::ArrayView< int >(),
        const ELayout layout = DEFAULT_ORDER,
        const Teuchos::RCP< Node > & node =
          Kokkos::DefaultNode::getDefaultNode());

  /** \brief Constructor with Teuchos::Comm and ParameterList
   *
   * \param teuchosComm [in] The Teuchos Communicator.  Note that an
   *        MDComm will be constructed from the information in plist.
   *
   * \param plist [in] ParameterList with construction information
   * \htmlonly
   * <iframe src="domi.xml" width="90%"height="400px"></iframe>
   * <hr />
   * \endhtmlonly
   *
   * \param node [in] the Kokkos node of the map
   */
  MDMap(const TeuchosCommRCP teuchosComm,
        Teuchos::ParameterList & plist,
        const Teuchos::RCP< Node > & node =
          Kokkos::DefaultNode::getDefaultNode());

  /** \brief Constructor with MDComm and ParameterList
   *
   * \param mdComm [in] an RCP of an MDComm (multi-dimensional
   *        communicator), on which this MDMap will be built.
   *
   * \param plist [in] ParameterList with construction information
   * \htmlonly
   * <iframe src="domi.xml" width="90%"height="400px"></iframe>
   * <hr />
   * \endhtmlonly
   *
   * \param node [in] the Kokkos node of the map
   */
  MDMap(const MDCommRCP mdComm,
        Teuchos::ParameterList & plist,
        const Teuchos::RCP< Node > & node =
          Kokkos::DefaultNode::getDefaultNode());

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
  MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
        int axis,
        GlobalOrd index);

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
  MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
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
  MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
        const Teuchos::ArrayView< Slice > & slices,
        const Teuchos::ArrayView< int > & bndryPad =
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
   * using the (parent,ordinal) or (parent,slices) constructor.  For a
   * full communicator, this method will always return true.
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
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
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
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is zero, then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is zero, then the returned lower
   * neighbor will be the highest axis rank processor along this axis.
   */
  int getLowerNeighbor(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is the highest axis rank along the axis,
   * then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is the highest axis rank along the axis,
   * then this method will return axis rank zero.
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
   * \param withBndryPad [in] specify whether the dimension should
   *        include boundary padding or not
   */
  GlobalOrd getGlobalDim(int axis, bool withBndryPad=false) const;

  /** \brief Get the bounds of the global problem
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the bounds should
   *        include boundary padding or not
   */
  Slice getGlobalBounds(int axis, bool withBndryPad=false) const;

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
  Slice getGlobalRankBounds(int axis, bool withBndryPad=false) const;

  /** \brief Get the local dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withPad [in] specify whether the dimension should include
   *        padding or not
   */
  LocalOrd getLocalDim(int axis, bool withPad=false) const;

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
  Slice getLocalBounds(int axis, bool withPad=false) const;

  //@}

  /** \name Storage order and communication and boundary padding */
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
  int getLowerPad(int axis) const;

  /** \brief Get the size of the upper padding along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that the returned padding can be either communication
   * padding or boundary padding as appropriate.
   */
  int getUpperPad(int axis) const;

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

  /** \brief Get the storage order
   */
  ELayout getLayout() const;

  //@}

  /** \name Conversions to other Maps */
  //@{

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
   *         MDMap
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
  getTpetraMap(bool withCommPad=true) const;

  /** \brief Return an RCP to a Tpetra::Map that represents the
   *         decomposition of this MDMap along the given axis
   *
   * \param axis [in] the requested axis
   *
   * \param withCommPad [in] flag whether to include the communication
   *        padding in the map
   *
   * Note that the boundary padding is always included in the map
   */
  Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
  getTpetraAxisMap(int axis,
                   bool withCommPad=true) const;

#endif

  //@}

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
  getLocalAxisIndex(LocalOrd localIndex) const;

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

  // Typedef for the exterior buffer padding type, a tuple that stores
  // lower and upper padding sizes, which can be different due to
  // processors on domain boundaries.
  typedef Teuchos::Tuple< int, 2 > padding;

  // A private method for computing the bounds and local dimensions,
  // after the global dimensions, communication and boundary padding
  // have been properly assigned.
  void computeBounds();

  // The underlying multi-dimensional communicator.
  MDCommRCP _mdComm;

  // The size of the global dimensions along each axis.  This includes
  // the values of the boundary padding.
  Teuchos::Array< GlobalOrd > _globalDims;

  // Store the start and stop indexes of the global problem along each
  // axis.  This includes the values of the boundary padding.
  Teuchos::Array< Slice > _globalBounds;

  // A double array of Slices that stores the starting and stopping
  // indexes for each axis processor along each axis.  This data
  // structure allows any processor to know the map bounds on any
  // other processor.  These bounds do NOT include the boundary
  // padding.  The outer array will have a size equal to the number of
  // dimensions.  Each inner array will have a size equal to the
  // number of axis processors along that axis.
  Teuchos::Array< Teuchos::Array< Slice > > _globalRankBounds;

  // The global stride between adjacent IDs.  This quantity is
  // "virtual" as it does not describe actual memory layout, but it is
  // useful for index conversion algorithms.
  Teuchos::Array< GlobalOrd > _globalStrides;

  // The minumum 1D index of the global data structure, including
  // boundary padding.  This will only be non-zero on a sub-map.
  GlobalOrd _globalMin;

  // The maximum 1D index of the global data structure, including
  // boundary padding.
  GlobalOrd _globalMax;

  // The size of the local dimensions along each axis.  This includes
  // the values of the communication padding.
  Teuchos::Array< LocalOrd > _localDims;

  // The local loop bounds along each axis, stored as an array of
  // Slices.  These bounds DO include the communication padding.
  Teuchos::Array< Slice > _localBounds;

  // The local stride between adjacent elements in memory.
  Teuchos::Array< LocalOrd > _localStrides;

  // The minimum 1D index of the local data structure, including
  // communication padding.
  LocalOrd _localMin;

  // The maximum 1D index of the local data structure, including
  // communnication padding.
  LocalOrd _localMax;

  // The communication padding that was specified at construction, one
  // value along each axis.
  Teuchos::Array< int > _commPadSizes;

  // The actual padding stored on this processor along each axis.  The
  // padding can be either communication padding or boundary padding
  // based on the processor position on the boundary of a domain.
  Teuchos::Array< padding > _pad;

  // The size of the boundary padding along each axis.
  Teuchos::Array< int > _bndryPadSizes;

  // The actual boundary padding stored along each axis.  For a full
  // communicator, the lower and upper boundary padding are both the
  // same as the corresponding value in _bndryPadSizes.  However, the
  // introduction of sub-maps creates the possibility of different
  // upper and lower boundary padding values.
  Teuchos::Array< padding > _bndryPad;

  // The storage order
  ELayout _layout;

  // The Kokkos node type
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
  // distribution of axis IDs along each axis, including communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array< Teuchos::RCP< const Epetra_Map > > _epetraAxisMaps;

  // An ArrayRCP that stores Epetra_Maps that represent the
  // distribution of axis IDs along each axis, excluding communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array< Teuchos::RCP< const Epetra_Map > > _epetraAxisOwnMaps;
#endif

#ifdef HAVE_TPETRA
  // An RCP pointing to a Tpetra::Map that is equivalent to this
  // MDMap, including communication padding.  It is mutable because we
  // do not construct it until it asked for by a get method that is
  // const.
  mutable Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
  _tpetraMap;

  // An RCP pointing to a Tpetra::Map that is equivalent to this
  // MDMap, excluding communication padding.  The returned Tpetra::Map
  // thus indicates what IDs are owned by this processor.  It is
  // mutable because we do not construct it until it asked for by a
  // get method that is const.
  mutable Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
  _tpetraOwnMap;

  // An ArrayRCP that stores Tpetra::Maps that represent the
  // distribution of axis IDs along each axis, including communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array<
    Teuchos::RCP<
      const Tpetra::Map< LocalOrd, GlobalOrd, Node > > > _tpetraAxisMaps;

  // An ArrayRCP that stores Tpetra::Maps that represent the
  // distribution of axis IDs along each axis, excluding communication
  // padding.  It is mutable because we do not construct it until it
  // asked for by a get method that is const.
  mutable Teuchos::Array<
    Teuchos::RCP<
      const Tpetra::Map< LocalOrd, GlobalOrd, Node > > > _tpetraAxisOwnMaps;
#endif

};

////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const MDCommRCP mdComm,
      const Teuchos::ArrayView< GlobalOrd > & dimensions,
      const Teuchos::ArrayView< int > & commPad,
      const Teuchos::ArrayView< int > & bndryPad,
      const ELayout layout,
      const Teuchos::RCP< Node > & node) :
  _mdComm(mdComm),
  _globalDims(mdComm->getNumDims()),
  _globalBounds(),
  _globalRankBounds(mdComm->getNumDims()),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(mdComm->getNumDims(), 0),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(),
  _commPadSizes(mdComm->getNumDims(), 0),
  _pad(),
  _bndryPadSizes(mdComm->getNumDims(), 0),
  _bndryPad(),
  _layout(layout),
  _node(node)
{
  // Temporarily store the number of dimensions
  int numDims = mdComm->getNumDims();

  // Check the global dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");

  // Copy the boundary padding sizes and compute the global dimensions
  // and bounds
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < bndryPad.size())
      _bndryPadSizes[axis] = bndryPad[axis];
    _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                       _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + 2*_bndryPadSizes[axis];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Copy the communication padding sizes and set the actual padding
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPad.size())
      _commPadSizes[axis] = commPad[axis];
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

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  computeBounds();
  _localMax = computeSize(_localDims());

  // Compute the strides
  _globalStrides = computeStrides(_globalDims, _layout);
  _localStrides  = computeStrides(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node>::
MDMap(TeuchosCommRCP teuchosComm,
      Teuchos::ParameterList & plist,
      const Teuchos::RCP< Node > & node) :
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
  _localMax(),
  _commPadSizes(),
  _pad(),
  _bndryPadSizes(),
  _bndryPad(),
  _layout(),
  _node(node)
{
  // Note that the call to the MDComm constructor in the constructor
  // initialization list will validate the ParameterList, so we don't
  // have to do it again here.

  // Temporarily store the number of dimensions
  int numDims = _mdComm->getNumDims();

  // Check the global dimensions
  Teuchos::Array< long long > dimensions =
    plist.get("dimensions", Teuchos::Array< long long >());
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");
  TEUCHOS_TEST_FOR_EXCEPTION(
    computeSize(dimensions) - 1 > std::numeric_limits< GlobalOrd >::max(),
    MapOrdinalError,
    "The maximum global ID of this MDMap is too large for GlobalOrd");

  // Copy the boundary padding sizes and compute the global dimensions
  // and bounds
  Teuchos::Array< int > bndryPad =
    plist.get("boundary pad", Teuchos::Array< int >());
  _bndryPadSizes.resize(numDims);
  _globalDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < bndryPad.size())
      _bndryPadSizes[axis] = bndryPad[axis];
    _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                       _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + 2*_bndryPadSizes[axis];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Copy the communication padding sizes and set the actual padding
  Teuchos::Array< int > commPad =
    plist.get("communication pad", Teuchos::Array< int >());
  _commPadSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPad.size())
      _commPadSizes[axis] = commPad[axis];
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

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  _globalRankBounds.resize(numDims);
  _localDims.resize(numDims);
  computeBounds();
  _localMax = computeSize(_localDims());

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
  _globalStrides = computeStrides(_globalDims, _layout);
  _localStrides  = computeStrides(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node>::
MDMap(MDCommRCP mdComm,
      Teuchos::ParameterList & plist,
      const Teuchos::RCP< Node > & node) :
  _mdComm(mdComm),
  _globalDims(mdComm->getNumDims()),
  _globalBounds(),
  _globalRankBounds(mdComm->getNumDims()),
  _globalStrides(),
  _globalMin(0),
  _globalMax(),
  _localDims(mdComm->getNumDims(), 0),
  _localBounds(),
  _localStrides(),
  _localMin(0),
  _localMax(),
  _commPadSizes(mdComm->getNumDims(), 0),
  _pad(),
  _bndryPadSizes(mdComm->getNumDims(), 0),
  _bndryPad(),
  _layout(),
  _node(node)
{
  // Validate the ParameterList
  plist.validateParameters(*getValidParameters());

  // Temporarily store the number of dimensions
  int numDims = _mdComm->getNumDims();

  // Check the global dimensions
  Teuchos::Array< long long > dimensions =
    plist.get("dimensions", Teuchos::Array< long long >());
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != dimensions.size(),
    InvalidArgument,
    "Size of dimensions does not match MDComm number of dimensions");
  TEUCHOS_TEST_FOR_EXCEPTION(
    computeSize(dimensions) - 1 > std::numeric_limits< GlobalOrd >::max(),
    MapOrdinalError,
    "The maximum global ID of this MDMap is too large for GlobalOrd");

  // Copy the boundary padding sizes and compute the global dimensions
  // and bounds
  Teuchos::Array< int > bndryPad =
    plist.get("boundary pad", Teuchos::Array< int >());
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < bndryPad.size())
      _bndryPadSizes[axis] = bndryPad[axis];
    _bndryPad.push_back(Teuchos::tuple(_bndryPadSizes[axis],
                                       _bndryPadSizes[axis]));
    _globalDims[axis] = dimensions[axis] + 2*_bndryPadSizes[axis];
    _globalBounds.push_back(ConcreteSlice(_globalDims[axis]));
  }

  // Compute the global size
  _globalMax = computeSize(_globalDims());

  // Copy the communication padding sizes and set the actual padding
  Teuchos::Array< int > commPad =
    plist.get("communication pad", Teuchos::Array< int >());
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < commPad.size())
      _commPadSizes[axis] = commPad[axis];
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

  // Compute _globalRankBounds, _localBounds, and _localDims.
  // Then compute the local size
  computeBounds();
  _localMax = computeSize(_localDims());

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
  _globalStrides = computeStrides(_globalDims, _layout);
  _localStrides  = computeStrides(_localDims , _layout);
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
      int axis,
      GlobalOrd index) :
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
  _layout(parent._layout),
  _node(parent._node)
{
  if (parent.onSubcommunicator())
  {
    int numDims = parent.getNumDims();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for communicator with " <<
        numDims << " dimensions");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ((index < parent.getGlobalBounds(axis,true).start()) ||
       (index >= parent.getGlobalBounds(axis,true).stop())),
      RangeError,
      "index = " << index  << " is invalid for MDMap axis " <<
      axis << " with bounds " << parent.getGlobalBounds(axis,true));

    // Determine the axis rank for the processor on which the index
    // lives, and construct the MDComm
    int thisAxisRank = -1;
    for (int axisRank = 0; axisRank < parent.getAxisCommSize(axis); ++axisRank)
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
    int numDims = parent.getNumDims();
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
      _localMax = _localMin;
      _commPadSizes.push_back(0);
      _pad.push_back(Teuchos::tuple(0,0));
      _bndryPadSizes.push_back(0);
      _bndryPad.push_back(Teuchos::tuple(0,0));
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
        }
        else
        {
          int axisRank = parent.getAxisRank(axis);
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

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
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
  _layout(parent._layout),
  _node(parent._node)
{
  if (parent.onSubcommunicator())
  {
    // Temporarily store the number of dimensions
    int numDims = parent.getNumDims();

    // Sanity check
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for MDMap with " <<
        numDims << " dimensions");

    // Copy the boundary padding sizes and set initial values for
    // _bndryPad
    _bndryPadSizes[axis] = bndryPad;
    _bndryPad[axis][0]   = bndryPad;
    _bndryPad[axis][1]   = bndryPad;

    // Convert the slice to concrete and check
    Slice bounds =
      slice.bounds(parent.getGlobalBounds(axis,true).stop());
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((bounds.start() < parent.getGlobalBounds(axis).start()) ||
       (bounds.stop() > parent.getGlobalBounds(axis).stop())),
      RangeError,
      "Slice along axis " << axis << " is " << bounds << " but must be within "
      << parent.getGlobalBounds(axis));
    // Adjust _globalRankBounds
    for (int axisRank = 0; axisRank < parent.getAxisCommSize(axis);
         ++axisRank)
    {
      typename Slice::Ordinal start = _globalRankBounds[axis][axisRank].start();
      typename Slice::Ordinal stop  = _globalRankBounds[axis][axisRank].stop();
      if (start < bounds.start()) start = bounds.start();
      if (stop  > bounds.stop() ) stop  = bounds.stop();
      _globalRankBounds[axis][axisRank] = ConcreteSlice(start, stop);
    }
    // Alter _bndryPad if necessary
    GlobalOrd start = bounds.start() - _bndryPadSizes[axis];
    if (start < 0)
    {
      _bndryPad[axis][0] = bounds.start();
      start = 0;
    }
    GlobalOrd stop = bounds.stop() + _bndryPadSizes[axis];
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

    // Build the slice for the MDComm sub-communicator constructor
    int pStart = -1;
    int pStop  = -1;
    for (int axisRank = 0; axisRank < parent.getAxisCommSize(axis);
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
      int axisRank = getAxisRank(axis);
      if (axisRank == 0)
        _pad[axis][0] = _bndryPad[axis][0];
      if (axisRank == _mdComm->getAxisCommSize(axis)-1)
        _pad[axis][1] = _bndryPad[axis][1];
      LocalOrd start = parent._localBounds[axis].start();
      if (_globalBounds[axis].start() >
          _globalRankBounds[axis][axisRank].start())
      {
        start = _globalBounds[axis].start() -
          _globalRankBounds[axis][axisRank].start();
      }
      else
      {
        if (_globalRankBounds[axis][axisRank].start() -
            _globalBounds[axis].start() < _pad[axis][0])
        {
          _pad[axis][0] = _globalRankBounds[axis][axisRank].start() -
            _globalBounds[axis].start();
        }
      }
      LocalOrd stop = parent._localBounds[axis].stop();
      if (_globalBounds[axis].stop() <
          _globalRankBounds[axis][axisRank].stop())
      {
        stop = _globalBounds[axis].stop() -
          _globalRankBounds[axis][axisRank].start();
      }
      else
      {
        if (_globalBounds[axis].stop() -
            _globalRankBounds[axis][axisRank].stop() < _pad[axis][1])
        {
          _pad[axis][1] = _globalBounds[axis].stop() -
            _globalRankBounds[axis][axisRank].stop();
        }
      }
      _localBounds[axis] = ConcreteSlice(start,stop);
      _localDims[axis]   = stop - start;
      _localMin         += start * _localStrides[axis];
      _localMax         -= (parent.getLocalBounds(axis,true).stop() - stop) *
                           _localStrides[axis];
      // The new sub-communicator may have fewer processors than the
      // parent communicator, so we need to fix _globalRankBounds
      Teuchos::Array< Slice > newRankBounds;
      for (int axisRank = 0; axisRank < parent.getAxisCommSize(axis);
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

template< class LocalOrd, class GlobalOrd, class Node >
MDMap< LocalOrd, GlobalOrd, Node >::
MDMap(const MDMap< LocalOrd, GlobalOrd, Node > & parent,
      const Teuchos::ArrayView< Slice > & slices,
      const Teuchos::ArrayView< int > & bndryPad)
{
  // Temporarily store the number of dimensions
  int numDims = parent.getNumDims();

  // Sanity check on dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    (slices.size() != numDims),
    InvalidArgument,
    "number of slices = " << slices.size() << " != parent MDMap number of "
    "dimensions = " << numDims);

  // Apply the single-Slice constructor to each axis in succession
  MDMap< LocalOrd, GlobalOrd, Node > tempMDMap1(parent);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int bndryPadding = (axis < bndryPad.size()) ? bndryPad[axis] : 0;
    MDMap< LocalOrd, GlobalOrd, Node > tempMDMap2(tempMDMap1,
                                                  axis,
                                                  slices[axis],
                                                  bndryPadding);
    tempMDMap1 = tempMDMap2;
  }
  *this = tempMDMap1;
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
             bool withBndryPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withBndryPad)
    return _globalDims[axis];
  else
    return _globalDims[axis] - _bndryPad[axis][0] - _bndryPad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalBounds(int axis,
                bool withBndryPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withBndryPad)
    return _globalBounds[axis];
  else
  {
    GlobalOrd start = _globalBounds[axis].start() + _bndryPad[axis][0];
    GlobalOrd stop  = _globalBounds[axis].stop()  - _bndryPad[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
LocalOrd
MDMap< LocalOrd, GlobalOrd, Node >::
getLocalDim(int axis,
            bool withPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withPad)
    return _localDims[axis];
  else
    return _localDims[axis] - _pad[axis][0] - _pad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Slice
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalRankBounds(int axis,
                    bool withBndryPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  int axisRank = getAxisRank(axis);
  if (withBndryPad)
  {
    GlobalOrd start = _globalRankBounds[axis][axisRank].start();
    GlobalOrd stop  = _globalRankBounds[axis][axisRank].stop();
    if (getAxisRank(axis) == 0)
      start -= _bndryPad[axis][0];
    if (getAxisRank(axis) == getAxisCommSize(axis)-1)
      stop += _bndryPad[axis][1];
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
               bool withPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (withPad)
    return _localBounds[axis];
  else
  {
    LocalOrd start = _localBounds[axis].start() + _pad[axis][0];
    LocalOrd stop  = _localBounds[axis].stop()  - _pad[axis][1];
    return ConcreteSlice(start, stop);
  }
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
bool
MDMap< LocalOrd, GlobalOrd, Node >::hasPadding() const
{
  bool result = false;
  for (int axis = 0; axis < getNumDims(); ++axis)
    if (_pad[axis][0] + _pad[axis][1]) result = true;
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getLowerPad(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _pad[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperPad(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _pad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getCommPadSize(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _commPadSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getLowerBndryPad(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _bndryPad[axis][0];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getUpperBndryPad(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _bndryPad[axis][1];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
int
MDMap< LocalOrd, GlobalOrd, Node >::getBndryPadSize(int axis) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _bndryPadSizes[axis];
}

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
ELayout
MDMap< LocalOrd, GlobalOrd, Node >::getLayout() const
{
  return _layout;
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::RCP< const Epetra_Map >
MDMap< LocalOrd, GlobalOrd, Node >::getEpetraMap(bool withCommPad) const
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
      int numDims = getNumDims();
      Teuchos::Array<size_type> localDims(numDims);
      for (int axis = 0; axis < numDims; ++axis)
        localDims[axis] = _localDims[axis];
      MDArray<int> myElements(localDims);
      Teuchos::Array<int> index(numDims);

      // Iterate over the local MDArray and assign global IDs
      for (MDArray<int>::iterator it = myElements.begin();
           it != myElements.end(); ++it)
      {
        int globalID = 0;
        for (int axis = 0; axis < numDims; ++axis)
        {
          int axisRank = getAxisRank(axis);
          int start    = _globalRankBounds[axis][axisRank].start() -
                         _pad[axis][0];
          globalID += (start + it.index(axis)) * _globalStrides[axis];
        }
        *it = globalID;
      }

      // Construct the Epetra_Map
      EpetraCommRCP epetraComm = _mdComm->getEpetraComm();
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
      int numDims = getNumDims();
      Teuchos::Array<int> index(numDims);
      Teuchos::Array<size_type> myDims(numDims);
      for (int axis = 0; axis < numDims; ++axis)
      {
        myDims[axis] = _localDims[axis] - _pad[axis][0] - _pad[axis][1];
        int axisRank = getAxisRank(axis);
        if (axisRank == 0)
          myDims[axis] += _bndryPad[axis][0];
        if (axisRank == getAxisCommSize(axis)-1)
          myDims[axis] += _bndryPad[axis][1];
      }
      MDArray<int> myElements(myDims());

      // Iterate over the local MDArray and assign global IDs
      for (MDArray<int>::iterator it = myElements.begin();
           it != myElements.end(); ++it)
      {
        int globalID = 0;
          for (int axis = 0; axis < numDims; ++axis)
          {
            int axisRank = getAxisRank(axis);
            int start    = _globalRankBounds[axis][axisRank].start();
            if (axisRank == 0)
              start -= _bndryPad[axis][0];
            if (axisRank == getAxisCommSize(axis)-1)
              start += _bndryPad[axis][1];
            globalID += (start + it.index(axis)) * _globalStrides[axis];
          }
      }

      // Construct the Epetra_Map
      EpetraCommRCP epetraComm = _mdComm->getEpetraComm();
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

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::RCP< const Epetra_Map >
MDMap< LocalOrd, GlobalOrd, Node >::
getEpetraAxisMap(int axis,
                 bool withCommPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if ((withCommPad     && (_epetraAxisMaps.size()    == 0)) ||
      (not withCommPad && (_epetraAxisOwnMaps.size() == 0)))
  {
    int numDims = getNumDims();
    EpetraCommRCP epetraComm = _mdComm->getEpetraComm();
    for (int axis=0; axis < numDims; ++axis)
    {
      Teuchos::Array<int> elements(getLocalDim(axis, withCommPad));
      int start = getGlobalRankBounds(axis,true).start();
      if (withCommPad && (getAxisRank(axis) != 0)) start -= _pad[axis][0];
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

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
MDMap< LocalOrd, GlobalOrd, Node >::getTpetraMap(bool withCommPad) const
{
  if (withCommPad)
  {
    if (_tpetraMap.is_null())
    {
      // Allocate the elementsMDArray and the index array
      int numDims = getNumDims();
      Teuchos::Array<size_type> localDims(numDims);
      for (int axis = 0; axis < numDims; ++axis)
        localDims[axis] = _localDims[axis];
      MDArray<GlobalOrd> elementMDArray(localDims);
      Teuchos::Array<LocalOrd> index(numDims);

      // Iterate over the local MDArray and assign global IDs
      for (typename MDArray<GlobalOrd>::iterator it = elementMDArray.begin();
           it != elementMDArray.end(); ++it)
      {
        GlobalOrd globalID = 0;
        for (int axis = 0; axis < numDims; ++axis)
        {
          int axisRank    = getAxisRank(axis);
          GlobalOrd start = _globalRankBounds[axis][axisRank].start() -
                            _pad[axis][0];
          globalID += (start + it.index(axis)) * _globalStrides[axis];
        }
        *it = globalID;
      }

      // Construct the Tpetra::Map
      const Teuchos::Array<GlobalOrd> & myElements = elementMDArray.array();
      TeuchosCommRCP teuchosComm = _mdComm->getTeuchosComm();
      _tpetraMap =
        Teuchos::rcp(new Tpetra::Map<LocalOrd,
                                     GlobalOrd,
                                     Node>(Teuchos::OrdinalTraits<
                                             Tpetra::global_size_t>::invalid(),
                                           myElements(),
                                           0,
                                           teuchosComm));
    }
    return _tpetraMap;
  }
  else
  {
    if (_tpetraOwnMap.is_null())
    {
      // Allocate the elementMDArray MDArray and the index array
      int numDims = getNumDims();
      Teuchos::Array<LocalOrd> index(numDims);
      Teuchos::Array<size_type> myDims(numDims);
      for (int axis = 0; axis < numDims; ++axis)
      {
        myDims[axis] =
          _localDims[axis] - _pad[axis][0] - _pad[axis][1];
        int axisRank = getAxisRank(axis);
        if (axisRank == 0)
          myDims[axis] += _bndryPad[axis][0];
        if (axisRank == getAxisCommSize(axis)-1)
          myDims[axis] += _bndryPad[axis][1];
      }
      MDArray<GlobalOrd> elementMDArray(myDims());

      // Iterate over the local MDArray and assign global IDs
      for (typename MDArray<GlobalOrd>::iterator it = elementMDArray.begin();
           it != elementMDArray.end(); ++it)
      {
        GlobalOrd globalID = 0;
          for (int axis = 0; axis < numDims; ++axis)
          {
            int axisRank    = getAxisRank(axis);
            GlobalOrd start = _globalRankBounds[axis][axisRank].start();
            if (axisRank == 0)
              start -= _bndryPad[axis][0];
            if (axisRank == getAxisCommSize(axis)-1)
              start += _bndryPad[axis][1];
            globalID += (start + it.index(axis)) * _globalStrides[axis];
          }
      }

      // Construct the Tpetra::Map
      const Teuchos::Array<GlobalOrd> & myElements = elementMDArray.array();
      TeuchosCommRCP teuchosComm = _mdComm->getTeuchosComm();
      _tpetraOwnMap =
        Teuchos::rcp(new Tpetra::Map<LocalOrd,
                                     GlobalOrd,
                                     Node>(Teuchos::OrdinalTraits<
                                             Tpetra::global_size_t>::invalid(),
                                           myElements(),
                                           0,
                                           teuchosComm));
    }
    return _tpetraOwnMap;
  }
}
#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::RCP< const Tpetra::Map< LocalOrd, GlobalOrd, Node > >
MDMap< LocalOrd, GlobalOrd, Node >::
getTpetraAxisMap(int axis,
                 bool withCommPad) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if ((withCommPad     && (_tpetraAxisMaps.size()    == 0)) ||
      (not withCommPad && (_tpetraAxisOwnMaps.size() == 0)))
  {
    int numDims = getNumDims();
    TeuchosCommRCP teuchosComm = _mdComm->getTeuchosComm();
    for (int axis=0; axis < numDims; ++axis)
    {
      Teuchos::Array<GlobalOrd> elements(getLocalDim(axis,withCommPad));
      GlobalOrd start = getGlobalRankBounds(axis,true).start();
      if (withCommPad && (getAxisRank(axis) != 0)) start -= _pad[axis][0];
      for (LocalOrd i = 0; i < elements.size(); ++i)
        elements[i] = i + start;
      if (withCommPad)
      {
        _tpetraAxisMaps.push_back(
          Teuchos::rcp(new Tpetra::Map<LocalOrd,
                                       GlobalOrd,
                                       Node>(Teuchos::OrdinalTraits<
                                               Tpetra::global_size_t>::invalid(),
                                             elements,
                                             0,
                                             teuchosComm)));
      }
      else
      {
        _tpetraAxisOwnMaps.push_back(
          Teuchos::rcp(new Tpetra::Map<LocalOrd,
                                       GlobalOrd,
                                       Node>(Teuchos::OrdinalTraits<
                                               Tpetra::global_size_t>::invalid(),
                                             elements,
                                             0,
                                             teuchosComm)));
      }
    }
  }
  if (withCommPad)
    return _tpetraAxisMaps[axis];
  else
    return _tpetraAxisOwnMaps[axis];
}
#endif

////////////////////////////////////////////////////////////////////////

template< class LocalOrd, class GlobalOrd, class Node >
Teuchos::Array< GlobalOrd >
MDMap< LocalOrd, GlobalOrd, Node >::
getGlobalAxisIndex(GlobalOrd globalIndex) const
{
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((globalIndex < _globalMin) || (globalIndex >= _globalMax)),
    RangeError,
    "invalid global index = " << globalIndex << " (should be between " <<
    _globalMin << " and " << _globalMax << ")");
#endif
  int numDims = getNumDims();
  Teuchos::Array< GlobalOrd > result(numDims);
  GlobalOrd index = globalIndex;
  if (_layout == LAST_INDEX_FASTEST)
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
    ((localIndex < _localMin) || (localIndex >= _localMax)),
    RangeError,
    "invalid local index = " << localIndex << " (should be between " <<
    _localMin << " and " << _localMax << ")");
#endif
  int numDims = getNumDims();
  Teuchos::Array< LocalOrd > result(numDims);
  LocalOrd index = localIndex;
  if (_layout == LAST_INDEX_FASTEST)
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
      _globalRankBounds[axis][getAxisRank(axis)].start() - _pad[axis][0];
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
    ((globalIndex < _globalMin) || (globalIndex >= _globalMax)),
    RangeError,
    "invalid global index = " << globalIndex << " (should be between " <<
    _globalMin << " and " << _globalMax << ")");
#endif
  Teuchos::Array< GlobalOrd > globalAxisIndex =
    getGlobalAxisIndex(globalIndex);
  LocalOrd result = 0;
  for (int axis = 0; axis < getNumDims(); ++axis)
  {
    LocalOrd localAxisIndex = globalAxisIndex[axis] -
      _globalRankBounds[axis][getAxisRank(axis)].start() + _pad[axis][0];
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
      // communication and boundary padding.
      LocalOrd  localDim  = (_globalDims[axis] - 2*_bndryPadSizes[axis]) /
                            axisCommSize;
      GlobalOrd axisStart = axisRank * localDim;

      // Adjustments for non-zero remainder.  Compute the remainder
      // using the mod operator.  If the remainder is > 0, then add an
      // element to the appropriate number of processors with the
      // highest axis ranks.  Note that this is the opposite of the
      // standard Tpetra::Map constructor (which adds an elements to
      // the lowest processor ranks), and provides better balance for
      // finite differencing systems with staggered data location.
      GlobalOrd remainder = (_globalDims[axis] - 2*_bndryPadSizes[axis]) %
                            axisCommSize;
      if (axisCommSize - axisRank - 1 < remainder)
      {
        ++localDim;
        axisStart += (remainder - axisCommSize + axisRank);
      }

      // Global adjustment for boundary padding
      axisStart += _bndryPadSizes[axis];

      // Compute and store the global axis bounds
      _globalRankBounds[axis].push_back(
        ConcreteSlice(axisStart, axisStart + localDim));

      // Set _localDims[axis] and _localBounds[axis] only if
      // axisRank equals the axis rank of this processor
      if (axisRank == getAxisRank(axis))
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
