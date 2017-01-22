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

#ifndef DOMI_MDVECTOR_HPP
#define DOMI_MDVECTOR_HPP

// #define DOMI_MDVECTOR_VERBOSE

// Standard includes
#include <ctime>

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_DefaultNode.hpp"
#include "Domi_MDMap.hpp"
#include "Domi_MDArrayRCP.hpp"

// Teuchos includes
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"

#ifdef HAVE_EPETRA
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#endif

#ifdef HAVE_TPETRA
#include "Tpetra_Vector.hpp"
#endif


#ifdef OPEN_MPI
// I provide dummy definitions for MPI structs so that
// typeid(MPI_Datatype) and typeid(MPI_Request) will compile.
// Defining empty structs like this should be fine since only the guts
// of OpenMPI will ever see the real definitions.  This is a silly game
// we are playing with the C++ type system here but it should work
// just fine.
struct ompi_datatype_t {};
//struct ompi_request_t {};
#endif

using std::cout;
using std::endl;

namespace Domi
{

/** \brief Multi-dimensional distributed vector
 *
 * The <tt>MDVector</tt> class is intended to perform the functions of
 * Epetra or Tpetra MultiVectors (or Vectors), with the additional
 * capability that the data can be multi-dimensional.  (Note that
 * Epetra and Tpetra MultiVectors are 2 dimensional and distributed in
 * an unstructured manner along only one of those dimensions.
 * <tt>MDVector</tt>s can be any number of dimensions, and can be
 * distributed in a structured fashion along all of them.)
 *
 * Where Epetra and Tpetra MultiVectors are built on top of
 * <tt>Epetra_Map</tt>s and <tt>Tpetra::Map</tt>s, respectively, the
 * <tt>MDVector</tt> is built on top of an <tt>MDMap</tt>.  As such,
 * they support boundary padding, communication padding, and
 * periodicity.
 *
 * The basic <tt>MDVector</tt> constructor takes an MDMap and an
 * optional boolean flag providing instruction whether to zero out the
 * data or not.  A slightly more advanced constructor takes an MDMap,
 * a leading dimension, an optional trailing dimension, and an
 * optional zero-out flag.  This allows the user to augment the given
 * <tt>MDMap</tt> with an undistributed dimension.  This can be
 * helpful when the user wants to store an array of quantities at each
 * grid point.  These quantities can be clustered or strided in memory
 * depending on whether a leading or trailing dimension is specified,
 * and what the data Layout is.
 *
 * There are additional constructors as well, based on copying data,
 * using <tt>Teuchos::ParameterList</tt>s to specify the
 * <tt>MDVector</tt> layout, and parent/slice constructors, where
 * slicing of the parent can be performed with an ordinal or a Slice.
 *
 * The ability to slice a parent <tt>MDVector</tt> and obtain a
 * sub-vector is an important feature of the <tt>MDVector</tt> class.
 * When you take a slice of a parent <tt>MDVector</tt>, the resulting
 * sub-vector may have a new, reduced communicator (both a new
 * <tt>MDComm</tt> and <tt>Teuchos::Comm</tt>).  When accessing a
 * sub-vector, you should always use the onSubcommunicator() method to
 * determine that the sub-vector exists on the current processor.
 *
 * Like the <tt>Teuchos::MultiVector</tt>, the user cannot index
 * directly into an <tt>MDVector</tt>.  Rather, he or she must use the
 * <tt>getData()</tt> and <tt>getDataNonConst()</tt> methods to obtain
 * an <tt>MDArrayView</tt> of the data for manipulation.  This is a
 * little different than the Tpetra behavior, where an
 * <tt>ArrayRCP</tt> is returned to the user.  For <tt>MDVector</tt>,
 * and <tt>MDArrayView</tt> is returned, as this supports the slicing
 * semantics of the class.  Note that each <tt>MDVector</tt> stores an
 * internal <tt>MDArrayRCP</tt> of either its original data or its
 * parent's data, and then it also stores an <tt>MDArrayView</tt> that
 * is a view into the <tt>MDArrayRCP</tt>.  So the
 * <tt>MDArrayView</tt> passed to the user will be valid as long as
 * the reference count of its <tt>MDArrayRCP</tt> is greater than
 * zero.  Looping bounds for this data can be obtained using the
 * <tt>getLocalBounds()</tt> method, which returns a concrete
 * <tt>Slice</tt> object, and takes a boolean that specifies whether
 * padding should be included in the start and stop indexes.
 *
 * There are a variety of methods for converting <tt>MDVector</tt>s to
 * Epetra or Tpetra <tt>MultiVector</tt>s or <tt>Vector</tt>s, using
 * either view or copy semantics.  The view semantic is only available
 * when the <tt>MDVector</tt> is contiguous in memory, meaning there
 * are no stride gaps due to slicing.  The copy semantic is always
 * available.  True Epetra or Tpetra <tt>MultiVector</tt>s, in which
 * the number of vectors is greater than 1, are only available when
 * there is a non-distributed leading dimension and the layout is
 * C-order, or there is a non-distributed trailing dimension and the
 * layout is Fortran-order.
 *
 * If the <tt>MDVector</tt> was built with communication padding, then
 * the user may use the <tt>updateCommPad()</tt> method to communicate
 * via message passing data from the owning processors into the
 * communication padding buffers.  If communication along a single
 * axis is desired, then the overloaded <tt>updateCommPad(int
 * axis)</tt> method can be called instead.  If asynchronous
 * communication is desired, where computations can be performed
 * between posting and receiving messages, then the
 * <tt>startUpdateCommPad(int axis)</tt> and <tt>endUpdateCommPad(int
 * axis)</tt> methods can be called.  Note that the message data
 * structures needed to coordinate these methods are stored
 * internally.
 */
template< class Scalar,
          class Node = DefaultNode::DefaultNodeType >
class MDVector : public Teuchos::Describable
{
public:

  /** \name Constructors and destructor */
  //@{

  /** \brief Main constructor
   *
   * \param mdMap [in] MDMap that describes the domain decomposition
   *        of this MDVector
   *
   * \param zeroOut [in] flag to initialize all data to zero
   */
  MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
           bool zeroOut = true);

  /** \brief Augmented constructor
   *
   * \param mdMap [in] MDMap that describes the domain decomposition
   *        of this MDVector
   *
   * \param leadingDim [in] leading augmentation dimension, typically
   *        used to provide additional degrees of freedom at each
   *        index defined by the MDMap
   *
   * \param trailingDim [in] trailing augmentation dimension, typically
   *        used to provide additional degrees of freedom at each
   *        index defined by the MDMap
   *
   * \param zeroOut [in] flag to initialize all data to zero
   *
   * If leadingDim or trailingDim is less than 2, then the MDMap will
   * not be augmented with a leading dimension or trailing dimension,
   * respectively.  Note that this constructor takes the given MDMap
   * and computes a new, augmented MDComm and a new MDMap for this
   * MDVector.
   */
  MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
           const dim_type leadingDim,
           const dim_type trailingDim = 0,
           bool zeroOut = true);

  /** \brief Constructor with initialization values (copy)
   *
   * \param mdMap [in] MDMap that describes the domain decomposition
   *        of this MDVector
   *
   * \param source [in] initialization values
   */
  MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
           const MDArrayView< Scalar > & source);

  /** \brief Constructor with managed array object (view)
   *
   * \param mdMap [in] MDMap that describes the domain decomposition
   *        of this MDVector
   *
   * \param source [in] memory-managed multi-dimensional array
   */
  MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
           const MDArrayRCP< Scalar > & source);

  /** \brief Copy constructor with view or copy semantics
   *
   * \param source [in] source MDVector
   *
   * \param access [in] enumeration specifying view or copy
   *        construction
   *
   * Create a new MDVector with the same structure of the source
   * MDVector and either perform a deep copy or take a view of the
   * source MDVector data for the new MDVector.  Default behavior is
   * to use view semantics.
   */
  MDVector(const MDVector< Scalar, Node > & source,
           Teuchos::DataAccess access = Teuchos::View);

  /** \brief Constructor with Teuchos::Comm and ParameterList
   *
   * \param teuchosComm [in] The Teuchos Communicator.  Note that an
   *        MDComm and MDMap will be constructed from the information
   *        in plist.
   *
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   */
  MDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
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
  MDVector(const Teuchos::RCP< const MDComm > mdComm,
           Teuchos::ParameterList & plist);

  /** \brief Parent/single global ordinal sub-vector constructor
   *
   * \param parent [in] an MDVector, from which this MDVector will be
   *        derived
   *
   * \param axis [in] the axis to which this index ordinal applies
   *
   * \param index [in] the global ordinal that defines the sub-vector
   */
  MDVector(const MDVector< Scalar, Node > & parent,
           int axis,
           dim_type index);

  /** \brief Parent/single slice sub-vector constructor
   *
   * \param parent [in] an MDVector, from which this MDVector will be
   *        derived
   *
   * \param axis [in] the axis to which this slice applies
   *
   * \param slice [in] a slice of global ordinals describing the
   *        sub-vector
   *
   * \param bndryPad [in] the boundary padding along the altered axis
   *        of the new sub-vector.  This may include indexes from the
   *        boundary padding of the parent MDVector, but it does not
   *        have to.
   */
  MDVector(const MDVector< Scalar, Node > & parent,
           int                              axis,
           const Slice &                    slice,
           int                              bndryPad = 0);

  /** \brief Parent/array of slices sub-vector constructor
   *
   * \param parent [in] an MDVector, from which this sub-vector will
   *        be derived.
   *
   * \param slices [in] an array of Slices of global axis indexes that
   *        defines the sub-vector.  These slices must not include
   *        indexes from the boundary padding along each axis.
   *
   * \param bndryPad [in] The boundary padding of the new sub-vector.
   *        These may include indexes from the boundary padding of the
   *        parent MDVector, but they do not have to.
   */
  MDVector(const MDVector< Scalar, Node > &    parent,
           const Teuchos::ArrayView< Slice > & slices,
           const Teuchos::ArrayView< int > &   bndryPad =
             Teuchos::ArrayView< int >());

  /** \brief Assignment operator
   *
   * \param source [in] source MDVector to be copied
   */
  MDVector< Scalar, Node > &
  operator=(const MDVector< Scalar, Node > & source);

  /** \brief Destructor
   */
  virtual ~MDVector();

  //@}

  /** \name MDMap methods */
  //@{

  /** \brief MDMap accessor method
   */
  inline const Teuchos::RCP< const MDMap< Node > >
  getMDMap() const;

  /** \brief Query whether this processor is on the sub-communicator
   *
   * Sub-communicators are formed when a parent MDVector is sliced by
   * using the (parent,ordinal) or (parent,slice) constructors.  For a
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
  inline int getCommIndex(int axis) const;

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

  /** \brief Get the global dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the dimension should
   *        include boundary padding or not
   */
  inline dim_type getGlobalDim(int axis, bool withBndryPad=false) const;

  /** \brief Get the bounds of the global problem
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPad [in] specify whether the bounds should
   *        include boundary padding or not
   */
  inline Slice getGlobalBounds(int axis, bool withBndryPadding=false) const;

  /** \brief Get the global loop bounds on this processor along the
   *         specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withBndryPadding [in] specify whether the dimension should
   *        include boundary padding or not
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.
   */
  inline Slice getGlobalRankBounds(int axis, bool withBndryPad=false) const;

  /** \brief Get the local dimension along the specified axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * \param withCommPad [in] specify whether the dimension should
   *        include communication padding or not
   */
  inline dim_type getLocalDim(int axis, bool withCommPad=false) const;

  /** \brief Get the local loop bounds along every axis
   *
   * The loop bounds are returned in the form of a <tt>Slice</tt>, in
   * which the <tt>start()</tt> method returns the loop begin value,
   * and the <tt>stop()</tt> method returns the non-inclusive end
   * value.  For this method, padding is included in the bounds.
   */
  inline Teuchos::ArrayView< const Slice > getLocalBounds() const;

  /** \brief Get the local looping bounds along the specified axis
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
  inline Slice getLocalBounds(int axis, bool withPad=false) const;

  /** \brief Get the local interior looping bounds along the specified
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
  inline Slice getLocalInteriorBounds(int axis) const;

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
  inline bool hasPadding() const;

  /** \brief Get the size of the lower padding along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that the returned padding can be either communication
   * padding or boundary padding as appropriate.
   */
  inline int getLowerPadSize(int axis) const;

  /** \brief Get the size of the upper padding along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * Note that the returned padding can be either communication
   * padding or boundary padding as appropriate.
   */
  inline int getUpperPadSize(int axis) const;

  /** \brief Get the communication padding size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the communication padding along the
   * given axis at the time of construction, regardless of the
   * processor's position relative to the domain boundary.
   */
  inline int getCommPadSize(int axis) const;

  /** \brief Get the size of the lower boundary padding along the
   *         given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  inline int getLowerBndryPad(int axis) const;

  /** \brief Get the size of the upper boundary padding along the
   *         given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  inline int getUpperBndryPad(int axis) const;

  /** \brief Get the boundary padding size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This returns the value of the boundary padding along the given
   * axis at the time of construction, regardless of whether a sub-map
   * reduced these values.
   */
  inline int getBndryPadSize(int axis) const;

  /** \brief Assign all elements of the lower pad a constant value
   *
   * \param axis [in] the axis along which the lower padding is chosen
   *
   * \param value [in] the value to be assigned to all elements of the
   *        lower padding
   */
  void setLowerPad(int axis, const Scalar value);

  /** \brief Assign all elements of the upper pad a constant value
   *
   * \param axis [in] the axis along which the upper padding is chosen
   *
   * \param value [in] the value to be assigned to all elements of the
   *        upper padding
   */
  void setUpperPad(int axis, const Scalar value);

  /** \brief Return whether the given axis supports a replicated boundary
   */
  inline bool isReplicatedBoundary(int axis) const;

  /** \brief Get the storage order
   */
  inline Layout getLayout() const;

  /** \brief True if there are no stride gaps on any processor
   *
   * An MDVector constructed from a communicator and dimensions will
   * always be contiguous.  An MDVector that is a slice of a parent
   * MDMVector will generally be non-contiguous, with some exceptions.
   * There are cases where some local data is contiguous and some is
   * not, but this method returns True only if all processes' local
   * data is contiguous.
   */
  inline bool isContiguous() const;

  //@}

  /** \name Conversion to other types of Vectors */
  //@{

#ifdef HAVE_EPETRA

  /** \brief Return a view of this MDVector as an Epetra_IntVector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as an Epetra_IntVector.  Currently, if the
   * MDVector is non-contiguous, a Domi::MDMapNoncontiguous exception
   * will be thrown, as Epetra_IntVectors are contiguous, and this is
   * a view transform.  Non-contiguous data is the result of slicing a
   * parent MDVector.  In this case, getEpetraIntVectorView() should
   * be called on the parent, or getEpetraIntVectorCopy() should be
   * called on this non-contiguous MDVector.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_IntVector, i.e. int, or a Domi::TypeError exception will
   * be thrown.
   */
  Teuchos::RCP< Epetra_IntVector > getEpetraIntVectorView() const;

  /** \brief Return a view of this MDVector as an Epetra_Vector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as an Epetra_Vector.  Currently, if the
   * MDVector is non-contiguous, a Domi::MDMapNoncontiguous exception
   * will be thrown, as Epetra_Vectors are contiguous, and this is a
   * view transform.  Non-contiguous data is the result of slicing a
   * parent MDVector.  In this case, getEpetraVectorView() should be
   * called on the parent, or getEpetraVectorCopy() should be called
   * on this non-contiguous MDVector.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_Vector, i.e. double, or a Domi::TypeError exception will
   * be thrown.
   */
  Teuchos::RCP< Epetra_Vector > getEpetraVectorView() const;

  /** \brief Return a view of this MDVector as an Epetra_MultiVector
   *
   * This method will attempt to utilize either the leading dimension
   * (if the layout is C-order) or the trailing dimension (if the
   * layout is Fortran-order) to serve as the Epetra_MultiVector
   * vector index.  This dimension must be non-distributed (i.e., the
   * corresponding commDim must be 1) in order for this mapping to
   * work.  There must also be no padding associated with this axis.
   * If these requirements are not met, then the entire MDVector will
   * be treated as a single vector, that is the number of vectors will
   * equal 1.
   *
   * Currently, if the MDVector is non-contiguous, a
   * Domi::MDMapNoncontiguous exception will be thrown, as
   * Epetra_MultiVectors are contiguous, and this is a view transform.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_MultiVector, i.e. double, or a Domi::TypeError exception
   * will be thrown.
   */
  Teuchos::RCP< Epetra_MultiVector > getEpetraMultiVectorView() const;

  /** \brief Return a copy of this MDVector as an Epetra_IntVector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as an Epetra_IntVector.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_IntVector, i.e. int, or a Domi::TypeError exception will
   * be thrown.
   */
  Teuchos::RCP< Epetra_IntVector > getEpetraIntVectorCopy() const;

  /** \brief Return a copy of this MDVector as an Epetra_Vector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as an Epetra_Vector.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_Vector, i.e. double, or a Domi::TypeError exception will
   * be thrown.
   */
  Teuchos::RCP< Epetra_Vector > getEpetraVectorCopy() const;

  /** \brief Return a copy of this MDVector as an Epetra_MultiVector
   *
   * This method will attempt to utilize either the leading dimension
   * (if the layout is C-order) or the trailing dimension (if the
   * layout is Fortran-order) to serve as the Epetra_MultiVector
   * vector index.  This dimension must be non-distributed (i.e., the
   * corresponding commDim must be 1) in order for this mapping to
   * work.  There must also be no padding associated with this axis.
   * If these requirements are not met, then the entire MDVector will
   * be treated as a single vector, that is the number of vectors will
   * equal 1.
   *
   * The MDVector Scalar template must be the same type as an
   * Epetra_MultiVector, i.e. double, or a Domi::TypeError exception
   * will be thrown.
   */
  Teuchos::RCP< Epetra_MultiVector > getEpetraMultiVectorCopy() const;

#endif

#ifdef HAVE_TPETRA

  /** \brief Return a view of this MDVector as a Tpetra::Vector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as a Tpetra::Vector.  Currently, if the
   * MDVector is non-contiguous, a Domi::MDMapNoncontiguous exception
   * will be thrown, as Tpetra::Vectors are contiguous, and this is a
   * view transform.  In this case, getTpetraVectorView() should be
   * called on the parent, or getTpetraVectorCopy() should be called
   * on this non-contiguous MDVector.
   */
  template< class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node2 = Node>
  Teuchos::RCP< Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node2 > >
  getTpetraVectorView() const;

  /** \brief Return a view of this MDVector as a Tpetra::MultiVector
   *
   * This method will attempt to utilize either the leading dimension
   * (if the layout is C-order) or the trailing dimension (if the
   * layout is Fortran-order) to serve as the Tpetra::MultiVector
   * vector index.  This dimension must be non-distributed (i.e., the
   * corresponding commDim must be 1) in order for this mapping to
   * work.  There must also be no padding associated with this axis.
   * If these requirements are not met, then the entire MDVector will
   * be treated as a single vector, that is the number of vectors will
   * equal 1.
   *
   * Currently, if the MDVector is non-contiguous, a
   * Domi::MDMapNoncontiguous exception will be thrown, as
   * Tpetra::MultiVectors are contiguous, and this is a view
   * transform.
   */
  template < class LocalOrdinal,
             class GlobalOrdinal = LocalOrdinal,
             class Node2 = Node >
  Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                     LocalOrdinal,
                                     GlobalOrdinal,
                                     Node2 > >
  getTpetraMultiVectorView() const;

  /** \brief Return a copy of this MDVector as a Tpetra::Vector
   *
   * The multiple dimensions of the MDVector will be flattened in
   * order to be expressed as a Tpetra::Vector.
   */
  template< class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node2 = Node >
  Teuchos::RCP< Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node2 > >
  getTpetraVectorCopy() const;

  /** \brief Return a copy of this MDVector as a Tpetra::MultiVector
   *
   * This method will attempt to utilize either the leading dimension
   * (if the layout is C-order) or the trailing dimension (if the
   * layout is Fortran-order) to serve as the Tpetra::MultiVector
   * vector index.  This dimension must be non-distributed (i.e., the
   * corresponding commDim must be 1) in order for this mapping to
   * work.  There must also be no padding associated with this axis.
   * If these requirements are not met, then the entire MDVector will
   * be treated as a single vector, that is the number of vectors will
   * equal 1.
   */
  template < class LocalOrdinal,
             class GlobalOrdinal = LocalOrdinal,
             class Node2 = Node >
  Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                     LocalOrdinal,
                                     GlobalOrdinal,
                                     Node2 > >
  getTpetraMultiVectorCopy() const;

#endif

  //@}

  /** \name Data extraction methods */
  //@{

  /** \brief Get a non-const view of the data as an MDArrayView
   *
   * \param includePadding [in] if true, include the boundary and
   *        communication padding in the returned MDArrayView
   */
  MDArrayView< Scalar > getDataNonConst(bool includePadding = true);

  /** \brief Get a const view of the data as an MDArrayView
   *
   * \param includePadding [in] if true, include the boundary and
   *        communication padding in the returned MDArrayView
   */
  MDArrayView< const Scalar > getData(bool includePadding = true) const;

  /** \brief Get a non-const view of the lower padding data along the
   *         given axis as an MDArrayView
   *
   * \param axis [in] the axis from which to extract the lower padding
   */
  MDArrayView< Scalar > getLowerPadDataNonConst(int axis);

  /** \brief Get a const view of the lower padding data along the
   *         given axis as an MDArrayView
   *
   * \param axis [in] the axis from which to extract the lower padding
   */
  MDArrayView< const Scalar > getLowerPadData(int axis) const;

  /** \brief Get a non-const view of the upper padding data along the
   *         given axis as an MDArrayView
   *
   * \param axis [in] the axis from which to extract the upper padding
   */
  MDArrayView< Scalar > getUpperPadDataNonConst(int axis);

  /** \brief Get a const view of the upper padding data along the
   *         given axis as an MDArrayView
   *
   * \param axis [in] the axis from which to extract the upper padding
   */
  MDArrayView< const Scalar > getUpperPadData(int axis) const;

  //@}

  /** \name Mathematical methods */
  //@{

  /** \brief Compute the dot product of this MDVector and MDVector a
   *
   * \param a [in] partner MDVector for performing dot product
   */
  Scalar
  dot(const MDVector< Scalar, Node > & a) const;

  /** \brief Compute the 1-norm of this MDVector
   */
  typename Teuchos::ScalarTraits< Scalar >::magnitudeType norm1() const;

  /** \brief Compute the 2-norm of this MDVector
   */
  typename Teuchos::ScalarTraits< Scalar >::magnitudeType norm2() const;

  /** \brief Compute the infinity-norm of this MDVector
   */
  typename Teuchos::ScalarTraits< Scalar >::magnitudeType normInf() const;

  /** \brief Compute the weighted norm of this
   *
   * \param weights [in] MDVector of weights for weighted norm
   */
  typename Teuchos::ScalarTraits< Scalar >::magnitudeType
  normWeighted(const MDVector< Scalar, Node > & weights) const;

  /** \brief Compute the mean (average) value of this MDVector
   */
  Scalar meanValue() const;

  //@}

  /** \name Implementation of the Teuchos::Describable interface */
  //@{

  /** \brief A simple one-line description of this MDVector
   */
  virtual std::string description() const;

  /** \brief Print the object with some verbosity level to a
   *         FancyOStream
   *
   * \param out [in] output stream
   *
   * \param verbLevel [in] verbosity level
   */
  virtual void
  describe(Teuchos::FancyOStream &out,
           const Teuchos::EVerbosityLevel verbLevel =
             Teuchos::Describable::verbLevel_default) const;

  //@}

  /** \name Global assignment methods */
  //@{

  /** \brief Set all values in the multivector with the given value.
   *
   * \param value [in] assignment value
   *
   * \param includePadding [in] if true, assign values to the boundary
   *        and communication padding as well
   */
  void putScalar(const Scalar & value,
                 bool includePadding = true);

  /** \brief Set all values in the multivector to pseudorandom numbers.
   */
  void randomize();

  //@}

  /** \name Global communication methods */
  //@{

  // /** \brief Sum values of a locally replicated multivector across all
  //  *         processes.
  //  */
  // void reduce();

  /** \brief The simplest method for updating the communication padding.
   *
   * This method will update the communication padding along all axes.
  */
  void updateCommPad();

  /* \brief Update the data in the communication padding along the
   *        specified axis
   *
   * \param axis [in] the axis along which communication will be
   *        performed
   */
  void updateCommPad(int axis);

  /** \brief Start an asyncronous update of the communication padding
   *
   * \param axis [in] the axis along which communication will be
   *        performed
   *
   * Post the non-blocking sends and receives for the communication
   * padding along the given axis
  */
  void startUpdateCommPad(int axis);

  /** \brief Complete an asyncronous update of the communication padding
   *
   * \param axis [in] the axis along which communication will be
   *        performed
   *
   * Wait for all of the non-blocking updates for the communication
   * padding along the given axis to complete
   */
  void endUpdateCommPad(int axis);

  //@}

  /** \name Sub-MDVector operators */
  //@{

  /** \brief Sub-vector access operator.
   *
   * \param index [in] index of the desired sub-vector.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDVector</tt>.
   *
   * The returned <tt>MDVector</tt> will have one fewer dimensions
   * than the calling <tt>MDVector</tt>.
   */
  MDVector< Scalar, Node >
  operator[](dim_type index) const;

  /** \brief Sub-vector access operator.
   *
   * \param slice [in] a Slice of global indexes that specifies the
   *        desired sub-vector.  Note that to obtain expected
   *        behavior, you should always chain together <tt>n</tt>
   *        square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDVector</tt>.
   *
   * The returned <tt>MDVector</tt> will have the same number of
   * dimensions as the calling <tt>MDVector</tt>.
   *
   * Note that if you wish to obtain a sub-vector that has boundary
   * padding along the axis being sliced, you will have to use the
   * MDVector constructor that takes a parent MDVector, an axis, a
   * Slice, and a boundary padding specification.
   */
  MDVector< Scalar, Node >
  operator[](Slice slice) const;

  //@}

  /** \name Input/Output */
  //@{

  /** \brief Write the MDVector to a binary file
   *
   * \param filename [in] name of the output file
   *
   * \param includeBndryPad [in] if true, include the boundary pad
   *        with the output data
   */
  void writeBinary(const std::string & filename,
                   bool includeBndryPad = false) const;

  /** \brief Read the MDVector from a binary file
   *
   * \param filename [in] name of the input file
   *
   * \param includeBndryPad [in] if true, include the boundary pad
   *        with the input data
   */
  void readBinary(const std::string & filename,
                  bool includeBndryPad = false);

  //@}

private:

  // The Teuchos communicator.  Note that this is always a reference
  // to the communicator of the _mdMap, and is stored only for
  // convenience
  Teuchos::RCP< const Teuchos::Comm< int > > _teuchosComm;

  // The MDMap that describes the domain decomposition of this
  // MDVector
  Teuchos::RCP< const MDMap< Node > > _mdMap;

  // The MDArrayRCP that stores the data of this MDVector
  MDArrayRCP< Scalar > _mdArrayRcp;

  // The MDArrayView of the (possibly whole) sub-view into this
  // MDVector's MDArrayRCP
  MDArrayView< Scalar > _mdArrayView;

  // The operator[](int) and operator[](Slice) methods are applied to
  // a specific axis, namely this internally stored and updated next
  // axis
  int _nextAxis;

  ///////////////////////////////////
  // *** Communication Support *** //
  ///////////////////////////////////

#ifdef HAVE_MPI
  // An array of MPI_Request objects for supporting non-blocking sends
  // and receives
  Teuchos::Array< MPI_Request > _requests;
#endif

  // Define a struct for storing all the information needed for a
  // single message: a pointer to the buffer, a reference to the
  // message's MPI data type, the rank of the communication partner,
  // and the axis alng which we are communicating
  struct MessageInfo
  {
    // Pointer to first element of data buffer
    void *buffer;
#ifdef HAVE_MPI
    // MPI data type (strided vector)
    Teuchos::RCP< MPI_Datatype > datatype;
#else
    // Teuchos ArrayView for periodic domains
    MDArrayView< Scalar > dataview;
#endif
    // Processor rank for communication partner
    int proc;
    // Communication is along this axis
    int axis;
  };

  // An array of MessageInfo objects representing all active send
  // buffers.  The outer array represents the axis and the inner
  // 2-Tuple represents the lower and upper boundaries.
  Teuchos::Array< Teuchos::Tuple< MessageInfo, 2 > > _sendMessages;

  // An array of MessageInfo objects representing all active receive
  // buffers.  The outer array represents the axis and the inner
  // 2-Tuple represents the lower and upper boundaries.
  Teuchos::Array< Teuchos::Tuple< MessageInfo, 2 > > _recvMessages;

  // A private method to initialize the _sendMessages and
  // _recvMessages arrays.
  void initializeMessages();

  //////////////////////////////////
  // *** Input/Output Support *** //
  //////////////////////////////////

  // Define a struct for storing all of the information needed to
  // write or read the MDVector to a file: arrays that store the file
  // shape, the local buffer shape, the local data shape, the file
  // starting coordinates, and the data starting coordinates.
  struct FileInfo
  {
    Teuchos::Array< dim_type > fileShape;
    Teuchos::Array< dim_type > bufferShape;
    Teuchos::Array< dim_type > dataShape;
    Teuchos::Array< dim_type > fileStart;
    Teuchos::Array< dim_type > dataStart;
#ifdef HAVE_MPI
    Teuchos::RCP< MPI_Datatype > filetype;
    Teuchos::RCP< MPI_Datatype > datatype;
#endif
  };

  // FileInfo struct for an input or output file that does not store
  // boundary padding.  This is mutable because it does not get set
  // until the first time the MDVector is read or written to a file,
  // and the writeBinary() method should logically be const.
  mutable Teuchos::RCP< FileInfo > _fileInfo;

  // FileInfo struct for an input or output file that does store
  // boundary padding.  This is mutable because it does not get set
  // until the first time the MDVector is read or written to a file,
  // and the writeBinary() method should logically be const.
  mutable Teuchos::RCP< FileInfo > _fileInfoWithBndry;

  // Compute either the _fileInfo or _fileInfoWithBndry data members.
  // This private method gets called by the writeBinary() and
  // readBinary() methods, and sets the requested fileInfo object,
  // unless it has already been set.  This is const so that it can be
  // called by writeBinary(), but its whole purpose is to change
  // mutable data members.
  Teuchos::RCP< FileInfo > & computeFileInfo(bool includeBndryPad) const;

};

/////////////////////////////
// *** Implementations *** //
/////////////////////////////

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
         bool zeroOut) :
  _teuchosComm(mdMap->getTeuchosComm()),
  _mdMap(mdMap),
  _mdArrayRcp(),
  _mdArrayView(),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  // Obtain the array of dimensions
  int numDims = _mdMap->numDims();
  Teuchos::Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = _mdMap->getLocalDim(axis,true);

  // Resize the MDArrayRCP and set the MDArrayView
  _mdArrayRcp.resize(dims);
  _mdArrayView = _mdArrayRcp();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
         const dim_type leadingDim,
         const dim_type trailingDim,
         bool zeroOut) :
  _teuchosComm(mdMap->getTeuchosComm()),
  _mdMap(mdMap->getAugmentedMDMap(leadingDim, trailingDim)),
  _mdArrayRcp(),
  _mdArrayView(),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  // Obtain the array of dimensions
  int numDims = _mdMap->numDims();
  Teuchos::Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = _mdMap->getLocalDim(axis,true);

  // Resize the MDArrayRCP and set the MDArrayView
  _mdArrayRcp.resize(dims);
  _mdArrayView = _mdArrayRcp();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
         const MDArrayView< Scalar > & source) :
  _mdMap(mdMap),
  _mdArrayRcp(source),
  _mdArrayView(_mdArrayRcp()),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");
  int numDims = _mdMap->numDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != _mdArrayRcp.numDims(),
    InvalidArgument,
    "MDMap and source array do not have the same number of dimensions");

  for (int axis = 0; axis < numDims; ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      _mdMap->getLocalDim(axis) != _mdArrayRcp.dimension(axis),
      InvalidArgument,
      "Axis " << axis << ": MDMap dimension = " << _mdMap->getLocalDim(axis)
      << ", MDArray dimension = " << _mdArrayRcp.dimension(axis));
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const MDMap< Node > > & mdMap,
         const MDArrayRCP< Scalar > & mdArrayRcp) :
  _mdMap(mdMap),
  _mdArrayRcp(mdArrayRcp),
  _mdArrayView(_mdArrayRcp()),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
#ifdef DOMI_MDVECTOR_VERBOSE
  std::cout << "_mdArrayRcp  = " << _mdArrayRcp  << std::endl;
  std::cout << "_mdArrayView.getRawPtr() = " << _mdArrayView.getRawPtr()
            << " (constructor)" << std::endl;
  std::cout << "_mdArrayView = " << _mdArrayView << std::endl;
#endif
  setObjectLabel("Domi::MDVector");
  int numDims = _mdMap->numDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != _mdArrayRcp.numDims(),
    InvalidArgument,
    "MDMap and source array do not have the same number of dimensions");

  for (int axis = 0; axis < numDims; ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      _mdMap->getLocalDim(axis) != _mdArrayRcp.dimension(axis),
      InvalidArgument,
      "Axis " << axis << ": MDMap dimension = " << _mdMap->getLocalDim(axis)
      << ", MDArray dimension = " << _mdArrayRcp.dimension(axis));
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const MDVector< Scalar, Node > & source,
         Teuchos::DataAccess access) :
  _teuchosComm(source.getMDMap()->getTeuchosComm()),
  _mdMap(source.getMDMap()),
  _mdArrayRcp(source._mdArrayRcp),
  _mdArrayView(source._mdArrayView),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  if (access == Teuchos::Copy)
  {
#ifdef DOMI_MDVECTOR_VERBOSE
    std::cout << "Inside MDVector copy constructor with copy access" << std::endl;
#endif
    // Obtain the array of dimensions
    int numDims = _mdMap->numDims();
    Teuchos::Array< dim_type > dims(numDims);
    for (int axis = 0; axis < numDims; ++axis)
      dims[axis] = _mdMap->getLocalDim(axis,true);

    // Reset the MDArrayRCP and set the MDArrayView
    _mdArrayRcp  = MDArrayRCP< Scalar >(dims, 0, source.getLayout());
    _mdArrayView = _mdArrayRcp();

    // Copy the source data to the new MDVector
    typedef typename MDArrayView< Scalar >::iterator iterator;
    typedef typename MDArrayView< Scalar >::const_iterator const_iterator;
    const_iterator src = source.getData().begin();
    for (iterator trg = _mdArrayView.begin(); trg != _mdArrayView.end(); ++trg)
    {
      *trg = *src;
      ++src;
    }
  }
#ifdef DOMI_MDVECTOR_VERBOSE
  else
  {
    std::cout << "Inside MDVector copy constructor with view access"
              << std::endl;
  }
#endif
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
         Teuchos::ParameterList & plist) :
  _teuchosComm(teuchosComm),
  _mdMap(),
  _mdArrayRcp(),
  _mdArrayView(),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  // Compute the MDComm and MDMap
  MDMap< Node > * myMdMap = new MDMap< Node >(teuchosComm, plist);
  dim_type leadingDim  = plist.get("leading dimension" , 0);
  dim_type trailingDim = plist.get("trailing dimension", 0);
  if (leadingDim + trailingDim > 0)
  {
    _mdMap = myMdMap->getAugmentedMDMap(leadingDim, trailingDim);
    delete myMdMap;
  }
  else
    _mdMap = Teuchos::rcp(myMdMap);

  // Obtain the array of dimensions
  int numDims = _mdMap->numDims();
  Teuchos::Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = _mdMap->getLocalDim(axis,true);

  // Resize the MDArrayRCP and set the MDArrayView
  _mdArrayRcp.resize(dims);
  _mdArrayView = _mdArrayRcp();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const Teuchos::RCP< const MDComm > mdComm,
         Teuchos::ParameterList & plist) :
  _teuchosComm(mdComm->getTeuchosComm()),
  _mdMap(Teuchos::rcp(new MDMap< Node >(mdComm, plist))),
  _mdArrayRcp(),
  _mdArrayView(),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  // Compute the MDMap
  MDMap< Node > * myMdMap = new MDMap< Node >(mdComm, plist);
  dim_type leadingDim  = plist.get("leading dimension" , 0);
  dim_type trailingDim = plist.get("trailing dimension", 0);
  if (leadingDim + trailingDim > 0)
  {
    _mdMap = myMdMap->getAugmentedMDMap(leadingDim, trailingDim);
    delete myMdMap;
  }
  else
    _mdMap = Teuchos::rcp(myMdMap);

  // Obtain the array of dimensions
  int numDims = _mdMap->numDims();
  Teuchos::Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = _mdMap->getLocalDim(axis,true);

  // Resize the MDArrayRCP and set the MDArrayView
  _mdArrayRcp.resize(dims);
  _mdArrayView = _mdArrayRcp();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const MDVector< Scalar, Node > & parent,
         int axis,
         dim_type globalIndex) :
  _teuchosComm(parent._teuchosComm),
  _mdMap(),
  _mdArrayRcp(parent._mdArrayRcp),
  _mdArrayView(parent._mdArrayView),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
  setObjectLabel("Domi::MDVector");

  // Obtain the parent MDMap
  Teuchos::RCP< const MDMap< Node > > parentMdMap = parent.getMDMap();

  // Obtain the new, sliced MDMap
  _mdMap = Teuchos::rcp(new MDMap< Node >(*parentMdMap,
                                          axis,
                                          globalIndex));

  // Check that we are on the new sub-communicator
  if (_mdMap->onSubcommunicator())
  {
    // Convert the index from global to local.  We start by
    // determining the starting global index on this processor along
    // the given axis, ignoring the boundary padding.  We then
    // subtract the lower padding, whether it is communication padding
    // or boundary padding.
    dim_type origin = parentMdMap->getGlobalRankBounds(axis,false).start() -
                      parentMdMap->getLowerPadSize(axis);

    // The local index along the given axis is the global axis minus
    // the starting index.  Since we are on the subcommunicator, this
    // should be valid.
    dim_type localIndex = globalIndex - origin;

    // Obtain the new MDArrayView using the local index
    MDArrayView< Scalar > newView(_mdArrayView, axis, localIndex);
    _mdArrayView = newView;
  }
  else
  {
    // We are not on the sub-communicator, so clear out the data
    // buffer and view
    _mdArrayRcp.clear();
    _mdArrayView = MDArrayView< Scalar >();
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const MDVector< Scalar, Node > & parent,
         int axis,
         const Slice & slice,
         int bndryPad) :
  _teuchosComm(),
  _mdMap(),
  _mdArrayRcp(parent._mdArrayRcp),
  _mdArrayView(parent._mdArrayView),
  _nextAxis(0),
#ifdef HAVE_MPI
  _requests(),
#endif
  _sendMessages(),
  _recvMessages()
{
#ifdef DOMI_MDVECTOR_VERBOSE
  std::cout << "slice axis " << axis << std::endl;
  std::cout << "  _mdArrayRcp  @ " << _mdArrayRcp.getRawPtr()  << std::endl;
  std::cout << "  _mdArrayView @ " << _mdArrayView.getRawPtr() << std::endl;
#endif

  setObjectLabel("Domi::MDVector");

  // Obtain the parent MDMap
  Teuchos::RCP< const MDMap< Node > > parentMdMap = parent.getMDMap();

  // Obtain the new, sliced MDMap
  _mdMap = Teuchos::rcp(new MDMap< Node >(*parentMdMap,
                                          axis,
                                          slice,
                                          bndryPad));
  _teuchosComm = _mdMap->getTeuchosComm();

  // Check that we are on the new sub-communicator
  if (_mdMap->onSubcommunicator())
  {
    // Get the concrete bounds
    Slice bounds = slice.bounds(parentMdMap->getGlobalDim(axis,true));

    // Convert the given Slice start index from global to local.  We
    // start by determining the starting global index on this
    // processor along the given axis, ignoring the boundary padding.
    // We then subtract the lower padding, whether it is communication
    // padding or boundary padding.
    dim_type origin = parentMdMap->getGlobalRankBounds(axis,false).start() -
                      parentMdMap->getLowerPadSize(axis);

    // Determine the starting index of our local slice.  This will be
    // the start of the given slice minus the starting global index on
    // this processor minus the given boundary pad.  If this is less
    // than zero, then the start is on a lower processor, so set the
    // local start to zero.
    dim_type start = std::max(0, bounds.start() - origin - bndryPad);

    // Now get the stop index of the local slice.  This will be the
    // stop of the given slice minus the starting global index on this
    // processor plus the given boundary pad.  If this is larger than
    // the local dimension, then set the local stop to the local
    // dimension.
    dim_type stop = std::min(bounds.stop() - origin + bndryPad,
                             parentMdMap->getLocalDim(axis,true));

    // Obtain the new MDArrayView using the local slice
    MDArrayView< Scalar > newView(_mdArrayView, axis, Slice(start,stop));
    _mdArrayView = newView;
  }
  else
  {
    // We are not on the sub-communicator, so clear out the data
    // buffer and view
    _mdArrayRcp.clear();
    _mdArrayView = MDArrayView< Scalar >();
  }
#ifdef DOMI_MDVECTOR_VERBOSE
  std::cout << "  _mdArrayView @ " << _mdArrayView.getRawPtr() << std::endl;
  std::cout << "  offset = " << int(_mdArrayView.getRawPtr() -
                                    _mdArrayRcp.getRawPtr()) << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::
MDVector(const MDVector< Scalar, Node > & parent,
         const Teuchos::ArrayView< Slice > & slices,
         const Teuchos::ArrayView< int > & bndryPad)
{
  setObjectLabel("Domi::MDVector");

  // Temporarily store the number of dimensions
  int numDims = parent.numDims();

  // Sanity check on dimensions
  TEUCHOS_TEST_FOR_EXCEPTION(
    (slices.size() != numDims),
    InvalidArgument,
    "number of slices = " << slices.size() << " != parent MDVector number of "
    "dimensions = " << numDims);

  // Apply the single-Slice constructor to each axis in succession
  MDVector< Scalar, Node > tempMDVector1(parent);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int bndryPadding = (axis < bndryPad.size()) ? bndryPad[axis] : 0;
    MDVector< Scalar, Node > tempMDVector2(tempMDVector1,
                                           axis,
                                           slices[axis],
                                           bndryPadding);
    tempMDVector1 = tempMDVector2;
  }
  *this = tempMDVector1;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node > &
MDVector< Scalar, Node >::
operator=(const MDVector< Scalar, Node > & source)
{
  _teuchosComm  = source._teuchosComm;
  _mdMap        = source._mdMap;
  _mdArrayRcp   = source._mdArrayRcp;
  _mdArrayView  = source._mdArrayView;
  _nextAxis     = source._nextAxis;
#ifdef HAVE_MPI
  _requests     = source._requests;
#endif
  _sendMessages = source._sendMessages;
  _recvMessages = source._recvMessages;
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >::~MDVector()
{
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
const Teuchos::RCP< const MDMap< Node > >
MDVector< Scalar, Node >::
getMDMap() const
{
  return _mdMap;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
bool
MDVector< Scalar, Node >::
onSubcommunicator() const
{
  return _mdMap->onSubcommunicator();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::RCP< const Teuchos::Comm< int > >
MDVector< Scalar, Node >::
getTeuchosComm() const
{
  return _mdMap->getTeuchosComm();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
numDims() const
{
  return _mdMap->numDims();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getCommDim(int axis) const
{
  return _mdMap->getCommDim(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
bool
MDVector< Scalar, Node >::
isPeriodic(int axis) const
{
  return _mdMap->isPeriodic(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getCommIndex(int axis) const
{
  return _mdMap->getCommIndex(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getLowerNeighbor(int axis) const
{
  return _mdMap->getLowerNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getUpperNeighbor(int axis) const
{
  return _mdMap->getUpperNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
dim_type
MDVector< Scalar, Node >::
getGlobalDim(int axis, bool withBndryPad) const
{
  return _mdMap->getGlobalDim(axis, withBndryPad);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Slice
MDVector< Scalar, Node >::
getGlobalBounds(int axis, bool withBndryPad) const
{
  return _mdMap->getGlobalBounds(axis, withBndryPad);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Slice
MDVector< Scalar, Node >::
getGlobalRankBounds(int axis, bool withBndryPad) const
{
  return _mdMap->getGlobalRankBounds(axis, withBndryPad);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
dim_type
MDVector< Scalar, Node >::
getLocalDim(int axis, bool withCommPad) const
{
  return _mdMap->getLocalDim(axis, withCommPad);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::ArrayView< const Slice >
MDVector< Scalar, Node >::
getLocalBounds() const
{
  return _mdMap->getLocalBounds();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Slice
MDVector< Scalar, Node >::
getLocalBounds(int axis, bool withCommPad) const
{
  return _mdMap->getLocalBounds(axis, withCommPad);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Slice
MDVector< Scalar, Node >::
getLocalInteriorBounds(int axis) const
{
  return _mdMap->getLocalInteriorBounds(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
bool
MDVector< Scalar, Node >::
hasPadding() const
{
  return _mdMap->hasPadding();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getLowerPadSize(int axis) const
{
  return _mdMap->getLowerPadSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getUpperPadSize(int axis) const
{
  return _mdMap->getUpperPadSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getCommPadSize(int axis) const
{
  return _mdMap->getCommPadSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getLowerBndryPad(int axis) const
{
  return _mdMap->getLowerBndryPad(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getUpperBndryPad(int axis) const
{
  return _mdMap->getUpperBndryPad(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
int
MDVector< Scalar, Node >::
getBndryPadSize(int axis) const
{
  return _mdMap->getBndryPadSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
setLowerPad(int axis,
            const Scalar value)
{
  MDArrayView< Scalar > lowerPad = getLowerPadDataNonConst(axis);
  lowerPad.assign(value);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
setUpperPad(int axis,
            const Scalar value)
{
  MDArrayView< Scalar > upperPad = getUpperPadDataNonConst(axis);
  upperPad.assign(value);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
bool
MDVector< Scalar, Node >::
isReplicatedBoundary(int axis) const
{
  return _mdMap->isReplicatedBoundary(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Layout
MDVector< Scalar, Node >::
getLayout() const
{
  return _mdMap->getLayout();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
bool
MDVector< Scalar, Node >::
isContiguous() const
{
  return _mdMap->isContiguous();
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA

// The getEpetraIntVectorView() method only makes sense for Scalar =
// int, because Epetra_IntVectors store data buffers of type int only.
// There is no convenient way to specialize just one (or a small
// handfull of) methods, instead you have to specialize the whole
// class.  So we allow getEpetraIntVectorView() to compile for any
// Scalar, but we will throw an exception if Scalar is not int.

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_IntVector >
MDVector< Scalar, Node >::
getEpetraIntVectorView() const
{
  // Throw an exception if Scalar is not int
  TEUCHOS_TEST_FOR_EXCEPTION(
    typeid(Scalar) != typeid(int),
    TypeError,
    "MDVector is of scalar type '" << typeid(Scalar).name() << "', but "
    "Epetra_IntVector requires scalar type 'int'");

  // Throw an exception if this MDVector's MDMap is not contiguous
  TEUCHOS_TEST_FOR_EXCEPTION(
    !isContiguous(),
    MDMapNoncontiguousError,
    "This MDVector's MDMap is non-contiguous.  This can happen when you take "
    "a slice of a parent MDVector.");

  // Get the Epetra_Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Epetra_Map > epetraMap = _mdMap->getEpetraMap(true);

  // Return the result.  We are changing a Scalar* to a double*, which
  // looks sketchy, but we have thrown an exception if Sca is not
  // double, so everything is kosher.
  void * buffer = (void*) _mdArrayView.getRawPtr();
  return Teuchos::rcp(new Epetra_IntVector(View,
                                           *epetraMap,
                                           (int*) buffer));
}

////////////////////////////////////////////////////////////////////////

// The getEpetraVectorView() method only makes sense for Scalar =
// double, because Epetra_Vectors store data buffers of type double
// only.  There is no convenient way to specialize just one (or a
// small handfull of) methods, instead you have to specialize the
// whole class.  So we allow getEpetraVectorView() to compile for any
// Scalar, but we will throw an exception if Scalar is not double.

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_Vector >
MDVector< Scalar, Node >::
getEpetraVectorView() const
{
  // Throw an exception if Scalar is not double
  TEUCHOS_TEST_FOR_EXCEPTION(
    typeid(Scalar) != typeid(double),
    TypeError,
    "MDVector is of scalar type '" << typeid(Scalar).name() << "', but "
    "Epetra_Vector requires scalar type 'double'");

  // Throw an exception if this MDVector's MDMap is not contiguous
  TEUCHOS_TEST_FOR_EXCEPTION(
    !isContiguous(),
    MDMapNoncontiguousError,
    "This MDVector's MDMap is non-contiguous.  This can happen when you take "
    "a slice of a parent MDVector.");

  // Get the Epetra_Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Epetra_Map > epetraMap = _mdMap->getEpetraMap(true);

  // Return the result.  We are changing a Scalar* to a double*, which
  // looks sketchy, but we have thrown an exception if Sca is not
  // double, so everything is kosher.
  void * buffer = (void*) _mdArrayView.getRawPtr();
  return Teuchos::rcp(new Epetra_Vector(View,
                                        *epetraMap,
                                        (double*) buffer));
}

////////////////////////////////////////////////////////////////////////

// The getEpetraMultiVectorView() method only makes sense for Scalar =
// double, because Epetra_MultiVectors store data buffers of type
// double only.  There is no convenient way to specialize just one (or
// a small handfull of) methods, instead you have to specialize the
// whole class.  So we allow getEpetraVectorView() to compile for any
// Scalar, but we will throw an exception if Scalar is not double.

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_MultiVector >
MDVector< Scalar, Node >::
getEpetraMultiVectorView() const
{
  // Throw an exception if Scalar is not double
  TEUCHOS_TEST_FOR_EXCEPTION(
    typeid(Scalar) != typeid(double),
    TypeError,
    "MDVector is of scalar type '" << typeid(Scalar).name() << "', but "
    "Epetra_Vector requires scalar type 'double'");

  // Determine the vector axis and related info
  int vectorAxis = (getLayout() == C_ORDER) ? 0 : numDims()-1;
  int padding    = getLowerPadSize(vectorAxis) + getUpperPadSize(vectorAxis);
  int commDim    = getCommDim(vectorAxis);
  int numVectors = getGlobalDim(vectorAxis);

  // Obtain the appropriate MDMap and check that it is contiguous
  Teuchos::RCP< const MDMap< Node > > newMdMap;
  if (padding == 0 && commDim == 1)
    newMdMap = Teuchos::rcp(new MDMap< Node >(*_mdMap, vectorAxis, 0));
  else
  {
    newMdMap = _mdMap;
    numVectors = 1;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! newMdMap->isContiguous(),
    MDMapNoncontiguousError,
    "This MDVector's MDMap is non-contiguous.  This can happen when you take "
    "a slice of a parent MDVector.");

  // Get the stride between vectors.  The MDMap strides are private,
  // but we know the new MDMap is contiguous, so we can calculate it
  // as the product of the new MDMap dimensions (including padding)
  size_type stride = newMdMap->getLocalDim(0,true);
  for (int axis = 0; axis < newMdMap->numDims(); ++axis)
    stride *= newMdMap->getLocalDim(axis,true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    stride*numVectors > Teuchos::OrdinalTraits<int>::max(),
    MapOrdinalError,
    "Buffer size " << stride*numVectors << " is too large for Epetra int "
    "ordinals");
  int lda = (int)stride;

  // Get the Epetra_Map that is equivalent to newMdMap
  Teuchos::RCP< const Epetra_Map > epetraMap = newMdMap->getEpetraMap(true);

  // Return the result.  We are changing a Scalar* to a double*, which
  // looks sketchy, but we have thrown an exception if Sca is not
  // double, so everything is kosher.
  void * buffer = (void*) _mdArrayView.getRawPtr();
  return Teuchos::rcp(new Epetra_MultiVector(View,
                                             *epetraMap,
                                             (double*) buffer,
                                             lda,
                                             numVectors));
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_IntVector >
MDVector< Scalar, Node >::
getEpetraIntVectorCopy() const
{
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  // Get the Epetra_Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Epetra_Map > epetraMap = _mdMap->getEpetraMap(true);

  // Construct the result
  Teuchos::RCP< Epetra_IntVector > result =
    Teuchos::rcp(new Epetra_IntVector(*epetraMap));

  // Copy the MDVector data buffer to the Epetra_IntVector, even if the
  // MDVector is non-contiguous
  int ii = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    result->operator[](ii++) = (int) *it;

  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_Vector >
MDVector< Scalar, Node >::
getEpetraVectorCopy() const
{
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  // Get the Epetra_Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Epetra_Map > epetraMap = _mdMap->getEpetraMap(true);

  // Construct the result
  Teuchos::RCP< Epetra_Vector > result =
    Teuchos::rcp(new Epetra_Vector(*epetraMap));

  // Copy the MDVector data buffer to the Epetra_Vector, even if the
  // MDVector is non-contiguous
  int ii = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    result->operator[](ii++) = (double) *it;

  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::RCP< Epetra_MultiVector >
MDVector< Scalar, Node >::
getEpetraMultiVectorCopy() const
{
  typedef typename MDArrayView< Scalar >::iterator iterator;
  typedef typename MDArrayView< const Scalar >::iterator citerator;

  // Determine the vector axis and related info
  int vectorAxis = (getLayout() == C_ORDER) ? 0 : numDims()-1;
  int padding    = getLowerPadSize(vectorAxis) + getUpperPadSize(vectorAxis);
  int commDim    = getCommDim(vectorAxis);
  int numVectors = getGlobalDim(vectorAxis);

  // Obtain the appropriate MDMap
  Teuchos::RCP< const MDMap< Node > > newMdMap;
  if (padding == 0 && commDim == 1)
    newMdMap = Teuchos::rcp(new MDMap< Node >(*_mdMap, vectorAxis, 0));
  else
  {
    newMdMap = _mdMap;
    numVectors = 1;
  }

  // Get the Epetra_Map that is equivalent to newMdMap
  Teuchos::RCP< const Epetra_Map > epetraMap = newMdMap->getEpetraMap(true);

  // Construct the result
  Teuchos::RCP< Epetra_MultiVector > result =
    Teuchos::rcp(new Epetra_MultiVector(*epetraMap, numVectors));

  // Copy the MDVector data to the Epetra_MultiVector, even if the
  // MDVector is non-contiguous
  int ii = 0;
  if (numVectors == 1)
  {
    for (citerator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
      result->operator[](0)[ii++] = (double) *it;
  }
  else
  {
    for (int iv = 0; iv < numVectors; ++iv)
    {
      ii = 0;
      MDArrayView< Scalar > subArray(_mdArrayView, vectorAxis, iv);
      for (iterator it = subArray.begin(); it != subArray.end(); ++it)
        result->operator[](iv)[ii++] = (double) *it;
    }
  }

  // Return the result
  return result;
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA

template< class Scalar, class Node >
template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node2 >
Teuchos::RCP< Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node2 > >
MDVector< Scalar, Node >::
getTpetraVectorView() const
{
  // Throw an exception if this MDVector's MDMap is not contiguous
  TEUCHOS_TEST_FOR_EXCEPTION(
    !isContiguous(),
    MDMapNoncontiguousError,
    "This MDVector's MDMap is non-contiguous.  This can happen when you take "
    "a slice of a parent MDVector.");

  // Get the Tpetra::Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2 > > tpetraMap =
    _mdMap->template getTpetraMap< LocalOrdinal, GlobalOrdinal, Node2 >(true);

  // Return the result
  return Teuchos::rcp(new Tpetra::Vector< Scalar,
                                          LocalOrdinal,
                                          GlobalOrdinal,
                                          Node2 >(tpetraMap,
                                                  _mdArrayView.arrayView()));
}

////////////////////////////////////////////////////////////////////////

template< class Scalar, class Node >
template< class LocalOrdinal,
          class GlobalOrdinal,
          class Node2 >
Teuchos::RCP< Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node2 > >
MDVector< Scalar, Node >::
getTpetraVectorCopy() const
{
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  // Get the Tpetra::Map that is equivalent to this MDVector's MDMap
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2 > > tpetraMap =
    _mdMap->template getTpetraMap< LocalOrdinal, GlobalOrdinal, Node2 >(true);

  // Construct the result
  Teuchos::RCP< Tpetra::Vector< Scalar,
                                LocalOrdinal,
                                GlobalOrdinal,
                                Node2 > > result =
    Teuchos::rcp(new Tpetra::Vector< Scalar,
                                     LocalOrdinal,
                                     GlobalOrdinal,
                                     Node2 >(tpetraMap) );

  // Copy the MDVector data to the Tpetra::Vector, even if the
  // MDVector is non-contiguous
  Teuchos::ArrayRCP< Scalar > tpetraVectorArray = result->getDataNonConst();
  int ii = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    tpetraVectorArray[ii++] = *it;

  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar, class Node >
template < class LocalOrdinal,
           class GlobalOrdinal,
           class Node2 >
Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                   LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2 > >
MDVector< Scalar, Node >::
getTpetraMultiVectorView() const
{
  // Determine the vector axis and related info
  int vectorAxis    = (getLayout() == C_ORDER) ? 0 : numDims()-1;
  int padding       = getLowerPadSize(vectorAxis) + getUpperPadSize(vectorAxis);
  int commDim       = getCommDim(vectorAxis);
  size_t numVectors = getGlobalDim(vectorAxis);

  // Obtain the appropriate MDMap and check that it is contiguous
  Teuchos::RCP< const MDMap< Node > > newMdMap;
  if (padding == 0 && commDim == 1)
    newMdMap = Teuchos::rcp(new MDMap< Node >(*_mdMap, vectorAxis, 0));
  else
  {
    newMdMap = _mdMap;
    numVectors = 1;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! newMdMap->isContiguous(),
    MDMapNoncontiguousError,
    "This MDVector's MDMap is non-contiguous.  This can happen when you take "
    "a slice of a parent MDVector.");

  // Get the stride between vectors.  The MDMap strides are private,
  // but we know the new MDMap is contiguous, so we can calculate it
  // as the product of the new MDMap dimensions (including padding)
  size_type stride = newMdMap->getLocalDim(0,true);
  for (int axis = 0; axis < newMdMap->numDims(); ++axis)
    stride *= newMdMap->getLocalDim(axis,true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    stride*numVectors > Teuchos::OrdinalTraits<GlobalOrdinal>::max(),
    MapOrdinalError,
    "Buffer size " << stride*numVectors << " is too large for Tpetra "
    "GlobalOrdinal = " << typeid(GlobalOrdinal).name() );
  size_t lda = (size_t)stride;

  // Get the Tpetra::Map that is equivalent to newMdMap
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2> > tpetraMap =
    newMdMap->template getTpetraMap< LocalOrdinal, GlobalOrdinal, Node2 >(true);

  // Return the result
  return Teuchos::rcp(new Tpetra::MultiVector< Scalar,
                                               LocalOrdinal,
                                               GlobalOrdinal,
                                               Node2 >(tpetraMap,
                                                       _mdArrayView.arrayView(),
                                                       lda,
                                                       numVectors));
}

////////////////////////////////////////////////////////////////////////

template< class Scalar, class Node >
template < class LocalOrdinal,
           class GlobalOrdinal,
           class Node2 >
Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                   LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2 > >
MDVector< Scalar, Node >::
getTpetraMultiVectorCopy() const
{
  typedef typename MDArrayView< Scalar >::iterator iterator;
  typedef typename MDArrayView< const Scalar >::iterator citerator;

  // Determine the vector axis and related info
  int vectorAxis = (getLayout() == C_ORDER) ? 0 : numDims()-1;
  int padding    = getLowerPadSize(vectorAxis) + getUpperPadSize(vectorAxis);
  int commDim    = getCommDim(vectorAxis);
  int numVectors = getGlobalDim(vectorAxis);

  // Obtain the appropriate MDMap
  Teuchos::RCP< const MDMap< Node > > newMdMap;
  if (padding == 0 && commDim == 1)
    newMdMap = Teuchos::rcp(new MDMap< Node >(*_mdMap, vectorAxis, 0));
  else
  {
    newMdMap = _mdMap;
    numVectors = 1;
  }

  // Get the Tpetra::Map that is equivalent to newMdMap
  Teuchos::RCP< const Tpetra::Map< LocalOrdinal,
                                   GlobalOrdinal,
                                   Node2 > > tpetraMap =
    newMdMap->template getTpetraMap< LocalOrdinal, GlobalOrdinal, Node2 >(true);

  // Construct the result
  Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                     LocalOrdinal,
                                     GlobalOrdinal,
                                     Node2 > > result =
    Teuchos::rcp(new Tpetra::MultiVector< Scalar,
                                          LocalOrdinal,
                                          GlobalOrdinal,
                                          Node2 >(tpetraMap, numVectors));

  // Copy the MDVector to the Tpetra::MultiVector, even if the
  // MDVector is non-contiguous
  int ii = 0;
  if (numVectors == 1)
  {
    Teuchos::ArrayRCP< Scalar > tpetraVectorArray = result->getDataNonConst(0);
    for (citerator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
      tpetraVectorArray[ii++] = (double) *it;
  }
  else
  {
    for (int iv = 0; iv < numVectors; ++iv)
    {
      ii = 0;
      Teuchos::ArrayRCP< Scalar > tpetraVectorArray =
        result->getDataNonConst(iv);
      MDArrayView< Scalar > subArray(_mdArrayView, vectorAxis, iv);
      for (iterator it = subArray.begin(); it != subArray.end(); ++it)
        tpetraVectorArray[ii++] = (double) *it;
    }
  }

  // Return the result
  return result;
}

#endif

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< Scalar >
MDVector< Scalar, Node >::
getDataNonConst(bool includePadding)
{
#ifdef DOMI_MDVECTOR_VERBOSE
  std::cout << "_mdArrayView.getRawPtr() = " << _mdArrayView.getRawPtr()
            << std::endl;
  std::cout << "_mdArrayView = " << _mdArrayView << std::endl;
#endif
  if (includePadding)
    return _mdArrayView;
  MDArrayView< Scalar > newArray(_mdArrayView);
  for (int axis = 0; axis < numDims(); ++axis)
  {
    int lo = getLowerPadSize(axis);
    int hi = getLocalDim(axis,true) - getUpperPadSize(axis);
    newArray = MDArrayView< Scalar >(newArray, axis, Slice(lo,hi));
  }
  return newArray;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< const Scalar >
MDVector< Scalar, Node >::
getData(bool includePadding) const
{
  if (includePadding)
    return _mdArrayView.getConst();
  MDArrayView< Scalar > newArray(_mdArrayView);
  for (int axis = 0; axis < numDims(); ++axis)
  {
    int lo = getLowerPadSize(axis);
    int hi = getLocalDim(axis,true) - getUpperPadSize(axis);
    newArray = MDArrayView< Scalar >(newArray, axis, Slice(lo,hi));
  }
  return newArray.getConst();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< Scalar >
MDVector< Scalar, Node >::
getLowerPadDataNonConst(int axis)
{
  MDArrayView< Scalar > newArrayView(_mdArrayView,
                                     axis,
                                     Slice(getLowerPadSize(axis)));
  return newArrayView;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< const Scalar >
MDVector< Scalar, Node >::
getLowerPadData(int axis) const
{
  MDArrayView< const Scalar > newArrayView(_mdArrayView.getConst(),
                                           axis,
                                           Slice(getLowerPadSize(axis)));
  return newArrayView;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< Scalar >
MDVector< Scalar, Node >::
getUpperPadDataNonConst(int axis)
{
  dim_type n  = getLocalDim(axis,true);
  int     pad = getUpperPadSize(axis);
  Slice   slice;
  if (pad) slice = Slice(n-pad,n);
  else     slice = Slice(n-1,n-1);
  MDArrayView< Scalar > newArrayView(_mdArrayView,
                                     axis,
                                     slice);
  return newArrayView;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDArrayView< const Scalar >
MDVector< Scalar, Node >::
getUpperPadData(int axis) const
{
  dim_type n  = getLocalDim(axis,true);
  int     pad = getUpperPadSize(axis);
  Slice   slice;
  if (pad) slice = Slice(n-pad,n);
  else     slice = Slice(n-1,n-1);
  MDArrayView< const Scalar > newArrayView(_mdArrayView.getConst(),
                                           axis,
                                           slice);
  return newArrayView;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Scalar
MDVector< Scalar, Node >::
dot(const MDVector< Scalar, Node > & a) const
{
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! _mdMap->isCompatible(*(a._mdMap)),
    MDMapError,
    "MDMap of calling MDVector and argument 'a' are incompatible");

  MDArrayView< const Scalar > aView = a.getData();
  Scalar local_dot = 0;
  iterator a_it = aView.begin();
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end();
       ++it, ++a_it)
    local_dot += *it * *a_it;
  Scalar global_dot = 0;
  Teuchos::reduceAll(*_teuchosComm,
                     Teuchos::REDUCE_SUM,
                     1,
                     &local_dot,
                     &global_dot);
  return global_dot;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
typename Teuchos::ScalarTraits< Scalar >::magnitudeType
MDVector< Scalar, Node >::
norm1() const
{
  typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType mag;
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  mag local_norm1 = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    local_norm1 += std::abs(*it);
  mag global_norm1 = 0;
  Teuchos::reduceAll(*_teuchosComm,
                     Teuchos::REDUCE_SUM,
                     1,
                     &local_norm1,
                     &global_norm1);
  return global_norm1;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
typename Teuchos::ScalarTraits< Scalar >::magnitudeType
MDVector< Scalar, Node >::
norm2() const
{
  typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType mag;
  mag norm2 = dot(*this);
  return Teuchos::ScalarTraits<mag>::squareroot(norm2);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
typename Teuchos::ScalarTraits< Scalar >::magnitudeType
MDVector< Scalar, Node >::
normInf() const
{
  typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType mag;
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  mag local_normInf = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    local_normInf = std::max(local_normInf, std::abs(*it));
  mag global_normInf = 0;
  Teuchos::reduceAll(*_teuchosComm,
                     Teuchos::REDUCE_MAX,
                     1,
                     &local_normInf,
                     &global_normInf);
  return global_normInf;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
typename Teuchos::ScalarTraits< Scalar >::magnitudeType
MDVector< Scalar, Node >::
normWeighted(const MDVector< Scalar, Node > & weights) const
{
  typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType mag;
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! _mdMap->isCompatible(*(weights._mdMap)),
    MDMapError,
    "MDMap of calling MDVector and argument 'weights' are incompatible");

  MDArrayView< const Scalar > wView = weights.getData();
  mag local_wNorm = 0;
  iterator w_it = wView.begin();
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end();
       ++it, ++w_it)
    local_wNorm += *it * *it * *w_it;
  mag global_wNorm = 0;
  Teuchos::reduceAll(*_teuchosComm,
                     Teuchos::REDUCE_SUM,
                     1,
                     &local_wNorm,
                     &global_wNorm);
  Teuchos::Array< dim_type > dimensions(numDims());
  for (int i = 0; i < numDims(); ++i)
    dimensions[i] = _mdMap->getGlobalDim(i);
  size_type n = computeSize(dimensions);
  if (n == 0) return 0;
  return Teuchos::ScalarTraits<mag>::squareroot(global_wNorm / n);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Scalar
MDVector< Scalar, Node >::
meanValue() const
{
  typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType mag;
  typedef typename MDArrayView< const Scalar >::iterator iterator;

  mag local_sum = 0;
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    local_sum += *it;
  mag global_sum = 0;
  Teuchos::reduceAll(*_teuchosComm,
                     Teuchos::REDUCE_SUM,
                     1,
                     &local_sum,
                     &global_sum);
  Teuchos::Array< dim_type > dimensions(numDims());
  for (int i = 0; i < numDims(); ++i)
    dimensions[i] = _mdMap->getGlobalDim(i);
  size_type n = computeSize(dimensions);
  if (n == 0) return 0;
  return global_sum / n;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
std::string
MDVector< Scalar, Node >::
description() const
{
  using Teuchos::TypeNameTraits;

  Teuchos::Array< dim_type > dims(numDims());
  for (int axis = 0; axis < numDims(); ++axis)
    dims[axis] = getGlobalDim(axis, true);

  std::ostringstream oss;
  oss << "\"Domi::MDVector\": {"
      << "Template parameters: {"
      << "Scalar: " << TypeNameTraits<Scalar>::name()
      << ", Node: " << TypeNameTraits< Node >::name()
      << "}";
  if (this->getObjectLabel() != "")
    oss << ", Label: \"" << this->getObjectLabel () << "\", ";
  oss << "Global dimensions: " << dims << " }";
  return oss.str();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
describe(Teuchos::FancyOStream &out,
         const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  const MDMap< Node > & mdMap = *(getMDMap());
  Teuchos::RCP< const Teuchos::Comm< int > > comm = mdMap.getTeuchosComm();
  const int myImageID = comm->getRank();
  const int numImages = comm->getSize();
  Teuchos::OSTab tab0(out);

  if (vl != VERB_NONE)
  {
    if (myImageID == 0)
    {
      out << "\"Domi::MDVector\":" << endl;
    }
    Teuchos::OSTab tab1(out);// applies to all processes
    if (myImageID == 0)
    {
      out << "Template parameters:";
      {
        Teuchos::OSTab tab2(out);
        out << "Scalar: " << TypeNameTraits<Scalar>::name() << endl
            << "Node: " << TypeNameTraits< Node >::name() << endl;
      }
      out << endl;
      if (this->getObjectLabel() != "")
      {
        out << "Label: \"" << getObjectLabel() << "\"" << endl;
      }
      Teuchos::Array< dim_type > globalDims(numDims());
      for (int axis = 0; axis < numDims(); ++axis)
        globalDims[axis] = getGlobalDim(axis, true);
      out << "Global dimensions: " << globalDims << endl;
    }
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr)
    {
      if (myImageID == imageCtr)
      {
        if (vl != VERB_LOW)
        {
          // VERB_MEDIUM and higher prints getLocalLength()
          out << "Process: " << myImageID << endl;
          Teuchos::OSTab tab2(out);
          Teuchos::Array< dim_type > localDims(numDims());
          for (int axis = 0; axis < numDims(); ++axis)
            localDims[axis] = getLocalDim(axis, true);
          out << "Local dimensions: " << localDims << endl;
        }
        std::flush(out); // give output time to complete
      }
      comm->barrier(); // give output time to complete
      comm->barrier();
      comm->barrier();
    }
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
putScalar(const Scalar & value,
          bool includePadding)
{
  typedef typename MDArrayView< Scalar >::iterator iterator;

  MDArrayView< Scalar > newArray = getDataNonConst(includePadding);
  for (iterator it = newArray.begin(); it != newArray.end(); ++it)
    *it = value;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
randomize()
{
  typedef typename MDArrayView< Scalar >::iterator iterator;

  Teuchos::ScalarTraits< Scalar >::seedrandom(time(NULL));
  for (iterator it = _mdArrayView.begin(); it != _mdArrayView.end(); ++it)
    *it = Teuchos::ScalarTraits< Scalar >::random();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
updateCommPad()
{
  for (int axis = 0; axis < numDims(); ++axis)
  {
    updateCommPad(axis);
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
updateCommPad(int axis)
{
  startUpdateCommPad(axis);
  endUpdateCommPad(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
startUpdateCommPad(int axis)
{
  // #define DOMI_MDVECTOR_OUTPUT_UPDATECOMMPAD

  // Initialize the _sendMessages and _recvMessages members on the
  // first call to startUpdateCommPad(int).
  if (_sendMessages.empty()) initializeMessages();

#ifdef HAVE_MPI
  int rank    = _teuchosComm->getRank();
  int numProc = _teuchosComm->getSize();
  int tag;
  // Since HAVE_MPI is defined, we know that _teuchosComm points to a
  // const Teuchos::MpiComm< int >.  We downcast, extract and
  // dereference so that we can get access to the MPI_Comm used to
  // construct it.
  Teuchos::RCP< const Teuchos::MpiComm< int > > mpiComm =
    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm< int > >(_teuchosComm);
  const Teuchos::OpaqueWrapper< MPI_Comm > & communicator =
    *(mpiComm->getRawMpiComm());

  // Post the non-blocking sends
  MPI_Request request;
  for (int boundary = 0; boundary < 2; ++boundary)
  {
    MessageInfo message = _sendMessages[axis][boundary];
    if (message.proc >= 0)
    {
      tag = 2 * (rank * numProc + message.proc) + boundary;

#ifdef DOMI_MDVECTOR_OUTPUT_UPDATECOMMPAD
      cout << rank << ": post send for axis " << axis << ", boundary "
           << boundary << ", buffer = " << message.buffer << ", proc = "
           << message.proc << ", tag = " << tag << endl;
#endif

      if (MPI_Isend(message.buffer,
                    1,
                    *(message.datatype),
                    message.proc,
                    tag,
                    communicator(),
                    &request))
        throw std::runtime_error("Domi::MDVector: Error in MPI_Isend");
      _requests.push_back(request);
    }
  }

  // Post the non-blocking receives
  for (int boundary = 0; boundary < 2; ++boundary)
  {
    MessageInfo message = _recvMessages[axis][boundary];
    if (message.proc >= 0)
    {
      tag = 2 * (message.proc * numProc + rank) + (1-boundary);

#ifdef DOMI_MDVECTOR_OUTPUT_UPDATECOMMPAD
      cout << rank << ": post recv for axis " << axis << ", boundary "
           << boundary << ", buffer = " << message.buffer << ", proc = "
           << message.proc << ", tag = " << tag << endl;
#endif

      if (MPI_Irecv(message.buffer,
                    1,
                    *(message.datatype),
                    message.proc,
                    tag,
                    communicator(),
                    &request))
        throw std::runtime_error("Domi::MDVector: Error in MPI_Irecv");
      _requests.push_back(request);
    }
  }
#else
  // HAVE_MPI is not defined, so we are on a single processor.
  // However, if the axis is periodic, we need to copy the appropriate
  // data to the communication padding.
  if (isPeriodic(axis))
  {
    for (int sendBndry = 0; sendBndry < 2; ++sendBndry)
    {
      int recvBndry = 1 - sendBndry;
      // Get the receive and send data views
      MDArrayView< Scalar > recvView = _recvMessages[axis][recvBndry].dataview;
      MDArrayView< Scalar > sendView = _sendMessages[axis][sendBndry].dataview;

      // Initialize the receive and send data view iterators
      typename MDArrayView< Scalar >::iterator it_recv = recvView.begin();
      typename MDArrayView< Scalar >::iterator it_send = sendView.begin();

      // Copy the send buffer to the receive buffer
      for ( ; it_recv != recvView.end(); ++it_recv, ++it_send)
        *it_recv = *it_send;
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
endUpdateCommPad(int axis)
{
#ifdef HAVE_MPI
  if (_requests.size() > 0)
  {
    Teuchos::Array< MPI_Status > status(_requests.size());
    if (MPI_Waitall(_requests.size(),
                    &(_requests[0]),
                    &(status[0]) ) )
      throw std::runtime_error("Domi::MDVector: Error in MPI_Waitall");
    _requests.clear();
  }
#endif
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >
MDVector< Scalar, Node >::
operator[](dim_type index) const
{
  MDVector< Scalar, Node > result(*this,
                                  _nextAxis,
                                  index);
  int newAxis = _nextAxis + 1;
  if (newAxis >= result.numDims()) newAxis = 0;
  result._nextAxis = newAxis;
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
MDVector< Scalar, Node >
MDVector< Scalar, Node >::
operator[](Slice slice) const
{
  MDVector< Scalar, Node > result(*this,
                                  _nextAxis,
                                  slice     );
  int newAxis = _nextAxis + 1;
  if (newAxis >= result.numDims()) newAxis = 0;
  result._nextAxis = newAxis;
  return result;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
initializeMessages()
{
  // #define DOMI_MDVECTOR_MESSAGE_INITIALIZE

  int ndims = numDims();
  Teuchos::Array<int> sizes(ndims);
  Teuchos::Array<int> subsizes(ndims);
  Teuchos::Array<int> starts(ndims);
  MessageInfo messageInfo;

  _sendMessages.resize(ndims);
  _recvMessages.resize(ndims);

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
  std::stringstream msg;
  int rank = _teuchosComm->getRank();
#endif

#ifdef HAVE_MPI
  int order = mpiOrder(getLayout());
  MPI_Datatype datatype = mpiType< Scalar >();
#endif

  // Loop over the axes we are going to send messages along
  for (int msgAxis = 0; msgAxis < ndims; ++msgAxis)
  {
    // Set the initial values for sizes, subsizes and starts
    for (int axis = 0; axis < ndims; ++axis)
    {
      sizes[axis]    = _mdArrayRcp.dimension(axis);
      subsizes[axis] = _mdArrayView.dimension(axis);
      starts[axis]   = 0;
    }

    //////////////////////////////////////////
    // Set the lower receive and send messages
    //////////////////////////////////////////

    int proc = getLowerNeighbor(msgAxis);

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
    msg << endl << "P" << rank << ": axis " << msgAxis << ", lower neighbor = "
        << proc << endl;
#endif

    // Fix the subsize along this message axis
    subsizes[msgAxis] = getLowerPadSize(msgAxis);
    // Fix the proc if the subsize is zero
    if (subsizes[msgAxis] == 0) proc = -1;
    // Assign the non-MPI members of messageInfo
    messageInfo.buffer = (void*) getData().getRawPtr();
    messageInfo.proc   = proc;
    messageInfo.axis   = msgAxis;

    if (proc >= 0)
    {
      ////////////////////////
      // Lower receive message
      ////////////////////////

#ifdef HAVE_MPI

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
      msg << endl << "P" << rank << ": axis " << msgAxis
          << ", Lower receive message:" << endl << "  ndims    = " << ndims
          << endl << "  sizes    = " << sizes << endl << "  subsizes = "
          << subsizes << endl << "  starts   = " << starts << endl
          << "  order    = " << order << endl;
#endif
      Teuchos::RCP< MPI_Datatype > commPad = Teuchos::rcp(new MPI_Datatype);
      MPI_Type_create_subarray(ndims,
                               &sizes[0],
                               &subsizes[0],
                               &starts[0],
                               order,
                               datatype,
                               commPad.get());
      MPI_Type_commit(commPad.get());
      messageInfo.datatype = commPad;
#else
      messageInfo.dataview = _mdArrayView;
      for (int axis = 0; axis < numDims(); ++axis)
      {
        Slice slice(starts[axis], starts[axis] + subsizes[axis]);
        messageInfo.dataview = MDArrayView< Scalar >(messageInfo.dataview,
                                                     axis,
                                                     slice);
      }
#endif

    }
    _recvMessages[msgAxis][0] = messageInfo;

    /////////////////////
    // Lower send message
    /////////////////////

    starts[msgAxis] = getLowerPadSize(msgAxis);
    if (isReplicatedBoundary(msgAxis) && getCommIndex(msgAxis) == 0)
      starts[msgAxis] += 1;
    if (proc >= 0)
    {

#ifdef HAVE_MPI

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
      msg << endl << "P" << rank << ": axis " << msgAxis
          << ", Lower send message:" << endl << "  ndims    = " << ndims
          << endl << "  sizes    = " << sizes << endl << "  subsizes = "
          << subsizes << endl << "  starts   = " << starts << endl
          << "  order    = " << order << endl;
#endif
      Teuchos::RCP< MPI_Datatype > commPad = Teuchos::rcp(new MPI_Datatype);
      MPI_Type_create_subarray(ndims,
                               &sizes[0],
                               &subsizes[0],
                               &starts[0],
                               order,
                               datatype,
                               commPad.get());
      MPI_Type_commit(commPad.get());
      messageInfo.datatype = commPad;
#else
      messageInfo.dataview = _mdArrayView;
      for (int axis = 0; axis < numDims(); ++axis)
      {
        Slice slice(starts[axis], starts[axis] + subsizes[axis]);
        messageInfo.dataview = MDArrayView< Scalar >(messageInfo.dataview,
                                                     axis,
                                                     slice);
      }
#endif

    }
    _sendMessages[msgAxis][0] = messageInfo;

    //////////////////////////////////////////
    // Set the upper receive and send messages
    //////////////////////////////////////////

    proc = getUpperNeighbor(msgAxis);

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
    msg << endl << "P" << rank << ": axis " << msgAxis << ", upper neighbor = "
        << proc << endl;
#endif

    subsizes[msgAxis] = getUpperPadSize(msgAxis);
    starts[msgAxis]  = _mdArrayView.dimension(msgAxis) -
                       getUpperPadSize(msgAxis);
    if (subsizes[msgAxis] == 0) proc = -1;
    messageInfo.proc = proc;
    if (proc >= 0)
    {
      ////////////////////////
      // Upper receive message
      ////////////////////////

#ifdef HAVE_MPI

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
      msg << endl << "P" << rank << ": axis " << msgAxis
          << ", Upper receive message:" << endl << "  ndims    = " << ndims
          << endl << "  sizes    = " << sizes << endl << "  subsizes = "
          << subsizes << endl << "  starts   = " << starts << endl
          << "  order    = " << order << endl;
#endif
      Teuchos::RCP< MPI_Datatype > commPad = Teuchos::rcp(new MPI_Datatype);
      MPI_Type_create_subarray(ndims,
                               &sizes[0],
                               &subsizes[0],
                               &starts[0],
                               order,
                               datatype,
                               commPad.get());
      MPI_Type_commit(commPad.get());
      messageInfo.datatype = commPad;
#else
      messageInfo.dataview = _mdArrayView;
      for (int axis = 0; axis < numDims(); ++axis)
      {
        Slice slice(starts[axis], starts[axis] + subsizes[axis]);
        messageInfo.dataview = MDArrayView< Scalar >(messageInfo.dataview,
                                                     axis,
                                                     slice);
      }
#endif
    }
    _recvMessages[msgAxis][1] = messageInfo;

    /////////////////////
    // Upper send message
    /////////////////////

    starts[msgAxis] -= getUpperPadSize(msgAxis);
    if (isReplicatedBoundary(msgAxis) &&
        getCommIndex(msgAxis) == getCommDim(msgAxis)-1)
      starts[msgAxis] -= 1;
    if (proc >= 0)
    {

#ifdef HAVE_MPI

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
      msg << endl << "P" << rank << ": axis " << msgAxis
          << ", Upper send message:" << endl << "  ndims    = " << ndims
          << endl << "  sizes    = " << sizes << endl << "  subsizes = "
          << subsizes << endl << "  starts   = " << starts << endl
          << "  order    = " << order << endl;
#endif
      Teuchos::RCP< MPI_Datatype > commPad = Teuchos::rcp(new MPI_Datatype);
      MPI_Type_create_subarray(ndims,
                               &sizes[0],
                               &subsizes[0],
                               &starts[0],
                               order,
                               datatype,
                               commPad.get());
      MPI_Type_commit(commPad.get());
      messageInfo.datatype = commPad;
#else
      messageInfo.dataview = _mdArrayView;
      for (int axis = 0; axis < numDims(); ++axis)
      {
        Slice slice(starts[axis], starts[axis] + subsizes[axis]);
        messageInfo.dataview = MDArrayView< Scalar >(messageInfo.dataview,
                                                     axis,
                                                     slice);
      }
#endif

    }
    _sendMessages[msgAxis][1] = messageInfo;
  }

#ifdef DOMI_MDVECTOR_MESSAGE_INITIALIZE
  for (int proc = 0; proc < _teuchosComm->getSize(); ++proc)
  {
    if (proc == rank)
    {
      cout << msg.str();
      msg.flush();
    }
    _teuchosComm->barrier();
  }
#endif

}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
writeBinary(const std::string & filename,
            bool includeBndryPad) const
{
  FILE * datafile;
  // If we are using MPI and overwriting an existing file, and the new
  // file is shorter than the old file, it appears that the new file
  // will retain the old file's length and trailing data.  To prevent
  // this behavior, we open and close the file to give it zero length.
  int pid = _teuchosComm->getRank();
  if (pid == 0)
  {
    datafile = fopen(filename.c_str(), "w");
    fclose(datafile);
  }
  _teuchosComm->barrier();

  // Get the pointer to this MDVector's MDArray, including all padding
  const Scalar * buffer = getData(true).getRawPtr();

  // Compute either _fileInfo or _fileInfoWithBndry, whichever is
  // appropriate, and return a reference to that fileInfo object
  Teuchos::RCP< FileInfo > & fileInfo = computeFileInfo(includeBndryPad);

  // Parallel output
#ifdef HAVE_MPI

  // Since HAVE_MPI is defined, we know that _teuchosComm points to a
  // const Teuchos::MpiComm< int >.  We downcast, extract and
  // dereference so that we can get access to the MPI_Comm used to
  // construct it.
  Teuchos::RCP< const Teuchos::MpiComm< int > > mpiComm =
    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm< int > >(_teuchosComm);
  const Teuchos::OpaqueWrapper< MPI_Comm > & communicator =
    *(mpiComm->getRawMpiComm());

  // Compute the access mode
  int access = MPI_MODE_WRONLY | MPI_MODE_CREATE;

  // I copy the filename C string, because the c_str() method returns
  // a const char*, and the MPI_File_open() function requires
  // (incorrectly) a non-const char*.
  char * cstr = new char[filename.size()+1];
  std::strcpy(cstr, filename.c_str());

  // Use MPI I/O to write the binary file
  MPI_File   mpiFile;
  MPI_Status status;
  char       datarep[7] = "native";
  MPI_File_open(communicator(), cstr, access, MPI_INFO_NULL, &mpiFile);
  MPI_File_set_view(mpiFile, 0, mpiType< Scalar >(),
                    *(fileInfo->filetype), datarep, MPI_INFO_NULL);
  MPI_File_write_all(mpiFile, (void*)buffer, 1, *(fileInfo->datatype),
                     &status);
  MPI_File_close(&mpiFile);

  // Delete the C string
  delete [] cstr;

  // Serial output
#else

  // Get the number of dimensions
  // int ndims = _mdMap->numDims();

  // Initialize the data file
  datafile = fopen(filename.c_str(), "w");

  // Obtain the data to write, including the boundary padding if requested
  MDArrayView< const Scalar > mdArrayView = getData(includeBndryPad);

  // Iterate over the data and write it to the data file
  typedef typename MDArrayView< Scalar >::const_iterator const_iterator;
  for (const_iterator it = mdArrayView.begin(); it != mdArrayView.end(); ++it)
  {
    fwrite((const void *) &(*it), sizeof(Scalar), 1, datafile);
  }

  // Close the data file
  fclose(datafile);

#endif

}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
void
MDVector< Scalar, Node >::
readBinary(const std::string & filename,
           bool includeBndryPad)
{
  // Get the pointer to this MDVector's MDArray, including all padding
  const Scalar * buffer = getDataNonConst(true).getRawPtr();

  // Compute either _fileInfo or _fileInfoWithBndry, whichever is
  // appropriate, and return a reference to that fileInfo object
  Teuchos::RCP< FileInfo > & fileInfo = computeFileInfo(includeBndryPad);

  // Parallel input
#ifdef HAVE_MPI

  // Since HAVE_MPI is defined, we know that _teuchosComm points to a
  // const Teuchos::MpiComm< int >.  We downcast, extract and
  // dereference so that we can get access to the MPI_Comm used to
  // construct it.
  Teuchos::RCP< const Teuchos::MpiComm< int > > mpiComm =
    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm< int > >(_teuchosComm);
  const Teuchos::OpaqueWrapper< MPI_Comm > & communicator =
    *(mpiComm->getRawMpiComm());

  // Compute the access mode
  int access = MPI_MODE_RDONLY;

  // I copy the filename C string, because the c_str() method returns
  // a const char*, and the MPI_File_open() function requires
  // (incorrectly) a non-const char*.
  char * cstr = new char[filename.size()+1];
  std::strcpy(cstr, filename.c_str());

  // Use MPI I/O to read the binary file
  MPI_File   mpiFile;
  MPI_Status status;
  char       datarep[7] = "native";
  MPI_File_open(communicator(), cstr, access, MPI_INFO_NULL, &mpiFile);
  MPI_File_set_view(mpiFile, 0, mpiType< Scalar >(),
                    *(fileInfo->filetype), datarep, MPI_INFO_NULL);
  MPI_File_read_all(mpiFile, (void*)buffer, 1, *(fileInfo->datatype),
                    &status);
  MPI_File_close(&mpiFile);

  // Delete the C string
  delete [] cstr;

  // Serial output
#else

  // Get the number of dimensions
  int ndims = _mdMap->numDims();

  // Initialize the data file
  FILE * datafile;
  datafile = fopen(filename.c_str(), "r");

  // Obtain the MDArrayView to read into, including the boundary
  // padding if requested
  MDArrayView< Scalar > mdArrayView = getDataNonConst(includeBndryPad);

  // Initialize the set of indexes
  Teuchos::Array< Ordinal > index(3);
  for (int axis = 0; axis < ndims; ++axis)
    index[axis] = fileInfo->dataStart[axis];

  // Iterate over the data and read it from the data file
  typedef typename MDArrayView< Scalar >::iterator iterator;
  for (iterator it = mdArrayView.begin(); it != mdArrayView.end(); ++it)
  {
    fread(&(*it), sizeof(Scalar), 1, datafile);
  }

  // Close the data file
  fclose(datafile);

#endif

}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class Node >
Teuchos::RCP< typename MDVector< Scalar, Node >::FileInfo > &
MDVector< Scalar, Node >::
computeFileInfo(bool includeBndryPad) const
{
  // Work with the appropriate FileInfo object.  By using a reference
  // here, we are working directly with the member data.
  Teuchos::RCP< MDVector< Scalar, Node >::FileInfo > & fileInfo =
    includeBndryPad ? _fileInfoWithBndry : _fileInfo;

  // If the fileInfo object already has been set, our work is done
  if (!fileInfo.is_null()) return fileInfo;

  // Initialize the new FileInfo object
  int ndims = _mdMap->numDims();
  fileInfo.reset(new MDVector< Scalar, Node >::FileInfo);
  fileInfo->fileShape.resize(ndims);
  fileInfo->bufferShape.resize(ndims);
  fileInfo->dataShape.resize(ndims);
  fileInfo->fileStart.resize(ndims);
  fileInfo->dataStart.resize(ndims);

  // Initialize the shapes and starts.
  for (int axis = 0; axis < ndims; ++axis)
  {
    // Initialize the FileInfo arrays using includeBndryPad where
    // appropriate
    fileInfo->fileShape[axis]   = getGlobalDim(axis,includeBndryPad);
    fileInfo->bufferShape[axis] = getLocalDim(axis,true );
    fileInfo->dataShape[axis]   = getLocalDim(axis,false);
    fileInfo->fileStart[axis]   = getGlobalRankBounds(axis,includeBndryPad).start();
    fileInfo->dataStart[axis]   = getLocalBounds(axis).start();
    // Modify dataShape and dataStart if boundary padding is included
    if (includeBndryPad)
    {
      int commIndex = _mdMap->getCommIndex(axis);
      if (commIndex == 0)
      {
        int pad = getLowerBndryPad(axis);
        fileInfo->dataShape[axis] += pad;
        fileInfo->dataStart[axis] -= pad;
      }
      if (commIndex == _mdMap->getCommDim(axis)-1)
      {
        fileInfo->dataShape[axis] += getUpperBndryPad(axis);
      }
    }
  }

#ifdef DOMI_MDVECTOR_DEBUG_IO
  cout << pid << ": fileShape   = " << fileInfo->fileShape()   << endl;
  cout << pid << ": bufferShape = " << fileInfo->bufferShape() << endl;
  cout << pid << ": dataShape   = " << fileInfo->dataShape()   << endl;
  cout << pid << ": fileStart   = " << fileInfo->fileStart()   << endl;
  cout << pid << ": dataStart   = " << fileInfo->dataStart()   << endl;
#endif

#ifdef HAVE_MPI
  int mpi_order = getLayout() == C_ORDER ? MPI_ORDER_C : MPI_ORDER_FORTRAN;
  // Build the MPI_Datatype for the file
  fileInfo->filetype = Teuchos::rcp(new MPI_Datatype);
  MPI_Type_create_subarray(ndims,
                           fileInfo->fileShape.getRawPtr(),
                           fileInfo->dataShape.getRawPtr(),
                           fileInfo->fileStart.getRawPtr(),
                           mpi_order,
                           mpiType< Scalar >(),
                           fileInfo->filetype.get());
  MPI_Type_commit(fileInfo->filetype.get());

  // Build the MPI_Datatype for the data
  fileInfo->datatype = Teuchos::rcp(new MPI_Datatype);
  MPI_Type_create_subarray(ndims,
                           fileInfo->bufferShape.getRawPtr(),
                           fileInfo->dataShape.getRawPtr(),
                           fileInfo->dataStart.getRawPtr(),
                           mpi_order,
                           mpiType< Scalar >(),
                           fileInfo->datatype.get());
  MPI_Type_commit(fileInfo->datatype.get());
#endif  // DGM_PARALLEL

  return fileInfo;
}

////////////////////////////////////////////////////////////////////////

}  // Namespace Domi

#endif
