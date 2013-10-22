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

#ifndef DOMI_MDVECTOR_HPP
#define DOMI_MDVECTOR_HPP

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_MDMap.hpp"

// Teuchos includes
#include "Teuchos_Describable.hpp"

namespace Domi
{

template< class Scalar,
          class LocalOrd,
          class GlobalOrd = LocalOrd,
          class Node = Kokkos::DefaultNode::DefaultNodeType >
class MDVector : Teuchos::Describable
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
  MDVector(const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > > & mdMap,
           bool zeroOut = true);

  /** \brief Constructor with initialization values (copy)
   *
   * \param mdMap [in] MDMap that describes the domain decomposition
   *        of this MDVector
   *
   * \param source [in] initialization values
   */
  MDVector(const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > > & mdMap,
           const MDArrayView< const Scalar > & source);

  /** \brief Copy constructor
   *
   * \param source [in] source MDVector
   */
  MDVector(const MDVector<Scalar, LocalOrd, GlobalOrd, Node > & source);

  /** \brief Parent/single global ordinal sub-vector constructor
   *
   * \param parent [in] an MDVector, from which this MDVector will be
   *        derived
   *
   * \param axis [in] the axis to which this index ordinal applies
   *
   * \param index [in] the global ordinal that defines the sub-vector
   */
  MDVector(const MDVector< Scalar, LocalOrd, GlobalOrd, Node > & parent,
           int axis,
           GlobalOrd index);

  /** \brief Parent/single slice sub-vector constructor
   *
   * \param parent [in] an MDVector, from which this MDVector will be
   *        derived
   *
   * \param axis [in] the axis to which this slice applies
   *
   * \param ghosts [in] the ghost region along the altered axis of the
   *        new sub-vector.  This may include indexes from the ghost
   *        region of the parent MDVector, but it does not have to.
   */
  MDVector(const MDVector< Scalar, LocalOrd, GlobalOrd, Node > & parent,
           int axis,
           const Slice & slice,
           int ghosts = 0);

  /** \brief Create a cloned MDVector for a different node type
   *
   * \param node2 [in] the new node
   */
  template< class Node2 >
  Teuchos::RCP< MDVector< Scalar, LocalOrd, GlobalOrd, Node2 >
  clone(const Teuchos::RCP< Node2 > & node2) const;

  /** \brief Destructor 
   */
  virtual ~MDVector();

  //@}

  /** \name MDMap methods */
  //@{

  /** \brief MDMap accessor method
   */
  const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > >
  getMdMap() const;

  /** \brief Query whether this processor is on the sub-communicator
   *
   * Sub-communicators are formed when a parent MDVector is sliced by
   * using the (parent,ordinal) or (parent,slice) constructors.  For a
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
   */
  LocalOrd getLocalDim(int axis, bool withHalos=false) const;

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

  /** \name Data extraction methods */
  //@{

  /** \brief Get a non-const view of the data as an MDArrayRCP
   */
  MDArrayRCP< Scalar > getDataNonConst();

  /** \brief Get a const view of the data as an MDArrayRCP
   */
  MDArrayRCP< const Scalar > getData() const;

  //@}

  /** \name Mathematical methods */
  //@{

  /** \brief Compute the dot product of this MDVector and MDVector a
   */
  Scalar
  dot(const MDVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > & a) const;

  /** \brief Compute the 1-norm of this MDVector
   */
  Teuchos::ScalarTraits< Scalar >::magnitudeType norm1() const;

  /** \brief Compute the 2-norm of this MDVector
   */
  Teuchos::ScalarTraits< Scalar >::magnitudeType norm2() const;

  /** \brief Compute the infinity-norm of this MDVector
   */
  Teuchos::ScalarTraits< Scalar >::magnitudeType normInf() const;

  /** \brief Compute the weighted norm of this 
   */
  Teuchos::ScalarTraits< Scalar >::magnitudeType
  normWeighted(const MDVector< Scalar,
                               LocalOrdinal,
                               GlobalOrdinal,
                               Node > & weights) const;

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
             FancyOStream
   */
  virtual void
  describe(Teuchos::FancyOStream &out,
           const Teuchos::EVerbosityLevel verbLevel =
             Teuchos::Describable::verbLevel_default) const;

  //@}

  /** \name Global assignment methods */
  //@{

  /// Set all values in the multivector with the given value. 
  void putScalar(const Scalar & value);

  /// Set all values in the multivector to pseudorandom numbers. 
  void randomize();

  //@}

  /** \name Global communication methods */
  //@{

  /// Sum values of a locally replicated multivector across all processes.
  void reduce();

  /// Update the halos along the specified axis
  void updateHalos(int axis);

  /// Update the halos along all axes
  void updateHalos();

  //@}

  /** \name Sub-MDVector operators */
  //@{

  /** \brief Sub-vector access operator.  The returned
   *  <tt>MDVector</tt> will have one fewer dimensions than the
   *  calling <tt>MDVector</tt>.
   *
   * \param index [in] index of the desired sub-vector.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDVector</tt>.
   */
  MDVector< Scalar, LocalOrd, GlobalOrd, Node >
  operator[](GlobalOrd index);

  /** \brief Sub-vector access operator.  The returned
   *  <tt>MDVector</tt> will have the same number of dimensions as the
   *  calling <tt>MDVector</tt>.
   *
   * \param slice [in] index of the desired sub-vector.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDVector</tt>.
   */
  MDVector< Scalar, LocalOrd, GlobalOrd, Node >
  operator[](Slice slice);

  //@}

private:

  // The MDMap that describes the domain decomposition of this
  // MDVector
  const RCP< const MDMap< LocalOrd, GlobalOrd, Node > > _mdMap;

  // The MDArrayRCP that stores the data of this MDVector
  const MDArrayRCP< Scalar > _mdArray;

};

/////////////////////
// Implementations //
/////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
MDVector(const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > > & mdMap,
         bool zeroOut) :
  Teuchos::Describable("Domi::MDVector"),
  _mdMap(mdMap),
  _mdArray()
{
  typedef Teuchos::ArrayView< Scalar >::size_type size_type;
  int numDims = _mdMap->getNumDims();
  Teuchos::Array< size_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = _mdMap->getLocalDim(axis);
  _mdArray.resize(dims);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
MDVector(const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > > & mdMap,
         const MDArrayView< const Scalar > & source) :
  Teuchos::Describable("Domi::MDVector"),
  _mdMap(mdMap),
  _mdArray(source)
{
  int numDims = _mdMap->getNumDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numDims != _mdArray.getNumDims(),
    InvalidArgument,
    "MDMap and source array do not have the same number of dimensions");

  for (int axis = 0; axis < numDims; ++axis)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      _mdMap->getLocalDim(axis) != _mdArray.getLocalDim(axis),
      InvalidArgument,
      "Axis " << axis << ": MDMap dimension = " << _mdMap->getLocalDim(axis)
      << ", MDArray dimension = " << _mdArray.getLocalDim(axis));
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
MDVector(const MDVector<Scalar, LocalOrd, GlobalOrd, Node > & source) :
  Teuchos::Describable("Domi::MDVector"),
  _mdMap(source.getMdMap()),
  _mdArray(source.getDataNonConst())
{
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
~MDVector();

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
const Teuchos::RCP< const MDMap< LocalOrd, GlobalOrd, Node > >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getMdMap() const
{
  return _mdMap;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
bool
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
onSubcommunicator() const
{
  return _mdMap->onSubcommunicator();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
TeuchosCommRCP 
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getTeuchosComm() const
{
  return _mdMap->getTeuchosComm();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
int
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getNumDims() const
{
  return _mdMap->getNumDims();
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
int
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getAxisCommSize(int axis) const
{
  return _mdMap->getAxisCommSize(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
bool
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
isPeriodic(int axis) const
{
  return _mdMap->isPeriodic(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
int
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getAxisRank(int axis) const
{
  return _mdMap->getAxisRank(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
int
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getLowerNeighbor(int axis) const
{
  return _mdMap->getLowerNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
int
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getUpperNeighbor(int axis) const
{
  return _mdMap->getUpperNeighbor(axis);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
GlobalOrd
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getGlobalDim(int axis, bool withGhosts=false) const
{
  return _mdMap->getGlobalDim(axis, withGhosts);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
Slice
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getGlobalBounds(int axis, bool withGhosts=false) const
{
  return _mdMap->getGlobalBounds(axis, withGhosts);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
Slice
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getGlobalRankBounds(int axis, bool withGhosts=false) const
{
  return _mdMap->getGlobalRankBounds(axis, withGhosts);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
LocalOrd
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getLocalDim(int axis, bool withHalos=false) const
{
  return _mdMap->getLocalDim(axis, withHalos);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
Slice
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getLocalBounds(int axis, bool withHalos=false) const
{
  return _mdMap->getLocalBounds(axis, withHalos);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDArrayRCP< Scalar >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getDataNonConst()
{
  return _mdArray;
}

////////////////////////////////////////////////////////////////////////

template< class Scalar,
          class LocalOrd,
          class GlobalOrd,
          class Node >
MDArrayRCP< const Scalar >
MDVector< Scalar, LocalOrd, GlobalOrd, Node >::
getData() const
{
  return _mdArray;
}

}

#endif
