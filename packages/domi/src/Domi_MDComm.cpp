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

// Domi includes
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"
#include "Domi_getValidParameters.hpp"
#include "Domi_MDComm.hpp"

// Teuchos includes
#include "Teuchos_TestForException.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Domi
{

MDComm::MDComm(const TeuchosCommRCP teuchosComm,
               const Teuchos::ArrayView< int > & axisCommSizes,
               const Teuchos::ArrayView< int > & periodic) :
  _teuchosComm(teuchosComm),
  _axisCommSizes(regularizeAxisSizes(teuchosComm->getSize(),
                                     axisCommSizes.size(),
                                     axisCommSizes)),
  _periodic(computePeriodic(axisCommSizes.size(), periodic)),
  _axisRanks(computeAxisRanks(teuchosComm->getRank(),
                              _axisCommSizes)),
  _axisStrides(computeStrides(_axisCommSizes,
                              DEFAULT_ORDER))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const TeuchosCommRCP teuchosComm,
               Teuchos::ParameterList & plist) :
  _teuchosComm(teuchosComm)
{
  // Validate the ParameterList
  plist.validateParameters(*getValidParameters());

  // Set the communicator sizes along each axis
  Teuchos::Array< int > axisCommSizes =
    plist.get("axis comm sizes", Teuchos::Array< int >());
  _axisCommSizes = regularizeAxisSizes(_teuchosComm->getSize(),
                                       axisCommSizes.size(),
                                       axisCommSizes());

  // Set the periodic flags along each axis
  Teuchos::Array< int > periodic =
    plist.get("periodic", Teuchos::Array< int >());
  _periodic = computePeriodic(_axisCommSizes.size(),
                              periodic);

  // Set the axis ranks for this processor
  _axisRanks = computeAxisRanks(_teuchosComm->getRank(),
                                _axisCommSizes);

  // Set the axis strides
  _axisStrides = computeStrides(_axisCommSizes,
                                DEFAULT_ORDER);
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const TeuchosCommRCP teuchosComm,
               int numDims) :
  _teuchosComm(teuchosComm),
  _axisCommSizes(regularizeAxisSizes(teuchosComm->getSize(),
                                     numDims,
                                     Teuchos::Array< int >())),
  _periodic(numDims, 0),
  _axisRanks(computeAxisRanks(teuchosComm->getRank(),
                              _axisCommSizes)),
  _axisStrides(computeStrides(_axisCommSizes,
                              DEFAULT_ORDER))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const TeuchosCommRCP teuchosComm,
               int numDims,
               const Teuchos::ArrayView< int > & axisCommSizes,
               const Teuchos::ArrayView< int > & periodic) :
  _teuchosComm(teuchosComm),
  _axisCommSizes(regularizeAxisSizes(teuchosComm->getSize(),
                                     numDims,
                                     axisCommSizes)),
  _periodic(computePeriodic(numDims, periodic)),
  _axisRanks(computeAxisRanks(teuchosComm->getRank(),
                              _axisCommSizes)),
  _axisStrides(computeStrides(_axisCommSizes,
                              DEFAULT_ORDER))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               int axis,
               int axisRank) :
  _axisCommSizes(),
  _periodic(),
  _axisRanks(),
  _axisStrides()
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
      ((axisRank < 0) || (axisRank >= parent.getAxisCommSize(axis))),
      RangeError,
      "axisRank = " << axisRank  << " is invalid for communicator axis " <<
        axis << " with " << parent.getAxisCommSize(axis) << " processors");

    if (numDims == 1)
    {
      // Set the communicator
      Teuchos::Array< int > ranks(1);
      ranks[0] = axisRank;
      _teuchosComm =
        parent.getTeuchosComm()->createSubcommunicator(ranks()).getConst();
      // If this processor is on the sub-communicator, set the MDComm
      // data members
      if (_teuchosComm.getRawPtr())
      {
        _axisCommSizes.push_back(1);
        _periodic.push_back(0);
        _axisRanks.push_back(0);
        _axisStrides.push_back(1);
      }
    }
    else
    {
      // Prepare to set the _teuchosComm data member.  First compute the looping
      // bounds.
      Teuchos::Array< Slice > bounds(numDims);
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
        if (myAxis == axis)
          bounds[myAxis] = Slice(axisRank, axisRank+1);
        else
          bounds[myAxis] = Slice(0, parent._axisCommSizes[myAxis]);

      // Now compute the size of the ranks array and the starting mdIndex
      size_type rankSize = 1;
      Teuchos::Array< int > mdIndex(numDims);
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        rankSize *= (bounds[myAxis].stop() - bounds[myAxis].start());
        mdIndex[myAxis] = bounds[myAxis].start();
      }
      Teuchos::Array< int > ranks(rankSize);

      // Compute the ranks array and set _teuchosComm
      bool done = false;
      for (int index = 0; index < rankSize; ++index)
      {
        // Compute ranks[index]
        int val = 0;
        for (int myAxis = 0; myAxis < numDims; ++myAxis)
          val += mdIndex[myAxis] * parent._axisStrides[myAxis];
        ranks[index] = val;
        // Increment mdIndex
        int myAxis = 0;
        done = false;
        while (not done)
        {
          ++mdIndex[myAxis];
          done = (mdIndex[myAxis] < bounds[myAxis].stop());
          if (not done)
          {
            mdIndex[myAxis] = bounds[myAxis].start();
            ++myAxis;
            if (myAxis >= numDims)
              done = true;
          }
        }
      }
      _teuchosComm =
        parent.getTeuchosComm()->createSubcommunicator(ranks()).getConst();

      // Set the remaining data members
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        if (myAxis != axis)
        {
          _axisCommSizes.push_back(parent._axisCommSizes[myAxis]);
          _periodic.push_back(parent._periodic[myAxis]);
          _axisRanks.push_back(parent._axisRanks[myAxis]);
        }
      }
      _axisStrides.push_back(1);
      for (int myAxis = 1; myAxis < _axisCommSizes.size(); ++myAxis)
        _axisStrides.push_back(_axisStrides[myAxis-1] *
                               _axisCommSizes[myAxis-1]);
    }
  }
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               int axis,
               const Slice & slice) :
  _axisCommSizes(parent._axisCommSizes),
  _periodic(parent._periodic),
  _axisRanks(parent._axisRanks),
  _axisStrides(parent.getNumDims())
{
  if (parent.onSubcommunicator())
  {
    // Sanity check
    size_type numDims = parent.getNumDims();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for communicator with " <<
        numDims << " dimensions");

    // Make sure the Slice we work with is concrete, and adjust
    // _axisCommSizes, _axisRanks, and _axisStrides
    Slice bounds = slice.bounds(parent.getAxisCommSize(axis));
    _axisCommSizes[axis] = (bounds.stop() - bounds.start());
    _axisRanks[axis] -= bounds.start();

    // Fix the periodic flag
    if ((bounds.start() == 0) &&
        (bounds.stop() == parent.getAxisCommSize(axis)))
      _periodic[axis] = parent._periodic[axis];
    else
      _periodic[axis] = 0;

    // Compute the new strides of the new sub-communicator
    _axisStrides[0] = 1;
    for (int myAxis = 1; myAxis < numDims; ++myAxis)
      _axisStrides[myAxis] = _axisStrides[myAxis-1] * _axisCommSizes[myAxis-1];

    // Compute the ranks of the subcommunicator
    size_type rankSize = 1;
    for (int myAxis = 0; myAxis < numDims; ++myAxis)
      rankSize *= _axisCommSizes[myAxis];
    Teuchos::Array< int > ranks(rankSize);
    Teuchos::Array< int > mdIndex(numDims, 0);
    bool done = false;
    for (int index = 0; index < rankSize; ++index)
    {
      // Compute ranks[index]
      int val = 0;
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        int start = (axis == myAxis) ? bounds.start() : 0;
        val += (mdIndex[myAxis] + start) * parent._axisStrides[myAxis];
      }
      ranks[index] = val;
      // Increment mdIndex
      int myAxis = 0;
      done = false;
      while (not done)
      {
        ++mdIndex[myAxis];
        done = (mdIndex[myAxis] < _axisCommSizes[myAxis]);
        if (not done)
        {
          mdIndex[myAxis] = 0;
          ++myAxis;
          if (myAxis >= numDims)
            done = true;
        }
      }
    }

    // Set the communicator
    _teuchosComm =
      parent.getTeuchosComm()->createSubcommunicator(ranks()).getConst();

    // On processors that are not a part of the subcommunicator, reset
    // the data attributes.
    if (_teuchosComm.getRawPtr() == 0)
    {
      _axisCommSizes.clear();
      _axisRanks.clear();
      _axisStrides.clear();
    }
  }
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               const Teuchos::ArrayView< Slice > & slices)
{
  // Sanity check
  size_type numDims = parent.getNumDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    slices.size() != numDims,
    InvalidArgument,
    "Length of array of slices does not match "
    "number of dimension of parent MDComm");

  // Apply the single-Slice constructor to each axis in succession
  MDComm tempMDComm1(parent);
  for (int axis = 0; axis < numDims; ++axis)
  {
    MDComm tempMDComm2 = MDComm(tempMDComm1, axis, slices[axis]);
    tempMDComm1 = tempMDComm2;
  }
  *this = tempMDComm1;
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & source) :
  _teuchosComm(source._teuchosComm),
#ifdef HAVE_EPETRA
  _epetraComm(source._epetraComm),
#endif
  _axisCommSizes(source._axisCommSizes),
  _periodic(source._periodic),
  _axisRanks(source._axisRanks),
  _axisStrides(source._axisStrides)
{
}

////////////////////////////////////////////////////////////////////////

MDComm::~MDComm()
{
}

////////////////////////////////////////////////////////////////////////

MDComm & MDComm::operator=(const MDComm & source)
{
  _teuchosComm   = source._teuchosComm;
#ifdef HAVE_EPETRA
  _epetraComm    = source._epetraComm;
#endif
  _axisCommSizes = source._axisCommSizes;
  _periodic      = source._periodic;
  _axisRanks     = source._axisRanks;
  _axisStrides   = source._axisStrides;
  return *this;
}

////////////////////////////////////////////////////////////////////////

bool
MDComm::onSubcommunicator() const
{
  return _teuchosComm.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

TeuchosCommRCP
MDComm::getTeuchosComm() const
{
  return _teuchosComm;
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA
EpetraCommRCP
MDComm::getEpetraComm() const
{
  if (_epetraComm.is_null() && not _teuchosComm.is_null())
  {
#ifdef HAVE_MPI
    Teuchos::RCP< const Teuchos::MpiComm<int> > teuchosMpiComm =
      Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >(_teuchosComm);
    Teuchos::RCP< const Teuchos::OpaqueWrapper< MPI_Comm > > opaqueMpiComm =
      teuchosMpiComm->getRawMpiComm();
    MPI_Comm mpiComm = opaqueMpiComm->operator()();
    _epetraComm =
      Teuchos::rcp< const Epetra_Comm >(new Epetra_MpiComm(mpiComm));
#else
    _epetraComm =
      Teuchos::rcp< const Epetra_Comm >(new Epetra_SerialComm());
#endif
  }

  return _epetraComm;
}
#endif

////////////////////////////////////////////////////////////////////////

int
MDComm::getNumDims() const
{
  return _axisCommSizes.size();
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getAxisCommSize(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getAxisCommSize()");
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _axisCommSizes[axis];
}

////////////////////////////////////////////////////////////////////////

bool
MDComm::isPeriodic(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::isPeriodic()");
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _periodic[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getAxisRank(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getAxisRank()");
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  return _axisRanks[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getLowerNeighbor(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getLowerNeighbor()");
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (_axisRanks[axis] == 0)
  {
    if (_periodic[axis])
    {
      return _teuchosComm->getRank() +
        (_axisCommSizes[axis]-1) * _axisStrides[axis];
    }
    else
    {
      return -1;
    }
  }
  return _teuchosComm->getRank() - _axisStrides[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getUpperNeighbor(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getUpperNeighbor()");
#if DOMI_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= getNumDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    getNumDims() << ")");
#endif
  if (_axisRanks[axis] == _axisCommSizes[axis] - 1)
  {
    if (_periodic[axis])
    {
      return _teuchosComm->getRank() -
        (_axisCommSizes[axis]-1) * _axisStrides[axis];
    }
    else
    {
      return -1;
    }
  }
  return _teuchosComm->getRank() + _axisStrides[axis];
}

}    // End namespace Domi
