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

// Teuchos includes
#include "Teuchos_TestForException.hpp"

// Domi includes
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDComm.hpp"

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

MDComm::MDComm(const Teuchos::RCP< const MDComm > parent,
               const Teuchos::ArrayView< Slice > & slices) :
  _axisCommSizes(parent->getNumDims()),  // just allocate the size
  _periodic(parent->getNumDims(), 0),    // will be adjusted below...
  _axisRanks(parent->_axisRanks),        // will be adjusted below...
  _axisStrides(parent->_axisStrides)     // copy constructor
{
  // Sanity check
  size_type numDims = parent->getNumDims();
  TEUCHOS_TEST_FOR_EXCEPTION(
    slices.size() != numDims,
    InvalidArgument,
    "Length of array of slices does not match "
    "number of dimension of parent MDComm");
  // Make sure the array of Slices we work with is concrete, compute
  // the _axisCommSizes and _axisRanks arrays, and compute the size of
  // the ranks array
  size_type rankSize = 1;
  Teuchos::Array< Slice > bounds;
  for (int axis = 0; axis < numDims; ++axis)
  {
    bounds.push_back(slices[axis].bounds(parent->getAxisCommSize(axis)));
    _axisCommSizes[axis] = (bounds[axis].stop() - bounds[axis].start());
    _axisRanks[axis] -= bounds[axis].start();
    rankSize *= _axisCommSizes[axis];
  }
  // Compute the ranks of the subcommunicator
  Teuchos::Array< int > ranks(rankSize);
  Teuchos::Array< int > mdIndex(numDims, 0);
  bool done = false;
  for (int index = 0; index < rankSize; ++index)
  {
    // Compute ranks[index]
    int val = 0;
    for (int axis = 0; axis < numDims; ++axis)
      val += (mdIndex[axis] + bounds[axis].start()) * _axisStrides[axis];
    ranks[index] = val;
    // Increment mdIndex
    int axis = 0;
    done = false;
    while (not done)
    {
      ++mdIndex[axis];
      done = (mdIndex[axis] < _axisCommSizes[axis]);
      if (not done)
      {
        mdIndex[axis] = 0;
        ++axis;
        if (axis >= numDims)
          done = true;
      }
    }
  }
  // Set the communicator
  _teuchosComm =
    parent->getTeuchosComm()->createSubcommunicator(ranks()).getConst();
  // On processors that are not a part of the subcommunicator, reset
  // the data attributes.
  if (_teuchosComm.getRawPtr() == 0)
  {
    _axisCommSizes.clear();
    _axisRanks.clear();
    _axisStrides.clear();
  }
  // Fix the periodic flags
  for (int axis = 0; axis < numDims; ++axis)
  {
    if ((bounds[axis].start() == 0) &&
        (bounds[axis].stop() == parent->getAxisCommSize(axis)))
      _periodic[axis] = parent->_periodic[axis];
  }
}

////////////////////////////////////////////////////////////////////////

MDComm::~MDComm()
{
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
  return _axisCommSizes[axis];
}

////////////////////////////////////////////////////////////////////////

bool
MDComm::isPeriodic(int axis) const
{
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
  if (_axisRanks[axis] == 0)
    if (_periodic[axis])
      return _teuchosComm->getRank() +
        (_axisCommSizes[axis]-1) * _axisStrides[axis];
    else
      return -1;
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
  if (_axisRanks[axis] == _axisCommSizes[axis] - 1)
    if (_periodic[axis])
      return _teuchosComm->getRank() -
        (_axisCommSizes[axis]-1) * _axisStrides[axis];
    else
      return -1;
  return _teuchosComm->getRank() + _axisStrides[axis];
}

}    // End namespace Domi
