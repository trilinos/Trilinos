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

// Domi includes
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"
#include "Domi_getValidParameters.hpp"
#include "Domi_MDComm.hpp"

// Teuchos includes
#include "Teuchos_TestForException.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Domi
{

////////////////////////////////////////////////////////////////////////

const Layout MDComm::commLayout = C_ORDER;

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const Teuchos::ArrayView< const int > & commDims,
               const Teuchos::ArrayView< const int > & periodic) :
  _teuchosComm(Teuchos::DefaultComm< int >::getComm()),
  _commDims(regularizeCommDims(_teuchosComm->getSize(),
                               commDims.size(),
                               commDims)),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(_teuchosComm->getRank(),
                                _commStrides)),
  _periodic(createArrayOfInts(commDims.size(), periodic))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               const Teuchos::ArrayView< const int > & commDims,
               const Teuchos::ArrayView< const int > & periodic) :
  _teuchosComm(teuchosComm),
  _commDims(regularizeCommDims(teuchosComm->getSize(),
                               commDims.size(),
                               commDims)),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(teuchosComm->getRank(),
                                _commStrides)),
  _periodic(createArrayOfInts(commDims.size(), periodic))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(Teuchos::ParameterList & plist) :
  _teuchosComm(Teuchos::DefaultComm< int >::getComm())
{
  // Validate the ParameterList
  plist.validateParameters(*getValidParameters());

  // Determine a first cut at the number of dimensions
  Teuchos::Array< int > dims =
    plist.get("dimensions", Teuchos::Array< int >());
  int numDims = dims.size();

  // Set the communicator sizes along each axis and reset the
  // dimensions if necessary
  Teuchos::Array< int > commDims =
    plist.get("comm dimensions", Teuchos::Array< int >());
  if (numDims == 0) numDims = commDims.size();
  _commDims = regularizeCommDims(_teuchosComm->getSize(),
                                 numDims,
                                 commDims());

  // Set the periodic flags along each axis
  Teuchos::Array< int > periodic =
    plist.get("periodic", Teuchos::Array< int >());
  _periodic = createArrayOfInts(numDims, periodic);

  // Set the axis strides
  _commStrides = computeStrides<int,int>(_commDims,
                                         commLayout);

  // Set the axis ranks for this processor
  _commIndex = computeCommIndexes(_teuchosComm->getRank(),
                                  _commStrides);
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               Teuchos::ParameterList & plist) :
  _teuchosComm(teuchosComm)
{
  // Validate the ParameterList
  plist.validateParameters(*getValidParameters());

  // Determine a first cut at the number of dimensions
  Teuchos::Array< int > dims =
    plist.get("dimensions", Teuchos::Array< int >());
  int numDims = dims.size();

  // Set the communicator sizes along each axis and reset the
  // dimensions if necessary
  Teuchos::Array< int > commDims =
    plist.get("comm dimensions", Teuchos::Array< int >());
  if (numDims == 0) numDims = commDims.size();
  _commDims = regularizeCommDims(_teuchosComm->getSize(),
                                 numDims,
                                 commDims());

  // Set the periodic flags along each axis
  Teuchos::Array< int > periodic =
    plist.get("periodic", Teuchos::Array< int >());
  _periodic = createArrayOfInts(numDims, periodic);

  // Set the axis strides
  _commStrides = computeStrides<int,int>(_commDims,
                                         commLayout);

  // Set the axis ranks for this processor
  _commIndex = computeCommIndexes(_teuchosComm->getRank(),
                                  _commStrides);
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(int numDims) :
  _teuchosComm(Teuchos::DefaultComm< int >::getComm()),
  _commDims(regularizeCommDims(_teuchosComm->getSize(),
                               numDims,
                               Teuchos::Array< int >())),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(_teuchosComm->getRank(),
                                _commStrides)),
  _periodic(numDims, 0)
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               int numDims) :
  _teuchosComm(teuchosComm),
  _commDims(regularizeCommDims(teuchosComm->getSize(),
                               numDims,
                               Teuchos::Array< int >())),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(teuchosComm->getRank(),
                                _commStrides)),
  _periodic(numDims, 0)
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(int numDims,
               const Teuchos::ArrayView< const int > & commDims,
               const Teuchos::ArrayView< const int > & periodic) :
  _teuchosComm(Teuchos::DefaultComm< int >::getComm()),
  _commDims(regularizeCommDims(_teuchosComm->getSize(),
                               numDims,
                               commDims)),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(_teuchosComm->getRank(),
                                _commStrides)),
  _periodic(createArrayOfInts(numDims, periodic))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               int numDims,
               const Teuchos::ArrayView< const int > & commDims,
               const Teuchos::ArrayView< const int > & periodic) :
  _teuchosComm(teuchosComm),
  _commDims(regularizeCommDims(teuchosComm->getSize(),
                               numDims,
                               commDims)),
  _commStrides(computeStrides<int,int>(_commDims,
                                       commLayout)),
  _commIndex(computeCommIndexes(teuchosComm->getRank(),
                                _commStrides)),
  _periodic(createArrayOfInts(numDims, periodic))
{
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               int axis,
               int axisRank) :
  _commDims(),
  _commStrides(),
  _commIndex(),
  _periodic()
{
  if (parent.onSubcommunicator())
  {
    int numDims = parent.numDims();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for communicator with " <<
        numDims << " dimensions");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axisRank < 0) || (axisRank >= parent.getCommDim(axis))),
      RangeError,
      "axisRank = " << axisRank  << " is invalid for communicator axis " <<
        axis << " with " << parent.getCommDim(axis) << " processors");

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
        _commDims.push_back(1);
        _periodic.push_back(0);
        _commIndex.push_back(0);
        _commStrides.push_back(1);
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
          bounds[myAxis] = Slice(0, parent._commDims[myAxis]);

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
          val += mdIndex[myAxis] * parent._commStrides[myAxis];
        ranks[index] = val;
        // Increment mdIndex
        int myAxis = numDims-1;
        done = false;
        while (not done)
        {
          ++mdIndex[myAxis];
          done = (mdIndex[myAxis] < bounds[myAxis].stop());
          if (not done)
          {
            mdIndex[myAxis] = bounds[myAxis].start();
            --myAxis;
            if (myAxis < 0)
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
          _commDims.push_back(parent._commDims[myAxis]);
          _periodic.push_back(parent._periodic[myAxis]);
          _commIndex.push_back(parent._commIndex[myAxis]);
        }
      }
      _commStrides = computeStrides<int,int>(_commDims, commLayout);
    }
  }
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               int axis,
               const Slice & slice) :
  _commDims(parent._commDims),
  _commStrides(parent.numDims()),
  _commIndex(parent._commIndex),
  _periodic(parent._periodic)
{
  if (parent.onSubcommunicator())
  {
    // Sanity check
    size_type numDims = parent.numDims();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ((axis < 0) || (axis >= numDims)),
      RangeError,
      "axis = " << axis  << " is invalid for communicator with " <<
        numDims << " dimensions");

    // Make sure the Slice we work with is concrete, and adjust
    // _commDims, _commIndex, and _commStrides
    Slice bounds = slice.bounds(parent.getCommDim(axis));
    _commDims[axis] = (bounds.stop() - bounds.start());
    _commIndex[axis] -= bounds.start();

    // Fix the periodic flag
    if ((bounds.start() == 0) &&
        (bounds.stop() == parent.getCommDim(axis)))
      _periodic[axis] = parent._periodic[axis];
    else
      _periodic[axis] = 0;

    // Compute the new strides of the new sub-communicator
    _commStrides = computeStrides<int,int>(_commDims, commLayout);

    // Compute the ranks of the subcommunicator
    size_type rankSize = 1;
    for (int myAxis = 0; myAxis < numDims; ++myAxis)
      rankSize *= _commDims[myAxis];
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
        val += (mdIndex[myAxis] + start) * parent._commStrides[myAxis];
      }
      ranks[index] = val;
      // Increment mdIndex
      int myAxis = numDims-1;
      done = false;
      while (not done)
      {
        ++mdIndex[myAxis];
        done = (mdIndex[myAxis] < _commDims[myAxis]);
        if (not done)
        {
          mdIndex[myAxis] = 0;
          --myAxis;
          if (myAxis < 0)
            done = true;
        }
      }
    }

    // Set the communicator
    _teuchosComm =
      parent.getTeuchosComm()->createSubcommunicator(ranks()).getConst();
  }

  // On processors that are not a part of the subcommunicator, reset
  // the data attributes.
  if (_teuchosComm.getRawPtr() == 0)
  {
    _commDims.clear();
    _commStrides.clear();
    _commIndex.clear();
    _periodic.clear();
  }
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & parent,
               const Teuchos::ArrayView< Slice > & slices)
{
  // Sanity check
  size_type numDims = parent.numDims();
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
  // std::cout << parent.getTeuchosComm()->getRank() << ": _commDims = " << _commDims
  //           << ", _commStrides = " << _commStrides << ", _commIndex = "
  //           << _commIndex << ", _periodic = " << _periodic << std::endl;
}

////////////////////////////////////////////////////////////////////////

MDComm::MDComm(const MDComm & source) :
  _teuchosComm(source._teuchosComm),
#ifdef HAVE_EPETRA
  _epetraComm(source._epetraComm),
#endif
  _commDims(source._commDims),
  _commStrides(source._commStrides),
  _commIndex(source._commIndex),
  _periodic(source._periodic)
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
  _commDims    = source._commDims;
  _commStrides = source._commStrides;
  _commIndex   = source._commIndex;
  _periodic    = source._periodic;
  return *this;
}

////////////////////////////////////////////////////////////////////////

bool
MDComm::onSubcommunicator() const
{
  return _teuchosComm.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Teuchos::Comm< int > >
MDComm::getTeuchosComm() const
{
  return _teuchosComm;
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA
Teuchos::RCP< const Epetra_Comm >
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
MDComm::numDims() const
{
  return _commDims.size();
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
MDComm::getCommDims() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getCommDim()");
  return _commDims;
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getCommDim(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getCommDim()");
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _commDims[axis];
}

////////////////////////////////////////////////////////////////////////

bool
MDComm::isPeriodic(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::isPeriodic()");
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _periodic[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getCommIndex(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getCommIndex()");
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  return _commIndex[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getLowerNeighbor(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getLowerNeighbor()");
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (_commIndex[axis] == 0)
  {
    if (_periodic[axis])
    {
      return _teuchosComm->getRank() +
        (_commDims[axis]-1) * _commStrides[axis];
    }
    else
    {
      return -1;
    }
  }
  return _teuchosComm->getRank() - _commStrides[axis];
}

////////////////////////////////////////////////////////////////////////

int
MDComm::getUpperNeighbor(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    not onSubcommunicator(),
    SubcommunicatorError,
    "MDComm::getUpperNeighbor()");
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    ((axis < 0) || (axis >= numDims())),
    RangeError,
    "invalid axis index = " << axis << " (number of dimensions = " <<
    numDims() << ")");
#endif
  if (_commIndex[axis] == _commDims[axis] - 1)
  {
    if (_periodic[axis])
    {
      return _teuchosComm->getRank() -
        (_commDims[axis]-1) * _commStrides[axis];
    }
    else
    {
      return -1;
    }
  }
  return _teuchosComm->getRank() + _commStrides[axis];
}

////////////////////////////////////////////////////////////////////////

Teuchos::ArrayView< Teuchos::RCP< const MDComm > >
MDComm::getAxisComms() const
{
  int nd = numDims();
  if (_axisComms.size() < nd)
  {
    for (int axis = 0; axis < nd; ++axis)
    {
      // This call will resize _axisComms if necessary, and then store
      // the generated axis comm in _axisComms
      getAxisComm(axis);
    }
  }
  return _axisComms();
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const MDComm >
MDComm::getAxisComm(int axis) const
{
  int nd = numDims();
  if (_axisComms.size() < nd) _axisComms.resize(nd);
  if (_axisComms[axis].is_null())
  {
    _axisComms[axis] = Teuchos::rcp(new MDComm(*this));
    int myAxis = 0;
    for (int i = 0; i < nd; ++i)
    {
      if (i != axis)
      {
        int commIndex = getCommIndex(i);
        Teuchos::RCP< MDComm > newMdComm =
          Teuchos::rcp(new MDComm(*_axisComms[axis],
                                  myAxis,
                                  commIndex));
        _axisComms[axis] = newMdComm;
      }
      else
      {
        ++myAxis;
      }
    }
  }
  return _axisComms[axis];
}

}    // End namespace Domi
