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
#include "Domi_MDMap.hpp"

namespace Domi
{

MDMap::
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

MDMap::
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

MDMap::
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

MDMap::
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

MDMap::
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

MDMap::
MDMap(const MDMap & source) :
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

MDMap::
MDMap(const MDMap & parent,
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

MDMap::
MDMap(const MDMap & parent,
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

MDMap::
MDMap(const MDMap & parent,
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
  MDMap tempMDMap1(parent);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int bndryPadding = (axis < bndryPad.size()) ? bndryPad[axis] : 0;
    MDMap tempMDMap2(tempMDMap1,
                             axis,
                             slices[axis],
                             bndryPadding);
    tempMDMap1 = tempMDMap2;
  }
  *this = tempMDMap1;
}

////////////////////////////////////////////////////////////////////////

MDMap::~MDMap()
{
}

////////////////////////////////////////////////////////////////////////

MDMap &
MDMap::operator=(const MDMap & source)
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

dim_type
MDMap::
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

Slice
MDMap::
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

dim_type
MDMap::
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

Teuchos::Array< dim_type >
MDMap::
getLocalDims() const
{
  return _localDims;
}

////////////////////////////////////////////////////////////////////////

Slice
MDMap::
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

Teuchos::ArrayView< const Slice >
MDMap::
getLocalBounds() const
{
  return _localBounds();
}

////////////////////////////////////////////////////////////////////////

Slice
MDMap::
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

Slice
MDMap::
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

bool
MDMap::hasPadding() const
{
  bool result = false;
  for (int axis = 0; axis < numDims(); ++axis)
    if (_pad[axis][0] + _pad[axis][1]) result = true;
  return result;
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getLowerPadSize(int axis) const
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

int
MDMap::getUpperPadSize(int axis) const
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

int
MDMap::getCommPadSize(int axis) const
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

int
MDMap::getLowerBndryPad(int axis) const
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

int
MDMap::getUpperBndryPad(int axis) const
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

Teuchos::Array< int >
MDMap::getBndryPadSizes() const
{
  return _bndryPadSizes;
}

////////////////////////////////////////////////////////////////////////

int
MDMap::getBndryPadSize(int axis) const
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

bool
MDMap::
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

bool
MDMap::
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

bool
MDMap::
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

bool
MDMap::isReplicatedBoundary(int axis) const
{
  return _mdComm->isPeriodic(axis) && bool(_replicatedBoundary[axis]);
}

////////////////////////////////////////////////////////////////////////

Layout
MDMap::getLayout() const
{
  return _layout;
}

////////////////////////////////////////////////////////////////////////

Teuchos::ArrayView< Teuchos::RCP< const MDMap > >
MDMap::getAxisMaps() const
{
  if (_axisMaps.size() == 0) _axisMaps.resize(numDims());
  for (int axis = 0; axis < numDims(); ++axis)
    if (_axisMaps[axis].is_null()) getAxisMap(axis);
  return _axisMaps();
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const MDMap >
MDMap::getAxisMap(int axis) const
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

Teuchos::RCP< const MDMap >
MDMap::getAugmentedMDMap(const dim_type leadingDim,
                                 const dim_type trailingDim) const
{
  // Construct the new MDMap
  MDMap * newMdMap = new MDMap(*this);

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

Teuchos::RCP< const Epetra_Map >
MDMap::getEpetraMap(bool withCommPad) const
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

Teuchos::RCP< const Epetra_Map >
MDMap::
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

Teuchos::Array< dim_type >
MDMap::
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

Teuchos::Array< dim_type >
MDMap::
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

size_type
MDMap::
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

size_type
MDMap::
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

size_type
MDMap::
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

size_type
MDMap::
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

bool
MDMap::isCompatible(const MDMap & mdMap) const
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

bool
MDMap::isSameAs(const MDMap & mdMap,
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

bool
MDMap::isContiguous() const
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

void
MDMap::computeBounds()
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
