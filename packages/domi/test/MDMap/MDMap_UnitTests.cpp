/*
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

// Teuchos includes
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Tuple.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDMap.hpp"

namespace
{

using std::string;
using Teuchos::Array;
using Teuchos::Tuple;
using Domi::TeuchosCommRCP;
using Domi::Slice;
using Domi::MDComm;
using Domi::MDCommRCP;
using Domi::MDMap;

int numDims = 2;
string axisCommSizesStr = "-1";
Array< int > axisCommSizes;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("numDims"  , &numDims     , "number of dimensions");
  clp.setOption("axisCommSizes", &axisCommSizesStr, "comma-separated list of "
                "number of processors along each axis");
}

//
// Templated Unit Tests
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, dimensionsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Check that the axisCommSizes are completely specified
  TEST_EQUALITY(axisCommSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisCommSizes[axis] > 0);
  }

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.getNumDims(), numDims);
  TEST_ASSERT(not mdMap.hasHalos());
  TEST_EQUALITY(mdMap.getStorageOrder(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis).stop(), dims[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis) , localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis);
    TEST_EQUALITY(globalRankBounds.start(), axisRank    *localDim);
    TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdMap.getLocalBounds(axis);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getLowerGhost(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperGhost(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), 0);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, halosConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct halos
  Array< int > halos(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    halos[axis] = axis+1;

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims(), halos());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.getNumDims(), numDims);
  TEST_EQUALITY(mdMap.hasHalos(), comm->getSize() > 1);
  TEST_EQUALITY(mdMap.getStorageOrder(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    int lowerHalo = 0;
    if (axisRank > 0) lowerHalo = halos[axis];
    int upperHalo = 0;
    if (axisRank < axisCommSizes[axis]-1) upperHalo = halos[axis];
    T myDim = localDim + lowerHalo + upperHalo;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis).stop(), dims[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis);
    TEST_EQUALITY(globalRankBounds.start(),  axisRank   *localDim);
    TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerHalo);
    TEST_EQUALITY(localBounds.stop(), lowerHalo+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), lowerHalo);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), upperHalo);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), halos[axis]);
    TEST_EQUALITY_CONST(mdMap.getLowerGhost(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperGhost(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), 0);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, ghostsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct halos
  Array< int > halos;

  // Construct ghosts
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    ghosts[axis] = axis+1;

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims(), halos(), ghosts());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.getNumDims(), numDims);
  TEST_ASSERT(mdMap.hasHalos());
  TEST_EQUALITY(mdMap.getStorageOrder(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    int lowerGhost = 0;
    if (axisRank == 0) lowerGhost = ghosts[axis];
    int upperGhost = 0;
    if (axisRank == axisCommSizes[axis]-1) upperGhost = ghosts[axis];
    T myDim = localDim + lowerGhost + upperGhost;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  ghosts[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    T myStart = ghosts[axis] +  axisRank    * localDim;
    T myStop  = ghosts[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= ghosts[axis];
    if (axisRank == axisCommSizes[axis]-1) myStop  += ghosts[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerGhost);
    TEST_EQUALITY(localBounds.stop(), lowerGhost+localDim);
    TEST_EQUALITY(mdMap.getLowerHalo(axis), lowerGhost);
    TEST_EQUALITY(mdMap.getUpperHalo(axis), upperGhost);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), 0);
    TEST_EQUALITY(mdMap.getLowerGhost(axis), ghosts[axis]);
    TEST_EQUALITY(mdMap.getUpperGhost(axis), ghosts[axis]);
    TEST_EQUALITY(mdMap.getGhostSize(axis), ghosts[axis]);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, halosAndGhostsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct halos and ghosts
  Array< int > halos(numDims);
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    halos[axis]  = axis+1;
    ghosts[axis] = axis+2;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims(), halos(), ghosts());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.getNumDims(), numDims);
  TEST_ASSERT(mdMap.hasHalos());
  TEST_EQUALITY(mdMap.getStorageOrder(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank  = mdMap.getAxisRank(axis);
    int lowerHalo = halos[axis];
    int upperHalo = halos[axis];
    if (axisRank == 0                ) lowerHalo = ghosts[axis];
    if (axisRank == axisCommSizes[axis]-1) upperHalo = ghosts[axis];
    T myDim = localDim + lowerHalo + upperHalo;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  ghosts[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    T myStart = ghosts[axis] +  axisRank    * localDim;
    T myStop  = ghosts[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= ghosts[axis];
    if (axisRank == axisCommSizes[axis]-1) myStop  += ghosts[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerHalo);
    TEST_EQUALITY(localBounds.stop(), lowerHalo+localDim);
    TEST_EQUALITY(mdMap.getLowerHalo(axis), lowerHalo);
    TEST_EQUALITY(mdMap.getUpperHalo(axis), upperHalo);
    TEST_EQUALITY(mdMap.getHaloSize(axis), halos[axis]);
    TEST_EQUALITY(mdMap.getLowerGhost(axis), ghosts[axis]);
    TEST_EQUALITY(mdMap.getUpperGhost(axis), ghosts[axis]);
    TEST_EQUALITY(mdMap.getGhostSize(axis), ghosts[axis]);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, periodic, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  MDCommRCP mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes, periodic));

  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm->isPeriodic(axis), (axis == 0));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, indexes, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct halos and ghosts
  Array< int > halos(numDims);
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    halos[axis]  = axis+1;
    ghosts[axis] = axis+2;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims(), halos(), ghosts());

  // Compute some quantities for testing
  Array< int > axisRanks(numDims);
  Array< int > lowerHalos(numDims);
  Array< int > upperHalos(numDims);
  Array< T >   myLocalDims(numDims);
  Array< T >   globalStrides(numDims, 1);
  Array< T >   localStrides(numDims, 1);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisRanks[axis]   = mdMap.getAxisRank(axis);
    lowerHalos[axis]  = mdMap.getLowerHalo(axis);
    upperHalos[axis]  = mdMap.getUpperHalo(axis);
    myLocalDims[axis] = localDim + lowerHalos[axis] + upperHalos[axis];
    if (axis > 0)
    {
      globalStrides[axis] = globalStrides[axis-1] *
        (dims[axis-1] + 2*ghosts[axis-1]);
      localStrides[axis] = localStrides[axis-1] * myLocalDims[axis-1];
    }
  }

  // Unit test globalAxisIndex <-> globalIndex
  Array< T > myGlobalAxisIndex(numDims);
  T myGlobalIndex = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    myGlobalAxisIndex[axis] = ghosts[axis] + dims[axis] / 2;
    myGlobalIndex += myGlobalAxisIndex[axis] * globalStrides[axis];
  }
  TEST_EQUALITY(mdMap.getGlobalAxisIndex(myGlobalIndex), myGlobalAxisIndex);
  TEST_EQUALITY(mdMap.getGlobalIndex(myGlobalAxisIndex()), myGlobalIndex);

  // Unit test localAxisIndex <-> localIndex
  Array< T > myLocalAxisIndex(numDims);
  T myLocalIndex = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    myLocalAxisIndex[axis] = lowerHalos[axis] + localDim / 2;
    myLocalIndex += myLocalAxisIndex[axis] * localStrides[axis];
  }
  TEST_EQUALITY(mdMap.getLocalAxisIndex(myLocalIndex), myLocalAxisIndex);
  TEST_EQUALITY(mdMap.getLocalIndex(myLocalAxisIndex()), myLocalIndex);

  // Test localIndex <-> globalIndex
  myGlobalIndex = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    myGlobalAxisIndex[axis] = myLocalAxisIndex[axis] - lowerHalos[axis] +
      axisRanks[axis] * localDim + ghosts[axis];
    myGlobalIndex += myGlobalAxisIndex[axis] * globalStrides[axis];
  }
  TEST_EQUALITY(mdMap.getGlobalIndex(myLocalIndex), myGlobalIndex);
  TEST_EQUALITY(mdMap.getLocalIndex(myGlobalIndex), myLocalIndex );
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, exceptions, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dims());

  // Unit test methods that should throw exceptions
#if DOMI_ENABLE_ABC
  TEST_THROW(mdMap.getAxisCommSize(    -1), Domi::RangeError);
  TEST_THROW(mdMap.isPeriodic(         -1), Domi::RangeError);
  TEST_THROW(mdMap.getAxisRank(        -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerNeighbor(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperNeighbor(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalDim(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalBounds(    -1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalDim(        -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalRankBounds(-1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalBounds(     -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerHalo(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperHalo(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getHaloSize(        -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerGhost(      -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperGhost(      -1), Domi::RangeError);
  TEST_THROW(mdMap.getGhostSize(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getAxisCommSize(    numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.isPeriodic(         numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getAxisRank(        numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerNeighbor(   numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperNeighbor(   numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalDim(       numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalBounds(    numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalDim(        numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalRankBounds(numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalBounds(     numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerHalo(       numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperHalo(       numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getHaloSize(        numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerGhost(      numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperGhost(      numDims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGhostSize(       numDims+1), Domi::RangeError);
#endif
}

//
// Instantiations
//

#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, dimensionsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, halosConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, ghostsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, halosAndGhostsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, indexes, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, exceptions, T )

UNIT_TEST_GROUP(int)
#if 1
UNIT_TEST_GROUP(long)
#endif

}  // namespace
