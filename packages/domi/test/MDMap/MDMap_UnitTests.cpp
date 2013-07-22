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
string axisSizesStr = "-1";
Array< int > axisSizes;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("numDims"  , &numDims     , "number of dimensions");
  clp.setOption("axisSizes", &axisSizesStr, "comma-separated list of number "
                "of processors along each axis");
}

//
// Templated Unit Tests
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, dimensionsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisSizes));

  // Check that the axisSizes are completely specified
  TEST_EQUALITY(axisSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisSizes[axis] > 0);
  }

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisSize(axis);

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
    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisSizes[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis) , localDim  );
    Slice globalBounds = mdMap.getGlobalAxisBounds(axis);
    TEST_EQUALITY(globalBounds.start(), axisRank    *localDim);
    TEST_EQUALITY(globalBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdMap.getLocalAxisBounds(axis);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop() , localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), 0);
    TEST_ASSERT(not mdMap.getPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, halosConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisSizes));

  // Check that the axisSizes are completely specified
  TEST_EQUALITY(axisSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisSizes[axis] > 0);
  }

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisSize(axis);

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
    if (axisRank < axisSizes[axis]-1) upperHalo = halos[axis];
    T myDim = localDim + lowerHalo + upperHalo;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisSizes[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalBounds = mdMap.getGlobalAxisBounds(axis);
    TEST_EQUALITY(globalBounds.start(),  axisRank   *localDim);
    TEST_EQUALITY(globalBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdMap.getLocalAxisBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop() , myDim);
    localBounds        = mdMap.getLocalAxisBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerHalo);
    TEST_EQUALITY(localBounds.stop(), lowerHalo+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), lowerHalo);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), upperHalo);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), halos[axis]);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), 0);
    TEST_ASSERT(not mdMap.getPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, ghostsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisSizes));

  // Check that the axisSizes are completely specified
  TEST_EQUALITY(axisSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisSizes[axis] > 0);
  }

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisSize(axis);

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
    if (axisRank == axisSizes[axis]-1) upperGhost = ghosts[axis];
    T myDim = localDim + lowerGhost + upperGhost;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisSizes[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalBounds = mdMap.getGlobalAxisBounds(axis,false);
    T myStart = ghosts[axis] +  axisRank    * localDim;
    T myStop  = ghosts[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalBounds.start(), myStart);
    TEST_EQUALITY(globalBounds.stop() , myStop );
    globalBounds = mdMap.getGlobalAxisBounds(axis,true);
    if (axisRank == 0                ) myStart -= ghosts[axis];
    if (axisRank == axisSizes[axis]-1) myStop  += ghosts[axis];
    TEST_EQUALITY(globalBounds.start(), myStart);
    TEST_EQUALITY(globalBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalAxisBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds        = mdMap.getLocalAxisBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerGhost);
    TEST_EQUALITY(localBounds.stop(), lowerGhost+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), lowerGhost);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), upperGhost);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), ghosts[axis]);
    TEST_ASSERT(not mdMap.getPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, halosAndGhostsConstructor, T )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisSizes));

  // Check that the axisSizes are completely specified
  TEST_EQUALITY(axisSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisSizes[axis] > 0);
  }

  // Construct dimensions
  T localDim = 10;
  Array< T > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisSize(axis);

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
    if (axisRank == axisSizes[axis]-1) upperHalo = ghosts[axis];
    T myDim = localDim + lowerHalo + upperHalo;

    TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisSizes[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*ghosts[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalBounds = mdMap.getGlobalAxisBounds(axis,false);
    T myStart = ghosts[axis] +  axisRank    * localDim;
    T myStop  = ghosts[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalBounds.start(), myStart);
    TEST_EQUALITY(globalBounds.stop() , myStop );
    globalBounds = mdMap.getGlobalAxisBounds(axis,true);
    if (axisRank == 0                ) myStart -= ghosts[axis];
    if (axisRank == axisSizes[axis]-1) myStop  += ghosts[axis];
    TEST_EQUALITY(globalBounds.start(), myStart);
    TEST_EQUALITY(globalBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalAxisBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds        = mdMap.getLocalAxisBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerHalo);
    TEST_EQUALITY(localBounds.stop(), lowerHalo+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerHalo(axis), lowerHalo);
    TEST_EQUALITY_CONST(mdMap.getUpperHalo(axis), upperHalo);
    TEST_EQUALITY_CONST(mdMap.getHaloSize(axis), halos[axis]);
    TEST_EQUALITY_CONST(mdMap.getGhostSize(axis), ghosts[axis]);
    TEST_ASSERT(not mdMap.getPeriodic(axis));
  }
}

//
// Instantiations
//

#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, dimensionsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, halosConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, ghostsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, halosAndGhostsConstructor, T )

UNIT_TEST_GROUP(int)
#if 0
UNIT_TEST_GROUP(long)
#endif

}  // namespace
