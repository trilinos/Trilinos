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
using Teuchos::ArrayView;
using Teuchos::Tuple;
typedef Domi::Ordinal Ordinal;
using Domi::TeuchosCommRCP;
using Domi::Slice;
const Ordinal & Default = Domi::Slice::Default;
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

#ifdef HAVE_EPETRA
  if (typeid(T) == typeid(int))
  {
    Teuchos::Array<int> gstrides(numDims);
    Teuchos::Array<int> lstrides(numDims);
    gstrides[0] = 1;
    lstrides[0] = 1;
    for (int axis=1; axis < numDims; ++axis)
    {
      gstrides[axis] = dims[axis-1] * gstrides[axis-1];
      lstrides[axis] = localDim     * lstrides[axis-1];
    }
    int gll = 0;
    int gur = 0;
    int lur = 0;
    for (int axis=0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      gll += ( axisRank   *localDim  ) * gstrides[axis];
      gur += ((axisRank+1)*localDim-1) * gstrides[axis];
      lur += (localDim-1) * lstrides[axis];
    }
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
    TEST_EQUALITY(epetraMap->GID(  0), gll);
    TEST_EQUALITY(epetraMap->GID(lur), gur);

    for (int axis = 0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
      TEST_EQUALITY(epetraMap->GID(         0),  axisRank   *localDim  );
      TEST_EQUALITY(epetraMap->GID(localDim-1), (axisRank+1)*localDim-1);
    }
  }
  else
  {
    TEST_THROW(mdMap.getEpetraMap(true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraMap(false), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,false), Domi::MapOrdinalError);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::Array<T> gstrides(numDims);
  Teuchos::Array<T> lstrides(numDims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < numDims; ++axis)
  {
    gstrides[axis] = dims[axis-1] * gstrides[axis-1];
    lstrides[axis] = localDim     * lstrides[axis-1];
  }
  T gll = 0;
  T gur = 0;
  T lll = 0;
  T lur = 0;
  for (int axis=0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    gll += ( axisRank   *localDim  ) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    lur += (localDim-1) * lstrides[axis];
  }
  Teuchos::RCP< const Tpetra::Map<T, T> > tpetraMap =
    mdMap.getTpetraMap();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    Teuchos::RCP< const Tpetra::Map<T,T> > tpetraMap =
      mdMap.getTpetraAxisMap(axis);
    TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim);
    TEST_EQUALITY(tpetraMap->getGlobalElement(localDim-1),
                  (axisRank+1)*localDim-1);
  }
#endif
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

#ifdef HAVE_EPETRA
  if (typeid(T) == typeid(int))
  {
    Teuchos::Array<int> gstrides(numDims);
    Teuchos::Array<int> lstrides(numDims);
    gstrides[0] = 1;
    lstrides[0] = 1;
    for (int axis=1; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis-1);
      int myHalos = 0;
      if (axisRank > 0) myHalos += halos[axis-1];
      if (axisRank != mdMap.getAxisCommSize(axis-1)-1) myHalos += halos[axis-1];
      gstrides[axis] = dims[axis-1]       * gstrides[axis-1];
      lstrides[axis] = (localDim+myHalos) * lstrides[axis-1];
    }
    int gll = 0;
    int gur = 0;
    int lll = 0;
    int lur = 0;
    for (int axis=0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      gll += ( axisRank   *localDim) * gstrides[axis];
      gur += ((axisRank+1)*localDim-1) * gstrides[axis];
      if (axisRank == 0)
        lur += (localDim-1) * lstrides[axis];
      else
      {
        lll += halos[axis] * lstrides[axis];
        lur += (localDim+halos[axis]-1) * lstrides[axis];
      }
    }
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
    TEST_EQUALITY(epetraMap->GID(lll), gll);
    TEST_EQUALITY(epetraMap->GID(lur), gur);

    for (int axis = 0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      int myHalo   = 0;
      if (axisRank > 0) myHalo = halos[axis];
      Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
      TEST_EQUALITY(epetraMap->GID(0), axisRank*localDim-myHalo);
      TEST_EQUALITY(epetraMap->GID(localDim+myHalo-1), (axisRank+1)*localDim-1);
    }
  }
  else
  {
    TEST_THROW(mdMap.getEpetraMap(true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraMap(false), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,false), Domi::MapOrdinalError);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::Array<T> gstrides(numDims);
  Teuchos::Array<T> lstrides(numDims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis-1);
    int myHalos = 0;
    if (axisRank > 0) myHalos += halos[axis-1];
    if (axisRank != mdMap.getAxisCommSize(axis-1)-1) myHalos += halos[axis-1];
    gstrides[axis] = dims[axis-1]       * gstrides[axis-1];
    lstrides[axis] = (localDim+myHalos) * lstrides[axis-1];
  }
  T gll = 0;
  T gur = 0;
  T lll = 0;
  T lur = 0;
  for (int axis=0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    gll += ( axisRank   *localDim) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    if (axisRank == 0)
      lur += (localDim-1) * lstrides[axis];
    else
    {
      lll += halos[axis] * lstrides[axis];
      lur += (localDim+halos[axis]-1) * lstrides[axis];
    }
  }
  Teuchos::RCP< const Tpetra::Map<T,T> > tpetraMap = mdMap.getTpetraMap();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    int myHalo   = 0;
    if (axisRank > 0) myHalo = halos[axis];
    Teuchos::RCP< const Tpetra::Map<T> > tpetraMap =
      mdMap.getTpetraAxisMap(axis);
    TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim-myHalo);
    TEST_EQUALITY(tpetraMap->getGlobalElement(localDim+myHalo-1),
                  (axisRank+1)*localDim-1);
  }
#endif
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

#ifdef HAVE_EPETRA
  if (typeid(T) == typeid(int))
  {
    Teuchos::Array<int> gstrides(numDims);
    Teuchos::Array<int> lstrides(numDims);
    gstrides[0] = 1;
    lstrides[0] = 1;
    for (int axis=1; axis < numDims; ++axis)
    {
      gstrides[axis] = (dims[axis-1]+2*ghosts[axis-1]) * gstrides[axis-1];
      int myGhosts = 0;
      int axisRank = mdMap.getAxisRank(axis-1);
      if (axisRank == 0) myGhosts += ghosts[axis-1];
      if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
        myGhosts += ghosts[axis-1];
      lstrides[axis] = (localDim+myGhosts) * lstrides[axis-1];
    }
    int gll = 0;
    int gur = 0;
    int lll = 0;
    int lur = 0;
    for (int axis=0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      gll += (axisRank*localDim+ghosts[axis]) * gstrides[axis];
      gur += ((axisRank+1)*localDim+ghosts[axis]-1) * gstrides[axis];
      if (axisRank == 0)
      {
        lll += ghosts[axis] * lstrides[axis];
        lur += (localDim+ghosts[axis]-1) * lstrides[axis];
      }
      else
        lur += (localDim-1) * lstrides[axis];
    }
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
    TEST_EQUALITY(epetraMap->GID(lll), gll);
    TEST_EQUALITY(epetraMap->GID(lur), gur);

    for (int axis = 0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      int myHalo   = 0;
      if (axisRank == 0) myHalo = ghosts[axis];
      Teuchos::RCP< const Epetra_Map > epetraMap =
        mdMap.getEpetraAxisMap(axis);
      TEST_EQUALITY(epetraMap->GID(myHalo), axisRank*localDim+ghosts[axis]);
      TEST_EQUALITY(epetraMap->GID(myHalo+localDim-1),
                    ghosts[axis]+(axisRank+1)*localDim-1);
    }
  }
  else
  {
    TEST_THROW(mdMap.getEpetraMap(true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraMap(false), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraAxisMap(0,false), Domi::MapOrdinalError);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::Array<T> gstrides(numDims);
  Teuchos::Array<T> lstrides(numDims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < numDims; ++axis)
  {
    gstrides[axis] = (dims[axis-1]+2*ghosts[axis-1]) * gstrides[axis-1];
    int myGhosts = 0;
    int axisRank = mdMap.getAxisRank(axis-1);
    if (axisRank == 0) myGhosts += ghosts[axis-1];
    if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
      myGhosts += ghosts[axis-1];
    lstrides[axis] = (localDim+myGhosts) * lstrides[axis-1];
  }
  T gll = 0;
  T gur = 0;
  T lll = 0;
  T lur = 0;
  for (int axis=0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    gll += (axisRank*localDim+ghosts[axis]) * gstrides[axis];
    gur += ((axisRank+1)*localDim+ghosts[axis]-1) * gstrides[axis];
    if (axisRank == 0)
    {
      lll += ghosts[axis] * lstrides[axis];
      lur += (localDim+ghosts[axis]-1) * lstrides[axis];
    }
    else
      lur += (localDim-1) * lstrides[axis];
  }
  Teuchos::RCP< const Tpetra::Map<T,T> > tpetraMap = mdMap.getTpetraMap();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdMap.getAxisRank(axis);
    int myHalo   = 0;
    if (axisRank == 0) myHalo = ghosts[axis];
    Teuchos::RCP< const Tpetra::Map<T> > tpetraMap =
      mdMap.getTpetraAxisMap(axis);
    TEST_EQUALITY(tpetraMap->getGlobalElement(myHalo),
                  axisRank*localDim+ghosts[axis]);
    TEST_EQUALITY(tpetraMap->getGlobalElement(myHalo+localDim-1),
                  ghosts[axis]+(axisRank+1)*localDim-1);
  }
#endif
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

#ifdef HAVE_EPETRA
  if (typeid(T) == typeid(int))
  {
    Teuchos::Array<int> gstrides(numDims);
    Teuchos::Array<int> lstrides(numDims);
    gstrides[0] = 1;
    lstrides[0] = 1;
    for (int axis=1; axis < numDims; ++axis)
    {
      gstrides[axis] = (dims[axis-1]+2*ghosts[axis-1]) * gstrides[axis-1];
      int myHalos = 0;
      int axisRank = mdMap.getAxisRank(axis-1);
      if (axisRank == 0)
        myHalos += ghosts[axis-1];
      else
        myHalos += halos[axis-1];
      if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
        myHalos += ghosts[axis-1];
      else
        myHalos += halos[axis-1];
      lstrides[axis] = (localDim+myHalos) * lstrides[axis-1];
    }
    int gll = 0;
    int gur = 0;
    int lll = 0;
    int lur = 0;
    for (int axis=0; axis < numDims; ++axis)
    {
      int axisRank = mdMap.getAxisRank(axis);
      gll += (axisRank*localDim+ghosts[axis]) * gstrides[axis];
      gur += ((axisRank+1)*localDim+ghosts[axis]-1) * gstrides[axis];
      if (axisRank == 0)
      {
        lll += ghosts[axis] * lstrides[axis];
        lur += (localDim+ghosts[axis]-1) * lstrides[axis];
      }
      else
      {
        lll += halos[axis] * lstrides[axis];
        lur += (localDim+halos[axis]-1) * lstrides[axis];
      }
    }
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
    TEST_EQUALITY(epetraMap->GID(lll), gll);
    TEST_EQUALITY(epetraMap->GID(lur), gur);
  }
  else
  {
    TEST_THROW(mdMap.getEpetraMap(true ), Domi::MapOrdinalError);
    TEST_THROW(mdMap.getEpetraMap(false), Domi::MapOrdinalError);
  }
#endif
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapLowerLeft, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions
  Array< T > dimensions(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dimensions[axis] = 10 * axisCommSizes[axis];

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions);
  
  // Figure out the lower left slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      T nd   = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdMap.getGlobalBounds(axis).start() >= newDims[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
      TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapLowerLeftWithHalo, T )
{
  // Construct the MDComm from command-line arguments and halos
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions and halos
  Array< T >   dimensions(numDims);
  Array< int > halos(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = 10 * axisCommSizes[axis];
    halos[axis]      = axis + 1;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions, halos);
  
  // Figure out the lower left slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      T nd   = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdMap.getGlobalBounds(axis).start() >= newDims[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
      TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (axisCommSizes[axis] > 1)
      {
        if (subMDMap.getAxisRank(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerHalo(axis), halos[axis]);
        }
        if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperHalo(axis), halos[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), 0);
        TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), 0);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapLowerRight, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions
  Array< T > dimensions(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dimensions[axis] = 10 * axisCommSizes[axis];

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions);
  
  // Figure out the lower right slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd,Default);
    }
    else if (axis == 1)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) >= newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis == 0)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapLowerRightWithGhost, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions and ghosts
  Array< T >   dimensions(numDims);
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = 10 * axisCommSizes[axis];
    ghosts[axis]     = axis + 2;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions, ArrayView<int>(), ghosts);
  
  // Figure out the lower right slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd+ghosts[axis],-ghosts[axis]);
    }
    else if (axis == 1)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(ghosts[axis],nd+ghosts[axis]);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice(ghosts[axis], -ghosts[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) >= newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis == 0)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + ghosts[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + ghosts[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapUpperLeft, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions
  Array< T > dimensions(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dimensions[axis] = 10 * axisCommSizes[axis];

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions);
  
  // Figure out the upper left slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd);
    }
    else if (axis == 1)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd,Default);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) >= newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis == 1)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapUpperLeftHaloGhost, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions, halos and ghosts
  Array< T > dimensions(numDims);
  Array< int > halos(numDims);
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = 10 * axisCommSizes[axis];
    halos[axis]      = axis + 1;
    ghosts[axis]     = axis + 2;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions, halos, ghosts);

  // Figure out the upper left slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      if (nc == 0) nc = 1;
      newCommSizes[axis] = nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(ghosts[axis],nd+ghosts[axis]);
    }
    else if (axis == 1)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd+ghosts[axis],-ghosts[axis]);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice(ghosts[axis], -ghosts[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) >= newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis == 1)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + ghosts[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + ghosts[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (axisCommSizes[axis] > 1)
      {
        if (subMDMap.getAxisRank(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerHalo(axis), halos[axis]);
        }
        if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperHalo(axis), halos[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), 0);
        TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), 0);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapUpperRight, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions
  Array< T > dimensions(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dimensions[axis] = 10 * axisCommSizes[axis];

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions);
  
  // Figure out the upper right slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - nc;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd,Default);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis < 2)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.getTeuchosComm().getRawPtr(), 0);
    TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapUpperRightNewGhosts, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the dimensions, halos and ghosts
  int mySize = 10;
  Array< T > dimensions(numDims);
  Array< int > halos(numDims);
  Array< int > ghosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = mySize * axisCommSizes[axis];
    halos[axis]      = axis + 1;
    ghosts[axis]     = axis + 2;
  }

  // Construct an MDMap
  MDMap< T > mdMap(mdComm, dimensions, halos, ghosts);

  // Compute the new ghost sizes and figure out the upper right slice
  Array< Slice > slices(numDims);
  Array< int >   newCommSizes(numDims);
  Array< T >     newDims(numDims);
  Array< int >   newGhosts(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Make the new ghosts one less than the old ghosts
    newGhosts[axis] = ghosts[axis] - 1;
    if (axis < 2)
    {
      T   nd = dimensions[axis] / 2;
      int nc = axisCommSizes[axis] / 2;
      newCommSizes[axis] = axisCommSizes[axis] - (nd-newGhosts[axis])/mySize;
      newDims[axis]      = nd;
      slices[axis]       = Slice(nd+ghosts[axis],-ghosts[axis]);
    }
    else
    {
      newCommSizes[axis] = axisCommSizes[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice(ghosts[axis], -ghosts[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap< T > subMDMap(mdMap, slices, newGhosts);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      if (axis < 2)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + ghosts[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      ghosts[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + ghosts[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (axisCommSizes[axis] > 1)
      {
        if (subMDMap.getAxisRank(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), newGhosts[axis]);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerHalo(axis), halos[axis]);
        }
        if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), newGhosts[axis]);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperHalo(axis), halos[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerHalo(axis), newGhosts[axis]);
        TEST_EQUALITY_CONST(subMDMap.getUpperHalo(axis), newGhosts[axis]);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
    TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapReduce, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct the parent MDMap
  Array< T > dimensions(numDims);
  Array<int> halos(numDims);
  Array<int> ghosts(numDims);
  int localSize = 10;
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = localSize*axisCommSizes[axis];
    halos[axis] = axis + 1;
    ghosts[axis] = axis + 2;
  }
  MDMap<T> mdMap(mdComm, dimensions, halos, ghosts);

  // We will reduce this parent MDMap several times by using the
  // single GlobalOrd constructor along each dimension
  for (int axis = 0; axis < numDims; ++axis)
  {
    T myOrd = ghosts[axis] + dimensions[axis] / 2;
    int newDims = (numDims > 1) ? numDims-1 : numDims;
    Slice bounds = mdMap.getGlobalRankBounds(axis);
    bool partOfSubMap = true;
    if ((myOrd < bounds.start()) || (bounds.stop() <= myOrd))
      partOfSubMap = false;
    MDMap<T> reducedMdMap(mdMap, axis, myOrd);
    if (partOfSubMap)
    {
      TEST_ASSERT(reducedMdMap.onSubcommunicator());
      TEST_EQUALITY(reducedMdMap.getNumDims(), newDims);
      for (int newAxis = 0; newAxis < newDims; ++newAxis)
      {
        int axisRank = reducedMdMap.getAxisRank(newAxis);
        if (numDims == 1)
        {
          TEST_EQUALITY_CONST(reducedMdMap.getGlobalDim(newAxis), 1);
          Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), myOrd  );
          TEST_EQUALITY(bounds.stop() , myOrd+1);
          bounds = reducedMdMap.getGlobalRankBounds(newAxis);
          TEST_EQUALITY(bounds.start(), myOrd  );
          TEST_EQUALITY(bounds.stop() , myOrd+1);
          TEST_EQUALITY_CONST(reducedMdMap.getLocalDim(newAxis), 1);
          bounds = reducedMdMap.getLocalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), 0);
          TEST_EQUALITY(bounds.stop(),  1);
          TEST_ASSERT(not reducedMdMap.hasHalos());
          TEST_EQUALITY_CONST(reducedMdMap.getLowerHalo(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getUpperHalo(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getHaloSize(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getLowerGhost(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getUpperGhost(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getGhostSize(newAxis), 0);
        }
        else if (newAxis < axis)
        {
          TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
                        dimensions[newAxis]);
          Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), ghosts[newAxis]                    );
          TEST_EQUALITY(bounds.stop() , ghosts[newAxis]+dimensions[newAxis]);
          bounds = reducedMdMap.getGlobalRankBounds(newAxis);
          TEST_EQUALITY(bounds.start(), ghosts[newAxis]+localSize* axisRank   );
          TEST_EQUALITY(bounds.stop() , ghosts[newAxis]+localSize*(axisRank+1));
          TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
          bounds = reducedMdMap.getLocalBounds(newAxis);
          TEST_ASSERT(reducedMdMap.hasHalos());
          if (reducedMdMap.getAxisRank(newAxis) == 0)
          {
            TEST_EQUALITY(bounds.start(), ghosts[newAxis]          );
            TEST_EQUALITY(bounds.stop() , ghosts[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis), ghosts[newAxis]);
          }
          else if (reducedMdMap.getAxisRank(newAxis) ==
                   reducedMdMap.getAxisCommSize(newAxis)-1)
          {
            TEST_EQUALITY(bounds.start(), halos[newAxis]          );
            TEST_EQUALITY(bounds.stop() , halos[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis), halos[newAxis]);
            TEST_EQUALITY(reducedMdMap.getUpperHalo(newAxis), ghosts[newAxis]);
          }
          else
          {
            TEST_EQUALITY(bounds.start(), halos[newAxis]          );
            TEST_EQUALITY(bounds.stop() , halos[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis), halos[newAxis]);
            TEST_EQUALITY(reducedMdMap.getUpperHalo(newAxis), halos[newAxis]);
          }
          TEST_EQUALITY(reducedMdMap.getLowerGhost(newAxis), ghosts[newAxis]);
          TEST_EQUALITY(reducedMdMap.getUpperGhost(newAxis), ghosts[newAxis]);
          TEST_EQUALITY(reducedMdMap.getGhostSize(newAxis) , ghosts[newAxis]);
        }
        else
        {
          TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
                        dimensions[newAxis+1]);
          Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), ghosts[newAxis+1]);
          TEST_EQUALITY(bounds.stop() , ghosts[newAxis+1]+
                        dimensions[newAxis+1]);
          bounds = reducedMdMap.getGlobalRankBounds(newAxis);
          TEST_EQUALITY(bounds.start(),
                        ghosts[newAxis+1]+localSize* axisRank   );
          TEST_EQUALITY(bounds.stop() ,
                        ghosts[newAxis+1]+localSize*(axisRank+1));
          TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
          bounds = reducedMdMap.getLocalBounds(newAxis);
          TEST_ASSERT(reducedMdMap.hasHalos());
          if (reducedMdMap.getAxisRank(newAxis) == 0)
          {
            TEST_EQUALITY(bounds.start(), ghosts[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , ghosts[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis),
                          ghosts[newAxis+1]);
          }
          else if (reducedMdMap.getAxisRank(newAxis) ==
                   reducedMdMap.getAxisCommSize(newAxis)-1)
          {
            TEST_EQUALITY(bounds.start(), halos[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , halos[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis),halos[newAxis+1]);
            TEST_EQUALITY(reducedMdMap.getUpperHalo(newAxis),ghosts[newAxis+1]);
          }
          else
          {
            TEST_EQUALITY(bounds.start(), halos[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , halos[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerHalo(newAxis), halos[newAxis+1]);
            TEST_EQUALITY(reducedMdMap.getUpperHalo(newAxis), halos[newAxis+1]);
          }
          TEST_EQUALITY(reducedMdMap.getLowerGhost(newAxis), ghosts[newAxis+1]);
          TEST_EQUALITY(reducedMdMap.getUpperGhost(newAxis), ghosts[newAxis+1]);
          TEST_EQUALITY(reducedMdMap.getGhostSize(newAxis) , ghosts[newAxis]+1);
        }
      }
    }
    else
    {
      TEST_ASSERT(not reducedMdMap.onSubcommunicator());
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDMap, subMapPeriodic, T )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  // Construct the periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  MDCommRCP mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes, periodic));
  
  // Figure out the lower slice
  Array< Slice > slices(numDims);  
  Array< T >     dimensions(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimensions[axis] = 10 * axisCommSizes[axis];
    if (axis == 1)
      slices[axis] = Slice(dimensions[axis]/2);
    else
      slices[axis] = Slice();
  }

  // Construct the MDMap and the sub-MDMap
  MDMap<T> mdMap(mdComm, dimensions);
  MDMap<T> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDMap?
  bool partOfSubComm = true;
  if (mdComm->getNumDims() > 1)
    if (mdComm->getAxisRank(1) > (mdComm->getAxisCommSize(1)-1)/2)
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDMap.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      if (axis == 1)
      {
        TEST_EQUALITY(subMDMap.getAxisCommSize(axis), (axisCommSizes[axis]+1)/2);
        TEST_EQUALITY(subMDMap.getGlobalDim(axis), 5*axisCommSizes[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getAxisCommSize(axis), axisCommSizes[axis]);
        TEST_EQUALITY(subMDMap.getGlobalDim(axis), 10*axisCommSizes[axis]);
      }
      TEST_EQUALITY(subMDMap.isPeriodic(axis), (axis == 0));
    }
  }
  else
  {
    TEST_EQUALITY(subMDMap.getNumDims(), 0);
  }
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
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, exceptions, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapLowerLeft, T )   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapLowerLeftWithHalo, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapLowerRight, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapLowerRightWithGhost, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapUpperLeft, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapUpperLeftHaloGhost, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapUpperRight, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapUpperRightNewGhosts, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapReduce, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDMap, subMapPeriodic, T )

UNIT_TEST_GROUP(int)
#if 0
UNIT_TEST_GROUP(long)
#endif

}  // namespace
