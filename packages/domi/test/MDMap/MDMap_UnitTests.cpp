/*
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
typedef Domi::size_type size_type;
typedef Domi::dim_type dim_type;
using Domi::Slice;
const dim_type & Default = Domi::Slice::Default;
using Domi::MDComm;
using Domi::MDMap;

int num_dims = 2;
string commDimsStr = "-1";
Array< int > commDims;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("numDims" , &num_dims   , "number of dimensions");
  clp.setOption("commDims", &commDimsStr, "comma-separated list of "
                "number of processors along each axis");
}

//
// Unit Tests
//

TEUCHOS_UNIT_TEST( MDMap, dimensionsConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);
  ArrayView< const dim_type > const_dims(dims());

  // Construct an MDMap
  MDMap<> mdMap(mdComm, const_dims);

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(    mdMap.onSubcommunicator());
  TEST_ASSERT(not mdMap.hasPadding()       );
  TEST_EQUALITY(mdMap.getTeuchosComm(), comm               );
  TEST_EQUALITY(mdMap.getMDComm()     , mdComm             );
  TEST_EQUALITY(mdMap.numDims()       , num_dims           );
  TEST_EQUALITY(mdMap.getLayout()     , Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int   axisRank            = mdMap.getCommIndex(axis);
    Slice globalBounds        = mdMap.getGlobalBounds(axis);
    Slice globalRankBounds    = mdMap.getGlobalRankBounds(axis);
    Slice localBounds         = mdMap.getLocalBounds(axis);
    Slice localInteriorBounds = mdMap.getLocalInteriorBounds(axis);
    int   intStart            = (mdMap.getLowerNeighbor(axis) == -1) ? 1 : 0;
    int   intStop             = (mdMap.getUpperNeighbor(axis) == -1) ?
                                localDim-1 : localDim;
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(      mdMap.getCommDim(axis)      , commDims[axis]       );
    TEST_EQUALITY(      mdMap.getGlobalDim(axis)    , dims[axis]           );
    TEST_EQUALITY_CONST(globalBounds.start()        , 0                    );
    TEST_EQUALITY(      globalBounds.stop()         , dims[axis]           );
    TEST_EQUALITY(      globalRankBounds.start()    , axisRank    *localDim);
    TEST_EQUALITY(      globalRankBounds.stop()     , (axisRank+1)*localDim);
    TEST_EQUALITY(      mdMap.getLocalDim(axis)     , localDim             );
    TEST_EQUALITY_CONST(localBounds.start()         , 0                    );
    TEST_EQUALITY(      localBounds.stop()          , localDim             );
    TEST_EQUALITY_CONST(localInteriorBounds.start() , intStart             );
    TEST_EQUALITY(      localInteriorBounds.stop()  , intStop              );
    TEST_EQUALITY_CONST(mdMap.getLowerPadSize(axis) , 0                    );
    TEST_EQUALITY_CONST(mdMap.getUpperPadSize(axis) , 0                    );
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis)  , 0                    );
    TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0                    );
    TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0                    );
    TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis) , 0                    );
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = dims[axis-1] * gstrides[axis-1];
    lstrides[axis] = localDim     * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += ( axisRank   *localDim  ) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    lur += (localDim-1) * lstrides[axis];
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(  0), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(         0),  axisRank   *localDim  );
    TEST_EQUALITY(epetraMap->GID(localDim-1), (axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(  0), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim);
    TEST_EQUALITY(tpetraMap->getGlobalElement(localDim-1),
                  (axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, pListDimensionsConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);

  // Construct dummy dimensions
  Array< dim_type > dims(num_dims, comm->getSize());

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("comm dimensions", commDims);
  plist.set("dimensions"     , dims    );

  // Get the actual communicator dimensions
  Teuchos::Array< int > actualCommDims =
    Domi::regularizeCommDims(comm->getSize(), plist);

  // Compute actual dimensions
  dim_type localDim = 10;
  for (int axis = 0; axis < num_dims; ++axis)
      dims[axis] = localDim * actualCommDims[axis];
  plist.set("dimensions", dims);

  // Construct an MDMap
  MDMap<> mdMap(comm, plist);

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  TEST_ASSERT(not mdMap.hasPadding());
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    TEST_EQUALITY(mdMap.getCommDim(axis), actualCommDims[axis]);
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
    TEST_EQUALITY_CONST(mdMap.getLowerPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = dims[axis-1] * gstrides[axis-1];
    lstrides[axis] = localDim     * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += ( axisRank   *localDim  ) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    lur += (localDim-1) * lstrides[axis];
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(  0), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(         0),  axisRank   *localDim  );
    TEST_EQUALITY(epetraMap->GID(localDim-1), (axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(  0), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim);
    TEST_EQUALITY(tpetraMap->getGlobalElement(localDim-1),
                  (axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, commPadConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication padding
  Array< int > commPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commPad[axis] = axis+1;

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dims(), commPad());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  TEST_EQUALITY(mdMap.hasPadding(), comm->getSize() > 1);
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    int lowerCommPad = 0;
    if (axisRank > 0) lowerCommPad = commPad[axis];
    int upperCommPad = 0;
    if (axisRank < commDims[axis]-1) upperCommPad = commPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), commDims[axis]);
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
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY_CONST(mdMap.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis-1);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad += commPad[axis-1];
    if (axisRank != mdMap.getCommDim(axis-1)-1)
      myCommPad += commPad[axis-1];
    gstrides[axis] = dims[axis-1]       * gstrides[axis-1];
    lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += ( axisRank   *localDim) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    if (axisRank == 0)
      lur += (localDim-1) * lstrides[axis];
    else
    {
      lll += commPad[axis] * lstrides[axis];
      lur += (localDim+commPad[axis]-1) * lstrides[axis];
    }
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad = commPad[axis];
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(0), axisRank*localDim-myCommPad);
    TEST_EQUALITY(epetraMap->GID(localDim+myCommPad-1), (axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad = commPad[axis];
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraAxisMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(0),
                  axisRank*localDim-myCommPad);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(localDim+myCommPad-1),
                  (axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, pListCommPadConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication padding
  Array< int > commPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"             , dims   );
  plist.set("communication pad sizes", commPad);

  // Construct an MDMap
  MDMap<> mdMap(mdComm, plist);

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  TEST_EQUALITY(mdMap.hasPadding(), comm->getSize() > 1);
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    int lowerCommPad = 0;
    if (axisRank > 0) lowerCommPad = commPad[axis];
    int upperCommPad = 0;
    if (axisRank < commDims[axis]-1) upperCommPad = commPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), commDims[axis]);
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
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY_CONST(mdMap.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY_CONST(mdMap.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis-1);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad += commPad[axis-1];
    if (axisRank != mdMap.getCommDim(axis-1)-1)
      myCommPad += commPad[axis-1];
    gstrides[axis] = dims[axis-1]         * gstrides[axis-1];
    lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += ( axisRank   *localDim) * gstrides[axis];
    gur += ((axisRank+1)*localDim-1) * gstrides[axis];
    if (axisRank == 0)
      lur += (localDim-1) * lstrides[axis];
    else
    {
      lll += commPad[axis] * lstrides[axis];
      lur += (localDim+commPad[axis]-1) * lstrides[axis];
    }
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad = commPad[axis];
    Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(0), axisRank*localDim-myCommPad);
    TEST_EQUALITY(epetraMap->GID(localDim+myCommPad-1), (axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank > 0) myCommPad = commPad[axis];
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraAxisMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(0),
                  axisRank*localDim-myCommPad);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(localDim+myCommPad-1),
                  (axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, bndryPadConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication padding
  Array< int > commPad;

  // Construct boundary padding
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    bndryPad[axis] = axis+1;

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dims(), commPad(), bndryPad());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  bool hasPadding = false;
  for (int axis = 0; axis < num_dims; axis++)
  {
    if (mdMap.getCommIndex(axis) == 0               ) hasPadding = true;
    if (mdMap.getCommIndex(axis) == commDims[axis]-1) hasPadding = true;
  }
  TEST_EQUALITY(mdMap.hasPadding(), hasPadding);
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    int lowerBndryPad = 0;
    if (axisRank == 0) lowerBndryPad = bndryPad[axis];
    int upperBndryPad = 0;
    if (axisRank == commDims[axis]-1) upperBndryPad = bndryPad[axis];
    dim_type myDim = localDim + lowerBndryPad + upperBndryPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), commDims[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]                 );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= bndryPad[axis];
    if (axisRank == commDims[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerBndryPad);
    TEST_EQUALITY(localBounds.stop(), lowerBndryPad+localDim);
    TEST_EQUALITY(mdMap.getLowerPadSize(axis), lowerBndryPad);
    TEST_EQUALITY(mdMap.getUpperPadSize(axis), upperBndryPad);
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
    TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
    int myBndryPad = 0;
    int axisRank = mdMap.getCommIndex(axis-1);
    if (axisRank == 0) myBndryPad += bndryPad[axis-1];
    if (axisRank == mdMap.getCommDim(axis-1)-1)
      myBndryPad += bndryPad[axis-1];
    lstrides[axis] = (localDim+myBndryPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
    gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
    if (axisRank == 0)
    {
      lll += bndryPad[axis] * lstrides[axis];
      lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
    }
    else
      lur += (localDim-1) * lstrides[axis];
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank == 0) myCommPad = bndryPad[axis];
    Teuchos::RCP< const Epetra_Map > epetraMap =
      mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(myCommPad), axisRank*localDim+bndryPad[axis]);
    TEST_EQUALITY(epetraMap->GID(myCommPad+localDim-1),
                  bndryPad[axis]+(axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal> > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank == 0) myCommPad = bndryPad[axis];
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraAxisMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(myCommPad),
                  axisRank*localDim+bndryPad[axis]);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(myCommPad+localDim-1),
                  bndryPad[axis]+(axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, pListBndryPadConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);

  // Get the actual communicator dimensions
  Teuchos::Array< int > actualCommDims =
    Domi::regularizeCommDims(comm->getSize(), num_dims, commDims);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * actualCommDims[axis];

  // Construct boundary padding
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    bndryPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("comm dimensions"   , commDims);
  plist.set("dimensions"        , dims    );
  plist.set("boundary pad sizes", bndryPad);

  // Construct an MDMap
  MDMap<> mdMap(comm, plist);

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  bool hasPadding = false;
  for (int axis = 0; axis < num_dims; axis++)
  {
    if (mdMap.getCommIndex(axis) == 0                     ) hasPadding = true;
    if (mdMap.getCommIndex(axis) == actualCommDims[axis]-1) hasPadding = true;
  }
  TEST_EQUALITY(mdMap.hasPadding(), hasPadding);
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    int lowerBndryPad = 0;
    if (axisRank == 0) lowerBndryPad = bndryPad[axis];
    int upperBndryPad = 0;
    if (axisRank == actualCommDims[axis]-1) upperBndryPad = bndryPad[axis];
    dim_type myDim = localDim + lowerBndryPad + upperBndryPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), actualCommDims[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]                 );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                     ) myStart -= bndryPad[axis];
    if (axisRank == actualCommDims[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerBndryPad);
    TEST_EQUALITY(localBounds.stop(), lowerBndryPad+localDim);
    TEST_EQUALITY(mdMap.getLowerPadSize(axis), lowerBndryPad);
    TEST_EQUALITY(mdMap.getUpperPadSize(axis), upperBndryPad);
    TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
    TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
    int myBndryPad = 0;
    int axisRank = mdMap.getCommIndex(axis-1);
    if (axisRank == 0) myBndryPad += bndryPad[axis-1];
    if (axisRank == mdMap.getCommDim(axis-1)-1)
      myBndryPad += bndryPad[axis-1];
    lstrides[axis] = (localDim+myBndryPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
    gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
    if (axisRank == 0)
    {
      lll += bndryPad[axis] * lstrides[axis];
      lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
    }
    else
      lur += (localDim-1) * lstrides[axis];
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank == 0) myCommPad = bndryPad[axis];
    Teuchos::RCP< const Epetra_Map > epetraMap =
      mdMap.getEpetraAxisMap(axis);
    TEST_EQUALITY(epetraMap->GID(myCommPad), axisRank*localDim+bndryPad[axis]);
    TEST_EQUALITY(epetraMap->GID(myCommPad+localDim-1),
                  bndryPad[axis]+(axisRank+1)*localDim-1);
  }
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank  = mdMap.getCommIndex(axis);
    int myCommPad = 0;
    if (axisRank == 0) myCommPad = bndryPad[axis];
    Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraAxisMap =
      mdMap.getTpetraAxisMap< Ordinal >(axis);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(myCommPad),
                  axisRank*localDim+bndryPad[axis]);
    TEST_EQUALITY(tpetraAxisMap->getGlobalElement(myCommPad+localDim-1),
                  bndryPad[axis]+(axisRank+1)*localDim-1);
  }
#endif
}

TEUCHOS_UNIT_TEST( MDMap, paddingConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication and boundary padding
  Array< int > commPad(num_dims);
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    commPad[axis]  = axis+1;
    bndryPad[axis] = axis+2;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dims(), commPad(), bndryPad());

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  TEST_ASSERT(mdMap.hasPadding());
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank     = mdMap.getCommIndex(axis);
    int lowerCommPad = commPad[axis];
    int upperCommPad = commPad[axis];
    if (axisRank == 0                ) lowerCommPad = bndryPad[axis];
    if (axisRank == commDims[axis]-1) upperCommPad = bndryPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), commDims[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]                 );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= bndryPad[axis];
    if (axisRank == commDims[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY(mdMap.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY(mdMap.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY(mdMap.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
    int myCommPad = 0;
    int axisRank  = mdMap.getCommIndex(axis-1);
    if (axisRank == 0)
      myCommPad += bndryPad[axis-1];
    else
      myCommPad += commPad[axis-1];
    if (axisRank == mdMap.getCommDim(axis-1)-1)
      myCommPad += bndryPad[axis-1];
    else
      myCommPad += commPad[axis-1];
    lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
    gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
    if (axisRank == 0)
    {
      lll += bndryPad[axis] * lstrides[axis];
      lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
    }
    else
    {
      lll += commPad[axis] * lstrides[axis];
      lur += (localDim+commPad[axis]-1) * lstrides[axis];
    }
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);
#endif
}

TEUCHOS_UNIT_TEST( MDMap, pListPaddingConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication and boundary padding
  Array< int > commPad(num_dims);
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    commPad[axis]  = axis+1;
    bndryPad[axis] = axis+2;
  }

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"             , dims    );
  plist.set("communication pad sizes", commPad );
  plist.set("boundary pad sizes"     , bndryPad);

  // Construct an MDMap
  MDMap<> mdMap(mdComm, plist);

  // Perform unit tests of MDMap as a whole
  TEST_ASSERT(mdMap.onSubcommunicator());
  TEST_EQUALITY(mdMap.numDims(), num_dims);
  TEST_ASSERT(mdMap.hasPadding());
  TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDMap axis quantities
  for (int axis = 0; axis < num_dims; ++axis)
  {
    int axisRank     = mdMap.getCommIndex(axis);
    int lowerCommPad = commPad[axis];
    int upperCommPad = commPad[axis];
    if (axisRank == 0               ) lowerCommPad = bndryPad[axis];
    if (axisRank == commDims[axis]-1) upperCommPad = bndryPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdMap.getCommDim(axis), commDims[axis]);
    TEST_ASSERT(not mdMap.isPeriodic(axis));
    TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]                 );
    TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= bndryPad[axis];
    if (axisRank == commDims[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdMap.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdMap.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY(mdMap.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY(mdMap.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY(mdMap.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
  }

  // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
  Teuchos::Array<int> gstrides(num_dims);
  Teuchos::Array<int> lstrides(num_dims);
  gstrides[0] = 1;
  lstrides[0] = 1;
  for (int axis=1; axis < num_dims; ++axis)
  {
    gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
    int myCommPad = 0;
    int axisRank  = mdMap.getCommIndex(axis-1);
    if (axisRank == 0)
      myCommPad += bndryPad[axis-1];
    else
      myCommPad += commPad[axis-1];
    if (axisRank == mdMap.getCommDim(axis-1)-1)
      myCommPad += bndryPad[axis-1];
    else
      myCommPad += commPad[axis-1];
    lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
  }
  int gll = 0;
  int gur = 0;
  int lll = 0;
  int lur = 0;
  for (int axis=0; axis < num_dims; ++axis)
  {
    int axisRank = mdMap.getCommIndex(axis);
    gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
    gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
    if (axisRank == 0)
    {
      lll += bndryPad[axis] * lstrides[axis];
      lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
    }
    else
    {
      lll += commPad[axis] * lstrides[axis];
      lur += (localDim+commPad[axis]-1) * lstrides[axis];
    }
  }

#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
  TEST_EQUALITY(epetraMap->GID(lll), gll);
  TEST_EQUALITY(epetraMap->GID(lur), gur);
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< const Tpetra::Map< Ordinal > > tpetraMap =
    mdMap.getTpetraMap< Ordinal >();
  TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
  TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);
#endif
}

TEUCHOS_UNIT_TEST( MDMap, indexes )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication and boundary padding
  Array< int > commPad(num_dims);
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    commPad[axis]  = axis+1;
    bndryPad[axis] = axis+2;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dims(), commPad(), bndryPad());

  // Compute some quantities for testing
  Array< int > commIndex(num_dims);
  Array< int > lowerPad(num_dims);
  Array< int > upperPad(num_dims);
  Array< dim_type > myLocalDims(num_dims);
  Array< size_type > globalStrides(num_dims, 1);
  Array< size_type > localStrides(num_dims, 1);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    commIndex[axis]   = mdMap.getCommIndex(axis);
    lowerPad[axis]    = mdMap.getLowerPadSize(axis);
    upperPad[axis]    = mdMap.getUpperPadSize(axis);
    myLocalDims[axis] = localDim + lowerPad[axis] + upperPad[axis];
    if (axis > 0)
    {
      globalStrides[axis] = globalStrides[axis-1] *
        (dims[axis-1] + 2*bndryPad[axis-1]);
      localStrides[axis] = localStrides[axis-1] * myLocalDims[axis-1];
    }
  }

  // Unit test globalIndex <-> globalID
  Array< dim_type > myGlobalIndex(num_dims);
  size_type myGlobalID = 0;
  for (int axis = 0; axis < num_dims; ++axis)
  {
    myGlobalIndex[axis] = bndryPad[axis] + dims[axis] / 2;
    myGlobalID += myGlobalIndex[axis] * globalStrides[axis];
  }
  TEST_EQUALITY(mdMap.getGlobalIndex(myGlobalID)  , myGlobalIndex);
  TEST_EQUALITY(mdMap.getGlobalID(myGlobalIndex()), myGlobalID   );

  // Unit test localIndex <-> localID
  Array< dim_type > myLocalIndex(num_dims);
  size_type myLocalID = 0;
  for (int axis = 0; axis < num_dims; ++axis)
  {
    myLocalIndex[axis] = lowerPad[axis] + localDim / 2;
    myLocalID += myLocalIndex[axis] * localStrides[axis];
  }
  TEST_EQUALITY(mdMap.getLocalIndex(myLocalID)  , myLocalIndex);
  TEST_EQUALITY(mdMap.getLocalID(myLocalIndex()), myLocalID   );

  // Test localID <-> globalID
  myGlobalID = 0;
  for (int axis = 0; axis < num_dims; ++axis)
  {
    myGlobalIndex[axis] = myLocalIndex[axis] - lowerPad[axis] +
      commIndex[axis] * localDim + bndryPad[axis];
    myGlobalID += myGlobalIndex[axis] * globalStrides[axis];
  }
  TEST_EQUALITY(mdMap.getGlobalID(myLocalID), myGlobalID);
  TEST_EQUALITY(mdMap.getLocalID(myGlobalID), myLocalID );
}

TEUCHOS_UNIT_TEST( MDMap, exceptions )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dims());

  // Unit test methods that should throw exceptions
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(mdMap.getCommDim(         -1), Domi::RangeError);
  TEST_THROW(mdMap.isPeriodic(         -1), Domi::RangeError);
  TEST_THROW(mdMap.getCommIndex(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerNeighbor(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperNeighbor(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalDim(       -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalBounds(    -1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalDim(        -1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalRankBounds(-1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalBounds(     -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerPadSize(    -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperPadSize(    -1), Domi::RangeError);
  TEST_THROW(mdMap.getCommPadSize(     -1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerBndryPad(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperBndryPad(   -1), Domi::RangeError);
  TEST_THROW(mdMap.getBndryPadSize(    -1), Domi::RangeError);
  TEST_THROW(mdMap.getCommDim(         num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.isPeriodic(         num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getCommIndex(       num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerNeighbor(   num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperNeighbor(   num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalDim(       num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalBounds(    num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalDim(        num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getGlobalRankBounds(num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLocalBounds(     num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerPadSize(    num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperPadSize(    num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getCommPadSize(     num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getLowerBndryPad(   num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getUpperBndryPad(   num_dims+1), Domi::RangeError);
  TEST_THROW(mdMap.getBndryPadSize(    num_dims+1), Domi::RangeError);
#endif
}

TEUCHOS_UNIT_TEST( MDMap, subMapLowerLeft )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions
  Array< dim_type > dimensions(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions);
  
  // Figure out the lower left slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis < 2)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < num_dims; ++axis)
    if (mdMap.getGlobalRankBounds(axis).start() >= newDims[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis), newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
      TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapLowerLeftWithCommPad )
{
  // Construct the MDComm from command-line arguments and communication padding
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions and communication padding
  Array< dim_type > dimensions(num_dims);
  Array< int > commPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = 10 * commDims[axis];
    commPad[axis]    = axis + 1;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions, commPad);
  
  // Figure out the lower left slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis < 2)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < num_dims; ++axis)
    if (mdMap.getGlobalRankBounds(axis).start() >= newDims[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis), newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
      TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
      TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (commDims[axis] > 1)
      {
        if (subMDMap.getCommIndex(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerPadSize(axis), commPad[axis]);
        }
        if (subMDMap.getCommIndex(axis) == subMDMap.getCommDim(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperPadSize(axis), commPad[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), 0);
        TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), 0);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapLowerRight )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions
  Array< dim_type > dimensions(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions);
  
  // Figure out the lower right slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis == 0)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd,Default);
    }
    else if (axis == 1)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) < commDims[0] - newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) >= newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]     );
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
    TEST_EQUALITY_CONST(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapLowerRightWithBndryPad )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions and boundary padding
  Array< dim_type > dimensions(num_dims);
  Array< int > bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = 10 * commDims[axis];
    bndryPad[axis]     = axis + 2;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions, ArrayView<int>(), bndryPad);
  
  // Figure out the lower right slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis == 0)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd+bndryPad[axis],-bndryPad[axis]);
    }
    else if (axis == 1)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(bndryPad[axis],nd+bndryPad[axis]);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]      = dimensions[axis];
      slices[axis]       = Slice(bndryPad[axis], -bndryPad[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) < commDims[0] - newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) >= newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]     );
      if (axis == 0)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + bndryPad[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + bndryPad[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapUpperLeft )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions
  Array< dim_type > dimensions(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions);
  
  // Figure out the upper left slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis == 0)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd);
    }
    else if (axis == 1)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd,Default);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) >= newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) < commDims[1] - newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]    );
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
    TEST_EQUALITY(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapUpperLeftPadding )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions, padding
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims);
  Array< int >      bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = 10 * commDims[axis];
    commPad[axis]    = axis + 1;
    bndryPad[axis]   = axis + 2;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Figure out the upper left slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis == 0)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      if (nc == 0) nc = 1;
      newCommDims[axis] = nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(bndryPad[axis],nd+bndryPad[axis]);
    }
    else if (axis == 1)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd+bndryPad[axis],-bndryPad[axis]);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice(bndryPad[axis], -bndryPad[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) >= newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) < commDims[1] - newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]    );
      if (axis == 1)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + bndryPad[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + bndryPad[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (commDims[axis] > 1)
      {
        if (subMDMap.getCommIndex(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerPadSize(axis), commPad[axis]);
        }
        if (subMDMap.getCommIndex(axis) == subMDMap.getCommDim(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), 0);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperPadSize(axis), commPad[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), 0);
        TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), 0);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapUpperRight )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions
  Array< dim_type > dimensions(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions);
  
  // Figure out the upper right slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    if (axis < 2)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - nc;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd,Default);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice();
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) < commDims[0] - newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) < commDims[1] - newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]    );
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
    TEST_EQUALITY_CONST(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapUpperRightNewBndryPad )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the dimensions and padding
  int mySize = 10;
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims);
  Array< int >      bndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = mySize * commDims[axis];
    commPad[axis]      = axis + 1;
    bndryPad[axis]     = axis + 2;
  }

  // Construct an MDMap
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Compute the new boundary padding sizes and figure out the upper
  // right slice
  Array< Slice >    slices(num_dims);
  Array< int >      newCommDims(num_dims);
  Array< dim_type > newDims(num_dims);
  Array< int >      newBndryPad(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    // Make the new boundary padding one less than the old boundary padding
    newBndryPad[axis] = bndryPad[axis] - 1;
    if (axis < 2)
    {
      dim_type nd = dimensions[axis] / 2;
      int      nc = commDims[axis] / 2;
      newCommDims[axis] = commDims[axis] - (nd-newBndryPad[axis])/mySize;
      newDims[axis]     = nd;
      slices[axis]      = Slice(nd+bndryPad[axis],-bndryPad[axis]);
    }
    else
    {
      newCommDims[axis] = commDims[axis];
      newDims[axis]     = dimensions[axis];
      slices[axis]      = Slice(bndryPad[axis], -bndryPad[axis]);
    }
  }

  // Construct the sub-MDMap
  MDMap<> subMDMap(mdMap, slices, newBndryPad);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm->getCommIndex(0) < commDims[0] - newCommDims[0])
    partOfSubComm = false;
  if (num_dims > 1)
    if (mdComm->getCommIndex(1) < commDims[1] - newCommDims[1])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_ASSERT(subMDMap.onSubcommunicator());
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      TEST_EQUALITY(subMDMap.getCommDim(axis)  , newCommDims[axis]);
      TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]    );
      if (axis < 2)
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      newDims[axis] + bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      dimensions[axis] + bndryPad[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
                      bndryPad[axis]);
        TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
                      newDims[axis] + bndryPad[axis]);
      }
      TEST_ASSERT(not subMDMap.isPeriodic(axis));
      if (commDims[axis] > 1)
      {
        if (subMDMap.getCommIndex(axis) == 0)
        {
          TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), newBndryPad[axis]);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getLowerPadSize(axis), commPad[axis]);
        }
        if (subMDMap.getCommIndex(axis) == subMDMap.getCommDim(axis)-1)
        {
          TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), newBndryPad[axis]);
        }
        else
        {
          TEST_EQUALITY(subMDMap.getUpperPadSize(axis), commPad[axis]);
        }
      }
      else
      {
        TEST_EQUALITY_CONST(subMDMap.getLowerPadSize(axis), newBndryPad[axis]);
        TEST_EQUALITY_CONST(subMDMap.getUpperPadSize(axis), newBndryPad[axis]);
      }
    }
  }
  else
  {
    TEST_ASSERT(not subMDMap.onSubcommunicator());
    TEST_EQUALITY_CONST(subMDMap.numDims(), 0);
    TEST_THROW(subMDMap.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapReduce )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the parent MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims);
  Array< int >      bndryPad(num_dims);
  int localSize = 10;
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = localSize*commDims[axis];
    commPad[axis] = axis + 1;
    bndryPad[axis] = axis + 2;
  }
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // We will reduce this parent MDMap several times by using the
  // single dim_type constructor along each dimension
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dim_type myOrd = bndryPad[axis] + dimensions[axis] / 2;
    int      newDims = (num_dims > 1) ? num_dims-1 : num_dims;
    Slice    bounds = mdMap.getGlobalRankBounds(axis);
    bool     partOfSubMap = true;
    if ((myOrd < bounds.start()) || (bounds.stop() <= myOrd))
      partOfSubMap = false;
    MDMap<> reducedMdMap(mdMap, axis, myOrd);
    if (partOfSubMap)
    {
      TEST_ASSERT(reducedMdMap.onSubcommunicator());
      TEST_EQUALITY(reducedMdMap.numDims(), newDims);
      for (int newAxis = 0; newAxis < newDims; ++newAxis)
      {
        int axisRank = reducedMdMap.getCommIndex(newAxis);
        if (num_dims == 1)
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
          TEST_ASSERT(not reducedMdMap.hasPadding());
          TEST_EQUALITY_CONST(reducedMdMap.getLowerPadSize(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getUpperPadSize(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getCommPadSize(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getLowerBndryPad(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getUpperBndryPad(newAxis), 0);
          TEST_EQUALITY_CONST(reducedMdMap.getBndryPadSize(newAxis), 0);
        }
        else if (newAxis < axis)
        {
          TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
                        dimensions[newAxis]);
          Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), bndryPad[newAxis]                    );
          TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+dimensions[newAxis]);
          bounds = reducedMdMap.getGlobalRankBounds(newAxis);
          TEST_EQUALITY(bounds.start(), bndryPad[newAxis]+localSize* axisRank   );
          TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+localSize*(axisRank+1));
          TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
          bounds = reducedMdMap.getLocalBounds(newAxis);
          TEST_ASSERT(reducedMdMap.hasPadding());
          if (reducedMdMap.getCommIndex(newAxis) == 0)
          {
            TEST_EQUALITY(bounds.start(), bndryPad[newAxis]          );
            TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis), bndryPad[newAxis]);
          }
          else if (reducedMdMap.getCommIndex(newAxis) ==
                   reducedMdMap.getCommDim(newAxis)-1)
          {
            TEST_EQUALITY(bounds.start(), commPad[newAxis]          );
            TEST_EQUALITY(bounds.stop() , commPad[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis), commPad[newAxis]);
            TEST_EQUALITY(reducedMdMap.getUpperPadSize(newAxis), bndryPad[newAxis]);
          }
          else
          {
            TEST_EQUALITY(bounds.start(), commPad[newAxis]          );
            TEST_EQUALITY(bounds.stop() , commPad[newAxis]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis), commPad[newAxis]);
            TEST_EQUALITY(reducedMdMap.getUpperPadSize(newAxis), commPad[newAxis]);
          }
          TEST_EQUALITY(reducedMdMap.getLowerBndryPad(newAxis), bndryPad[newAxis]);
          TEST_EQUALITY(reducedMdMap.getUpperBndryPad(newAxis), bndryPad[newAxis]);
          TEST_EQUALITY(reducedMdMap.getBndryPadSize(newAxis) , bndryPad[newAxis]);
        }
        else
        {
          TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
                        dimensions[newAxis+1]);
          Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
          TEST_EQUALITY(bounds.start(), bndryPad[newAxis+1]);
          TEST_EQUALITY(bounds.stop() , bndryPad[newAxis+1]+
                        dimensions[newAxis+1]);
          bounds = reducedMdMap.getGlobalRankBounds(newAxis);
          TEST_EQUALITY(bounds.start(),
                        bndryPad[newAxis+1]+localSize* axisRank   );
          TEST_EQUALITY(bounds.stop() ,
                        bndryPad[newAxis+1]+localSize*(axisRank+1));
          TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
          bounds = reducedMdMap.getLocalBounds(newAxis);
          TEST_ASSERT(reducedMdMap.hasPadding());
          if (reducedMdMap.getCommIndex(newAxis) == 0)
          {
            TEST_EQUALITY(bounds.start(), bndryPad[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , bndryPad[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis),
                          bndryPad[newAxis+1]);
          }
          else if (reducedMdMap.getCommIndex(newAxis) ==
                   reducedMdMap.getCommDim(newAxis)-1)
          {
            TEST_EQUALITY(bounds.start(), commPad[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , commPad[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis),commPad[newAxis+1]);
            TEST_EQUALITY(reducedMdMap.getUpperPadSize(newAxis),bndryPad[newAxis+1]);
          }
          else
          {
            TEST_EQUALITY(bounds.start(), commPad[newAxis+1]          );
            TEST_EQUALITY(bounds.stop() , commPad[newAxis+1]+localSize);
            TEST_EQUALITY(reducedMdMap.getLowerPadSize(newAxis), commPad[newAxis+1]);
            TEST_EQUALITY(reducedMdMap.getUpperPadSize(newAxis), commPad[newAxis+1]);
          }
          TEST_EQUALITY(reducedMdMap.getLowerBndryPad(newAxis), bndryPad[newAxis+1]);
          TEST_EQUALITY(reducedMdMap.getUpperBndryPad(newAxis), bndryPad[newAxis+1]);
          TEST_EQUALITY(reducedMdMap.getBndryPadSize(newAxis) , bndryPad[newAxis]+1);
        }
      }
    }
    else
    {
      TEST_ASSERT(not reducedMdMap.onSubcommunicator());
    }
  }
}

TEUCHOS_UNIT_TEST( MDMap, subMapPeriodic )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  // Construct the periodic flags
  Array< int > periodic(num_dims, 0);
  periodic[0] = 1;
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims, periodic));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Figure out the lower slice
  Array< Slice >    slices(num_dims);  
  Array< dim_type > dimensions(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
  {
    dimensions[axis] = 10 * commDims[axis];
    if (axis == 1)
      slices[axis] = Slice(dimensions[axis]/2);
    else
      slices[axis] = Slice();
  }

  // Construct the MDMap and the sub-MDMap
  MDMap<> mdMap(mdComm, dimensions);
  MDMap<> subMDMap(mdMap, slices);

  // Should this processor be a part of the sub-MDMap?
  bool partOfSubComm = true;
  if (mdComm->numDims() > 1)
    if (mdComm->getCommIndex(1) > (mdComm->getCommDim(1)-1)/2)
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDMap.numDims(), num_dims);
    for (int axis = 0; axis < num_dims; ++axis)
    {
      if (axis == 1)
      {
        TEST_EQUALITY(subMDMap.getCommDim(axis), (commDims[axis]+1)/2);
        TEST_EQUALITY(subMDMap.getGlobalDim(axis), 5*commDims[axis]);
      }
      else
      {
        TEST_EQUALITY(subMDMap.getCommDim(axis), commDims[axis]);
        TEST_EQUALITY(subMDMap.getGlobalDim(axis), 10*commDims[axis]);
      }
      TEST_EQUALITY(subMDMap.isPeriodic(axis), (axis == 0));
    }
  }
  else
  {
    TEST_EQUALITY(subMDMap.numDims(), 0);
  }
}

TEUCHOS_UNIT_TEST( MDMap, isPad )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims,1);
  Array< int >      bndryPad(num_dims,2);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // I want to implement a test that is relatively simple and works
  // with different numbers of dimensions (1D, 2D, 3D).  The MDArray
  // class has a serial iterator for arbitrary number of dimensions.
  // So I allocate an MDArray using the MDMap local dimensions
  // (including padding) and iterate over the MDArray as a proxy for
  // iterating over the MDMap indexes.
  Array< dim_type > localDims(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    localDims[axis] = mdMap.getLocalDim(axis,true);
  Domi::MDArray< int > proxy(localDims);
  Array< dim_type > index(num_dims);
  typedef Domi::MDArray< int >::iterator iterator;
  for (iterator it = proxy.begin(); it != proxy.end(); ++it)
  {
    bool isPad      = false;
    bool isCommPad  = false;
    bool isBndryPad = false;
    for (int axis = 0; axis < num_dims; ++axis)
    {
      index[axis] = it.index(axis);
      // Check the lower boundary
      if (mdComm->getLowerNeighbor(axis) == -1)
      {
        if (index[axis] < bndryPad[axis])
        {
          isBndryPad = true;
        }
      }
      else
      {
        if (index[axis] < commPad[axis])
        {
          isCommPad = true;
        }
      }
      // Check the upper boundary
      if (mdComm->getUpperNeighbor(axis) == -1)
      {
        if (localDims[axis] - index[axis] <= bndryPad[axis])
        {
          isBndryPad = true;
        }
      }
      else
      {
        if (localDims[axis] - index[axis] <= commPad[axis])
        {
          isCommPad = true;
        }
      }
    }
    isPad = isCommPad || isBndryPad;
    TEST_EQUALITY(mdMap.isPad(index)     , isPad     );
    TEST_EQUALITY(mdMap.isCommPad(index) , isCommPad );
    TEST_EQUALITY(mdMap.isBndryPad(index), isBndryPad);
  }
}

TEUCHOS_UNIT_TEST( MDMap, augmentedLeading )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims,1);
  Array< int >      bndryPad(num_dims,2);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Compute the augmented MDMap with leading dimension
  Teuchos::RCP< const MDMap<> > newMdMap = mdMap.getAugmentedMDMap(3);

  // Check global MDMap attributes
  TEST_ASSERT(newMdMap->onSubcommunicator());
  TEST_EQUALITY(newMdMap->numDims(), num_dims+1);

  // Check the MDMap attributes for axis 0
  TEST_ASSERT(! newMdMap->isPeriodic(0));
  TEST_EQUALITY_CONST(newMdMap->getCommDim(0)                 , 1);
  TEST_EQUALITY_CONST(newMdMap->getGlobalDim(0)               , 3);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(0).start()    , 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(0).stop()     , 3);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(0).start(), 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(0).stop() , 3);
  TEST_EQUALITY_CONST(newMdMap->getLocalDim(0)                , 3);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(0).start()     , 0);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(0).stop()      , 3);
  TEST_EQUALITY_CONST(newMdMap->getLowerPadSize(0)            , 0);
  TEST_EQUALITY_CONST(newMdMap->getUpperPadSize(0)            , 0);

  // Check the MDMap attributes for the remaining axes
  for (int dof = 0; dof < 3; ++dof)
  {
    MDMap<> oldMdMap(*newMdMap, 0, dof);
    TEST_ASSERT(mdMap.isSameAs(oldMdMap));
  }
}

TEUCHOS_UNIT_TEST( MDMap, augmentedTrailing )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims,1);
  Array< int >      bndryPad(num_dims,2);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Compute the augmented MDMap with trailing dimension
  Teuchos::RCP< const MDMap<> > newMdMap = mdMap.getAugmentedMDMap(0,5);

  // Check global MDMap attributes
  TEST_ASSERT(newMdMap->onSubcommunicator());
  TEST_EQUALITY(newMdMap->numDims(), num_dims+1);

  // Check the MDMap attributes for axis num_dims
  TEST_ASSERT(! newMdMap->isPeriodic(num_dims));
  TEST_EQUALITY_CONST(newMdMap->getCommDim(num_dims)                 , 1);
  TEST_EQUALITY_CONST(newMdMap->getGlobalDim(num_dims)               , 5);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(num_dims).start()    , 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(num_dims).stop()     , 5);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(num_dims).start(), 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(num_dims).stop() , 5);
  TEST_EQUALITY_CONST(newMdMap->getLocalDim(num_dims)                , 5);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(num_dims).start()     , 0);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(num_dims).stop()      , 5);
  TEST_EQUALITY_CONST(newMdMap->getLowerPadSize(num_dims)            , 0);
  TEST_EQUALITY_CONST(newMdMap->getUpperPadSize(num_dims)            , 0);

  // Check the MDMap attributes for the remaining axes
  for (int dof = 0; dof < 5; ++dof)
  {
    MDMap<> oldMdMap(*newMdMap, num_dims, dof);
    TEST_ASSERT(mdMap.isSameAs(oldMdMap));
  }
}

TEUCHOS_UNIT_TEST( MDMap, augmentedBoth )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims,1);
  Array< int >      bndryPad(num_dims,2);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 10 * commDims[axis];
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Compute the augmented MDMap with trailing dimension
  Teuchos::RCP< const MDMap<> > newMdMap = mdMap.getAugmentedMDMap(2,3);

  // Check global MDMap attributes
  TEST_ASSERT(newMdMap->onSubcommunicator());
  TEST_EQUALITY(newMdMap->numDims(), num_dims+2);

  // Check the MDMap attributes for axis 0
  TEST_ASSERT(! newMdMap->isPeriodic(0));
  TEST_EQUALITY_CONST(newMdMap->getCommDim(0)                 , 1);
  TEST_EQUALITY_CONST(newMdMap->getGlobalDim(0)               , 2);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(0).start()    , 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(0).stop()     , 2);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(0).start(), 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(0).stop() , 2);
  TEST_EQUALITY_CONST(newMdMap->getLocalDim(0)                , 2);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(0).start()     , 0);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(0).stop()      , 2);
  TEST_EQUALITY_CONST(newMdMap->getLowerPadSize(0)            , 0);
  TEST_EQUALITY_CONST(newMdMap->getUpperPadSize(0)            , 0);

  // Check the MDMap attributes for axis num_dims+1
  TEST_ASSERT(! newMdMap->isPeriodic(num_dims+1));
  TEST_EQUALITY_CONST(newMdMap->getCommDim(num_dims+1)                 , 1);
  TEST_EQUALITY_CONST(newMdMap->getGlobalDim(num_dims+1)               , 3);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(num_dims+1).start()    , 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalBounds(num_dims+1).stop()     , 3);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(num_dims+1).start(), 0);
  TEST_EQUALITY_CONST(newMdMap->getGlobalRankBounds(num_dims+1).stop() , 3);
  TEST_EQUALITY_CONST(newMdMap->getLocalDim(num_dims+1)                , 3);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(num_dims+1).start()     , 0);
  TEST_EQUALITY_CONST(newMdMap->getLocalBounds(num_dims+1).stop()      , 3);
  TEST_EQUALITY_CONST(newMdMap->getLowerPadSize(num_dims+1)            , 0);
  TEST_EQUALITY_CONST(newMdMap->getUpperPadSize(num_dims+1)            , 0);

  // Check the MDMap attributes for the remaining axes
  for (int ldof = 0; ldof < 2; ++ldof)
  {
    MDMap<> tempMdMap(*newMdMap, 0, ldof);
    for (int tdof = 0; tdof < 3; ++tdof)
    {
      MDMap<> oldMdMap(tempMdMap, num_dims, tdof);
      TEST_ASSERT(mdMap.isSameAs(oldMdMap));
    }
  }
}

TEUCHOS_UNIT_TEST( MDMap, contiguous )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, num_dims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(num_dims);
  for (int axis = 0; axis < num_dims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct the MDMap
  Array< dim_type > dimensions(num_dims);
  Array< int >      commPad(num_dims,2);
  Array< int >      bndryPad(num_dims,3);
  for (int axis = 0; axis < num_dims; ++axis)
    dimensions[axis] = 8 * commDims[axis];
  MDMap<> mdMap(mdComm, dimensions, commPad, bndryPad);

  // Initial MDMap should be contiguous
  TEST_ASSERT(mdMap.isContiguous());

  // Take a slice of the MDMap
  MDMap<> slicedMap(mdMap, 0, 4);

  // New MDMap should not be contiguous, unless the original MDMap was 1D
  if (slicedMap.onSubcommunicator())
    TEST_EQUALITY(slicedMap.isContiguous(), (num_dims==1));
}

}  // namespace
