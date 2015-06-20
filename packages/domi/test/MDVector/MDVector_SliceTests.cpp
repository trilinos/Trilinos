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

// System include
#include <cstdlib>
#include <algorithm>

// Teuchos includes
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Tuple.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDVector.hpp"

namespace
{

using std::string;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Tuple;
using Teuchos::rcp;
typedef Domi::Ordinal Ordinal;
typedef Domi::size_type size_type;
typedef Domi::dim_type dim_type;
using Domi::MDArray;
using Domi::MDArrayView;
using Domi::Slice;
const Ordinal & Default = Domi::Slice::Default;
using Domi::MDComm;
using Domi::MDMap;
using Domi::MDVector;

int numDims = 2;
string commDimsStr = "-1";
Array< int > commDims;

////////////////////////////////////////////////////////////////////////

/** \brief Assign each element of an MDVector a unique global ID
 *
 * \param mdv [in] an MDVector, whose values will be re-assigned
 *
 */
template< typename Sca, class Node >
void assignGlobalIDsToMDVector(MDVector< Sca, Node > & mdv)
{
  // Give the iterator a simpler name
  typedef typename MDArrayView< Sca >::iterator iterator;
  // Get a pointer to the underlying MDMap
  const Teuchos::RCP< const MDMap< Node > > mdMap = mdv.getMDMap();
  // Loop over the underlying MDArrayView
  MDArrayView< Sca > mdav = mdv.getDataNonConst();
  for (iterator it = mdav.begin(); it != mdav.end(); ++it)
  {
    // Obtain the local index
    Array< dim_type > localIndex(mdv.numDims());
    for (int axis = 0; axis < mdv.numDims(); ++axis)
      localIndex[axis] = it.index(axis);
    // Convert the local index to a local ID, and then convert to a
    // global ID and store
    size_type localID = mdMap->getLocalID(localIndex);
    *it = (Sca) mdMap->getGlobalID(localID);
  }
}

////////////////////////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("numDims" , &numDims    , "number of dimensions");
  clp.setOption("commDims", &commDimsStr, "comma-separated list of "
                "number of processors along each axis");
}

////////////////////////////////////////////////////////////////////////

//
// Templated Unit Tests
//

////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, SliceLow, Sca )
{
  // Initialize the MDComm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct global dimensions and strides.  This test will have
  // boundary padding only
  dim_type localDim = 10;
  int commSize = 0;
  int bndrySize = 2;
  Array< dim_type > dims(numDims);
  Array< dim_type > strides(numDims);
  Array< int > commPad(numDims);
  Array< int > bndryPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[axis] = localDim * mdComm->getCommDim(axis);
    commPad[axis]  = commSize;
    bndryPad[axis] = bndrySize;
    if (axis == 0) strides[0] = 1;
    else strides[axis] = strides[axis-1] * (2*bndryPad[axis] + dims[axis-1]);
  }

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims(), commPad(), bndryPad()));
  MDVector< Sca > mdVector(mdMap);
  assignGlobalIDsToMDVector(mdVector);

  // Perform tests along each axis
  dim_type width = 2 * localDim / 3;
  Slice globalSlice(bndrySize, width+bndrySize);
  Slice localSlice(0, width);
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Construct the sub-vector
    MDVector< Sca > subVector(mdVector, axis, globalSlice);

    // Check whether on-processor or not
    if (commDims[axis] > 1 && mdComm->getCommIndex(axis) > 0)
    {
      TEST_ASSERT(!subVector.onSubcommunicator());
    }
    else
    {
      // Compute the sub-vector statistics
      bool contig = (axis == numDims-1);
      Sca  begin  = 0;
      Sca  end    = 0;
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        int myCommIndex = mdMap->getCommIndex(myAxis);
        dim_type min    =  myCommIndex    * localDim;
        dim_type max    = (myCommIndex+1) * localDim + bndrySize;
        if (myCommIndex  > 0                 ) min += bndrySize;
        if (myCommIndex == commDims[myAxis]-1) max += bndrySize;

        if (myAxis == axis)
        {
          begin +=  bndrySize          * strides[myAxis];
          end   += (bndrySize+width-1) * strides[myAxis];
        }
        else
        {
          begin +=  min    * strides[myAxis];
          end   += (max-1) * strides[myAxis];
        }
      }
      MDArrayView< const Sca > subArray = subVector.getData();

      // Perform the unit tests
      TEST_ASSERT(subVector.onSubcommunicator());
      TEST_EQUALITY(subVector.isContiguous()           , contig     );
      TEST_EQUALITY(subVector.numDims()                , numDims    );
      TEST_EQUALITY(subVector.getGlobalDim(axis)       , width      );
      TEST_EQUALITY(subVector.getGlobalBounds(axis)    , globalSlice);
      TEST_EQUALITY(subVector.getGlobalRankBounds(axis), globalSlice);
      TEST_EQUALITY(subVector.getLocalDim(axis)        , width      );
      TEST_EQUALITY(subVector.getLocalBounds(axis)     , localSlice );
      TEST_EQUALITY(subVector.getLowerBndryPad(axis)   , 0          );
      TEST_EQUALITY(subVector.getUpperBndryPad(axis)   , 0          );
      TEST_EQUALITY(subVector.getBndryPadSize(axis)    , 0          );
      TEST_EQUALITY(subVector.getCommPadSize(axis)     , 0          );
      TEST_EQUALITY(subVector.getLowerPadSize(axis)    , 0          );
      TEST_EQUALITY(subVector.getUpperPadSize(axis)    , 0          );
      TEST_EQUALITY(*(subArray.begin() ), begin);
      TEST_EQUALITY(*(subArray.rbegin()), end  );
    }
  }
}

////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, SliceMid, Sca )
{
  // Initialize the MDComm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct global dimensions and strides.  This test will have
  // communication padding only.
  dim_type localDim = 10;
  int      commSize = 1;
  Array< dim_type > dims(numDims);
  Array< dim_type > strides(numDims);
  Array< int >      commPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[axis] = localDim * mdComm->getCommDim(axis);
    if (axis == 0) strides[0] = 1;
    else strides[axis] = strides[axis-1] * dims[axis-1];
    commPad[axis] = commSize;
  }

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims(), commPad()));
  MDVector< Sca > mdVector(mdMap);
  assignGlobalIDsToMDVector(mdVector);

  // Perform tests along each axis
  dim_type start = localDim / 3;
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Construct the sub-vector
    dim_type stop  = dims[axis] - localDim / 3;
    dim_type width = stop - start;
    Slice slice(start, stop);
    MDVector< Sca > subVector(mdVector, axis, slice);

    // Compute sub-vector stats
    bool contig = (axis == numDims-1);
    dim_type globalRankStart = start;
    dim_type globalRankStop  = stop;
    dim_type myStart = 0;
    dim_type myStop  = localDim;
    int commIndex = mdComm->getCommIndex(axis);
    int loPad = 0;
    int hiPad = 0;
    if (commIndex > 0               ) loPad = commSize;
    if (commIndex < commDims[axis]-1) hiPad = commSize;
    int pad = loPad + hiPad;
    if (commIndex == 0)
      myStart = localDim / 3;
    else
    {
      globalRankStart = localDim * commIndex - commSize;
      myStop += commSize;
    }
    if (commIndex == commDims[axis]-1)
      myStop -= localDim / 3;
    else
    {
      globalRankStop  = localDim * (commIndex+1) + commSize;
      myStop += commSize;
    }
    dim_type myWidth = myStop - myStart - pad;
    Slice grs   = Slice(globalRankStart+loPad, globalRankStop-hiPad);
    Sca   begin = 0;
    Sca   end   = 0;
    for (int myAxis = 0; myAxis < numDims; ++myAxis)
    {
      if (myAxis == axis)
      {
        begin +=  globalRankStart   * strides[myAxis];
        end   += (globalRankStop-1) * strides[myAxis];
      }
      else
      {
        int myCommIndex = mdMap->getCommIndex(myAxis);
        int myLoPad = 0;
        int myHiPad = 0;
        if (myCommIndex > 0                 ) myLoPad = commSize;
        if (myCommIndex < commDims[myAxis]-1) myHiPad = commSize;
        dim_type min =  myCommIndex    * localDim - myLoPad;
        dim_type max = (myCommIndex+1) * localDim + myHiPad;
        begin       +=  min    * strides[myAxis];
        end         += (max-1) * strides[myAxis];
      }
    }
    MDArrayView< const Sca > subArray = subVector.getData();

    // Perform the unit tests
    Slice mySlice(loPad, loPad+myWidth);
    TEST_ASSERT(subVector.onSubcommunicator());
    TEST_EQUALITY(subVector.isContiguous()           , contig     );
    TEST_EQUALITY(subVector.numDims()                , numDims    );
    TEST_EQUALITY(subVector.getGlobalDim(axis)       , width      );
    TEST_EQUALITY(subVector.getGlobalBounds(axis)    , slice      );
    TEST_EQUALITY(subVector.getGlobalRankBounds(axis), grs        );
    TEST_EQUALITY(subVector.getLocalDim(axis)        , myWidth    );
    TEST_EQUALITY(subVector.getLocalDim(axis,true)   , myWidth+pad);
    TEST_EQUALITY(subVector.getLocalBounds(axis)     , mySlice    );
    TEST_EQUALITY(subVector.getLowerPadSize(axis)    , loPad      );
    TEST_EQUALITY(subVector.getUpperPadSize(axis)    , hiPad      );
    TEST_EQUALITY(subVector.getCommPadSize(axis)     , commSize   );
    TEST_EQUALITY(subVector.getLowerBndryPad(axis)   , 0          );
    TEST_EQUALITY(subVector.getUpperBndryPad(axis)   , 0          );
    TEST_EQUALITY(subVector.getBndryPadSize(axis)    , 0          );
    TEST_EQUALITY(*(subArray.begin() ), begin);
    TEST_EQUALITY(*(subArray.rbegin()), end  );
  }
}

////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, SliceHi, Sca )
{
  // Initialize the MDComm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct global dimensions and strides.  This test will have
  // both boundary and communication padding, set to the same size to
  // keep things simple(r).
  dim_type localDim = 10;
  int padSize = 2;
  Array< dim_type > dims(numDims);
  Array< dim_type > strides(numDims);
  Array< int >      padding(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[axis] = localDim * mdComm->getCommDim(axis);
    if (axis == 0) strides[0] = 1;
    else strides[axis] = strides[axis-1] * (dims[axis-1] + 2*padSize);
    padding[axis] = padSize;
  }

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims(), padding(), padding()));
  MDVector< Sca > mdVector(mdMap);
  assignGlobalIDsToMDVector(mdVector);

  // Perform tests along each axis
  dim_type width = 2 * localDim / 3;
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Construct the sub-vector
    Slice slice(padSize+dims[axis]-width, padSize+dims[axis]);
    Slice fullSlice(dims[axis]-width, dims[axis]+2*padSize);
    MDVector< Sca > subVector(mdVector, axis, slice, padSize);

    // Check whether on-processor or not
    if (commDims[axis] > 1 &&
        mdComm->getCommIndex(axis) < mdComm->getCommDim(axis)-1)
    {
      TEST_ASSERT(!subVector.onSubcommunicator());
    }
    else
    {
      // Compute the sub-vector statistics
      bool contig = (axis == numDims-1);
      Slice local(padSize, width+padSize);
      Sca  begin  = 0;
      Sca  end    = 0;
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        if (myAxis == axis)
        {
          begin += (dims[axis]-width      ) * strides[myAxis];
          end   += (2*padSize+dims[axis]-1) * strides[myAxis];
        }
        else
        {
          dim_type min =  mdMap->getCommIndex(myAxis)    * localDim;
          dim_type max = (mdMap->getCommIndex(myAxis)+1) * localDim + 2*padSize;
          begin +=  min    * strides[myAxis];
          end   += (max-1) * strides[myAxis];
        }
      }
      MDArrayView< const Sca > subArray = subVector.getData();

      // Perform the unit tests
      TEST_ASSERT(subVector.onSubcommunicator());
      TEST_EQUALITY(subVector.isContiguous()                , contig         );
      TEST_EQUALITY(subVector.numDims()                     , numDims        );
      TEST_EQUALITY(subVector.getGlobalDim(axis)            , width          );
      TEST_EQUALITY(subVector.getGlobalDim(axis,true)       , width+2*padSize);
      TEST_EQUALITY(subVector.getGlobalBounds(axis)         , slice          );
      TEST_EQUALITY(subVector.getGlobalBounds(axis,true)    , fullSlice      );
      TEST_EQUALITY(subVector.getGlobalRankBounds(axis)     , slice          );
      TEST_EQUALITY(subVector.getGlobalRankBounds(axis,true), fullSlice      );
      TEST_EQUALITY(subVector.getLocalDim(axis)             , width          );
      TEST_EQUALITY(subVector.getLocalBounds(axis)          , local          );
      TEST_EQUALITY(subVector.getLowerPadSize(axis)         , padSize        );
      TEST_EQUALITY(subVector.getUpperPadSize(axis)         , padSize        );
      TEST_EQUALITY(subVector.getCommPadSize(axis)          , padSize        );
      TEST_EQUALITY(subVector.getLowerBndryPad(axis)        , padSize        );
      TEST_EQUALITY(subVector.getUpperBndryPad(axis)        , padSize        );
      TEST_EQUALITY(subVector.getBndryPadSize(axis)         , padSize        );
      TEST_EQUALITY(*(subArray.begin() ), begin);
      TEST_EQUALITY(*(subArray.rbegin()), end  );
    }
  }
}

////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, SliceIndex, Sca )
{
  // Initialize the MDComm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct global dimensions and strides.
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  Array< dim_type > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[axis] = localDim * mdComm->getCommDim(axis);
    if (axis == 0) strides[0] = 1;
    else strides[axis] = strides[axis-1] * dims[axis-1];
  }

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims()));
  MDVector< Sca > mdVector(mdMap);
  assignGlobalIDsToMDVector(mdVector);

  // Perform tests along each axis
  dim_type loc = 2;
  for (int axis = 0; axis < numDims; ++axis)
  {
    MDVector< Sca > subVector(mdVector, axis, loc);

    // Check whether on-processor or not
    if (commDims[axis] > 1 && mdComm->getCommIndex(axis) > 0)
    {
      TEST_ASSERT(!subVector.onSubcommunicator());
    }
    else
    {
      int  newDims = std::max(numDims-1,1);
      bool contig  = (axis == numDims-1 || numDims == 1);
      Sca  begin   = 0;
      Sca  end     = 0;
      for (int myAxis = 0; myAxis < numDims; ++myAxis)
      {
        if (myAxis == axis)
        {
          begin += loc * strides[myAxis];
          end   += loc * strides[myAxis];
        }
        else
        {
          dim_type min =  mdMap->getCommIndex(myAxis)    * localDim;
          dim_type max = (mdMap->getCommIndex(myAxis)+1) * localDim;
          begin +=  min    * strides[myAxis];
          end   += (max-1) * strides[myAxis];
        }
      }
      MDArrayView< const Sca > subArray = subVector.getData();

      // Perform the unit tests
      TEST_ASSERT(subVector.onSubcommunicator());
      TEST_EQUALITY(subVector.isContiguous(), contig );
      TEST_EQUALITY(subVector.numDims()     , newDims);
      TEST_EQUALITY(*(subArray.begin() )    , begin  );
      TEST_EQUALITY(*(subArray.rbegin())    , end    );
    }
  }
}

////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, CompoundSlice, Sca )
{
  // Initialize the MDComm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct global dimensions and strides.
  dim_type localDim  = 10;
  int      commSize  = 1;
  int      bndrySize = 2;
  Array< dim_type > dims(numDims);
  Array< dim_type > strides(numDims);
  Array< int >      commPad(numDims);
  Array< int >      bndryPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[axis] = localDim * mdComm->getCommDim(axis);
    if (axis == 0) strides[0] = 1;
    else strides[axis] = strides[axis-1] * (dims[axis-1] + 2*bndrySize);
    commPad[axis]  = commSize;
    bndryPad[axis] = bndrySize;
  }

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims, commPad, bndryPad));
  MDVector< Sca > mdVector(mdMap);
  assignGlobalIDsToMDVector(mdVector);

  // Take a slice that shaves off one element from the lower end and
  // one element from the upper end, and do this along each axis
  MDVector< Sca > subVector(mdVector);
  Slice slice(bndrySize+1, -bndrySize-1);
  for (int axis = 0; axis < numDims; ++axis)
    subVector = MDVector< Sca >(subVector, axis, slice, bndrySize);

  // Compute the begin and end global indexes
  Sca begin = 0;
  Sca end   = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    dim_type start =  mdComm->getCommIndex(axis)    * localDim + bndrySize;
    dim_type stop  = (mdComm->getCommIndex(axis)+1) * localDim + bndrySize;
    if (mdComm->getCommIndex(axis) == 0)
      start = 1;
    else
      start -= commSize;
    if (mdComm->getCommIndex(axis) == mdComm->getCommDim(axis)-1)
      stop += bndrySize - 1;
    else
      stop += commSize;
    begin +=  start   * strides[axis];
    end   += (stop-1) * strides[axis];
  }

  // Test the new sub-MDVector
  bool contig  = (numDims == 1);
  TEST_ASSERT(subVector.onSubcommunicator());
  TEST_EQUALITY(subVector.isContiguous(), contig );
  TEST_EQUALITY(subVector.numDims()     , numDims);

  MDArrayView< const Sca > subArray = subVector.getData();
  TEST_EQUALITY(*(subArray.begin() ), begin);
  TEST_EQUALITY(*(subArray.rbegin()), end  );

  for (int axis = 0; axis < numDims; ++axis)
  {
    int loOff = 0;
    int hiOff = 0;
    if (mdComm->getCommIndex(axis) == 0) loOff = 1;
    if (mdComm->getCommIndex(axis) == mdComm->getCommDim(axis)-1) hiOff = 1;
    // Perform the unit tests
    TEST_EQUALITY(subVector.getGlobalDim(axis),
                  mdVector.getGlobalDim(axis)-2);
    TEST_EQUALITY(subVector.getGlobalDim(axis,true),
                  mdVector.getGlobalDim(axis,true)-2);
    TEST_EQUALITY(subVector.getGlobalBounds(axis).start(),
                  mdVector.getGlobalBounds(axis).start()+1);
    TEST_EQUALITY(subVector.getGlobalBounds(axis).stop(),
                  mdVector.getGlobalBounds(axis).stop()-1);
    TEST_EQUALITY(subVector.getGlobalBounds(axis,true).start(),
                  mdVector.getGlobalBounds(axis,true).start()+1);
    TEST_EQUALITY(subVector.getGlobalBounds(axis,true).stop(),
                  mdVector.getGlobalBounds(axis,true).stop()-1);
    TEST_EQUALITY(subVector.getGlobalRankBounds(axis).start(),
                  mdVector.getGlobalRankBounds(axis).start()+loOff);
    TEST_EQUALITY(subVector.getGlobalRankBounds(axis).stop(),
                  mdVector.getGlobalRankBounds(axis).stop()-hiOff);
    TEST_EQUALITY(subVector.getGlobalRankBounds(axis,true).start(),
                  mdVector.getGlobalRankBounds(axis,true).start()+loOff);
    TEST_EQUALITY(subVector.getGlobalRankBounds(axis,true).stop(),
                  mdVector.getGlobalRankBounds(axis,true).stop()-hiOff);
    TEST_EQUALITY(subVector.getLocalDim(axis),
                  mdVector.getLocalDim(axis)-loOff-hiOff);
    TEST_EQUALITY(subVector.getLocalBounds(axis).start(),
                  mdVector.getLocalBounds(axis).start());
    TEST_EQUALITY(subVector.getLocalBounds(axis).stop(),
                  mdVector.getLocalBounds(axis).stop()-loOff-hiOff);
    TEST_EQUALITY(subVector.getLowerPadSize(axis),
                  mdVector.getLowerPadSize(axis));
    TEST_EQUALITY(subVector.getUpperPadSize(axis),
                  mdVector.getUpperPadSize(axis));
    TEST_EQUALITY(subVector.getCommPadSize(axis),
                  mdVector.getCommPadSize(axis));
    TEST_EQUALITY(subVector.getLowerBndryPad(axis),
                  mdVector.getLowerBndryPad(axis));
    TEST_EQUALITY(subVector.getUpperBndryPad(axis),
                  mdVector.getUpperBndryPad(axis));
    TEST_EQUALITY(subVector.getBndryPadSize(axis),
                  mdVector.getBndryPadSize(axis));
  }
}

////////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, SliceLow     , Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, SliceMid     , Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, SliceHi      , Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, SliceIndex   , Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, CompoundSlice, Sca )

UNIT_TEST_GROUP(double)
#if 1
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(long)
#endif

}  // namespace
