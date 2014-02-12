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
#include "Domi_MDVector.hpp"

namespace
{

using std::string;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Tuple;
using Teuchos::rcp;
typedef Domi::Ordinal Ordinal;
using Domi::TeuchosCommRCP;
using Domi::MDArrayView;
using Domi::Slice;
const Ordinal & Default = Domi::Slice::Default;
using Domi::MDComm;
using Domi::MDCommRCP;
using Domi::MDMap;
using Domi::MDVector;

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

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, dimensionsConstructor, Sca, Ord )
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
  Ord localDim = 10;
  Array< Ord > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap< Ord > > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap< Ord >(mdComm, dims()));
  MDVector< Sca, Ord > mdVector(mdMap);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.getNumDims(), numDims);
  TEST_ASSERT(not mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getAxisRank(axis);
    TEST_EQUALITY(mdVector.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdVector.isPeriodic(axis));
    TEST_EQUALITY(mdVector.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY_CONST(mdVector.getGlobalBounds(axis).start(), 0);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis).stop(), dims[axis]);
    TEST_EQUALITY(mdVector.getLocalDim(axis) , localDim);
    Slice globalRankBounds = mdVector.getGlobalRankBounds(axis);
    TEST_EQUALITY(globalRankBounds.start(), axisRank    *localDim);
    TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdVector.getLocalBounds(axis);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), localDim);
    TEST_EQUALITY_CONST(mdVector.getLowerPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getBndryPadSize(axis), 0);
  }

  // Perform tests of MDVector data
  Sca scalar = 5;
  mdVector.putScalar(scalar);
  MDArrayView< const Sca > cData = mdVector.getData();
  MDArrayView< Sca > data = mdVector.getDataNonConst();
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(cData.dimension(axis), localDim);
    TEST_EQUALITY(data.dimension(axis), localDim);
  }
  for (typename MDArrayView< Sca >::const_iterator it = cData.cbegin();
       it != cData.cend(); ++it)
    TEST_EQUALITY(*it, scalar);
  for (typename MDArrayView< Sca >::iterator it = data.begin();
       it != data.end(); ++it)
    TEST_EQUALITY(*it, scalar);

  // Test dot product and norms
  Sca norm1 = scalar * Domi::computeSize(dims);
  Sca dot   = scalar * norm1 ;
  Sca norm2 = Teuchos::ScalarTraits< Sca >::squareroot(dot);
  Sca tolerance = 1.0e-14;
  TEST_COMPARE(std::abs(mdVector.dot(mdVector) - dot   ), <, tolerance);
  TEST_COMPARE(std::abs(mdVector.norm1()       - norm1 ), <, tolerance);
  TEST_COMPARE(std::abs(mdVector.norm2()       - norm2 ), <, tolerance);
  TEST_COMPARE(std::abs(mdVector.normInf()     - scalar), <, tolerance);
  TEST_COMPARE(std::abs(mdVector.normWeighted(mdVector) -
                        std::pow(scalar, 3./2.)        ), <, tolerance);
  TEST_COMPARE(std::abs(mdVector.meanValue()   - scalar), <, tolerance);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, randomize, Sca, Ord )
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
  Ord localDim = 10;
  Array< Ord > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap< Ord > > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap< Ord >(mdComm, dims()));
  MDVector< Sca, Ord > mdVector(mdMap, true);

  // Test that all values are zero
  typedef typename MDArrayView< const Sca >::const_iterator const_iterator;
  MDArrayView< const Sca > view = mdVector.getData();
  for (const_iterator it = view.cbegin(); it != view.cend(); ++it)
    TEST_EQUALITY_CONST(*it, 0);

  // Test the randomize method
  mdVector.randomize();
  int count = 0;
  for (const_iterator it = view.cbegin(); it != view.cend(); ++it)
  {
    if (*it != 0) ++count;
    TEST_COMPARE(*it, >=, -1);
    TEST_COMPARE(*it, <=,  1);
  }
  TEST_COMPARE(count, >, 0);
  
}

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListDimensionsConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);

//   // Check that the axisCommSizes are completely specified
//   TEST_EQUALITY(axisCommSizes.size(), numDims)
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     TEST_ASSERT(axisCommSizes[axis] > 0);
//   }

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< long long > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * axisCommSizes[axis];

//   // Construct a ParameterList
//   Teuchos::ParameterList plist;
//   plist.set("axis comm sizes", axisCommSizes);
//   plist.set("dimensions"     , dims         );

//   // Construct an MDMap
//   MDMap< Ord > mdMap(comm, plist);

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_ASSERT(not mdMap.hasPadding());
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis).stop(), dims[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis) , localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis);
//     TEST_EQUALITY(globalRankBounds.start(), axisRank    *localDim);
//     TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
//     Slice localBounds  = mdMap.getLocalBounds(axis);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), localDim);
//     TEST_EQUALITY_CONST(mdMap.getLowerPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getUpperPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     gstrides[axis] = dims[axis-1] * gstrides[axis-1];
//     lstrides[axis] = localDim     * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += ( axisRank   *localDim  ) * gstrides[axis];
//     gur += ((axisRank+1)*localDim-1) * gstrides[axis];
//     lur += (localDim-1) * lstrides[axis];
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(  0), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
//     TEST_EQUALITY(epetraMap->GID(         0),  axisRank   *localDim  );
//     TEST_EQUALITY(epetraMap->GID(localDim-1), (axisRank+1)*localDim-1);
//   }
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord, Ord> > tpetraMap =
//     mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(  0), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap =
//       mdMap.getTpetraAxisMap(axis);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(localDim-1),
//                   (axisRank+1)*localDim-1);
//   }
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, commPadConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< Ord > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication padding
//   Array< int > commPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     commPad[axis] = axis+1;

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dims(), commPad());

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_EQUALITY(mdMap.hasPadding(), comm->getSize() > 1);
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     int lowerCommPad = 0;
//     if (axisRank > 0) lowerCommPad = commPad[axis];
//     int upperCommPad = 0;
//     if (axisRank < axisCommSizes[axis]-1) upperCommPad = commPad[axis];
//     T myDim = localDim + lowerCommPad + upperCommPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis).stop(), dims[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis);
//     TEST_EQUALITY(globalRankBounds.start(),  axisRank   *localDim);
//     TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
//     TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
//     TEST_EQUALITY_CONST(mdMap.getLowerPad(axis), lowerCommPad);
//     TEST_EQUALITY_CONST(mdMap.getUpperPad(axis), upperCommPad);
//     TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), commPad[axis]);
//     TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis-1);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad += commPad[axis-1];
//     if (axisRank != mdMap.getAxisCommSize(axis-1)-1)
//       myCommPad += commPad[axis-1];
//     gstrides[axis] = dims[axis-1]       * gstrides[axis-1];
//     lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += ( axisRank   *localDim) * gstrides[axis];
//     gur += ((axisRank+1)*localDim-1) * gstrides[axis];
//     if (axisRank == 0)
//       lur += (localDim-1) * lstrides[axis];
//     else
//     {
//       lll += commPad[axis] * lstrides[axis];
//       lur += (localDim+commPad[axis]-1) * lstrides[axis];
//     }
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad = commPad[axis];
//     Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
//     TEST_EQUALITY(epetraMap->GID(0), axisRank*localDim-myCommPad);
//     TEST_EQUALITY(epetraMap->GID(localDim+myCommPad-1), (axisRank+1)*localDim-1);
//   }
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad = commPad[axis];
//     Teuchos::RCP< const Tpetra::Map<T> > tpetraMap =
//       mdMap.getTpetraAxisMap(axis);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim-myCommPad);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(localDim+myCommPad-1),
//                   (axisRank+1)*localDim-1);
//   }
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListCommPadConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< long long > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication padding
//   Array< int > commPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     commPad[axis] = axis+1;

//   // Construct a ParameterList
//   Teuchos::ParameterList plist;
//   plist.set("dimensions"       , dims   );
//   plist.set("communication pad", commPad);

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, plist);

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_EQUALITY(mdMap.hasPadding(), comm->getSize() > 1);
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     int lowerCommPad = 0;
//     if (axisRank > 0) lowerCommPad = commPad[axis];
//     int upperCommPad = 0;
//     if (axisRank < axisCommSizes[axis]-1) upperCommPad = commPad[axis];
//     T myDim = localDim + lowerCommPad + upperCommPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis), dims[axis]);
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis).stop(), dims[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis);
//     TEST_EQUALITY(globalRankBounds.start(),  axisRank   *localDim);
//     TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
//     TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
//     TEST_EQUALITY_CONST(mdMap.getLowerPad(axis), lowerCommPad);
//     TEST_EQUALITY_CONST(mdMap.getUpperPad(axis), upperCommPad);
//     TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), commPad[axis]);
//     TEST_EQUALITY_CONST(mdMap.getLowerBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getUpperBndryPad(axis), 0);
//     TEST_EQUALITY_CONST(mdMap.getBndryPadSize(axis), 0);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis-1);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad += commPad[axis-1];
//     if (axisRank != mdMap.getAxisCommSize(axis-1)-1)
//       myCommPad += commPad[axis-1];
//     gstrides[axis] = dims[axis-1]       * gstrides[axis-1];
//     lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += ( axisRank   *localDim) * gstrides[axis];
//     gur += ((axisRank+1)*localDim-1) * gstrides[axis];
//     if (axisRank == 0)
//       lur += (localDim-1) * lstrides[axis];
//     else
//     {
//       lll += commPad[axis] * lstrides[axis];
//       lur += (localDim+commPad[axis]-1) * lstrides[axis];
//     }
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad = commPad[axis];
//     Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraAxisMap(axis);
//     TEST_EQUALITY(epetraMap->GID(0), axisRank*localDim-myCommPad);
//     TEST_EQUALITY(epetraMap->GID(localDim+myCommPad-1), (axisRank+1)*localDim-1);
//   }
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank > 0) myCommPad = commPad[axis];
//     Teuchos::RCP< const Tpetra::Map<T> > tpetraMap =
//       mdMap.getTpetraAxisMap(axis);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(0), axisRank*localDim-myCommPad);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(localDim+myCommPad-1),
//                   (axisRank+1)*localDim-1);
//   }
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, bndryPadConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< Ord > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication padding
//   Array< int > commPad;

//   // Construct boundary padding
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     bndryPad[axis] = axis+1;

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dims(), commPad(), bndryPad());

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_ASSERT(mdMap.hasPadding());
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     int lowerBndryPad = 0;
//     if (axisRank == 0) lowerBndryPad = bndryPad[axis];
//     int upperBndryPad = 0;
//     if (axisRank == axisCommSizes[axis]-1) upperBndryPad = bndryPad[axis];
//     T myDim = localDim + lowerBndryPad + upperBndryPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
//                   2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
//                   bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
//     T myStart = bndryPad[axis] +  axisRank    * localDim;
//     T myStop  = bndryPad[axis] + (axisRank+1) * localDim;
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop() , myStop );
//     globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
//     if (axisRank == 0                    ) myStart -= bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop( ), myStop );
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerBndryPad);
//     TEST_EQUALITY(localBounds.stop(), lowerBndryPad+localDim);
//     TEST_EQUALITY(mdMap.getLowerPad(axis), lowerBndryPad);
//     TEST_EQUALITY(mdMap.getUpperPad(axis), upperBndryPad);
//     TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
//     TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
//     int myBndryPad = 0;
//     int axisRank = mdMap.getAxisRank(axis-1);
//     if (axisRank == 0) myBndryPad += bndryPad[axis-1];
//     if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
//       myBndryPad += bndryPad[axis-1];
//     lstrides[axis] = (localDim+myBndryPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
//     gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
//     if (axisRank == 0)
//     {
//       lll += bndryPad[axis] * lstrides[axis];
//       lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
//     }
//     else
//       lur += (localDim-1) * lstrides[axis];
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank == 0) myCommPad = bndryPad[axis];
//     Teuchos::RCP< const Epetra_Map > epetraMap =
//       mdMap.getEpetraAxisMap(axis);
//     TEST_EQUALITY(epetraMap->GID(myCommPad), axisRank*localDim+bndryPad[axis]);
//     TEST_EQUALITY(epetraMap->GID(myCommPad+localDim-1),
//                   bndryPad[axis]+(axisRank+1)*localDim-1);
//   }
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank == 0) myCommPad = bndryPad[axis];
//     Teuchos::RCP< const Tpetra::Map<T> > tpetraMap =
//       mdMap.getTpetraAxisMap(axis);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(myCommPad),
//                   axisRank*localDim+bndryPad[axis]);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(myCommPad+localDim-1),
//                   bndryPad[axis]+(axisRank+1)*localDim-1);
//   }
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListBndryPadConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< long long > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * axisCommSizes[axis];

//   // Construct boundary padding
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     bndryPad[axis] = axis+1;

//   // Construct a ParameterList
//   Teuchos::ParameterList plist;
//   plist.set("axis comm sizes", axisCommSizes);
//   plist.set("dimensions"     , dims         );
//   plist.set("boundary pad"   , bndryPad     );

//   // Construct an MDMap
//   MDMap< Ord > mdMap(comm, plist);

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_ASSERT(mdMap.hasPadding());
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     int lowerBndryPad = 0;
//     if (axisRank == 0) lowerBndryPad = bndryPad[axis];
//     int upperBndryPad = 0;
//     if (axisRank == axisCommSizes[axis]-1) upperBndryPad = bndryPad[axis];
//     T myDim = localDim + lowerBndryPad + upperBndryPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
//                   2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
//                   bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
//     T myStart = bndryPad[axis] +  axisRank    * localDim;
//     T myStop  = bndryPad[axis] + (axisRank+1) * localDim;
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop() , myStop );
//     globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
//     if (axisRank == 0                    ) myStart -= bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop( ), myStop );
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerBndryPad);
//     TEST_EQUALITY(localBounds.stop(), lowerBndryPad+localDim);
//     TEST_EQUALITY(mdMap.getLowerPad(axis), lowerBndryPad);
//     TEST_EQUALITY(mdMap.getUpperPad(axis), upperBndryPad);
//     TEST_EQUALITY_CONST(mdMap.getCommPadSize(axis), 0);
//     TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
//     int myBndryPad = 0;
//     int axisRank = mdMap.getAxisRank(axis-1);
//     if (axisRank == 0) myBndryPad += bndryPad[axis-1];
//     if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
//       myBndryPad += bndryPad[axis-1];
//     lstrides[axis] = (localDim+myBndryPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
//     gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
//     if (axisRank == 0)
//     {
//       lll += bndryPad[axis] * lstrides[axis];
//       lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
//     }
//     else
//       lur += (localDim-1) * lstrides[axis];
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank == 0) myCommPad = bndryPad[axis];
//     Teuchos::RCP< const Epetra_Map > epetraMap =
//       mdMap.getEpetraAxisMap(axis);
//     TEST_EQUALITY(epetraMap->GID(myCommPad), axisRank*localDim+bndryPad[axis]);
//     TEST_EQUALITY(epetraMap->GID(myCommPad+localDim-1),
//                   bndryPad[axis]+(axisRank+1)*localDim-1);
//   }
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);

//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank  = mdMap.getAxisRank(axis);
//     int myCommPad = 0;
//     if (axisRank == 0) myCommPad = bndryPad[axis];
//     Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap =
//       mdMap.getTpetraAxisMap(axis);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(myCommPad),
//                   axisRank*localDim+bndryPad[axis]);
//     TEST_EQUALITY(tpetraMap->getGlobalElement(myCommPad+localDim-1),
//                   bndryPad[axis]+(axisRank+1)*localDim-1);
//   }
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, paddingConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< Ord > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication and boundary padding
//   Array< int > commPad(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     commPad[axis]  = axis+1;
//     bndryPad[axis] = axis+2;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dims(), commPad(), bndryPad());

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_ASSERT(mdMap.hasPadding());
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank     = mdMap.getAxisRank(axis);
//     int lowerCommPad = commPad[axis];
//     int upperCommPad = commPad[axis];
//     if (axisRank == 0                ) lowerCommPad = bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) upperCommPad = bndryPad[axis];
//     T myDim = localDim + lowerCommPad + upperCommPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
//                   2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
//                   bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
//     T myStart = bndryPad[axis] +  axisRank    * localDim;
//     T myStop  = bndryPad[axis] + (axisRank+1) * localDim;
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop() , myStop );
//     globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
//     if (axisRank == 0                    ) myStart -= bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop( ), myStop );
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
//     TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
//     TEST_EQUALITY(mdMap.getLowerPad(axis), lowerCommPad);
//     TEST_EQUALITY(mdMap.getUpperPad(axis), upperCommPad);
//     TEST_EQUALITY(mdMap.getCommPadSize(axis), commPad[axis]);
//     TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
//     int myCommPad = 0;
//     int axisRank  = mdMap.getAxisRank(axis-1);
//     if (axisRank == 0)
//       myCommPad += bndryPad[axis-1];
//     else
//       myCommPad += commPad[axis-1];
//     if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
//       myCommPad += bndryPad[axis-1];
//     else
//       myCommPad += commPad[axis-1];
//     lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
//     gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
//     if (axisRank == 0)
//     {
//       lll += bndryPad[axis] * lstrides[axis];
//       lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
//     }
//     else
//     {
//       lll += commPad[axis] * lstrides[axis];
//       lur += (localDim+commPad[axis]-1) * lstrides[axis];
//     }
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListPaddingConstructor, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< long long > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication and boundary padding
//   Array< int > commPad(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     commPad[axis]  = axis+1;
//     bndryPad[axis] = axis+2;
//   }

//   // Construct a ParameterList
//   Teuchos::ParameterList plist;
//   plist.set("dimensions"       , dims    );
//   plist.set("communication pad", commPad );
//   plist.set("boundary pad"     , bndryPad);

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, plist);

//   // Perform unit tests of MDMap as a whole
//   TEST_ASSERT(mdMap.onSubcommunicator());
//   TEST_EQUALITY(mdMap.getNumDims(), numDims);
//   TEST_ASSERT(mdMap.hasPadding());
//   TEST_EQUALITY(mdMap.getLayout(), Domi::DEFAULT_ORDER);

//   // Perform unit tests of MDMap axis quantities
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     int axisRank     = mdMap.getAxisRank(axis);
//     int lowerCommPad = commPad[axis];
//     int upperCommPad = commPad[axis];
//     if (axisRank == 0                ) lowerCommPad = bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) upperCommPad = bndryPad[axis];
//     T myDim = localDim + lowerCommPad + upperCommPad;

//     TEST_EQUALITY(mdMap.getAxisCommSize(axis), axisCommSizes[axis]);
//     TEST_ASSERT(not mdMap.isPeriodic(axis));
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalDim(axis,false), dims[axis]               );
//     TEST_EQUALITY_CONST(mdMap.getGlobalBounds(axis,true ).start(), 0);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).start(), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,true ).stop(), dims[axis]+
//                   2*bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getGlobalBounds(axis,false).stop(), dims[axis]+
//                   bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getLocalDim(axis,true ), myDim   );
//     TEST_EQUALITY(mdMap.getLocalDim(axis,false), localDim);
//     Slice globalRankBounds = mdMap.getGlobalRankBounds(axis,false);
//     T myStart = bndryPad[axis] +  axisRank    * localDim;
//     T myStop  = bndryPad[axis] + (axisRank+1) * localDim;
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop() , myStop );
//     globalRankBounds = mdMap.getGlobalRankBounds(axis,true);
//     if (axisRank == 0                    ) myStart -= bndryPad[axis];
//     if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
//     TEST_EQUALITY(globalRankBounds.start(), myStart);
//     TEST_EQUALITY(globalRankBounds.stop( ), myStop );
//     Slice localBounds  = mdMap.getLocalBounds(axis,true);
//     TEST_EQUALITY_CONST(localBounds.start(), 0);
//     TEST_EQUALITY(localBounds.stop(), myDim);
//     localBounds = mdMap.getLocalBounds(axis,false);
//     TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
//     TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
//     TEST_EQUALITY(mdMap.getLowerPad(axis), lowerCommPad);
//     TEST_EQUALITY(mdMap.getUpperPad(axis), upperCommPad);
//     TEST_EQUALITY(mdMap.getCommPadSize(axis), commPad[axis]);
//     TEST_EQUALITY(mdMap.getLowerBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getUpperBndryPad(axis), bndryPad[axis]);
//     TEST_EQUALITY(mdMap.getBndryPadSize(axis), bndryPad[axis]);
//   }

//   // This data will only be used if HAVE_EPETRA or HAVE_TPETRA is defined
//   Teuchos::Array<int> gstrides(numDims);
//   Teuchos::Array<int> lstrides(numDims);
//   gstrides[0] = 1;
//   lstrides[0] = 1;
//   for (int axis=1; axis < numDims; ++axis)
//   {
//     gstrides[axis] = (dims[axis-1]+2*bndryPad[axis-1]) * gstrides[axis-1];
//     int myCommPad = 0;
//     int axisRank  = mdMap.getAxisRank(axis-1);
//     if (axisRank == 0)
//       myCommPad += bndryPad[axis-1];
//     else
//       myCommPad += commPad[axis-1];
//     if (axisRank == mdMap.getAxisCommSize(axis-1)-1)
//       myCommPad += bndryPad[axis-1];
//     else
//       myCommPad += commPad[axis-1];
//     lstrides[axis] = (localDim+myCommPad) * lstrides[axis-1];
//   }
//   int gll = 0;
//   int gur = 0;
//   int lll = 0;
//   int lur = 0;
//   for (int axis=0; axis < numDims; ++axis)
//   {
//     int axisRank = mdMap.getAxisRank(axis);
//     gll += (axisRank*localDim+bndryPad[axis]) * gstrides[axis];
//     gur += ((axisRank+1)*localDim+bndryPad[axis]-1) * gstrides[axis];
//     if (axisRank == 0)
//     {
//       lll += bndryPad[axis] * lstrides[axis];
//       lur += (localDim+bndryPad[axis]-1) * lstrides[axis];
//     }
//     else
//     {
//       lll += commPad[axis] * lstrides[axis];
//       lur += (localDim+commPad[axis]-1) * lstrides[axis];
//     }
//   }

// #ifdef HAVE_EPETRA
//   Teuchos::RCP< const Epetra_Map > epetraMap = mdMap.getEpetraMap();
//   TEST_EQUALITY(epetraMap->GID(lll), gll);
//   TEST_EQUALITY(epetraMap->GID(lur), gur);
// #endif

// #ifdef HAVE_TPETRA
//   Teuchos::RCP< const Tpetra::Map<Ord,Ord> > tpetraMap = mdMap.getTpetraMap();
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lll), gll);
//   TEST_EQUALITY(tpetraMap->getGlobalElement(lur), gur);
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, indexes, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< Ord > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct communication and boundary padding
//   Array< int > commPad(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     commPad[axis]  = axis+1;
//     bndryPad[axis] = axis+2;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dims(), commPad(), bndryPad());

//   // Compute some quantities for testing
//   Array< int > axisRanks(numDims);
//   Array< int > lowerCommPad(numDims);
//   Array< int > upperCommPad(numDims);
//   Array< Ord >   myLocalDims(numDims);
//   Array< Ord >   globalStrides(numDims, 1);
//   Array< Ord >   localStrides(numDims, 1);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     axisRanks[axis]    = mdMap.getAxisRank(axis);
//     lowerCommPad[axis] = mdMap.getLowerPad(axis);
//     upperCommPad[axis] = mdMap.getUpperPad(axis);
//     myLocalDims[axis]  = localDim + lowerCommPad[axis] + upperCommPad[axis];
//     if (axis > 0)
//     {
//       globalStrides[axis] = globalStrides[axis-1] *
//         (dims[axis-1] + 2*bndryPad[axis-1]);
//       localStrides[axis] = localStrides[axis-1] * myLocalDims[axis-1];
//     }
//   }

//   // Unit test globalAxisIndex <-> globalIndex
//   Array< Ord > myGlobalAxisIndex(numDims);
//   T myGlobalIndex = 0;
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     myGlobalAxisIndex[axis] = bndryPad[axis] + dims[axis] / 2;
//     myGlobalIndex += myGlobalAxisIndex[axis] * globalStrides[axis];
//   }
//   TEST_EQUALITY(mdMap.getGlobalAxisIndex(myGlobalIndex), myGlobalAxisIndex);
//   TEST_EQUALITY(mdMap.getGlobalIndex(myGlobalAxisIndex()), myGlobalIndex);

//   // Unit test localAxisIndex <-> localIndex
//   Array< Ord > myLocalAxisIndex(numDims);
//   T myLocalIndex = 0;
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     myLocalAxisIndex[axis] = lowerCommPad[axis] + localDim / 2;
//     myLocalIndex += myLocalAxisIndex[axis] * localStrides[axis];
//   }
//   TEST_EQUALITY(mdMap.getLocalAxisIndex(myLocalIndex), myLocalAxisIndex);
//   TEST_EQUALITY(mdMap.getLocalIndex(myLocalAxisIndex()), myLocalIndex);

//   // Test localIndex <-> globalIndex
//   myGlobalIndex = 0;
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     myGlobalAxisIndex[axis] = myLocalAxisIndex[axis] - lowerCommPad[axis] +
//       axisRanks[axis] * localDim + bndryPad[axis];
//     myGlobalIndex += myGlobalAxisIndex[axis] * globalStrides[axis];
//   }
//   TEST_EQUALITY(mdMap.getGlobalIndex(myLocalIndex), myGlobalIndex);
//   TEST_EQUALITY(mdMap.getLocalIndex(myGlobalIndex), myLocalIndex );
// }
// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, exceptions, Sca, Ord )
// {
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   // Note: axisCommSizes from command line should be fully specified
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct dimensions
//   Ord localDim = 10;
//   Array< Ord > dims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dims[axis] = localDim * mdComm->getAxisCommSize(axis);

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dims());

//   // Unit test methods that should throw exceptions
// #if DOMI_ENABLE_ABC
//   TEST_THROW(mdMap.getAxisCommSize(    -1), Domi::RangeError);
//   TEST_THROW(mdMap.isPeriodic(         -1), Domi::RangeError);
//   TEST_THROW(mdMap.getAxisRank(        -1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerNeighbor(   -1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperNeighbor(   -1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalDim(       -1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalBounds(    -1), Domi::RangeError);
//   TEST_THROW(mdMap.getLocalDim(        -1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalRankBounds(-1), Domi::RangeError);
//   TEST_THROW(mdMap.getLocalBounds(     -1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerPad(        -1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperPad(        -1), Domi::RangeError);
//   TEST_THROW(mdMap.getCommPadSize(     -1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerBndryPad(   -1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperBndryPad(   -1), Domi::RangeError);
//   TEST_THROW(mdMap.getBndryPadSize(    -1), Domi::RangeError);
//   TEST_THROW(mdMap.getAxisCommSize(    numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.isPeriodic(         numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getAxisRank(        numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerNeighbor(   numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperNeighbor(   numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalDim(       numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalBounds(    numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getLocalDim(        numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getGlobalRankBounds(numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getLocalBounds(     numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerPad(        numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperPad(        numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getCommPadSize(     numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getLowerBndryPad(   numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getUpperBndryPad(   numDims+1), Domi::RangeError);
//   TEST_THROW(mdMap.getBndryPadSize(    numDims+1), Domi::RangeError);
// #endif
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapLowerLeft, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions
//   Array< Ord > dimensions(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dimensions[axis] = 10 * axisCommSizes[axis];

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions);
  
//   // Figure out the lower left slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis < 2)
//     {
//       T nd   = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice();
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   for (int axis = 0; axis < numDims; ++axis)
//     if (mdMap.getGlobalBounds(axis).start() >= newDims[axis])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
//       TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapLowerLeftWithCommPad, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments and communication padding
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions and communication padding
//   Array< Ord >   dimensions(numDims);
//   Array< int > commPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = 10 * axisCommSizes[axis];
//     commPad[axis]    = axis + 1;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions, commPad);
  
//   // Figure out the lower left slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis < 2)
//     {
//       T nd   = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice();
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   for (int axis = 0; axis < numDims; ++axis)
//     if (mdMap.getGlobalBounds(axis).start() >= newDims[axis])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
//       TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//       if (axisCommSizes[axis] > 1)
//       {
//         if (subMDMap.getAxisRank(axis) == 0)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), 0);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getLowerPad(axis), commPad[axis]);
//         }
//         if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), 0);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getUpperPad(axis), commPad[axis]);
//         }
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), 0);
//         TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), 0);
//       }
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapLowerRight, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions
//   Array< Ord > dimensions(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dimensions[axis] = 10 * axisCommSizes[axis];

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions);
  
//   // Figure out the lower right slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis == 0)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd,Default);
//     }
//     else if (axis == 1)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice();
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) >= newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis == 0)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapLowerRightWithBndryPad, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions and boundary padding
//   Array< Ord >   dimensions(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = 10 * axisCommSizes[axis];
//     bndryPad[axis]     = axis + 2;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions, ArrayView<int>(), bndryPad);
  
//   // Figure out the lower right slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis == 0)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd+bndryPad[axis],-bndryPad[axis]);
//     }
//     else if (axis == 1)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(bndryPad[axis],nd+bndryPad[axis]);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice(bndryPad[axis], -bndryPad[axis]);
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) >= newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis == 0)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       newDims[axis] + bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       dimensions[axis] + bndryPad[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       newDims[axis] + bndryPad[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapUpperLeft, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions
//   Array< Ord > dimensions(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dimensions[axis] = 10 * axisCommSizes[axis];

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions);
  
//   // Figure out the upper left slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis == 0)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd);
//     }
//     else if (axis == 1)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd,Default);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice();
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) >= newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis == 1)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapUpperLeftPadding, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions, padding
//   Array< Ord > dimensions(numDims);
//   Array< int > commPad(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = 10 * axisCommSizes[axis];
//     commPad[axis]      = axis + 1;
//     bndryPad[axis]     = axis + 2;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions, commPad, bndryPad);

//   // Figure out the upper left slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis == 0)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       if (nc == 0) nc = 1;
//       newCommSizes[axis] = nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(bndryPad[axis],nd+bndryPad[axis]);
//     }
//     else if (axis == 1)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd+bndryPad[axis],-bndryPad[axis]);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice(bndryPad[axis], -bndryPad[axis]);
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) >= newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis == 1)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       newDims[axis] + bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       dimensions[axis] + bndryPad[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       newDims[axis] + bndryPad[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//       if (axisCommSizes[axis] > 1)
//       {
//         if (subMDMap.getAxisRank(axis) == 0)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), 0);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getLowerPad(axis), commPad[axis]);
//         }
//         if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), 0);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getUpperPad(axis), commPad[axis]);
//         }
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), 0);
//         TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), 0);
//       }
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapUpperRight, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions
//   Array< Ord > dimensions(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//     dimensions[axis] = 10 * axisCommSizes[axis];

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions);
  
//   // Figure out the upper right slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     if (axis < 2)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - nc;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd,Default);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice();
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis < 2)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(), newDims[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), dimensions[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getGlobalBounds(axis).start(), 0);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(), newDims[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY_CONST(subMDMap.getTeuchosComm().getRawPtr(), 0);
//     TEST_EQUALITY_CONST(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapUpperRightNewBndryPad, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the dimensions and padding
//   int mySize = 10;
//   Array< Ord > dimensions(numDims);
//   Array< int > commPad(numDims);
//   Array< int > bndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = mySize * axisCommSizes[axis];
//     commPad[axis]      = axis + 1;
//     bndryPad[axis]     = axis + 2;
//   }

//   // Construct an MDMap
//   MDMap< Ord > mdMap(mdComm, dimensions, commPad, bndryPad);

//   // Compute the new boundary padding sizes and figure out the upper
//   // right slice
//   Array< Slice > slices(numDims);
//   Array< int >   newCommSizes(numDims);
//   Array< Ord >     newDims(numDims);
//   Array< int >   newBndryPad(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     // Make the new boundary padding one less than the old boundary padding
//     newBndryPad[axis] = bndryPad[axis] - 1;
//     if (axis < 2)
//     {
//       T   nd = dimensions[axis] / 2;
//       int nc = axisCommSizes[axis] / 2;
//       newCommSizes[axis] = axisCommSizes[axis] - (nd-newBndryPad[axis])/mySize;
//       newDims[axis]      = nd;
//       slices[axis]       = Slice(nd+bndryPad[axis],-bndryPad[axis]);
//     }
//     else
//     {
//       newCommSizes[axis] = axisCommSizes[axis];
//       newDims[axis]      = dimensions[axis];
//       slices[axis]       = Slice(bndryPad[axis], -bndryPad[axis]);
//     }
//   }

//   // Construct the sub-MDMap
//   MDMap< Ord > subMDMap(mdMap, slices, newBndryPad);

//   // Should this processor be a part of the sub-MDComm?
//   bool partOfSubComm = true;
//   if (mdComm->getAxisRank(0) < axisCommSizes[0] - newCommSizes[0])
//     partOfSubComm = false;
//   if (numDims > 1)
//     if (mdComm->getAxisRank(1) < axisCommSizes[1] - newCommSizes[1])
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_ASSERT(subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       TEST_EQUALITY(subMDMap.getAxisCommSize(axis), newCommSizes[axis]);
//       TEST_EQUALITY(subMDMap.getGlobalDim(axis), newDims[axis]);
//       if (axis < 2)
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       newDims[axis] + bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       dimensions[axis] + bndryPad[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).start(),
//                       bndryPad[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalBounds(axis).stop(),
//                       newDims[axis] + bndryPad[axis]);
//       }
//       TEST_ASSERT(not subMDMap.isPeriodic(axis));
//       if (axisCommSizes[axis] > 1)
//       {
//         if (subMDMap.getAxisRank(axis) == 0)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), newBndryPad[axis]);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getLowerPad(axis), commPad[axis]);
//         }
//         if (subMDMap.getAxisRank(axis) == subMDMap.getAxisCommSize(axis)-1)
//         {
//           TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), newBndryPad[axis]);
//         }
//         else
//         {
//           TEST_EQUALITY(subMDMap.getUpperPad(axis), commPad[axis]);
//         }
//       }
//       else
//       {
//         TEST_EQUALITY_CONST(subMDMap.getLowerPad(axis), newBndryPad[axis]);
//         TEST_EQUALITY_CONST(subMDMap.getUpperPad(axis), newBndryPad[axis]);
//       }
//     }
//   }
//   else
//   {
//     TEST_ASSERT(not subMDMap.onSubcommunicator());
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//     TEST_THROW(subMDMap.getAxisCommSize(0) , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.isPeriodic(0)      , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getAxisRank(0)     , Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getLowerNeighbor(0), Domi::SubcommunicatorError);
//     TEST_THROW(subMDMap.getUpperNeighbor(0), Domi::SubcommunicatorError);
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapReduce, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

//   // Construct the parent MDMap
//   Array< Ord > dimensions(numDims);
//   Array<int> commPad(numDims);
//   Array<int> bndryPad(numDims);
//   int localSize = 10;
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = localSize*axisCommSizes[axis];
//     commPad[axis] = axis + 1;
//     bndryPad[axis] = axis + 2;
//   }
//   MDMap<T> mdMap(mdComm, dimensions, commPad, bndryPad);

//   // We will reduce this parent MDMap several times by using the
//   // single GlobalOrd constructor along each dimension
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     T myOrd = bndryPad[axis] + dimensions[axis] / 2;
//     int newDims = (numDims > 1) ? numDims-1 : numDims;
//     Slice bounds = mdMap.getGlobalRankBounds(axis);
//     bool partOfSubMap = true;
//     if ((myOrd < bounds.start()) || (bounds.stop() <= myOrd))
//       partOfSubMap = false;
//     MDMap<T> reducedMdMap(mdMap, axis, myOrd);
//     if (partOfSubMap)
//     {
//       TEST_ASSERT(reducedMdMap.onSubcommunicator());
//       TEST_EQUALITY(reducedMdMap.getNumDims(), newDims);
//       for (int newAxis = 0; newAxis < newDims; ++newAxis)
//       {
//         int axisRank = reducedMdMap.getAxisRank(newAxis);
//         if (numDims == 1)
//         {
//           TEST_EQUALITY_CONST(reducedMdMap.getGlobalDim(newAxis), 1);
//           Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), myOrd  );
//           TEST_EQUALITY(bounds.stop() , myOrd+1);
//           bounds = reducedMdMap.getGlobalRankBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), myOrd  );
//           TEST_EQUALITY(bounds.stop() , myOrd+1);
//           TEST_EQUALITY_CONST(reducedMdMap.getLocalDim(newAxis), 1);
//           bounds = reducedMdMap.getLocalBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), 0);
//           TEST_EQUALITY(bounds.stop(),  1);
//           TEST_ASSERT(not reducedMdMap.hasPadding());
//           TEST_EQUALITY_CONST(reducedMdMap.getLowerPad(newAxis), 0);
//           TEST_EQUALITY_CONST(reducedMdMap.getUpperPad(newAxis), 0);
//           TEST_EQUALITY_CONST(reducedMdMap.getCommPadSize(newAxis), 0);
//           TEST_EQUALITY_CONST(reducedMdMap.getLowerBndryPad(newAxis), 0);
//           TEST_EQUALITY_CONST(reducedMdMap.getUpperBndryPad(newAxis), 0);
//           TEST_EQUALITY_CONST(reducedMdMap.getBndryPadSize(newAxis), 0);
//         }
//         else if (newAxis < axis)
//         {
//           TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
//                         dimensions[newAxis]);
//           Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), bndryPad[newAxis]                    );
//           TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+dimensions[newAxis]);
//           bounds = reducedMdMap.getGlobalRankBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), bndryPad[newAxis]+localSize* axisRank   );
//           TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+localSize*(axisRank+1));
//           TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
//           bounds = reducedMdMap.getLocalBounds(newAxis);
//           TEST_ASSERT(reducedMdMap.hasPadding());
//           if (reducedMdMap.getAxisRank(newAxis) == 0)
//           {
//             TEST_EQUALITY(bounds.start(), bndryPad[newAxis]          );
//             TEST_EQUALITY(bounds.stop() , bndryPad[newAxis]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis), bndryPad[newAxis]);
//           }
//           else if (reducedMdMap.getAxisRank(newAxis) ==
//                    reducedMdMap.getAxisCommSize(newAxis)-1)
//           {
//             TEST_EQUALITY(bounds.start(), commPad[newAxis]          );
//             TEST_EQUALITY(bounds.stop() , commPad[newAxis]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis), commPad[newAxis]);
//             TEST_EQUALITY(reducedMdMap.getUpperPad(newAxis), bndryPad[newAxis]);
//           }
//           else
//           {
//             TEST_EQUALITY(bounds.start(), commPad[newAxis]          );
//             TEST_EQUALITY(bounds.stop() , commPad[newAxis]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis), commPad[newAxis]);
//             TEST_EQUALITY(reducedMdMap.getUpperPad(newAxis), commPad[newAxis]);
//           }
//           TEST_EQUALITY(reducedMdMap.getLowerBndryPad(newAxis), bndryPad[newAxis]);
//           TEST_EQUALITY(reducedMdMap.getUpperBndryPad(newAxis), bndryPad[newAxis]);
//           TEST_EQUALITY(reducedMdMap.getBndryPadSize(newAxis) , bndryPad[newAxis]);
//         }
//         else
//         {
//           TEST_EQUALITY(reducedMdMap.getGlobalDim(newAxis),
//                         dimensions[newAxis+1]);
//           Slice bounds = reducedMdMap.getGlobalBounds(newAxis);
//           TEST_EQUALITY(bounds.start(), bndryPad[newAxis+1]);
//           TEST_EQUALITY(bounds.stop() , bndryPad[newAxis+1]+
//                         dimensions[newAxis+1]);
//           bounds = reducedMdMap.getGlobalRankBounds(newAxis);
//           TEST_EQUALITY(bounds.start(),
//                         bndryPad[newAxis+1]+localSize* axisRank   );
//           TEST_EQUALITY(bounds.stop() ,
//                         bndryPad[newAxis+1]+localSize*(axisRank+1));
//           TEST_EQUALITY(reducedMdMap.getLocalDim(newAxis), localSize);
//           bounds = reducedMdMap.getLocalBounds(newAxis);
//           TEST_ASSERT(reducedMdMap.hasPadding());
//           if (reducedMdMap.getAxisRank(newAxis) == 0)
//           {
//             TEST_EQUALITY(bounds.start(), bndryPad[newAxis+1]          );
//             TEST_EQUALITY(bounds.stop() , bndryPad[newAxis+1]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis),
//                           bndryPad[newAxis+1]);
//           }
//           else if (reducedMdMap.getAxisRank(newAxis) ==
//                    reducedMdMap.getAxisCommSize(newAxis)-1)
//           {
//             TEST_EQUALITY(bounds.start(), commPad[newAxis+1]          );
//             TEST_EQUALITY(bounds.stop() , commPad[newAxis+1]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis),commPad[newAxis+1]);
//             TEST_EQUALITY(reducedMdMap.getUpperPad(newAxis),bndryPad[newAxis+1]);
//           }
//           else
//           {
//             TEST_EQUALITY(bounds.start(), commPad[newAxis+1]          );
//             TEST_EQUALITY(bounds.stop() , commPad[newAxis+1]+localSize);
//             TEST_EQUALITY(reducedMdMap.getLowerPad(newAxis), commPad[newAxis+1]);
//             TEST_EQUALITY(reducedMdMap.getUpperPad(newAxis), commPad[newAxis+1]);
//           }
//           TEST_EQUALITY(reducedMdMap.getLowerBndryPad(newAxis), bndryPad[newAxis+1]);
//           TEST_EQUALITY(reducedMdMap.getUpperBndryPad(newAxis), bndryPad[newAxis+1]);
//           TEST_EQUALITY(reducedMdMap.getBndryPadSize(newAxis) , bndryPad[newAxis]+1);
//         }
//       }
//     }
//     else
//     {
//       TEST_ASSERT(not reducedMdMap.onSubcommunicator());
//     }
//   }
// }

// TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, subMapPeriodic, Sca, Ord )
// {
//   // Construct the MDComm from command-line arguments
//   TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
//   Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
//   // Construct the periodic flags
//   Array< int > periodic(numDims, 0);
//   periodic[0] = 1;
//   MDCommRCP mdComm =
//     Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes, periodic));
  
//   // Figure out the lower slice
//   Array< Slice > slices(numDims);  
//   Array< Ord >     dimensions(numDims);
//   for (int axis = 0; axis < numDims; ++axis)
//   {
//     dimensions[axis] = 10 * axisCommSizes[axis];
//     if (axis == 1)
//       slices[axis] = Slice(dimensions[axis]/2);
//     else
//       slices[axis] = Slice();
//   }

//   // Construct the MDMap and the sub-MDMap
//   MDMap<T> mdMap(mdComm, dimensions);
//   MDMap<T> subMDMap(mdMap, slices);

//   // Should this processor be a part of the sub-MDMap?
//   bool partOfSubComm = true;
//   if (mdComm->getNumDims() > 1)
//     if (mdComm->getAxisRank(1) > (mdComm->getAxisCommSize(1)-1)/2)
//       partOfSubComm = false;

//   // Do some unit tests
//   if (partOfSubComm)
//   {
//     TEST_EQUALITY(subMDMap.getNumDims(), numDims);
//     for (int axis = 0; axis < numDims; ++axis)
//     {
//       if (axis == 1)
//       {
//         TEST_EQUALITY(subMDMap.getAxisCommSize(axis), (axisCommSizes[axis]+1)/2);
//         TEST_EQUALITY(subMDMap.getGlobalDim(axis), 5*axisCommSizes[axis]);
//       }
//       else
//       {
//         TEST_EQUALITY(subMDMap.getAxisCommSize(axis), axisCommSizes[axis]);
//         TEST_EQUALITY(subMDMap.getGlobalDim(axis), 10*axisCommSizes[axis]);
//       }
//       TEST_EQUALITY(subMDMap.isPeriodic(axis), (axis == 0));
//     }
//   }
//   else
//   {
//     TEST_EQUALITY(subMDMap.getNumDims(), 0);
//   }
// }

//
// Instantiations
//

#define UNIT_TEST_GROUP( Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, dimensionsConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, randomize, Sca, Ord ) // \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListDimensionsConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, commPadConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListCommPadConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, bndryPadConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListBndryPadConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, paddingConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListPaddingConstructor, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, indexes, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, exceptions, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapLowerLeft, Sca, Ord )   \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapLowerLeftWithCommPad, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapLowerRight, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapLowerRightWithBndryPad, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapUpperLeft, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapUpperLeftPadding, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapUpperRight, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapUpperRightNewBndryPad, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapReduce, Sca, Ord ) \
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, subMapPeriodic, Sca, Ord )

UNIT_TEST_GROUP(double, int)
#if 0
UNIT_TEST_GROUP(float, int)
UNIT_TEST_GROUP(float, long)
UNIT_TEST_GROUP(double, long)
#endif

}  // namespace
