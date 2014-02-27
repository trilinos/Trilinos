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
using Domi::MDArray;
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
    TEST_EQUALITY_CONST(mdVector.getLowerPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperPadSize(axis), 0);
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
  Sca interiorScalar = 9;
  mdVector.putScalar(interiorScalar, false);  // false = exclude padding
  cData = mdVector.getData(false);
  for (typename MDArrayView< Sca >::const_iterator it = cData.cbegin();
       it != cData.cend(); ++it)
    TEST_EQUALITY(*it, interiorScalar);
  for (int axis = 0; axis < numDims; ++axis)
  {
    cData = mdVector.getLowerPadData(axis);
    for (typename MDArrayView< Sca >::const_iterator it = cData.cbegin();
         it != cData.cend(); ++it)
      TEST_EQUALITY(*it, scalar);
    cData = mdVector.getUpperPadData(axis);
    for (typename MDArrayView< Sca >::const_iterator it = cData.cbegin();
         it != cData.cend(); ++it)
      TEST_EQUALITY(*it, scalar);
  }

  // Test dot product and norms
  mdVector.putScalar(scalar);
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

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, initializationConstructor, Sca, Ord )
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

  // Construct an MDMap
  typedef Teuchos::RCP< MDMap< Ord > > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap< Ord >(mdComm, dims()));

  // Construct an MDArrayView with initialization values
  typedef typename MDArray< Sca >::size_type size_type;
  Array< size_type > initDims;
  Sca scalar = 2;
  for (int axis = 0; axis < numDims; ++axis)
    initDims.push_back(mdMap->getLocalDim(axis,true));
  MDArray< Sca > initVals(initDims(), scalar);

  // Construct an MDVector using the initialization values
  MDVector< Sca, Ord > mdVector(mdMap, initVals);

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
    TEST_EQUALITY_CONST(mdVector.getLowerPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getBndryPadSize(axis), 0);
  }

  // Perform tests of MDVector data
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
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, copyConstructor, Sca, Ord )
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
  Ord localDim = 9;
  Array< Ord > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct an MDMap and a source MDVector
  typedef Teuchos::RCP< MDMap< Ord > > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap< Ord >(mdComm, dims()));
  MDVector< Sca, Ord > source(mdMap);
  Sca scalar = 27;
  source.putScalar(scalar);

  // Construct an MDVector using the copy constructor
  MDVector< Sca, Ord > mdVector(source);

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
    TEST_EQUALITY_CONST(mdVector.getLowerPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getBndryPadSize(axis), 0);
  }

  // Perform tests of MDVector data
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
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListDimensionsConstructor, Sca, Ord )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);

  // Check that the axisCommSizes are completely specified
  TEST_EQUALITY(axisCommSizes.size(), numDims)
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_ASSERT(axisCommSizes[axis] > 0);
  }

  // Construct dimensions
  Ord localDim = 10;
  Array< long int > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * axisCommSizes[axis];

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("axis comm sizes", axisCommSizes);
  plist.set("dimensions"     , dims         );

  // Construct an MDVector
  MDVector< Sca, Ord > mdVector(comm, plist);

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
    TEST_EQUALITY_CONST(mdVector.getLowerPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getBndryPadSize(axis), 0);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListCommPadConstructor, Sca, Ord )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  Ord localDim = 10;
  Array< long int > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct communication padding
  Array< int > commPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"       , dims   );
  plist.set("communication pad", commPad);

  // Construct an MDVector
  MDVector< Sca, Ord > mdVector(mdComm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.getNumDims(), numDims);
  TEST_EQUALITY(mdVector.hasPadding(), comm->getSize() > 1);
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getAxisRank(axis);
    int lowerCommPad = 0;
    if (axisRank > 0) lowerCommPad = commPad[axis];
    int upperCommPad = 0;
    if (axisRank < axisCommSizes[axis]-1) upperCommPad = commPad[axis];
    Ord myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdVector.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdVector.isPeriodic(axis));
    TEST_EQUALITY(mdVector.getGlobalDim(axis), dims[axis]);
    TEST_EQUALITY_CONST(mdVector.getGlobalBounds(axis).start(), 0);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis).stop(), dims[axis]);
    TEST_EQUALITY(mdVector.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdVector.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdVector.getGlobalRankBounds(axis);
    TEST_EQUALITY(globalRankBounds.start(),  axisRank   *localDim);
    TEST_EQUALITY(globalRankBounds.stop() , (axisRank+1)*localDim);
    Slice localBounds  = mdVector.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdVector.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY_CONST(mdVector.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY_CONST(mdVector.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY_CONST(mdVector.getLowerBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getUpperBndryPad(axis), 0);
    TEST_EQUALITY_CONST(mdVector.getBndryPadSize(axis), 0);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListBndryPadConstructor, Sca, Ord )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);

  // Construct dimensions
  Ord localDim = 10;
  Array< long int > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * axisCommSizes[axis];

  // Construct boundary padding
  Array< int > bndryPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    bndryPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("axis comm sizes", axisCommSizes);
  plist.set("dimensions"     , dims         );
  plist.set("boundary pad"   , bndryPad     );

  // Construct an MDVector
  MDVector< Sca, Ord > mdVector(comm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.getNumDims(), numDims);
  TEST_ASSERT(mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getAxisRank(axis);
    int lowerBndryPad = 0;
    if (axisRank == 0) lowerBndryPad = bndryPad[axis];
    int upperBndryPad = 0;
    if (axisRank == axisCommSizes[axis]-1) upperBndryPad = bndryPad[axis];
    Ord myDim = localDim + lowerBndryPad + upperBndryPad;

    TEST_EQUALITY(mdVector.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdVector.isPeriodic(axis));
    TEST_EQUALITY(mdVector.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY_CONST(mdVector.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdVector.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdVector.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdVector.getGlobalRankBounds(axis,false);
    Ord myStart = bndryPad[axis] +  axisRank    * localDim;
    Ord myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdVector.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= bndryPad[axis];
    if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdVector.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdVector.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerBndryPad);
    TEST_EQUALITY(localBounds.stop(), lowerBndryPad+localDim);
    TEST_EQUALITY(mdVector.getLowerPadSize(axis), lowerBndryPad);
    TEST_EQUALITY(mdVector.getUpperPadSize(axis), upperBndryPad);
    TEST_EQUALITY_CONST(mdVector.getCommPadSize(axis), 0);
    TEST_EQUALITY(mdVector.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getBndryPadSize(axis), bndryPad[axis]);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MDVector, pListPaddingConstructor, Sca, Ord )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Note: axisCommSizes from command line should be fully specified
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Construct dimensions
  Ord localDim = 10;
  Array< long int > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getAxisCommSize(axis);

  // Construct communication and boundary padding
  Array< int > commPad(numDims);
  Array< int > bndryPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    commPad[axis]  = axis+1;
    bndryPad[axis] = axis+2;
  }

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"       , dims    );
  plist.set("communication pad", commPad );
  plist.set("boundary pad"     , bndryPad);

  // Construct an MDVector
  MDVector< Sca, Ord > mdVector(mdComm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.getNumDims(), numDims);
  TEST_ASSERT(mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank     = mdVector.getAxisRank(axis);
    int lowerCommPad = commPad[axis];
    int upperCommPad = commPad[axis];
    if (axisRank == 0                ) lowerCommPad = bndryPad[axis];
    if (axisRank == axisCommSizes[axis]-1) upperCommPad = bndryPad[axis];
    Ord myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdVector.getAxisCommSize(axis), axisCommSizes[axis]);
    TEST_ASSERT(not mdVector.isPeriodic(axis));
    TEST_EQUALITY(mdVector.getGlobalDim(axis,true ), dims[axis]+2*bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalDim(axis,false), dims[axis]               );
    TEST_EQUALITY_CONST(mdVector.getGlobalBounds(axis,true ).start(), 0);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,false).start(), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,true ).stop(), dims[axis]+
                  2*bndryPad[axis]);
    TEST_EQUALITY(mdVector.getGlobalBounds(axis,false).stop(), dims[axis]+
                  bndryPad[axis]);
    TEST_EQUALITY(mdVector.getLocalDim(axis,true ), myDim   );
    TEST_EQUALITY(mdVector.getLocalDim(axis,false), localDim);
    Slice globalRankBounds = mdVector.getGlobalRankBounds(axis,false);
    Ord myStart = bndryPad[axis] +  axisRank    * localDim;
    Ord myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdVector.getGlobalRankBounds(axis,true);
    if (axisRank == 0                    ) myStart -= bndryPad[axis];
    if (axisRank == axisCommSizes[axis]-1) myStop  += bndryPad[axis];
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop( ), myStop );
    Slice localBounds  = mdVector.getLocalBounds(axis,true);
    TEST_EQUALITY_CONST(localBounds.start(), 0);
    TEST_EQUALITY(localBounds.stop(), myDim);
    localBounds = mdVector.getLocalBounds(axis,false);
    TEST_EQUALITY_CONST(localBounds.start(), lowerCommPad);
    TEST_EQUALITY(localBounds.stop(), lowerCommPad+localDim);
    TEST_EQUALITY(mdVector.getLowerPadSize(axis), lowerCommPad);
    TEST_EQUALITY(mdVector.getUpperPadSize(axis), upperCommPad);
    TEST_EQUALITY(mdVector.getCommPadSize(axis), commPad[axis]);
    TEST_EQUALITY(mdVector.getLowerBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getUpperBndryPad(axis), bndryPad[axis]);
    TEST_EQUALITY(mdVector.getBndryPadSize(axis), bndryPad[axis]);
  }
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

////////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, dimensionsConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, initializationConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, copyConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListDimensionsConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListCommPadConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListBndryPadConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, pListPaddingConstructor, Sca, Ord ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MDVector, randomize, Sca, Ord )

UNIT_TEST_GROUP(double, int)
#if 0
UNIT_TEST_GROUP(float, int)
UNIT_TEST_GROUP(float, long)
UNIT_TEST_GROUP(double, long)
#endif

}  // namespace
