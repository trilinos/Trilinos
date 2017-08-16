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

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("numDims" , &numDims    , "number of dimensions");
  clp.setOption("commDims", &commDimsStr, "comma-separated list of "
                "number of processors along each axis");
}

//
// Templated Unit Tests
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, dimensionsConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 4;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims()));
  MDVector< Sca > mdVector(mdMap);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_ASSERT(not mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    TEST_EQUALITY(mdVector.getCommDim(axis), commDims[axis]);
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

  // Test dot product and norms for non-ordinal types
  if (! Teuchos::ScalarTraits< Sca >::isOrdinal)
  {
    mdVector.putScalar(scalar);
    Sca norm1 = scalar * Domi::computeSize(dims);
    Sca dot   = scalar * norm1 ;
    typename Teuchos::ScalarTraits< Sca >::magnitudeType norm2 =
      Teuchos::ScalarTraits< Sca >::squareroot(dot);
    typename Teuchos::ScalarTraits< Sca >::magnitudeType tolerance =
      Sca(1.0e-7);
    TEST_COMPARE(std::abs(mdVector.dot(mdVector) - dot   ), <, tolerance);
    TEST_COMPARE(std::abs(mdVector.norm1()       - norm1 ), <, tolerance);
    TEST_COMPARE(std::abs(mdVector.norm2()       - norm2 ), <, tolerance);
    TEST_COMPARE(std::abs(mdVector.normInf()     - scalar), <, tolerance);
    TEST_COMPARE(std::abs(mdVector.normWeighted(mdVector) -
                          std::pow(scalar, 3./2.)        ), <, tolerance);
    TEST_COMPARE(std::abs(mdVector.meanValue()   - scalar), <, tolerance);
  }

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, initializationConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims()));

  // Construct an MDArrayView with initialization values
  typedef typename Domi::size_type size_type;
  Array< dim_type > initDims;
  Sca scalar = 2;
  for (int axis = 0; axis < numDims; ++axis)
    initDims.push_back(mdMap->getLocalDim(axis,true));
  MDArray< Sca > initVals(initDims(), scalar);

  // Construct an MDVector using the initialization values
  MDVector< Sca > mdVector(mdMap, initVals);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_ASSERT(not mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    TEST_EQUALITY(mdVector.getCommDim(axis), commDims[axis]);
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, copyConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 9;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap and a source MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims()));
  MDVector< Sca > source(mdMap);
  Sca scalar = 27;
  source.putScalar(scalar);

  // Construct an MDVector using the copy constructor
  MDVector< Sca > mdVector(source);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_ASSERT(not mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    TEST_EQUALITY(mdVector.getCommDim(axis), commDims[axis]);
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, pListDimensionsConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);

  // Construct dummy dimensions
  Array< dim_type > dims(numDims, comm->getSize());

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("comm dimensions", commDims);
  plist.set("dimensions"     , dims    );

  // Get the actual communicator dimensions
  Teuchos::Array< int > actualCommDims =
    Domi::regularizeCommDims(comm->getSize(), plist);

  // Compute actual dimensions
  dim_type localDim = 10;
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * actualCommDims[axis];
  plist.set("dimensions", dims);

  // Construct an MDVector
  MDVector< Sca > mdVector(comm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_ASSERT(not mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    TEST_EQUALITY(mdVector.getCommDim(axis), actualCommDims[axis]);
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, pListCommPadConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct communication padding
  Array< int > commPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"             , dims   );
  plist.set("communication pad sizes", commPad);

  // Construct an MDVector
  MDVector< Sca > mdVector(mdComm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_EQUALITY(mdVector.hasPadding(), comm->getSize() > 1);
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    int lowerCommPad = 0;
    if (axisRank > 0) lowerCommPad = commPad[axis];
    int upperCommPad = 0;
    if (axisRank < commDims[axis]-1) upperCommPad = commPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdVector.getCommDim(axis), commDims[axis]);
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, pListBndryPadConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);

  // Get the actual communicator dimensions
  Teuchos::Array< int > actualCommDims =
    Domi::regularizeCommDims(comm->getSize(), numDims, commDims);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * actualCommDims[axis];

  // Construct boundary padding
  Array< int > bndryPad(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    bndryPad[axis] = axis+1;

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("comm dimensions"   , commDims);
  plist.set("dimensions"        , dims    );
  plist.set("boundary pad sizes", bndryPad);

  // Construct an MDVector
  MDVector< Sca > mdVector(comm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  bool hasPadding = false;
  for (int axis = 0; axis < numDims; axis++)
  {
    if (mdVector.getCommIndex(axis)==0                     ) hasPadding = true;
    if (mdVector.getCommIndex(axis)==actualCommDims[axis]-1) hasPadding = true;
  }
  TEST_EQUALITY(mdVector.hasPadding(), hasPadding);
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank = mdVector.getCommIndex(axis);
    int lowerBndryPad = 0;
    if (axisRank == 0) lowerBndryPad = bndryPad[axis];
    int upperBndryPad = 0;
    if (axisRank == actualCommDims[axis]-1) upperBndryPad = bndryPad[axis];
    dim_type myDim = localDim + lowerBndryPad + upperBndryPad;

    TEST_EQUALITY(mdVector.getCommDim(axis), actualCommDims[axis]);
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
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdVector.getGlobalRankBounds(axis,true);
    if (axisRank == 0                     ) myStart -= bndryPad[axis];
    if (axisRank == actualCommDims[axis]-1) myStop  += bndryPad[axis];
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

    Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
      mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, pListPaddingConstructor, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

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
  plist.set("dimensions"             , dims    );
  plist.set("communication pad sizes", commPad );
  plist.set("boundary pad sizes"     , bndryPad);

  // Construct an MDVector
  MDVector< Sca > mdVector(mdComm, plist);

  // Perform unit tests of MDVector as a whole
  TEST_ASSERT(mdVector.onSubcommunicator());
  TEST_EQUALITY(mdVector.numDims(), numDims);
  TEST_ASSERT(mdVector.hasPadding());
  TEST_EQUALITY(mdVector.getLayout(), Domi::DEFAULT_ORDER);

  // Perform unit tests of MDVector axis quantities
  for (int axis = 0; axis < numDims; ++axis)
  {
    int axisRank     = mdVector.getCommIndex(axis);
    int lowerCommPad = commPad[axis];
    int upperCommPad = commPad[axis];
    if (axisRank == 0               ) lowerCommPad = bndryPad[axis];
    if (axisRank == commDims[axis]-1) upperCommPad = bndryPad[axis];
    dim_type myDim = localDim + lowerCommPad + upperCommPad;

    TEST_EQUALITY(mdVector.getCommDim(axis), commDims[axis]);
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
    dim_type myStart = bndryPad[axis] +  axisRank    * localDim;
    dim_type myStop  = bndryPad[axis] + (axisRank+1) * localDim;
    TEST_EQUALITY(globalRankBounds.start(), myStart);
    TEST_EQUALITY(globalRankBounds.stop() , myStop );
    globalRankBounds = mdVector.getGlobalRankBounds(axis,true);
    if (axisRank == 0               ) myStart -= bndryPad[axis];
    if (axisRank == commDims[axis]-1) myStop  += bndryPad[axis];
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

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector > epetraVectorView =
      mdVector.getEpetraVectorView();
  }
  else
    TEST_THROW(mdVector.getEpetraVectorView(), Domi::TypeError);

  Teuchos::RCP< Epetra_Vector > epetraVectorCopy =
    mdVector.getEpetraVectorCopy();
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorView =
    mdVector.template getTpetraVectorView< int, int >();
  Teuchos::RCP< Tpetra::Vector< Sca, int, int > > tpetraVectorCopy =
    mdVector.template getTpetraVectorCopy< int, int >();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, augmentedConstruction, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap
  Teuchos::RCP< const MDMap<> > mdMap = Teuchos::rcp(new MDMap<>(mdComm, dims));

  // Construct and test an MDVector with a leading dimension
  MDVector< Sca > mdv1(mdMap, 3);
  TEST_EQUALITY(mdv1.numDims()      , numDims+1);
  TEST_EQUALITY(mdv1.getCommDim(0)  , 1        );
  TEST_EQUALITY(mdv1.getGlobalDim(0), 3        );

  // Construct and test an MDVector with a trailing dimension
  MDVector< Sca > mdv2(mdMap, 0, 2);
  TEST_EQUALITY(mdv2.numDims()            , numDims+1);
  TEST_EQUALITY(mdv2.getCommDim(numDims)  , 1        );
  TEST_EQUALITY(mdv2.getGlobalDim(numDims), 2        );

  // Construct a ParameterList
  Teuchos::ParameterList plist;
  plist.set("dimensions"        , dims);
  plist.set("trailing dimension", 5   );

  // Construct and test an MDVector with a trailing dimension
  MDVector< Sca > mdv3(mdComm, plist);
  TEST_EQUALITY(mdv3.numDims()            , numDims+1);
  TEST_EQUALITY(mdv3.getCommDim(numDims)  , 1        );
  TEST_EQUALITY(mdv3.getGlobalDim(numDims), 5        );

#ifdef HAVE_EPETRA
  if (typeid(Sca) == typeid(double))
  {
    Teuchos::RCP< Epetra_Vector      > ev1  = mdv1.getEpetraVectorView();
    Teuchos::RCP< Epetra_MultiVector > emv2 = mdv2.getEpetraMultiVectorView();
    Teuchos::RCP< Epetra_MultiVector > emv3 = mdv3.getEpetraMultiVectorView();
    TEST_EQUALITY_CONST(emv2->NumVectors(), 2);
    TEST_EQUALITY_CONST(emv3->NumVectors(), 5);
  }
  else
  {
    TEST_THROW(mdv1.getEpetraVectorView()     , Domi::TypeError);
    TEST_THROW(mdv2.getEpetraMultiVectorView(), Domi::TypeError);
    TEST_THROW(mdv2.getEpetraMultiVectorView(), Domi::TypeError);
  }

  Teuchos::RCP< Epetra_Vector      > ev1c  = mdv1.getEpetraVectorCopy();
  Teuchos::RCP< Epetra_MultiVector > emv2c = mdv2.getEpetraMultiVectorCopy();
  Teuchos::RCP< Epetra_MultiVector > emv3c = mdv3.getEpetraMultiVectorCopy();
  TEST_EQUALITY_CONST(emv2c->NumVectors(), 2);
  TEST_EQUALITY_CONST(emv3c->NumVectors(), 5);
#endif

#ifdef HAVE_TPETRA
  Teuchos::RCP< Tpetra::Vector< Sca, int >      > tv1  =
    mdv1.template getTpetraVectorView< int >();
  Teuchos::RCP< Tpetra::MultiVector< Sca, int > > tmv2 =
    mdv2.template getTpetraMultiVectorView< int >();
  Teuchos::RCP< Tpetra::MultiVector< Sca, int > > tmv3 =
    mdv3.template getTpetraMultiVectorView< int >();
  TEST_EQUALITY_CONST(tmv2->getNumVectors(), 2);
  TEST_EQUALITY_CONST(tmv3->getNumVectors(), 5);

  Teuchos::RCP< Tpetra::Vector< Sca, int >      > tv1c  =
    mdv1.template getTpetraVectorCopy< int >();
  Teuchos::RCP< Tpetra::MultiVector< Sca, int > > tmv2c =
    mdv2.template getTpetraMultiVectorCopy< int >();
  Teuchos::RCP< Tpetra::MultiVector< Sca, int > > tmv3c =
    mdv3.template getTpetraMultiVectorCopy< int >();
  TEST_EQUALITY_CONST(tmv2c->getNumVectors(), 2);
  TEST_EQUALITY_CONST(tmv3c->getNumVectors(), 5);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, randomize, Sca )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, commDims));

  // Ensure that the commDims are completely specified
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm->getCommDim(axis);

  // Construct dimensions
  dim_type localDim = 10;
  Array< dim_type > dims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    dims[axis] = localDim * mdComm->getCommDim(axis);

  // Construct an MDMap and MDVector
  typedef Teuchos::RCP< MDMap<> > MDMapRCP;
  MDMapRCP mdMap = rcp(new MDMap<>(mdComm, dims()));
  MDVector< Sca > mdVector(mdMap, true);

  // Test that all values are zero
  typedef typename MDArrayView< const Sca >::const_iterator const_iterator;
  MDArrayView< const Sca > view = mdVector.getData();
  for (const_iterator it = view.cbegin(); it != view.cend(); ++it)
    TEST_EQUALITY_CONST(*it, 0);

  // Test the randomize method
  mdVector.randomize();
  Sca rand_min = -1;
  Sca rand_max =  1;
  if (Teuchos::ScalarTraits< Sca >::isOrdinal)
  {
    rand_min = 0;
    rand_max = RAND_MAX;
  }
  int count = 0;
  for (const_iterator it = view.cbegin(); it != view.cend(); ++it)
  {
    if (*it != 0) ++count;
    TEST_COMPARE(*it, >=, rand_min);
    TEST_COMPARE(*it, <=, rand_max);
  }
  TEST_COMPARE(count, >, 0);
  
}

////////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, dimensionsConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, initializationConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, copyConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, pListDimensionsConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, pListCommPadConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, pListBndryPadConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, pListPaddingConstructor, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, augmentedConstruction, Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, randomize, Sca )

UNIT_TEST_GROUP(double)
#if 1
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(long)
#endif

}  // namespace
