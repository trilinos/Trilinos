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

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDComm.hpp"
#include "Domi_Exceptions.hpp"

// Teuchos includes
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

namespace
{

using std::string;
using Teuchos::Array;
using Domi::TeuchosCommRCP;
using Domi::MDComm;
using Domi::MDCommRCP;
using Domi::Slice;
typedef Domi::Ordinal Ordinal;
const Ordinal & Default = Domi::Slice::Default;

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

#if 1
TEUCHOS_UNIT_TEST( MDComm, regularizeAxisSizes )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Array< int > axisCommSizesActual =
    Domi::regularizeAxisSizes(comm->getSize(), numDims, axisCommSizes);
  
  TEST_EQUALITY(axisCommSizesActual.size(), numDims);

  int numProcs = 1;
  for (int axis = 0; axis < numDims; ++axis)
    numProcs *= axisCommSizesActual[axis];

  TEST_EQUALITY(numProcs, comm->getSize());
}

TEUCHOS_UNIT_TEST( MDComm, axisCommSizesConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // If axisCommSizes is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = axisCommSizes.size(); axis < numDims; ++axis)
  {
    axisCommSizes.push_back(-1);
  }
  MDComm mdComm(comm, axisCommSizes);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Comm > epetraComm = mdComm.getEpetraComm();
  TEST_EQUALITY(epetraComm->NumProc(), comm->getSize());
  TEST_EQUALITY(epetraComm->MyPID(), comm->getRank());
#endif
}

TEUCHOS_UNIT_TEST( MDComm, axisCommSizesPeriodicConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // If axisCommSizes is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = axisCommSizes.size(); axis < numDims; ++axis)
  {
    axisCommSizes.push_back(-1);
  }
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;

  // Construct an MDComm
  MDComm mdComm(comm, axisCommSizes, periodic);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==0));
  }
}
#endif

TEUCHOS_UNIT_TEST( MDComm, pListConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // If axisCommSizes is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = axisCommSizes.size(); axis < numDims; ++axis)
  {
    axisCommSizes.push_back(-1);
  }

  Teuchos::ParameterList plist;
  plist.set("axis comm sizes", axisCommSizes);

  MDComm mdComm(comm, plist);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Comm > epetraComm = mdComm.getEpetraComm();
  TEST_EQUALITY(epetraComm->NumProc(), comm->getSize());
  TEST_EQUALITY(epetraComm->MyPID(), comm->getRank());
#endif
}

TEUCHOS_UNIT_TEST( MDComm, pListConstructorPeriodic )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // If axisCommSizes is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = axisCommSizes.size(); axis < numDims; ++axis)
  {
    axisCommSizes.push_back(-1);
  }
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;

  Teuchos::ParameterList plist;
  plist.set("axis comm sizes", axisCommSizes);
  plist.set("periodic"       , periodic     );

  // Construct an MDComm
  MDComm mdComm(comm, plist);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==0));
  }
}

TEUCHOS_UNIT_TEST( MDComm, pListConstructorBad )
{
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Teuchos::ParameterList plist;

  // Test for illegal parameter name
  std::string illegalName("illegal parameter name");
  plist.set(illegalName, 0);
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterName);

  // Test for illegal parameter type
  plist.remove(illegalName);
  plist.set("axis comm sizes", "illegal parameter type");
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterType);

  // Test for illegal parameter value
  Teuchos::Array< int > axisCommSizes(3);
  axisCommSizes[0] = -2;
  axisCommSizes[1] = -1;
  axisCommSizes[2] = 0;
  plist.set("axis comm sizes", axisCommSizes);
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterValue);
}

TEUCHOS_UNIT_TEST( MDComm, numDimsConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  TEST_EQUALITY(mdComm.getAxisCommSize(0), comm->getSize());
  TEST_ASSERT(not mdComm.isPeriodic(0));
  for (int axis = 1; axis < numDims; ++axis)
  {
    TEST_EQUALITY_CONST(mdComm.getAxisCommSize(axis), 1);
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsAxisSizesConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisCommSizes);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims && axis < axisCommSizes.size(); ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsAxisSizesPeriodicConstructor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[numDims-1] = 1;
  // Construct the MDComm
  MDComm mdComm(comm, numDims, axisCommSizes, periodic);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims && axis < axisCommSizes.size(); ++axis)
  {
    if (axisCommSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisCommSize(axis), axisCommSizes[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==numDims-1));
  }
}

TEUCHOS_UNIT_TEST( MDComm, getTeuchosComm )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisCommSizes);
  TeuchosCommRCP newComm = mdComm.getTeuchosComm();

  TEST_EQUALITY_CONST(newComm.shares_resource(comm), true);

}

TEUCHOS_UNIT_TEST( MDComm, getAxisRank )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisCommSizes);
  // Get the final axisCommSizes and compute the strides
  axisCommSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisCommSizes[axis-1];
  }
  // Compute what the axis ranks should be for this processor
  Array< int > axisRanks(numDims);
  int rank = comm->getRank();
  for (int axis = numDims-1; axis > 0; --axis)
  {
    axisRanks[axis] = rank / strides[axis];
    rank            = rank % strides[axis];
  }
  axisRanks[0] = rank;

  // Test the getAxisRank() method
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getAxisRank(axis), axisRanks[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, getLowerNeighbor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisCommSizes);
  // Get the final axisCommSizes and compute the strides
  axisCommSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisCommSizes[axis-1];
  }
  // Get the axis ranks for this processor
  Array< int > axisRanks(numDims);
  for (int axis = numDims-1; axis > 0; --axis)
    axisRanks[axis] = mdComm.getAxisRank(axis);
  // Compute the ranks of the lower neighbors
  Array< int > lowerNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    lowerNeighborRanks[axis] = comm->getRank() - strides[axis];

  // Test the lower neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisRanks[axis] == 0)
    {
      TEST_EQUALITY(mdComm.getLowerNeighbor(axis), -1);
    }
    else
    {
      TEST_EQUALITY(mdComm.getLowerNeighbor(axis), lowerNeighborRanks[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, getLowerNeighborPeriodic )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  // Construct MDComm
  MDComm mdComm(comm, numDims, axisCommSizes, periodic);
  // Get the final axisCommSizes and compute the strides
  axisCommSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisCommSizes[axis-1];
  }
  // Get the axis ranks for this processor
  Array< int > axisRanks(numDims);
  for (int axis = numDims-1; axis > 0; --axis)
    axisRanks[axis] = mdComm.getAxisRank(axis);
  // Compute the ranks of the lower neighbors
  Array< int > lowerNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    if ((axis == 0) && (axisRanks[axis] == 0))
      lowerNeighborRanks[axis] = comm->getRank() +
        (axisCommSizes[axis]-1)*strides[axis];
    else
      lowerNeighborRanks[axis] = comm->getRank() - strides[axis];

  // Test the lower neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    if ((axis != 0) && (axisRanks[axis] == 0))
    {
      TEST_EQUALITY(mdComm.getLowerNeighbor(axis), -1);
    }
    else
    {
      TEST_EQUALITY(mdComm.getLowerNeighbor(axis), lowerNeighborRanks[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, getUpperNeighbor )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisCommSizes);
  // Get the final axisCommSizes and compute the strides
  axisCommSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisCommSizes[axis-1];
  }
  // Get the axis ranks for this processor
  Array< int > axisRanks(numDims);
  for (int axis = numDims-1; axis > 0; --axis)
    axisRanks[axis] = mdComm.getAxisRank(axis);
  // Compute the ranks of the upper neighbors
  Array< int > upperNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    upperNeighborRanks[axis] = comm->getRank() + strides[axis];

  // Test the upper neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisRanks[axis] == mdComm.getAxisCommSize(axis)-1)
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), -1);
    }
    else
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), upperNeighborRanks[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, getUpperNeighborPeriodic )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[numDims-1] = 1;
  // Construct MDComm
  MDComm mdComm(comm, numDims, axisCommSizes, periodic);
  // Get the final axisCommSizes and compute the strides
  axisCommSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisCommSizes[axis-1];
  }
  // Get the axis ranks for this processor
  Array< int > axisRanks(numDims);
  for (int axis = numDims-1; axis > 0; --axis)
    axisRanks[axis] = mdComm.getAxisRank(axis);
  // Compute the ranks of the upper neighbors
  Array< int > upperNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    if ((axis == numDims-1) && (axisRanks[axis] == axisCommSizes[axis]-1))
      upperNeighborRanks[axis] = comm->getRank() -
        (axisCommSizes[axis]-1)*strides[axis];
    else
      upperNeighborRanks[axis] = comm->getRank() + strides[axis];

  // Test the upper neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    if ((axis != numDims-1) && (axisRanks[axis] == axisCommSizes[axis]-1))
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), -1);
    }
    else
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), upperNeighborRanks[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, exceptions )
{
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // Construct MDComm
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Unit test methods that should throw exceptions
#if DOMI_ENABLE_ABC
  TEST_THROW(mdComm.getAxisCommSize( -1), Domi::RangeError);
  TEST_THROW(mdComm.isPeriodic(      -1), Domi::RangeError);
  TEST_THROW(mdComm.getAxisRank(     -1), Domi::RangeError);
  TEST_THROW(mdComm.getLowerNeighbor(-1), Domi::RangeError);
  TEST_THROW(mdComm.getUpperNeighbor(-1), Domi::RangeError);
  TEST_THROW(mdComm.getAxisCommSize( numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.isPeriodic(      numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getAxisRank(     numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getLowerNeighbor(numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getUpperNeighbor(numDims+1), Domi::RangeError);
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommLowerLeft )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
  
  // Figure out the lower left slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      int n = axisCommSizes[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = axisCommSizes[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm.getAxisRank(axis) >= newSizes[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY(subMDComm.numDims(), 0);
  }
#ifdef HAVE_EPETRA
  if (partOfSubComm)
  {
    Teuchos::RCP< const Epetra_Comm > epetraComm = subMDComm.getEpetraComm();
    TeuchosCommRCP teuchosComm = subMDComm.getTeuchosComm();
    TEST_EQUALITY(epetraComm->NumProc(), teuchosComm->getSize());
    TEST_EQUALITY(epetraComm->MyPID(), teuchosComm->getRank());
  }
  else
  {
    TEST_ASSERT(subMDComm.getEpetraComm().is_null());
  }
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommLowerRight )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
  
  // Figure out the lower right slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      int n = axisCommSizes[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = axisCommSizes[axis] - n;
    }
    else if (axis == 1)
    {
      int n = axisCommSizes[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = axisCommSizes[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getAxisRank(0) < axisCommSizes[0] - newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getAxisRank(1) >= newSizes[1])
      partOfSubComm = false;

#if 0
  if (partOfSubComm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.numDims(), 0);
  }
#ifdef HAVE_EPETRA
  if (partOfSubComm)
  {
    Teuchos::RCP< const Epetra_Comm > epetraComm = subMDComm.getEpetraComm();
    TeuchosCommRCP teuchosComm = subMDComm.getTeuchosComm();
    TEST_EQUALITY(epetraComm->NumProc(), teuchosComm->getSize());
    TEST_EQUALITY(epetraComm->MyPID(), teuchosComm->getRank());
  }
  else
  {
    TEST_ASSERT(subMDComm.getEpetraComm().is_null());
  }
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommUpperLeft )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
  
  // Figure out the upper left slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      int n = axisCommSizes[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else if (axis == 1)
    {
      int n = axisCommSizes[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = axisCommSizes[axis] - n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = axisCommSizes[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getAxisRank(0) >= newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getAxisRank(1) < axisCommSizes[1] - newSizes[1])
      partOfSubComm = false;

#if 0
  if (partOfSubComm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.numDims(), 0);
  }
#ifdef HAVE_EPETRA
  if (partOfSubComm)
  {
    Teuchos::RCP< const Epetra_Comm > epetraComm = subMDComm.getEpetraComm();
    TeuchosCommRCP teuchosComm = subMDComm.getTeuchosComm();
    TEST_EQUALITY(epetraComm->NumProc(), teuchosComm->getSize());
    TEST_EQUALITY(epetraComm->MyPID(), teuchosComm->getRank());
  }
  else
  {
    TEST_ASSERT(subMDComm.getEpetraComm().is_null());
  }
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommUpperRight )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
  
  // Figure out the upper right slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      int n = axisCommSizes[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = axisCommSizes[axis] - n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = axisCommSizes[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getAxisRank(0) < axisCommSizes[0] - newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getAxisRank(1) < axisCommSizes[1] - newSizes[1])
      partOfSubComm = false;

#if 0
  if (partOfSubComm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY_CONST(subMDComm.onSubcommunicator(), true);
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.onSubcommunicator()         , false);
    TEST_EQUALITY_CONST(subMDComm.getTeuchosComm().getRawPtr(), 0    );
    TEST_EQUALITY_CONST(subMDComm.numDims()                , 0    );

    TEST_THROW(subMDComm.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
#ifdef HAVE_EPETRA
  if (partOfSubComm)
  {
    Teuchos::RCP< const Epetra_Comm > epetraComm = subMDComm.getEpetraComm();
    TeuchosCommRCP teuchosComm = subMDComm.getTeuchosComm();
    TEST_EQUALITY(epetraComm->NumProc(), teuchosComm->getSize());
    TEST_EQUALITY(epetraComm->MyPID(), teuchosComm->getRank());
  }
  else
  {
    TEST_ASSERT(subMDComm.getEpetraComm().is_null());
  }
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommReduce )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  MDComm mdComm(comm, numDims, axisCommSizes);

  // Get the actual axisCommSizes, in case the command-line
  // specification was incomplete
  Array< int > finalCommSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    finalCommSizes[axis] = mdComm.getAxisCommSize(axis);

  // Compute the axis strides
  Array< int > axisStrides(numDims);
  Array< int > newAxisStrides;
  axisStrides[0] = 1;
  for (int axis = 1; axis < numDims; ++axis)
    axisStrides[axis] = axisStrides[axis-1] * finalCommSizes[axis-1];

  // We will reduce this MDComm several times by using the axis-rank
  // constructor along each dimension
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Compute the new axis strides
    newAxisStrides.clear();
    for (int newAxis = 0; newAxis < numDims; ++newAxis)
      if (newAxis != axis)
        newAxisStrides.push_back(axisStrides[newAxis]);

    int redAxisRank = finalCommSizes[axis] / 2;
    int newDims = (numDims > 1) ? numDims - 1 : numDims;
    bool partOfSubComm = true;
    if (mdComm.getAxisRank(axis) != redAxisRank)
      partOfSubComm = false;
    MDComm reducedMdComm(mdComm, axis, redAxisRank);
    if (partOfSubComm)
    {
      TEST_ASSERT(reducedMdComm.onSubcommunicator());
      TEST_EQUALITY(reducedMdComm.numDims(), newDims);
      for (int newAxis = 0; newAxis < newDims; ++newAxis)
      {
        int axisRank = reducedMdComm.getAxisRank(newAxis);
        if (numDims == 1)
        {
          TEST_EQUALITY_CONST(reducedMdComm.getAxisCommSize(newAxis), 1);
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY_CONST(reducedMdComm.getAxisRank(newAxis),0);
          TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
        }
        else if (newAxis < axis)
        {
          TEST_EQUALITY(reducedMdComm.getAxisCommSize(newAxis),
                        finalCommSizes[newAxis]);
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY(reducedMdComm.getAxisRank(newAxis),
                        mdComm.getAxisRank(newAxis));
          if (axisRank == 0)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          }
          else
          {
            int lowerNeighbor = 0;
            for (int i = 0; i < newDims; ++i)
              lowerNeighbor += reducedMdComm.getAxisRank(i) * newAxisStrides[i];
            lowerNeighbor -= newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getLowerNeighbor(newAxis),
                          lowerNeighbor);
          }
          if (axisRank == reducedMdComm.getAxisCommSize(newAxis)-1)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
          }
          else
          {
            int upperNeighbor = 0;
            for (int i = 0; i < newDims; ++i)
              upperNeighbor += reducedMdComm.getAxisRank(i) * newAxisStrides[i];
            upperNeighbor += newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getUpperNeighbor(newAxis),
                          upperNeighbor);
          }
        }
        else
        {
          TEST_EQUALITY(reducedMdComm.getAxisCommSize(newAxis),
                        finalCommSizes[newAxis+1])
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY(reducedMdComm.getAxisRank(newAxis),
                        mdComm.getAxisRank(newAxis+1));
          if (axisRank == 0)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          }
          else
          {
            int lowerNeighbor = 0;
            for (int i = 0; i < newDims; ++i)
              lowerNeighbor += reducedMdComm.getAxisRank(i) * newAxisStrides[i];
            lowerNeighbor -= newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getLowerNeighbor(newAxis),
                          lowerNeighbor);
          }
          if (axisRank == reducedMdComm.getAxisCommSize(newAxis)-1)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
          }
          else
          {
            int upperNeighbor = 0;
            for (int i = 0; i < newDims; ++i)
              upperNeighbor += reducedMdComm.getAxisRank(i) * newAxisStrides[i];
            upperNeighbor += newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getUpperNeighbor(newAxis),
                          upperNeighbor);
          }
        }
      }
    }
    else
    {
      TEST_ASSERT(not reducedMdComm.onSubcommunicator());
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, subCommPeriodic )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  axisCommSizes = Domi::splitStringOfIntsWithCommas(axisCommSizesStr);
  // Construct the periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  MDComm mdComm(comm, numDims, axisCommSizes, periodic);

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm.getAxisCommSize(axis);
  
  // Figure out the lower slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 1)
    {
      int n = axisCommSizes[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = axisCommSizes[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm.getAxisRank(axis) >= newSizes[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_EQUALITY(subMDComm.isPeriodic(axis), (axis == 0));
    }
  }
  else
  {
    TEST_EQUALITY(subMDComm.numDims(), 0);
  }
}

}  // namespace
