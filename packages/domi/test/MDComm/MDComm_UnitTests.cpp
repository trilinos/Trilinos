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

// Domi includes
#include "Domi_Utils.hpp"
#include "Domi_MDComm.hpp"

namespace
{

using std::string;
using Teuchos::Array;
using Domi::TeuchosCommRCP;
using Domi::MDComm;

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

TEUCHOS_UNIT_TEST( MDComm, regularizeAxisSizes )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Array< int > axisSizesActual =
    Domi::regularizeAxisSizes(comm->getSize(), numDims, axisSizes);
  
  TEST_EQUALITY(axisSizesActual.size(), numDims);

  int numProcs = 1;
  for (int axis = 0; axis < numDims; ++axis)
    numProcs *= axisSizesActual[axis];

  TEST_EQUALITY(numProcs, comm->getSize());
}

TEUCHOS_UNIT_TEST( MDComm, axisSizesConstructor )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  // If axisSizes is shorter than numDims, pad it with -1 values at
  // the end
  for (int axis = axisSizes.size(); axis < numDims; ++axis)
  {
    axisSizes.push_back(-1);
  }
  MDComm mdComm(comm, axisSizes);

  TEST_EQUALITY(mdComm.getNumDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axisSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisSize(axis), axisSizes[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsConstructor )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims);

  TEST_EQUALITY(mdComm.getNumDims(), numDims);

  TEST_EQUALITY(mdComm.getAxisSize(0), comm->getSize());
  for (int axis = 1; axis < numDims; ++axis)
  {
    TEST_EQUALITY_CONST(mdComm.getAxisSize(axis), 1);
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsAxisSizesConstructor )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisSizes);

  TEST_EQUALITY(mdComm.getNumDims(), numDims);

  for (int axis = 0; axis < numDims && axis < axisSizes.size(); ++axis)
  {
    if (axisSizes[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getAxisSize(axis), axisSizes[axis]);
    }
  }
}

TEUCHOS_UNIT_TEST( MDComm, getTeuchosComm )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisSizes);
  TeuchosCommRCP newComm = mdComm.getTeuchosComm();

  TEST_EQUALITY_CONST(newComm.shares_resource(comm), true);

}

TEUCHOS_UNIT_TEST( MDComm, getAxisRank )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisSizes);
  // Get the final axisSizes and compute the strides
  axisSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisSizes[axis] = mdComm.getAxisSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisSizes[axis-1];
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
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisSizes);
  // Get the final axisSizes and compute the strides
  axisSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisSizes[axis] = mdComm.getAxisSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisSizes[axis-1];
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

TEUCHOS_UNIT_TEST( MDComm, getUpperNeighbor )
{
  Domi::splitStringOfIntsWithCommas(axisSizesStr, axisSizes);
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, axisSizes);
  // Get the final axisSizes and compute the strides
  axisSizes.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    axisSizes[axis] = mdComm.getAxisSize(axis);
    if (axis == 0)
      strides[axis] = 1;
    else
      strides[axis] = strides[axis-1] * axisSizes[axis-1];
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
    if (axisRanks[axis] == mdComm.getAxisSize(axis)-1)
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), -1);
    }
    else
    {
      TEST_EQUALITY(mdComm.getUpperNeighbor(axis), upperNeighborRanks[axis]);
    }
  }
}

}  // namespace
