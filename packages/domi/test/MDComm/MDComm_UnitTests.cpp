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
using Domi::MDComm;
using Domi::Slice;
typedef Domi::Ordinal Ordinal;
const Ordinal & Default = Domi::Slice::Default;

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

#if 1
TEUCHOS_UNIT_TEST( MDComm, regularizeCommDims )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  Array< int > commDimsActual =
    Domi::regularizeCommDims(comm->getSize(), numDims, commDims);
  
  TEST_EQUALITY(commDimsActual.size(), numDims);

  int numProcs = 1;
  for (int axis = 0; axis < numDims; ++axis)
    numProcs *= commDimsActual[axis];

  TEST_EQUALITY(numProcs, comm->getSize());
}

TEUCHOS_UNIT_TEST( MDComm, commDimsConstructor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // If commDims is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = commDims.size(); axis < numDims; ++axis)
  {
    commDims.push_back(-1);
  }
  MDComm mdComm(comm, commDims);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
    }
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
#ifdef HAVE_EPETRA
  Teuchos::RCP< const Epetra_Comm > epetraComm = mdComm.getEpetraComm();
  TEST_EQUALITY(epetraComm->NumProc(), comm->getSize());
  TEST_EQUALITY(epetraComm->MyPID(), comm->getRank());
#endif
}

TEUCHOS_UNIT_TEST( MDComm, commDimsPeriodicConstructor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // If commDims is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = commDims.size(); axis < numDims; ++axis)
  {
    commDims.push_back(-1);
  }
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;

  // Construct an MDComm
  MDComm mdComm(comm, commDims, periodic);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==0));
  }
}
#endif

TEUCHOS_UNIT_TEST( MDComm, pListConstructor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // If commDims is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = commDims.size(); axis < numDims; ++axis)
  {
    commDims.push_back(-1);
  }

  Teuchos::ParameterList plist;
  plist.set("comm dimensions", commDims);

  MDComm mdComm(comm, plist);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
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
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // If commDims is shorter than numDims, pad it with -1 values
  // at the end
  for (int axis = commDims.size(); axis < numDims; ++axis)
  {
    commDims.push_back(-1);
  }
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;

  Teuchos::ParameterList plist;
  plist.set("comm dimensions", commDims);
  plist.set("periodic"       , periodic     );

  // Construct an MDComm
  MDComm mdComm(comm, plist);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims; ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==0));
  }
}

TEUCHOS_UNIT_TEST( MDComm, pListConstructorBad )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  Teuchos::ParameterList plist;

  // Test for illegal parameter name
  std::string illegalName("illegal parameter name");
  plist.set(illegalName, 0);
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterName);

  // Test for illegal parameter type
  plist.remove(illegalName);
  plist.set("comm dimensions", "illegal parameter type");
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterType);

  // Test for illegal parameter value
  Teuchos::Array< int > commDims(3);
  commDims[0] = -2;
  commDims[1] = -1;
  commDims[2] = 0;
  plist.set("comm dimensions", commDims);
  TEST_THROW(MDComm(comm, plist), Teuchos::Exceptions::InvalidParameterValue);
}

TEUCHOS_UNIT_TEST( MDComm, numDimsConstructor )
{
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  TEST_EQUALITY(mdComm.getCommDim(0), comm->getSize());
  TEST_ASSERT(not mdComm.isPeriodic(0));
  for (int axis = 1; axis < numDims; ++axis)
  {
    TEST_EQUALITY_CONST(mdComm.getCommDim(axis), 1);
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsCommDimsConstructor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, commDims);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims && axis < commDims.size(); ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
    }
    TEST_ASSERT(not mdComm.isPeriodic(axis));
  }
}

TEUCHOS_UNIT_TEST( MDComm, numDimsCommDimsPeriodicConstructor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[numDims-1] = 1;
  // Construct the MDComm
  MDComm mdComm(comm, numDims, commDims, periodic);

  TEST_EQUALITY(mdComm.numDims(), numDims);

  for (int axis = 0; axis < numDims && axis < commDims.size(); ++axis)
  {
    if (commDims[axis] > 0)
    {
      TEST_EQUALITY(mdComm.getCommDim(axis), commDims[axis]);
    }
    TEST_EQUALITY(mdComm.isPeriodic(axis), (axis==numDims-1));
  }
}

TEUCHOS_UNIT_TEST( MDComm, getTeuchosComm )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, commDims);
  Teuchos::RCP< const Teuchos::Comm< int > > newComm = mdComm.getTeuchosComm();

  TEST_EQUALITY_CONST(newComm.shares_resource(comm), true);

}

TEUCHOS_UNIT_TEST( MDComm, getCommIndex )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, commDims);
  // Get the final commDims and compute the strides
  commDims.resize(numDims);
  Array< int > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int index1 = numDims-axis-1;
    int index2 = numDims-axis;
    if (axis == 0)
      strides[index1] = 1;
    else
      strides[index1] = strides[index2] * commDims[index2];
  }
  // Compute what the axis ranks should be for this processor
  Array< int > commIndex(numDims);
  int rank = comm->getRank();
  for (int axis = 0; axis < numDims; ++axis)
  {
    commIndex[axis] = rank / strides[axis];
    rank            = rank % strides[axis];
  }

  // Test the getCommIndex() method
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getCommIndex(axis), commIndex[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, getLowerNeighbor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims, compute the commStrides, commIndex and
  // lower neighbor ranks
  commDims.resize(numDims);
  Array< int > commStrides(numDims);
  Array< int > commIndex(numDims);
  Array< int > lowerNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int index1 = numDims-axis-1;
    int index2 = numDims-axis;
    if (axis == 0)
      commStrides[index1] = 1;
    else
      commStrides[index1] = commStrides[index2] * commDims[index2];
  }
  for (int axis = 0; axis < numDims; ++axis)
  {
    commIndex[axis] = mdComm.getCommIndex(axis);
    if (commIndex[axis] == 0)
      lowerNeighborRanks[axis] = -1;
    else
      lowerNeighborRanks[axis] = comm->getRank() - commStrides[axis];
  }

  // Test the lower neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getLowerNeighbor(axis), lowerNeighborRanks[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, getLowerNeighborPeriodic )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();

  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;

  // Construct MDComm
  MDComm mdComm(comm, numDims, commDims, periodic);

  // Get the final commDims, compute the commStrides, commIndex, and
  // lower neighbor ranks
  commDims.resize(numDims);
  Array< int > commStrides(numDims);
  Array< int > commIndex(numDims);
  Array< int > lowerNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int index1 = numDims-axis-1;
    int index2 = numDims-axis;
    if (axis == 0)
      commStrides[index1] = 1;
    else
      commStrides[index1] = commStrides[index2] * commDims[index2];
  }
  for (int axis = 0; axis < numDims; ++axis)
  {
    commIndex[axis] = mdComm.getCommIndex(axis);
    if (periodic[axis] && commIndex[axis] == 0)
      lowerNeighborRanks[axis] = comm->getRank() + (commDims[axis]-1) *
        commStrides[axis];
    else if (commIndex[axis] == 0)
      lowerNeighborRanks[axis] = -1;
    else
      lowerNeighborRanks[axis] = comm->getRank() - commStrides[axis];
  }

  // Test the lower neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getLowerNeighbor(axis), lowerNeighborRanks[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, getUpperNeighbor )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims, compute the commStrides, commIndex and
  // upper neighbor ranks
  commDims.resize(numDims);
  Array< int > commStrides(numDims);
  Array< int > commIndex(numDims);
  Array< int > upperNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int index1 = numDims-axis-1;
    int index2 = numDims-axis;
    if (axis == 0)
      commStrides[index1] = 1;
    else
      commStrides[index1] = commStrides[index2] * commDims[index2];
  }
  for (int axis = 0; axis < numDims; ++axis)
  {
    commIndex[axis] = mdComm.getCommIndex(axis);
    if (commIndex[axis] == commDims[axis]-1)
      upperNeighborRanks[axis] = -1;
    else
      upperNeighborRanks[axis] = comm->getRank() + commStrides[axis];
  }

  // Test the upper neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getUpperNeighbor(axis), upperNeighborRanks[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, getUpperNeighborPeriodic )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();

  // Construct periodic flags
  Array< int > periodic(numDims, 0);
  periodic[numDims-1] = 1;

  // Construct MDComm
  MDComm mdComm(comm, numDims, commDims, periodic);

  // Get the final commDims, and compute the commStrides, commIndex,
  // and upper neighbor ranks
  commDims.resize(numDims);
  Array< int > commStrides(numDims);
  Array< int > commIndex(numDims);
  Array< int > upperNeighborRanks(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  for (int axis = 0; axis < numDims; ++axis)
  {
    int index1 = numDims-axis-1;
    int index2 = numDims-axis;
    if (axis == 0)
      commStrides[index1] = 1;
    else
      commStrides[index1] = commStrides[index2] * commDims[index2];
  }
  for (int axis = 0; axis < numDims; ++axis)
  {
    commIndex[axis] = mdComm.getCommIndex(axis);
    if (periodic[axis] && commIndex[axis] == commDims[axis]-1)
      upperNeighborRanks[axis] = comm->getRank() - (commDims[axis]-1) *
        commStrides[axis];
    else if (commIndex[axis] == commDims[axis]-1)
      upperNeighborRanks[axis] = -1;
    else
      upperNeighborRanks[axis] = comm->getRank() + commStrides[axis];
  }

  // Test the upper neighbor
  for (int axis = 0; axis < numDims; ++axis)
  {
    TEST_EQUALITY(mdComm.getUpperNeighbor(axis), upperNeighborRanks[axis]);
  }
}

TEUCHOS_UNIT_TEST( MDComm, exceptions )
{
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  // Construct MDComm
  MDComm mdComm(comm, numDims, commDims);

  // Unit test methods that should throw exceptions
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(mdComm.getCommDim(      -1), Domi::RangeError);
  TEST_THROW(mdComm.isPeriodic(      -1), Domi::RangeError);
  TEST_THROW(mdComm.getCommIndex(    -1), Domi::RangeError);
  TEST_THROW(mdComm.getLowerNeighbor(-1), Domi::RangeError);
  TEST_THROW(mdComm.getUpperNeighbor(-1), Domi::RangeError);
  TEST_THROW(mdComm.getCommDim(      numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.isPeriodic(      numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getCommIndex(    numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getLowerNeighbor(numDims+1), Domi::RangeError);
  TEST_THROW(mdComm.getUpperNeighbor(numDims+1), Domi::RangeError);
#endif
}

TEUCHOS_UNIT_TEST( MDComm, subCommLowerLeft )
{
  // Construct the MDComm from command-line arguments
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  
  // Figure out the lower left slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      int n = commDims[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = commDims[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm.getCommIndex(axis) >= newSizes[axis])
      partOfSubComm = false;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMDComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getCommDim(axis), newSizes[axis]);
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
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      subMDComm.getTeuchosComm();
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
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  
  // Figure out the lower right slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      int n = commDims[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = commDims[axis] - n;
    }
    else if (axis == 1)
    {
      int n = commDims[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = commDims[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getCommIndex(0) < commDims[0] - newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getCommIndex(1) >= newSizes[1])
      partOfSubComm = false;

#if 0
  std::cout << "P" << comm->getRank() << ": commIndex = ("
            << mdComm.getCommIndex(0) << "," << mdComm.getCommIndex(1) << ")"
            << std::endl;
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
      TEST_EQUALITY(subMDComm.getCommDim(axis), newSizes[axis]);
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
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      subMDComm.getTeuchosComm();
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
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  
  // Figure out the upper left slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 0)
    {
      int n = commDims[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else if (axis == 1)
    {
      int n = commDims[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = commDims[axis] - n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = commDims[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getCommIndex(0) >= newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getCommIndex(1) < commDims[1] - newSizes[1])
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
      TEST_EQUALITY(subMDComm.getCommDim(axis), newSizes[axis]);
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
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      subMDComm.getTeuchosComm();
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
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  MDComm mdComm(comm, numDims, commDims);

  // Get the final commDims
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  
  // Figure out the upper right slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis < 2)
    {
      int n = commDims[axis] / 2;
      slices[axis] = Slice(n,Default);
      newSizes[axis] = commDims[axis] - n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = commDims[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMDComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  if (mdComm.getCommIndex(0) < commDims[0] - newSizes[0])
    partOfSubComm = false;
  if (numDims > 1)
    if (mdComm.getCommIndex(1) < commDims[1] - newSizes[1])
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
      TEST_EQUALITY(subMDComm.getCommDim(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.onSubcommunicator()         , false);
    TEST_EQUALITY_CONST(subMDComm.getTeuchosComm().getRawPtr(), 0    );
    TEST_EQUALITY_CONST(subMDComm.numDims()                   , 0    );

    TEST_THROW(subMDComm.getCommDim(0)      , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getCommIndex(0)    , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
#ifdef HAVE_EPETRA
  if (partOfSubComm)
  {
    Teuchos::RCP< const Epetra_Comm > epetraComm = subMDComm.getEpetraComm();
    Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm =
      subMDComm.getTeuchosComm();
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
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);
  MDComm mdComm(comm, numDims, commDims);

  // Get the actual commDims, in case the command-line specification
  // was incomplete
  Array< int > finalCommDims(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    finalCommDims[axis] = mdComm.getCommDim(axis);

  // Compute the axis commStrides
  Array< int > commStrides =
    Domi::computeStrides<int,int>(finalCommDims,Domi::MDComm::commLayout);

  // We will reduce this MDComm several times by using the comm index
  // constructor along each dimension
  Array< int > newCommDims;
  Array< int > newAxisStrides;
  for (int axis = 0; axis < numDims; ++axis)
  {
    // Compute the new comm dimensions and comm strides
    int newNumDims = (numDims > 1) ? numDims - 1 : numDims;
    newCommDims.clear();
    newAxisStrides.clear();
    for (int newAxis = 0; newAxis < numDims; ++newAxis)
    {
      if (newAxis != axis)
        newCommDims.push_back(finalCommDims[newAxis]);
    }
    if (numDims == 1) newCommDims.push_back(finalCommDims[0]);
    newAxisStrides =
      Domi::computeStrides<int,int>(newCommDims, Domi::MDComm::commLayout);

    int redAxisRank = finalCommDims[axis] / 2;
    bool partOfSubComm = true;
    if (mdComm.getCommIndex(axis) != redAxisRank)
      partOfSubComm = false;

    // Construct the reduced MDComm
    MDComm reducedMdComm(mdComm, axis, redAxisRank);

    // Test the reduced MDComm
    if (partOfSubComm)
    {
      TEST_ASSERT(reducedMdComm.onSubcommunicator());
      TEST_EQUALITY(reducedMdComm.numDims(), newNumDims);
      for (int newAxis = 0; newAxis < newNumDims; ++newAxis)
      {
        int axisRank = reducedMdComm.getCommIndex(newAxis);
        if (numDims == 1)
        {
          TEST_EQUALITY_CONST(reducedMdComm.getCommDim(newAxis), 1);
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY_CONST(reducedMdComm.getCommIndex(newAxis),0);
          TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
        }
        else if (newAxis < axis)
        {
          TEST_EQUALITY(reducedMdComm.getCommDim(newAxis),
                        finalCommDims[newAxis]);
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY(reducedMdComm.getCommIndex(newAxis),
                        mdComm.getCommIndex(newAxis));
          if (axisRank == 0)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          }
          else
          {
            int lowerNeighbor = 0;
            for (int i = 0; i < newNumDims; ++i)
              lowerNeighbor += reducedMdComm.getCommIndex(i) * newAxisStrides[i];
            lowerNeighbor -= newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getLowerNeighbor(newAxis),
                          lowerNeighbor);
          }
          if (axisRank == reducedMdComm.getCommDim(newAxis)-1)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
          }
          else
          {
            int upperNeighbor = 0;
            for (int i = 0; i < newNumDims; ++i)
              upperNeighbor += reducedMdComm.getCommIndex(i) * newAxisStrides[i];
            upperNeighbor += newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getUpperNeighbor(newAxis),
                          upperNeighbor);
          }
        }
        else
        {
          TEST_EQUALITY(reducedMdComm.getCommDim(newAxis),
                        finalCommDims[newAxis+1])
          TEST_ASSERT(not reducedMdComm.isPeriodic(newAxis));
          TEST_EQUALITY(reducedMdComm.getCommIndex(newAxis),
                        mdComm.getCommIndex(newAxis+1));
          if (axisRank == 0)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getLowerNeighbor(newAxis),-1);
          }
          else
          {
            int lowerNeighbor = 0;
            for (int i = 0; i < newNumDims; ++i)
              lowerNeighbor += reducedMdComm.getCommIndex(i) * newAxisStrides[i];
            lowerNeighbor -= newAxisStrides[newAxis];
            TEST_EQUALITY(reducedMdComm.getLowerNeighbor(newAxis),
                          lowerNeighbor);
          }
          if (axisRank == reducedMdComm.getCommDim(newAxis)-1)
          {
            TEST_EQUALITY_CONST(reducedMdComm.getUpperNeighbor(newAxis),-1);
          }
          else
          {
            int upperNeighbor = 0;
            for (int i = 0; i < newNumDims; ++i)
              upperNeighbor += reducedMdComm.getCommIndex(i) * newAxisStrides[i];
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
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  commDims = Domi::splitStringOfIntsWithCommas(commDimsStr);

  // Construct the periodic flags and the MDComm
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  MDComm mdComm(comm, numDims, commDims, periodic);

  // Get the final commDims
  commDims.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    commDims[axis] = mdComm.getCommDim(axis);
  
  // Figure out the lower slice
  Array< Slice > slices(numDims);
  Array< int >   newSizes(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (axis == 1)
    {
      int n = commDims[axis] / 2;
      if (n == 0) n = 1;
      slices[axis] = Slice(n);
      newSizes[axis] = n;
    }
    else
    {
      slices[axis] = Slice();
      newSizes[axis] = commDims[axis];
    }
  }

  // Construct the sub-MDComm
  MDComm subMdComm(mdComm, slices);

  // Should this processor be a part of the sub-MDComm?
  bool partOfSubComm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm.getCommIndex(axis) >= newSizes[axis])
      partOfSubComm = false;

  // std::cout << comm->getRank() << ": slices = " << slices << ", newSizes = "
  //           << newSizes << ", part of sub-communicator = " << partOfSubComm
  //           << std::endl;

  // Do some unit tests
  if (partOfSubComm)
  {
    TEST_EQUALITY(subMdComm.numDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMdComm.getCommDim(axis), newSizes[axis]);
      TEST_EQUALITY(subMdComm.isPeriodic(axis), (axis == 0));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMdComm.numDims(), 0);
  }
}

}  // namespace
