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
#include "Domi_Exceptions.hpp"

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

TEUCHOS_UNIT_TEST( MDComm_Subcomm, lowerLeft )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm->getAxisCommSize(axis);
  
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
  bool partOfSubcomm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm->getAxisRank(axis) >= newSizes[axis])
      partOfSubcomm = false;

  // Do some unit tests
  if (partOfSubcomm)
  {
    TEST_EQUALITY(subMDComm.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY(subMDComm.getNumDims(), 0);
  }
}

TEUCHOS_UNIT_TEST( MDComm_Subcomm, lowerRight )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm->getAxisCommSize(axis);
  
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
  bool partOfSubcomm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newSizes[0])
    partOfSubcomm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) >= newSizes[1])
      partOfSubcomm = false;

#if 0
  if (partOfSubcomm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubcomm)
  {
    TEST_EQUALITY(subMDComm.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.getNumDims(), 0);
  }
}

TEUCHOS_UNIT_TEST( MDComm_Subcomm, upperLeft )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm->getAxisCommSize(axis);
  
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
  bool partOfSubcomm = true;
  if (mdComm->getAxisRank(0) >= newSizes[0])
    partOfSubcomm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newSizes[1])
      partOfSubcomm = false;

#if 0
  if (partOfSubcomm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubcomm)
  {
    TEST_EQUALITY(subMDComm.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_ASSERT(not subMDComm.isPeriodic(axis));
    }
  }
  else
  {
    TEST_EQUALITY_CONST(subMDComm.getNumDims(), 0);
  }
}

TEUCHOS_UNIT_TEST( MDComm_Subcomm, upperRight )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  MDCommRCP mdComm = Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes));

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm->getAxisCommSize(axis);
  
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
  bool partOfSubcomm = true;
  if (mdComm->getAxisRank(0) < axisCommSizes[0] - newSizes[0])
    partOfSubcomm = false;
  if (numDims > 1)
    if (mdComm->getAxisRank(1) < axisCommSizes[1] - newSizes[1])
      partOfSubcomm = false;

#if 0
  if (partOfSubcomm)
    std::cout << "P" << comm->getRank() << ": IS part of sub-comm "
              << newSizes << std::endl;
  else
    std::cout << "P" << comm->getRank() << ": is NOT part of sub-comm"
              << std::endl;
#endif

  // Do some unit tests
  if (partOfSubcomm)
  {
    TEST_EQUALITY_CONST(subMDComm.onSubcommunicator(), true);
    TEST_EQUALITY(subMDComm.getNumDims(), numDims);
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
    TEST_EQUALITY_CONST(subMDComm.getNumDims()                , 0    );

    TEST_THROW(subMDComm.getAxisCommSize(0) , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getAxisRank(0)     , Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getLowerNeighbor(0), Domi::SubcommunicatorError);
    TEST_THROW(subMDComm.getUpperNeighbor(0), Domi::SubcommunicatorError);
  }
}

TEUCHOS_UNIT_TEST( MDComm_Subcomm, periodic )
{
  // Construct the MDComm from command-line arguments
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  Domi::splitStringOfIntsWithCommas(axisCommSizesStr, axisCommSizes);
  // Construct the periodic flags
  Array< int > periodic(numDims, 0);
  periodic[0] = 1;
  MDCommRCP mdComm =
    Teuchos::rcp(new MDComm(comm, numDims, axisCommSizes, periodic));

  // Get the final axisCommSizes
  axisCommSizes.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    axisCommSizes[axis] = mdComm->getAxisCommSize(axis);
  
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
  bool partOfSubcomm = true;
  for (int axis = 0; axis < numDims; ++axis)
    if (mdComm->getAxisRank(axis) >= newSizes[axis])
      partOfSubcomm = false;

  // Do some unit tests
  if (partOfSubcomm)
  {
    TEST_EQUALITY(subMDComm.getNumDims(), numDims);
    for (int axis = 0; axis < numDims; ++axis)
    {
      TEST_EQUALITY(subMDComm.getAxisCommSize(axis), newSizes[axis]);
      TEST_EQUALITY(subMDComm.isPeriodic(axis), (axis == 0));
    }
  }
  else
  {
    TEST_EQUALITY(subMDComm.getNumDims(), 0);
  }
}

}
