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

// System includes
#include <string>
#include <stdlib.h>

// Domi includes
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"
#include "Domi_getValidParameters.hpp"

namespace Domi
{

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
factor(int n)
{
  Teuchos::Array< int > factors;

  if (n < 2)
  {
    factors.push_back(n);
    n = 0;
  }

  int z = 2;
  while (z*z <= n)
  {
    if (n % z == 0)
    {
      factors.push_back(z);
      n /= z;
    }
    else
    {
      z += 1;
    }
  }

  if (n > 1) factors.push_back(n);

  return factors;
}

////////////////////////////////////////////////////////////////////////

int indexOfMax(const Teuchos::ArrayView< const float > & seq)
{
  int result = 0;
  for (int i = 1; i < seq.size(); ++i)
    if (seq[i] > seq[result])
      result = i;
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
decomposeProcs(int nprocs,
               const Teuchos::ArrayView< dim_type > & dimensions)
{
  // Compute a float version of dimensions
  Teuchos::Array< float > dims;
  for (int i = 0; i < dimensions.size(); ++i)
    dims.push_back(float(dimensions[i]));
  Teuchos::Array< int > ifact = factor(nprocs);
  // Compute the result
  Teuchos::Array< int > result(dimensions.size(), 1);
  for (int i = ifact.size()-1; i >= 0; --i)
  {
    int imax      = indexOfMax(dims);
    dims[imax]   /= ifact[i];
    result[imax] *= ifact[i];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
regularizeCommDims(int numProcs,
                   int numDims,
                   const Teuchos::ArrayView< const int > & commDims)
{
  // Allocate the return array, filled with the value -1
  Teuchos::Array< int > result(numDims, -1);
  // Copy the candidate array into the return array
  for (int axis = 0; axis < numDims && axis < commDims.size(); ++axis)
    result[axis] = commDims[axis];
  // Compute the block of processors accounted for, and the number of
  // unspecified axis sizes
  int block       = 1;
  int unspecified = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (result[axis] <= 0)
      ++unspecified;
    else
      block *= result[axis];
  }
  // If all processor counts are specified, check the processor block
  // against the total number of processors and return
  if (unspecified == 0)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (block != numProcs),
      InvalidArgument,
      "Product of axis processor sizes (" << block << ") does not "
      "equal total number of processors (" << numProcs << ")");
    return result;
  }
  // For underspecified processor partitions, give the remainder to
  // the first unspecified axis and set all the rest to 1
  TEUCHOS_TEST_FOR_EXCEPTION(
    (numProcs % block),
    InvalidArgument,
    "Number of processors (" << numProcs << ") do not divide evenly by "
    << block);
  int quotient = numProcs / block;
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (result[axis] <= 0)
    {
      result[axis] = quotient;
      quotient = 1;
    }
  }
  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
regularizeCommDims(int numProcs,
                   Teuchos::ParameterList & plist)
{
  // Obtain the commDims and dimensions from the ParameterList, and
  // determine the number of dimensions
  Teuchos::Array< int > commDims = plist.get("comm dimensions",
                                             Teuchos::Array< int >());
  Teuchos::Array< dim_type > dims = plist.get("dimensions",
                                              Teuchos::Array< dim_type >());
  int numDims = dims.size();
  if (numDims == 0) numDims = commDims.size();
  if (numDims == 0)
  {
    numDims = 1;
    commDims.push_back(numProcs);
  }
  // Allocate the return array, filled with the value -1
  Teuchos::Array< int > result(numDims, -1);
  // Copy the candidate array into the return array
  for (int axis = 0; axis < numDims && axis < commDims.size(); ++axis)
    result[axis] = commDims[axis];
  // Compute the block of processors accounted for, and the number of
  // unspecified axis sizes
  int block       = 1;
  int unspecified = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (result[axis] <= 0)
      ++unspecified;
    else
      block *= result[axis];
  }
  // If all processor counts are specified, check the processor block
  // against the total number of processors and return
  if (unspecified == 0)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (block != numProcs),
      InvalidArgument,
      "Product of axis processor sizes (" << block << ") does not "
      "equal total number of processors (" << numProcs << ")");
    return result;
  }
  // Make sure that the number of processors divides equally into the
  // product of all of the specified axis processors
  TEUCHOS_TEST_FOR_EXCEPTION(
    (numProcs % block),
    InvalidArgument,
    "Number of processors (" << numProcs << ") do not divide evenly by "
    << block);
  // Create an array of dimensions with entries for every commDim that
  // is not specified
  Teuchos::Array< dim_type > myDims;
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (result[axis] <= 0)
    {
      if (axis < dims.size())
        myDims.push_back(dims[axis]);
      else
        // If "dimensions" was not specified, then myDims will be an
        // array of equal values, namely the value of numProcs
        myDims.push_back(numProcs);
    }
  }
  // For underspecified processor partitions, we'll use the
  // decomposeProcs() function
  Teuchos::Array< int > myCommDims = decomposeProcs(numProcs/block, myDims);
  // Copy the results into the result array and return
  int ii = 0;
  for (int axis = 0; axis < numDims; ++axis)
  {
    if (result[axis] <= 0)
      result[axis] = myCommDims[ii++];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
computeCommIndexes(int rank,
                   const Teuchos::ArrayView< int > & commStrides)
{
  Teuchos::Array< int > result(commStrides.size());
  for (int axis = 0; axis < commStrides.size(); ++axis)
  {
    result[axis] = rank / commStrides[axis];
    rank         = rank % commStrides[axis];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
createArrayOfInts(int numDims,
                  const Teuchos::ArrayView< const int > & source)
{
  Teuchos::Array< int > result(numDims, 0);
  for (int axis = 0; axis < numDims && axis < source.size(); ++axis)
    result[axis] = source[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int >
splitStringOfIntsWithCommas(std::string data)
{
  Teuchos::Array< int > result;
  size_t current = 0;
  while (current < data.size())
  {
    size_t next = data.find(",", current);
    if (next == std::string::npos) next = data.size();
    result.push_back(atoi(data.substr(current, next-current).c_str()));
    current = next + 1;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

// Specializations of mpiType<T> for all basic data types supported by
// MPI

template<>
MPI_Datatype mpiType<char>()
{
  return MPI_CHAR;
}

template<>
MPI_Datatype mpiType<signed char>()
{
  return MPI_SIGNED_CHAR;
}

template<>
MPI_Datatype mpiType<unsigned char>()
{
  return MPI_UNSIGNED_CHAR;
}

template<>
MPI_Datatype mpiType<short>()
{
  return MPI_SHORT;
}

template<>
MPI_Datatype mpiType<unsigned short>()
{
  return MPI_UNSIGNED_SHORT;
}

template<>
MPI_Datatype mpiType<int>()
{
  return MPI_INT;
}

template<>
MPI_Datatype mpiType<unsigned int>()
{
  return MPI_UNSIGNED;
}

template<>
MPI_Datatype mpiType<long>()
{
  return MPI_LONG;
}

template<>
MPI_Datatype mpiType<unsigned long>()
{
  return MPI_UNSIGNED_LONG;
}

template<>
MPI_Datatype mpiType<long long>()
{
  return MPI_LONG_LONG;
}

template<>
MPI_Datatype mpiType<unsigned long long>()
{
  return MPI_UNSIGNED_LONG_LONG;
}

template<>
MPI_Datatype mpiType<float>()
{
  return MPI_FLOAT;
}

template<>
MPI_Datatype mpiType<double>()
{
  return MPI_DOUBLE;
}

template<>
MPI_Datatype mpiType<long double>()
{
  return MPI_LONG_DOUBLE;
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

int mpiOrder(Layout layout)
{
  int result = MPI_ORDER_FORTRAN;
  if (layout == C_ORDER) result = MPI_ORDER_C;
  return result;
}

#endif

} // end namespace Domi
