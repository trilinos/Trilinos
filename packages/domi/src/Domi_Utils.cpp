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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Domi_Utils.hpp"

namespace Domi
{

Teuchos::Array< int > computeAxisRanks(int rank,
                                       Teuchos::ArrayView< int > axisSizes)
{
  Teuchos::Array< int > result(axisSizes.size());
  int relRank = rank;
  int stride = 1;
  for (int axis = 0; axis < axisSizes.size()-1; ++axis)
    stride *= axisSizes[axis];
  for (int axis = axisSizes.size()-1; axis > 0; --axis)
  {
    result[axis] = relRank / stride;
    relRank      = relRank % stride;
    stride       = stride / axisSizes[axis-1];
  }
  result[0] = relRank;
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::Array< int > computeAxisRanks(int rank,
                                       int offset,
                                       Teuchos::ArrayView< int > axisStrides)
{
  Teuchos::Array< int > result(axisStrides.size());
  int relRank = rank - offset;
  for (int axis = axisStrides.size()-1; axis >= 0; --axis)
  {
    result[axis] = relRank / axisStrides[axis];
    relRank      = relRank % axisStrides[axis];
  }
  return result;
}

} // end namespace Domi
