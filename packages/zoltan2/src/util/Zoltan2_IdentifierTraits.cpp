// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_IdentifierTraits.cpp
   \brief Defines basic traits for user global identifiers.
*/

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_AlltoAll.hpp>
#include <Zoltan2_IdentifierTraits.hpp>

#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>
#include <Teuchos_ReductionOp.hpp>

#include <utility>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <cstdlib>


namespace Zoltan2
{

/*! \brief helper to hash values larger than int to an int.
 *  Hash values do not need to be unique, but there should be
 *  as few overlaps as possible.
 */
int getHashCode(const unsigned char *a, size_t len)
{
  int total=0;
  unsigned char *to = reinterpret_cast<unsigned char *>(&total);
  int c=0;
  for (size_t i=0; i < len; i++){
    to[c++] += a[i];
    if (c == sizeof(int))
      c = 0;
  }
  if (total < 0)
    total *= -1;
  return total;
}

} // namespace Z2

