// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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

#include "Teuchos_EReductionType.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

#ifdef HAVE_TEUCHOS_MPI
namespace Details {

MPI_Op
getMpiOpForEReductionType (const enum EReductionType reductionType)
{
  switch (reductionType) {
  case REDUCE_SUM: return MPI_SUM;
  case REDUCE_MIN: return MPI_MIN;
  case REDUCE_MAX: return MPI_MAX;
  case REDUCE_AND: return MPI_LAND; // logical AND, not bitwise AND
  case REDUCE_BOR: return MPI_BOR; // bitwise OR
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      "The given EReductionType value is invalid.");
  }
}

} // namespace Details
#endif // HAVE_TEUCHOS_MPI

const char*
toString (const EReductionType reductType)
{
  switch (reductType) {
  case REDUCE_SUM: return "REDUCE_SUM";
  case REDUCE_MIN: return "REDUCE_MIN";
  case REDUCE_MAX: return "REDUCE_MAX";
  case REDUCE_AND: return "REDUCE_AND";
  case REDUCE_BOR: return "REDUCE_BOR";
  default:
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "Teuchos::toString(EReductionType): "
       "Invalid EReductionType value " << reductType << ".  Valid values "
       "are REDUCE_SUM = " << REDUCE_SUM << ", REDUCE_MIN = " << REDUCE_MIN
       << ", REDUCE_MAX = " << REDUCE_MIN << ", REDUCE_AND = " << REDUCE_AND
       << ", and REDUCE_BOR = " << REDUCE_BOR
       << ".");
  }
}

} // namespace Teuchos
