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

// Teuchos includes
#include "Teuchos_OrdinalTraits.hpp"

// Domi includes
#include "Domi_Slice.hpp"

namespace Domi
{

// Default value
const Ordinal Slice::Default(Teuchos::OrdinalTraits<Ordinal>::max());

////////////////////////////////////////////////////////////////////////

Slice Slice::bounds(Ordinal size) const
{
  // If this Slice is already bounded and the size is large enough,
  // just return a copy of this Slice.  This is done for performance
  // reasons to avoid, when we can, all of the conditionals below.
  if (_bounded_pos && size >= _stop ) return *this;
  if (_bounded_neg && size >  _start) return *this;
  // Initialize the lo and hi indexes
  Ordinal lo = _start;
  Ordinal hi = _stop;
  // Convert Default values to concrete indexes
  if (lo == Default) lo = (_step < 0) ? size-1 :  0;
  if (hi == Default) hi = (_step > 0) ? size   : -1;
  // Convert negative values to nonnegative indexes
  while (lo < 0) lo += size;
  while (hi < 0) hi += size;
  // Clip too-large values
  if ((_step < 0) && (lo > size)) lo = size;
  if ((_step > 0) && (hi > size)) hi = size;
  // Fine-tune hi to be a precise stopping index
  Ordinal numSteps = (hi - lo) / _step;
  if ((hi - lo) % _step) numSteps += 1;
  hi = lo + _step * numSteps;
  // Return the new slice representing the true bounds.  This is
  // returned as a ConcreteSlice, which has a optimally efficient
  // bounds() method for future calls.
  return ConcreteSlice(lo, hi, _step);
}

////////////////////////////////////////////////////////////////////////

std::string Slice::toString() const
{
  std::stringstream ss;
  ss << "[";
  if (_step > 0 && _start != 0      ) ss << _start;
  if (_step < 0 && _start != Default) ss << _start;
  ss << ":";
  if (_step > 0 && _stop != Default) ss << _stop;
  if (_step < 0 && _stop != 0      ) ss << _stop;
  if (_step != 1)                    ss << ":" << _step;
  ss << "]";
  return ss.str();
}

////////////////////////////////////////////////////////////////////////

std::ostream & operator<<(std::ostream & os,
                          const Slice & slice)
{
  return os << slice.toString();
}

////////////////////////////////////////////////////////////////////////

ConcreteSlice::ConcreteSlice(Ordinal stopVal) :
  Slice(0,stopVal,1)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (stop() < 0), InvalidArgument,
    "ConcreteSlice stop value cannot be negative"
    );
}

////////////////////////////////////////////////////////////////////////

ConcreteSlice::ConcreteSlice(Ordinal startVal,
                             Ordinal stopVal,
                             Ordinal stepVal) :
  Slice(startVal,stopVal,stepVal)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (start() < 0), InvalidArgument,
    "ConcreteSlice start value cannot be negative"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    (stop() < 0), InvalidArgument,
    "ConcreteSlice stop value cannot be negative"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    (step() == 0), InvalidArgument,
    "ConcreteSlice step interval cannot be zero"
    );
}

}  // namespace Domi
