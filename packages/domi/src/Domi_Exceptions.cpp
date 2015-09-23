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

#include "Domi_Exceptions.hpp"

namespace Domi
{

////////////////////////////////////////////////////////////////////////

InvalidArgument::InvalidArgument(std::string msg) :
  std::invalid_argument(msg)
{
}

////////////////////////////////////////////////////////////////////////

RangeError::RangeError(std::string msg) :
  std::range_error(msg)
{
}

////////////////////////////////////////////////////////////////////////

SubcommunicatorError::SubcommunicatorError(std::string op) :
  std::domain_error(op + " attempted on a processor that does "
                    "not belong to the operative subcommunicator")
{
}

////////////////////////////////////////////////////////////////////////

MapOrdinalError::MapOrdinalError(std::string msg) :
  std::runtime_error(msg)
{
}

////////////////////////////////////////////////////////////////////////

MDMapError::MDMapError(std::string msg) :
  std::runtime_error(msg)
{
}

////////////////////////////////////////////////////////////////////////

MDMapNoncontiguousError::MDMapNoncontiguousError(std::string msg) :
  std::runtime_error(msg)
{
}

////////////////////////////////////////////////////////////////////////

TypeError::TypeError(std::string msg) :
  std::runtime_error(msg)
{
}

////////////////////////////////////////////////////////////////////////

BoundsError::BoundsError(std::string msg) :
  std::runtime_error(msg)
{
}

}
