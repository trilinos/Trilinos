//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_ApplyType.hpp>
#include <stdexcept>


namespace TSQR {
  ApplyType::ApplyType (const std::string& op) :
    type_ (decide_apply_type (op)),
    lapackString_ (enumToLapackString (type_))
  {}

  ApplyType::ApplyType (const ApplyType& rhs) :
    type_ (rhs.type_),
    lapackString_ (rhs.lapackString_)
  {}

  ApplyType& ApplyType::operator= (const ApplyType& rhs) {
    type_ = rhs.type_;
    lapackString_ = rhs.lapackString_;
    return *this;
  }

  const ApplyType ApplyType::NoTranspose = ApplyType ("N");
  const ApplyType ApplyType::Transpose = ApplyType ("T");
  const ApplyType ApplyType::ConjugateTranspose = ApplyType ("C");

  std::string 
  ApplyType::enumToLapackString (const ApplyType::ApplyType_ theType)
  {
    if (theType == NoTranspose_)
      return std::string("N");
    else if (theType == Transpose_)
      return std::string("T");
    else if (theType == ConjugateTranspose_)
      return std::string("C");
    else
      throw std::logic_error("Invalid ApplyType: should never get here");
  }

  bool 
  ApplyType::decide_transposed (const std::string& op) const 
  {
    if (op[0] == 'N' || op[0] == 'n')
      return false;
    else
      {
	const char validTransposeLetters[] = "TtCcHh";
	const int numValidTransposeLetters = 6;

	for (int k = 0; k < numValidTransposeLetters; ++k)
	  if (op[0] == validTransposeLetters[k])
	    return true;

	throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
      }
  }

  ApplyType::ApplyType_
  ApplyType::decide_apply_type (const std::string& op) const 
  {
    if (op[0] == 'T' || op[0] == 't')
      return Transpose_;
    else if (op[0] == 'N' || op[0] == 'n')
      return NoTranspose_;
    else if (op[0] == 'C' || op[0] == 'c' || op[0] == 'H' || op[0] == 'h')
      return ConjugateTranspose_;
    else
      throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
  }
}

