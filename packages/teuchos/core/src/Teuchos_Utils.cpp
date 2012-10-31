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

#include "Teuchos_Utils.hpp"
#include "Teuchos_GlobalMPISession.hpp"

namespace Teuchos {

double Utils::chopVal_ = 1.0e-16;

double Utils::chop(const double& x) 
{
	if (std::fabs(x) < chopVal_) return 0;
	return x;
}

std::string Utils::trimWhiteSpace( const std::string& str )
{
  typedef std::string::size_type size_type;
  const size_type len = str.length();
  if (len==0) {
      return str;
    }
  size_type first_non_white = 0;
  for(
    first_non_white = 0 ;
    isWhiteSpace(str[first_non_white]) && first_non_white < len ;
    ++first_non_white
    );
  // Above, if only whitespace is found, then first_non_white==len on
  // termination of the loop!
  size_type last_non_white = 0;
  for(
    last_non_white = len-1 ;
    isWhiteSpace(str[last_non_white]) && (last_non_white != 0);
    --last_non_white
    );
  // Above, if only whitespace is found, last_non_white==0 on termination of
  // the loop!
  if( first_non_white > last_non_white )
    return std::string(""); // The std::string is all whitespace!
  return str.substr(first_non_white,last_non_white-first_non_white+1);
}

std::string Utils::toString(const int& x)
{
	char s[100];
	std::sprintf(s, "%d", x);
	return std::string(s);
}

std::string Utils::toString(const unsigned int& x)
{
	char s[100];
	std::sprintf(s, "%d", x);
	return std::string(s);
}

std::string Utils::toString(const double& x)
{
	char s[100];
	std::sprintf(s, "%g", x);
	return std::string(s);
}

std::string Utils::getParallelExtension(
  int         procRank_in
  ,int        numProcs_in
  )
{

  int        procRank = -1;
  int        numProcs = -1;
  if( numProcs_in > 0 ) {
    procRank = procRank_in;
    numProcs = numProcs_in;
  }
  else {
    procRank = Teuchos::GlobalMPISession::getRank();
    numProcs = Teuchos::GlobalMPISession::getNProc();
  }

  int maxProcOrder = 1;
  double tmp = numProcs;
  for( int i = 0; i < 10; ++i, tmp *= 0.1 ) {
    if(tmp >= 1.0)
      ++maxProcOrder;
    else
      break;
  }

  std::ostringstream parallelExtension;
  parallelExtension
    << std::setfill('0')
    << std::right << std::setw(maxProcOrder)
    << numProcs
    << "."
    << std::setfill('0')
    << std::right << std::setw(maxProcOrder)
    << procRank;
  return parallelExtension.str();
}

} // end namespace Teuchos
