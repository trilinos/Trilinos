// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
	std::snprintf(s, sizeof(s), "%d", x);
	return std::string(s);
}

std::string Utils::toString(const long long& x)
{
	char s[100];
	std::snprintf(s, sizeof(s), "%lld", x);
	return std::string(s);
}

std::string Utils::toString(const unsigned int& x)
{
	char s[100];
	std::snprintf(s, sizeof(s), "%d", x);
	return std::string(s);
}

std::string Utils::toString(const double& x)
{
	char s[100];
	std::snprintf(s, sizeof(s), "%g", x);
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
