//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

#include "MatrixMarket_generic.hpp"
#include "MatrixMarket_split.hpp"
#include <algorithm>

namespace Tpetra {
  namespace MatrixMarket {

    int maxLineLength() { return 1024; }

    bool 
    checkCommentLine (const std::string& line, 
		      size_t& start, 
		      size_t& size,
		      const size_t lineNumber,
		      const bool tolerant)
    {
      // In tolerant mode, empty lines are considered comment lines.
      if (line.empty()) 
	{
	  if (tolerant)
	    return true;
	  else
	    {
	      std::ostringstream os;
	      os << "Line " << lineNumber << " contains no characters";
	      throw std::invalid_argument (os.str());
	    }
	}
      // The line of comments or data "starts" after any whitespace
      // characters.  Whitespace-only lines are considered "empty."
      start = line.find_first_not_of (" \t");
      if (start == std::string::npos)
	{ // It's a whitespace-only line
	  if (tolerant)
	    return true;
	  else
	    {
	      std::ostringstream os;
	      os << "Line " << lineNumber << " contains only whitespace";
	      throw std::invalid_argument (os.str());
	    }
	}
      // Position of the first comment character, if any.
      const size_t commentPos = line.find_first_of("%#", start);
      if (commentPos == std::string::npos)
	{ // There are no comment characters in the line.
	  // line.substr(start,npos) gives the substring of line
	  // containing valid data.
	  size = std::string::npos;
	  return false;
	}
      else 
	{ // [start, start+size-1] is the (inclusive) range of
	  // characters (if any) between the first nonwhitespace
	  // character, and the first comment character.
	  size = commentPos - start;
	  // It's not a "comment line," because there could be valid
	  // data before the first comment character.
	  return false;
	}
    }

  } // namespace MatrixMarket
} // namespace Tpetra
