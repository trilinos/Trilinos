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

#include "Teuchos_FileInputStream.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

FileInputStream::FileInputStream(const std::string& filename)
	: XMLInputStream(), file_(std::fopen(filename.c_str(), "rb"))
{
  TEST_FOR_EXCEPTION(file_ == NULL,
                     std::runtime_error,
                     "FileInputStream ctor failed to open file: " 
                     << filename);
}

unsigned int FileInputStream::readBytes(unsigned char* const toFill, 
																				const unsigned int maxToRead)
{
	if (
#if defined(ICL) || defined(__sgi)
    feof(file_)
#else
    std::feof(file_)
#endif
    )
    return (size_t)0;
	int n = std::fread((void*) toFill, sizeof(char), maxToRead, file_);
  if (n==0) return (size_t)0; 

  const bool
    is_eof
#if defined(ICL) || defined(__sgi)
    = feof(file_)
#else
    = std::feof(file_)
#endif
    ;

	TEST_FOR_EXCEPTION(
    n < 0 || (n<(int) maxToRead && !is_eof),
    std::runtime_error,
    "FileInputStream::readBytes error"
    );
	
	return (size_t) n;
}

// 2007/07/08: rabartl: Above, for some reason the header <cstdio> on the
// Intel C++ compiler 8.0 for Windows XP just does not have the function
// feof(...) in the correct std namespace.  On this compiler, you just have to
// access it from the global namespace.
