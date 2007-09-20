// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
