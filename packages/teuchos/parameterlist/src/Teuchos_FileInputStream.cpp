// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FileInputStream.hpp"
#include "Teuchos_Assert.hpp"

using namespace Teuchos;

FileInputStream::FileInputStream(const std::string& filename)
	: XMLInputStream(), file_(std::fopen(filename.c_str(), "rb"))
{
  TEUCHOS_TEST_FOR_EXCEPTION(file_ == NULL,
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

	TEUCHOS_TEST_FOR_EXCEPTION(
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
