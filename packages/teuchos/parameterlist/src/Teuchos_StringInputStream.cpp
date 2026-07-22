// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StringInputStream.hpp"

using namespace Teuchos;


unsigned int StringInputStream::readBytes(unsigned char* const toFill,
																					const unsigned int maxToRead)
{
	if (pos_ == text_.length()) return 0;
	
	int toRead = static_cast<int>(text_.length() - pos_);
	if ((int) maxToRead < toRead) toRead = maxToRead;

  std::strncpy((char*) toFill, text_.c_str()+pos_, toRead);

	pos_ += toRead;
	
	return (size_t) toRead;
}

