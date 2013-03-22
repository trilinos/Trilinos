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

#ifndef TEUCHOS_FILESTREAM_H
#define TEUCHOS_FILESTREAM_H

//! \file Teuchos_FILEstream.hpp

#include <streambuf>

namespace Teuchos
{

  //! Teuchos::FILEstream: Combined C FILE and C++ stream

  /*! Teuchos::FILEstream is a class that defines an object that is
      simultaneously a C FILE object and a C++ stream object.  The
      utility of this class is in connecting existing C++ code that
      uses streams and C code that uses FILEs.  An important example
      of this situation is the python wrappers for Trilinos packages.
      Trilinos is of course written primarily in C++, but the python
      wrappers must interface to the python C API.  Wrappers for
      Trilinos methods or operators that expect a stream can be given
      a Teuchos::FILEstream, which then behaves as a FILE within the
      python C API.  This is a low-level object that should not be
      needed at the user level.
   */

  class FILEstream : public std::streambuf {

  public:

    //! Constructor

    /*! The only constructor for Teuchos:FILEstream, and it requires a
        pointer to a C FILE struct.
     */
    FILEstream(std::FILE* file): self_file(file) {}

  protected:

    std::streambuf::int_type overflow(std::streambuf::int_type c) {
      return std::fputc(c, self_file) == EOF?
	std::streambuf::traits_type::eof(): c;
    }

    FILE* self_file;
  };
}

#endif
