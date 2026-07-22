// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
