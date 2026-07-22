// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FILEINPUTSTREAM_H
#define ROL_FILEINPUTSTREAM_H

/*! \file ROL_FileInputStream.hpp
    \brief Definition of an input stream derived class for reading XML from a file
*/

#include "Teuchos_ConfigDefs.hpp"
#include <cstdio>

namespace ROL
{
  using std::string;

  class FileInputStream
    {
    public:
      /** \brief Construct with a filename */
      FileInputStream(const std::string& filename);

      /** \brief Destructor */
      virtual ~FileInputStream() { if( file_ ) { std::fclose(file_); } }

      /** \brief Read up to maxToRead bytes */
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead);

    private:
      std::FILE* file_;
    };
}
#endif

