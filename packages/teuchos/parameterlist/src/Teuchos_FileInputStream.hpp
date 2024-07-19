// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_FILEINPUTSTREAM_H
#define TEUCHOS_FILEINPUTSTREAM_H

/*! \file Teuchos_FileInputStream.hpp
    \brief Definition of XMLInputStream derived class for reading XML from a file
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputStream.hpp"

namespace Teuchos
{
  using std::string;

  /** \ingroup XML
   * \brief Instantiation of XMLInputStream class for reading an entire document from a file.
   * \note This is a low-level object and should not be needed at the user level.
   * FileInputSource is the user-level object.
   */
  class FileInputStream : public XMLInputStream
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

