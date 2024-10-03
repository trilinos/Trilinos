// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STRINGINPUTSTREAM_H
#define TEUCHOS_STRINGINPUTSTREAM_H

/*! \file Teuchos_StringInputStream.hpp
    \brief Definition of XMLInputStream derived class for reading XML from a std::string
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputStream.hpp"


namespace Teuchos
{
  using std::string;

  /**
   * \brief Instantiation of XMLInputStream for reading an entire document from a std::string
   *
   * This is a low-level object and should not be needed at the user level.
   * FileInputSource is the user-level object.
   */

  class StringInputStream : public XMLInputStream
    {
    public:

      //! Construct with the std::string from which data will be read
      StringInputStream(const std::string& text)
        : XMLInputStream(), text_(text), pos_(0) {;}

      //! Destructor
      virtual ~StringInputStream() {;}

      //! Read up to maxToRead bytes
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead);

    private:
      std::string text_;
      unsigned int pos_;
    };

}
#endif

