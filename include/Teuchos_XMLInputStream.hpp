// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XMLINPUTSTREAM_H
#define TEUCHOS_XMLINPUTSTREAM_H

/*! \file Teuchos_XMLInputStream.hpp
    \brief A base class for defining a XML input stream
*/

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos
{
  /** \brief XMLInputStream represents an XML input stream
  *	that can be used by a XMLInputSource
  */
  class XMLInputStream
    {
    public:
      /** \brief Constructor */
      XMLInputStream(){;}

      /** \brief Destructor */
      inline virtual ~XMLInputStream(){;}

      /** \brief Read up to maxToRead bytes from the stream */
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead) = 0 ;

      /** \brief Identify current position */
      virtual unsigned int curPos() const ;

    };
}
#endif

