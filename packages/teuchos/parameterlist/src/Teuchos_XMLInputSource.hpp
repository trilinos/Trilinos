// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XMLINPUTSOURCE_H
#define TEUCHOS_XMLINPUTSOURCE_H

/*! \file Teuchos_XMLInputSource.hpp
    \brief A base class for defining a source of XML input.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLInputStream.hpp"

namespace Teuchos
{
  /**
   * \brief XMLInputSource represents a source of XML input that can be parsed
   * to produce an XMLObject.
   *
   * \note
   *	<ul>
   *	<li>The source might be a file, a socket, a
   * std::string. The XMLObject is created with a call to the getObject() method.
   *
   *    <li> The source gets its data from a XMLInputStream object that is
   * created (internally) to work with this source.
   *
   *    <li> getObject() is implemented with EXPAT if Teuchos is configured with
   * <tt>--enable-expat</tt>.
   *	</ul>
   */
  class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLInputSource
    {
    public:
      /** \brief Empty constructor */
      XMLInputSource(){;}

      /** \brief Destructor */
      virtual ~XMLInputSource(){;}

      /** \brief Virtual input source interface */
      virtual RCP<XMLInputStream> stream() const = 0 ;

      /** \brief Get an object by invoking the TreeBuildingXMLHandler on the
       * input data */
      XMLObject getObject() const ;

    };

}
#endif

