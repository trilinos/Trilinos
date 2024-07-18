// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_LIBXML2HANDLERADAPTER_H
#define TEUCHOS_LIBXML2HANDLERADAPTER_H

/*! \file Teuchos_Libxml2HandlerAdapter.hpp
    \brief libxml2 adapter for the TreeBuildingXMLHandler
*/

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOS_LIBXML2
#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_RCP.hpp"

#include <libxml/parser.h>

extern "C"
{
  /** \ingroup libXML2 callback for start of an XML element. */
  void xmlSAX2StartElement(void* context,
                           const xmlChar* name,
                           const xmlChar** attr);

  /** \ingroup libXML2 callback for end of an XML element. */
  void xmlSAX2EndElement(void* context,
                         const xmlChar* name);

  /** \ingroup libXML2 callback for character data. */
  void xmlSAX2Characters(void* context,
                         const xmlChar* s,
                         int len);
};

#endif


#endif
