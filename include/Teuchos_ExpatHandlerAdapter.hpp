// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_EXPATHANDLERADAPTER_H
#define TEUCHOS_EXPATHANDLERADAPTER_H

/*! \file Teuchos_ExpatHandlerAdapter.hpp
    \brief Expat adapter for the TreeBuildingXMLHandler
*/

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPAT

#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_RCP.hpp"

#include "expat.h"

extern "C"
{
  /** \ingroup XML Expat callback for start of an XML element. */
  void expatStartElementHandler(void* handler,
                                const XML_Char* name,
                                const XML_Char** attr);

  /** \ingroup XML Expat callback for end of an XML element. */
  void expatEndElementHandler(void* handler,
                              const XML_Char* name);

  /** \ingroup XML Expat callback for character data. */
  void expatCharacterDataHandler(void* handler,
                                 const XML_Char* s,
                                 int len);
}

#endif // HAVE_TEUCHOS_EXPAT

#endif // TEUCHOS_EXPATHANDLERADAPTER_H
