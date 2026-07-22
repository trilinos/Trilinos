// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Libxml2HandlerAdapter.hpp"

#ifdef HAVE_TEUCHOS_LIBXML2

#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace Teuchos;

void xmlSAX2StartElement(void* handler,
                         const xmlChar* name,
                         const xmlChar** attr)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	std::string tag = (const char*) name;
	Teuchos::map<std::string, std::string> attributes;
	
	/* the attribute data is stored in a C array of C strings, in order
	 * {key1, val1, key2, val2, ...}. */

	for (int i=0; attr[i] != 0; i+=2)
		{
			std::string key = (const char*) attr[i];
			std::string val = (const char*) attr[i+1];
			attributes[key] = val;
		}

	h->startElement(tag, attributes);
}

void xmlSAX2EndElement(void* handler,
                       const xmlChar* name)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	std::string tag = (const char*) name;

	h->endElement(tag);
}

void xmlSAX2CharacterData(void* handler,
                          const xmlChar* s,
                          int len)
{
	char* str = new char[len+1];
	strncpy(str, (const char*) s, len);


  str[len] = '\0';
	std::string chars = str;

	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	h->characters(chars);
  delete [] str;
}

#endif
