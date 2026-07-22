// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ExpatHandlerAdapter.hpp"

#ifdef HAVE_TEUCHOS_EXPAT

#include "Teuchos_TreeBuildingXMLHandler.hpp"

using namespace Teuchos;

void expatStartElementHandler(void* handler,
															const XML_Char* name,
															const XML_Char** attr)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	std::string tag = name;
	Teuchos::map<std::string, std::string> attributes;
	
	/* the attribute data is stored in a C array of C strings, in order
	 * {key1, val1, key2, val2, ...}. */

	for (int i=0; attr[i] != 0; i+=2)
		{
			std::string key = attr[i];
			std::string val = attr[i+1];
			attributes[key] = val;
		}

	h->startElement(tag, attributes);
}

void expatEndElementHandler(void* handler,
														const XML_Char* name)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	std::string tag = name;

	h->endElement(tag);
}

void expatCharacterDataHandler(void* handler,
															 const XML_Char* s,
															 int len)
{
	char* str = new char[len+1];
	strncpy(str, s, len);


  str[len] = '\0';
	std::string chars = str;

	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	h->characters(chars);
  delete [] str;
}

#endif // HAVE_TEUCHOS_EXPAT
