// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_TreeBuildingXMLHandler.hpp"


#ifdef HAVE_TEUCHOS_EXPAT
#  include "Teuchos_ExpatHandlerAdapter.hpp"
#  define EXPAT_BUFSIZE 8192
#else
#  include "Teuchos_XMLParser.hpp"
#endif

using namespace Teuchos;

XMLObject XMLInputSource::getObject() const
{

#ifdef HAVE_TEUCHOS_EXPAT

	RCP<TreeBuildingXMLHandler> handler = rcp(new TreeBuildingXMLHandler());

	XML_Parser parser = XML_ParserCreate(NULL);

	XML_SetElementHandler(parser, expatStartElementHandler,
												expatEndElementHandler);

	XML_SetCharacterDataHandler(parser, expatCharacterDataHandler);

	XML_SetUserData(parser, (void*) &(*handler));

	RCP<XMLInputStream> s = stream();

	bool done = false;
	unsigned int bufsize = EXPAT_BUFSIZE;
	unsigned char buf[EXPAT_BUFSIZE];

	while (!done)
		{
			unsigned int nRead = s->readBytes(buf, bufsize);
			if (nRead < bufsize)
				{
					done = true;
				}
			XML_Parse(parser, (char*) buf, bufsize, done);
		}

	return handler->getObject();

#else

  XMLParser parser(stream());

  return parser.parse();

#endif

}
			
			
