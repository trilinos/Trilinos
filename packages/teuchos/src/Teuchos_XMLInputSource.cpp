// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_TreeBuildingXMLHandler.hpp"


#ifdef HAVE_EXPAT
#include "Teuchos_ExpatHandlerAdapter.hpp"
#define EXPAT_BUFSIZE 8192
#endif

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#endif

using namespace Teuchos;

XMLObject XMLInputSource::getObject() const
{
#ifdef HAVE_EXPAT

	RefCountPtr<TreeBuildingXMLHandler> handler = rcp(new TreeBuildingXMLHandler());

	XML_Parser parser = XML_ParserCreate(NULL);

	XML_SetElementHandler(parser, expatStartElementHandler, 
												expatEndElementHandler);

	XML_SetCharacterDataHandler(parser, expatCharacterDataHandler);

	XML_SetUserData(parser, (void*) &(*handler));

	RefCountPtr<XMLInputStream> s = stream();

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

  TEST_FOR_EXCEPTION(true, logic_error, "XMLInputSource::getObject() - no XML parser installed");
  return XMLObject(); // -Wall
#endif

}
			
			
