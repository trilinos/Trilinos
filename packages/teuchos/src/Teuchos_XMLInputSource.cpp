// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
			
			
