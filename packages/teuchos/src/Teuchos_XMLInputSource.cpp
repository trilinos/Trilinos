#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_TreeBuildingXMLHandler.hpp"


#ifdef HAVE_XERCES
#include "XercesHandlerAdapter.h"
#include "XercesInputSourceAdapter.h"
#include "XercesInputStreamAdapter.h"
#include <util/PlatformUtils.hpp>
#endif

#ifdef HAVE_EXPAT
#include "Teuchos_ExpatHandlerAdapter.hpp"
#define EXPAT_BUFSIZE 8192
#endif


using namespace Teuchos;

XMLObject XMLInputSource::getObject() const
{
#ifdef HAVE_XERCES

	static bool first = true;
	if (first)
		{
			XMLPlatformUtils::Initialize();
			first = false;
		}

	SAXParser parser;
	XercesHandlerAdapter handler(new TreeBuildingXMLHandler());
	XercesInputSourceAdapter inputSource(this);

	parser.setDocumentHandler(&handler);

	try 
		{
			parser.parse(inputSource);
		}
	catch(exception& e)
		{
			TEST_FOR_EXCEPTION(true, runtime_error, 
                         "exception detected in SAX parsing: " << e.what());
		}
	return handler.getObject();
#endif

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
			
			
