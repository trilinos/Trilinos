#include "Teuchos_ExpatHandlerAdapter.hpp"

#ifdef HAVE_EXPAT

#include "Teuchos_TreeBuildingXMLHandler.hpp"

using namespace Teuchos;

void expatStartElementHandler(void* handler, 
															const XML_Char* name, 
															const XML_Char** attr)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	string tag = name;
	Hashtable<string, string> attributes;
	
	/* the attribute data is stored in a C array of C strings, in order 
	 * {key1, val1, key2, val2, ...}. */

	for (int i=0; attr[i] != 0; i+=2)
		{
			string key = attr[i];
			string val = attr[i+1];
			attributes.put(key, val);
		}

	h->startElement(tag, attributes);
}

void expatEndElementHandler(void* handler, 
														const XML_Char* name)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	string tag = name;

	h->endElement(tag);
}

void expatCharacterDataHandler(void* handler, 
															 const XML_Char* s,
															 int len)
{
	char* str = new char[len+1];
	strncpy(str, s, len);


  str[len] = '\0';
	string chars = str;

	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	h->characters(chars, chars.length());
  delete [] str;
}

#endif
