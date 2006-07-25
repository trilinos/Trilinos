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

#include "Teuchos_ExpatHandlerAdapter.hpp"

#ifdef HAVE_TEUCHOS_EXPAT

#include "Teuchos_TreeBuildingXMLHandler.hpp"

using namespace Teuchos;

void expatStartElementHandler(void* handler, 
															const XML_Char* name, 
															const XML_Char** attr)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	string tag = name;
	Teuchos::map<string, string> attributes;
	
	/* the attribute data is stored in a C array of C strings, in order 
	 * {key1, val1, key2, val2, ...}. */

	for (int i=0; attr[i] != 0; i+=2)
		{
			string key = attr[i];
			string val = attr[i+1];
			attributes[key] = val;
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
	h->characters(chars);
  delete [] str;
}

#endif // HAVE_TEUCHOS_EXPAT
