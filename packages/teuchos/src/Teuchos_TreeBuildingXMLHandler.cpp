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

#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

TreeBuildingXMLHandler::TreeBuildingXMLHandler()
	: root_(), current_(), path_()
{ 
	current_ = root_;
}

void TreeBuildingXMLHandler::characters(const string& chars, 
																				const unsigned int length)
{
  if (StrUtils::isWhite(chars)) return;
  TEST_FOR_EXCEPTION(current_.isEmpty(), std::logic_error,
                     "TreeBuildingXMLHandler::trying to add content to an empty node");
                        
  current_.addContent(StrUtils::fixUnprintableCharacters(chars));
}

void TreeBuildingXMLHandler::startElement(const string& tag, 
																					const Map<string,string>& attributes)
{
  XMLObject parent;
  
  if (current_.isEmpty()) 
    {
      root_ = XMLObject("root");
      current_ = root_;
    }
  parent = current_;
  path_.push(current_);
  current_ = XMLObject(tag);
  parent.addChild(current_);

  for (Map::const_iterator i=attributes.begin(); i != attributes.end(); ++i)
    {
      const string& key = (*i).first;
      const string& val = (*i).second;
      current_.addAttribute(key, val);
    }
}

void TreeBuildingXMLHandler::endElement(const string& tag)
{
  if (path_.size() > 0)
    {
      current_ = path_.top();
      path_.pop();
    }
}

