// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_Assert.hpp"

using namespace Teuchos;

TreeBuildingXMLHandler::TreeBuildingXMLHandler()
	: root_(), current_(), path_()
{
	current_ = root_;
}

void TreeBuildingXMLHandler::characters(const std::string& chars)
{
  TEUCHOS_TEST_FOR_EXCEPTION(current_.isEmpty(), std::logic_error,
                     "TreeBuildingXMLHandler::trying to add content to an empty node");

  size_t k = current_.numContentLines();
  // a new line indicates the start of a new ContentLine. Only add it if current line is not empty.
  if(chars.compare("\n")==0) {
    if((k>0) && (current_.getContentLine(k-1).length()>0))
      current_.addContent("");
  } else {
    // If no contentLine exists to append create one, else add to the last contentLine
    // Do not add white characters at the start of a new line
    if(k==0) {
      if(!StrUtils::isWhite(chars))
        current_.addContent(StrUtils::fixUnprintableCharacters(chars));
    } else {
      if((!StrUtils::isWhite(chars)) || (current_.getContentLine(k-1).length()>0))
        current_.appendContentLine(k-1,StrUtils::fixUnprintableCharacters(chars));
    }
  }
}

void TreeBuildingXMLHandler::startElement(const std::string& tag,
																					const Map& attributes)
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
      const std::string& key = (*i).first;
      const std::string& val = (*i).second;
      current_.addAttribute(key, val);
    }
}

int TreeBuildingXMLHandler::endElement(const std::string& tag)
{
  int error = 0;
  if (path_.size() > 0)
    {
      if (current_.getTag() != tag) {
        error = 1; // error: tags must be balanced
      }
      // remove empty contentLines at the end.
      size_t k = current_.numContentLines();
      while( (k>0) && (current_.getContentLine(--k).length() == 0))
        current_.removeContentLine(k);

      current_ = path_.top();
      path_.pop();
    }
  else
    {
      error = 1; // error: cannot end element that wasn't started
    }
  return error;
}

