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
																					const Hashtable<string,string>& attributes)
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
  
  Array<string> keys;
  Array<string> vals;
  attributes.arrayify(keys, vals);
  for (int i=0; i<keys.length(); i++)
    {
      current_.addAttribute(keys[i], vals[i]);
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

