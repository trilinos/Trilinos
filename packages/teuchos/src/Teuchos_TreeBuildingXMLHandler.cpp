#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;

TreeBuildingXMLHandler::TreeBuildingXMLHandler()
	: root_(), current_(), path_()
{ 
	current_ = root_;
}

void TreeBuildingXMLHandler::characters(const string& chars, 
																				const unsigned int length)
{
	try
		{
			if (StrUtils::isWhite(chars)) return;
			if (!current_.isEmpty())
				{
					current_.addContent(StrUtils::fixUnprintableCharacters(chars));
				}
			else
				{
					Error::raise("trying to add content to an empty node");
				}
		}
	catch(exception& e)
		{
			Error::trace(e, "in TreeBuildingXMLHandler::characters");
		}
}

void TreeBuildingXMLHandler::startElement(const string& tag, 
																					const Hashtable<string,string>& attributes)
{
	try
		{
			XMLObject parent;

			try
				{
					if (current_.isEmpty()) 
						{
							root_ = XMLObject("root");
							current_ = root_;
						}
					parent = current_;
					path_.push(current_);
					current_ = XMLObject(tag);
					parent.addChild(current_);
				}
			catch(exception& e)
				{
					Error::trace(e, "in creating current and parent");
				}

			try
				{
					Array<string> keys;
					Array<string> vals;
					attributes.arrayify(keys, vals);
					for (int i=0; i<keys.length(); i++)
						{
							current_.addAttribute(keys[i], vals[i]);
						}
				}
			catch(exception& e)
				{
					Error::trace(e, "in assigning attributes");
				}
		}
	catch(exception& e)
		{
			Error::trace(e, "in TreeBuildingXMLHandler::startElement(tag=" + tag + ")");
		}
}

void TreeBuildingXMLHandler::endElement(const string& tag)
{
	try
		{
			if (path_.size() > 0)
				{
					current_ = path_.pop();
				}
		}
	catch(exception& e)
		{
			Error::trace(e, "in TreeBuildingXMLHandler::endElement(tag=" + tag + ")");
		}
}

