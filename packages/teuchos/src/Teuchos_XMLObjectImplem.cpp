#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;


XMLObjectImplem::XMLObjectImplem(const string& tag)
	: tag_(tag), attributes_(), children_(0), content_(0)
{;}

XMLObjectImplem* XMLObjectImplem::deepCopy() const 
{
	XMLObjectImplem* rtn = new XMLObjectImplem(tag_);
	if (rtn==0) Error::raise("XMLObjectImplem::deepCopy()");
	rtn->attributes_ = attributes_;
	rtn->content_ = content_;
	
	for (int i=0; i<children_.length(); i++)
		{
			rtn->addChild(children_[i].deepCopy());
		}

	return rtn;
}

int XMLObjectImplem::numChildren() const {return children_.length();}

void XMLObjectImplem::addAttribute(const string& name, const string& value)
{
	try
		{
			attributes_.put(name, value);
		}
	catch(exception& e)
		{
			Error::trace(e,  "XMLObjectImplem::addAttribute(), adding attribute " + name);
		}
}

void XMLObjectImplem::addChild(const XMLObject& child)
{
	try
		{
			children_.append(child);
		}
	catch(exception& e)
		{
			Error::trace(e,  "XMLObjectImplem::addChild(), adding child " 
                   + child.getTag());
		}
}

void XMLObjectImplem::addContent(const string& contentLine)
{
	try
		{
			content_.append(contentLine);
		}
	catch(exception& e)
		{
			Error::trace(e,  
                   "in XMLObjectImplem::addContent(), adding content line " 
                   + contentLine);
		}
}

const XMLObject& XMLObjectImplem::getChild(int i) const 
{
	return children_[i];
}

string XMLObjectImplem::header() const
{
	string rtn = "<" + tag_;
      
	Array<string> names;
	Array<string> values;
	attributes_.arrayify(names, values);

	for (int i=0; i<names.length(); i++)
		{
			rtn += " " + names[i] + "=\"" + values[i] + "\"";
		}
	rtn += ">";
	return rtn;
}

string XMLObjectImplem::toString() const
{
  try
    {
      string rtn = "<" + tag_;
      
      Array<string> names;
      Array<string> values;
      attributes_.arrayify(names, values);
      int i = 0;
      for (i=0; i<names.length(); i++)
				{
					rtn += " " + names[i] + "=\"" + values[i] + "\"";
				}
      if (content_.length()==0 && children_.length()==0) 
				{
					rtn += "/>\n" ;
				}
      else
				{
					rtn += ">\n";
          bool allBlankContent = true;
          for (i=0; i<content_.length(); i++)
						{
              if (!StrUtils::isWhite(content_[i])) 
                {
                  allBlankContent=false;
                  break;
                }
						}
          if (allBlankContent)
            {
              for (i=0; i<content_.length(); i++)
                {
                  rtn += content_[i] + "\n";
                }
            }
					for (i=0; i<children_.length(); i++)
						{
							rtn += children_[i].toString();
						}
					rtn += "</" + tag_ + ">\n";
				}
      return rtn;
    }
  catch(exception& e)
    {
      Error::trace(e,  "XMLObjectImplem::toString(), element tag=" + tag_);
    }
  return ""; // -Wall
}

