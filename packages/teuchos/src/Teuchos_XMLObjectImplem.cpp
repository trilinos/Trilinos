#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;


XMLObjectImplem::XMLObjectImplem(const string& tag)
	: tag_(tag), attributes_(), children_(0), content_(0)
{;}

XMLObjectImplem* XMLObjectImplem::deepCopy() const 
{
	XMLObjectImplem* rtn = new XMLObjectImplem(tag_);
	TEST_FOR_EXCEPTION(rtn==0, runtime_error, "XMLObjectImplem::deepCopy()");
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
  attributes_.put(name, value);
}

void XMLObjectImplem::addChild(const XMLObject& child)
{
  children_.append(child);
}

void XMLObjectImplem::addContent(const string& contentLine)
{
  content_.append(contentLine);
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

