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
  attributes_[name] = value;
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
      
  for (Map::const_iterator i=attributes_.begin(); i!=attributes_.end(); ++i)
		{
			rtn += " " + (*i).first + "=\"" + (*i).second + "\"";
		}

  if (content_.length()==0 && children_.length()==0) 
    {
      rtn += "/>";
    }
  else
    {
      rtn += ">";
    }
	return rtn;
}

string XMLObjectImplem::toString() const
{
  string rtn = header() + "\n";
      
  
  if (content_.length()==0 && children_.length()==0) 
    {
      return rtn;
    }
  else
    {
      bool allBlankContent = true;
      for (int i=0; i<content_.length(); i++)
        {
          if (!StrUtils::isWhite(content_[i])) 
            {
              allBlankContent=false;
              break;
            }
        }
      if (allBlankContent)
        {
          for (int i=0; i<content_.length(); i++)
            {
              rtn += content_[i] + "\n";
            }
        }
      for (int i=0; i<children_.length(); i++)
        {
          rtn += children_[i].toString();
        }
      rtn += "</" + tag_ + ">\n";
    }
  return rtn;
}

void XMLObjectImplem::print(ostream& os, int indent) const
{
  for (int i=0; i<indent; i++) os << " ";
  
  os << header() << endl;
  
  if (content_.length()==0 && children_.length()==0) 
    {
      return;
    }
  else
    {
      printContent(os, indent+2);

      for (int i=0; i<children_.length(); i++)
        {
          children_[i].print(os, indent+2);
        }
      for (int i=0; i<indent; i++) os << " ";
      os << "</" << tag_ << ">\n";
    }
}

void XMLObjectImplem::printContent(ostream& os, int indent) const 
{
  string space = "";
  for (int i=0; i<indent; i++) space += " ";

  bool allBlankContent = true;
  for (int i=0; i<content_.length(); i++)
    {
      if (!StrUtils::isWhite(content_[i])) 
        {
          allBlankContent=false;
          break;
        }
    }
  
  if (allBlankContent)
    {
      for (int i=0; i<content_.length(); i++)
        {
          os << space << content_[i];
        }
    }
}

