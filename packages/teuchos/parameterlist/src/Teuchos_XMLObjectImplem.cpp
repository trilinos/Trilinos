// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StrUtils.hpp"
#include <cstring>

using namespace Teuchos;


XMLObjectImplem::XMLObjectImplem(const std::string& tag)
  : tag_(tag), attributes_(), children_(0), content_(0)
{;}

XMLObjectImplem* XMLObjectImplem::deepCopy() const
{
  XMLObjectImplem* rtn = new XMLObjectImplem(tag_);
  TEUCHOS_TEST_FOR_EXCEPTION(rtn==0, std::runtime_error, "XMLObjectImplem::deepCopy()");
  rtn->attributes_ = attributes_;
  rtn->content_ = content_;
	
  for (int i=0; i<children_.length(); i++)
  {
    rtn->addChild(children_[i].deepCopy());
  }

  return rtn;
}

int XMLObjectImplem::numChildren() const
{
  return children_.length();
}

void XMLObjectImplem::addAttribute(const std::string& name,
				   const std::string& value)
{
  attributes_[name] = value;
}

void XMLObjectImplem::addChild(const XMLObject& child)
{
  children_.append(child);
}

void XMLObjectImplem::addContent(const std::string& contentLine)
{
  content_.append(contentLine);
}

void XMLObjectImplem::removeContentLine(const size_t& i)
{
  Array<std::string>::iterator pos = content_.begin()+i;
  // does bound checking within content_.erase if BoundaryChecks are enabled
  content_.erase(pos);
}

const XMLObject& XMLObjectImplem::getChild(int i) const
{
  return children_[i];
}

std::string XMLObjectImplem::header(bool strictXML) const
{
  std::string rtn = "<" + tag_;
  for (Map::const_iterator i=attributes_.begin(); i!=attributes_.end(); ++i)
  {
    if (strictXML)
    {
      rtn += " "
	+ (*i).first
	+ "="
	+ XMLifyAttVal((*i).second);
    }
    else
    {
      rtn += " " + (*i).first + "=\"" + (*i).second + "\"";
    }
  }

  rtn += ">";
  return rtn;
}

std::string XMLObjectImplem::XMLifyAttVal(const std::string &attval) {
  std::string ret;
  bool hasQuot, hasApos;
  char delim;

  if (attval.find("\"") == std::string::npos)
  {
    hasQuot = false;
  }
  else
  {
    hasQuot = true;
  }

  if (attval.find("\'") == std::string::npos)
  {
    hasApos = false;
  }
  else
  {
    hasApos = true;
  }

  if (!hasQuot || hasApos)
  {
    delim = '\"'; // wrap the attribute value in "
  }
  else
  {
    delim = '\''; // wrap the attribute value in '
  }

  // Rules:
  // "-wrapped std::string cannot contain a literal "
  // '-wrapped std::string cannot contain a literal '
  // attribute value cannot contain a literal <
  // attribute value cannot contain a literal &
  ret.push_back(delim);
  for (std::string::const_iterator i=attval.begin(); i != attval.end(); i++)
  {
    if (*i == delim)
    {
      if (delim == '\'') ret.append("&apos;");
      else if (delim == '\"') ret.append("&quot;");
    }
    else if (*i == '&')
    {
      ret.append("&amp;");
    }
    else if (*i == '<')
    {
      ret.append("&lt;");
    }
    else
    {
      ret.push_back(*i);
    }
  }
  ret.push_back(delim);

  return ret;
}

std::string XMLObjectImplem::terminatedHeader(bool strictXML) const
{
  std::string rtn = "<" + tag_;
  for (Map::const_iterator i=attributes_.begin(); i!=attributes_.end(); ++i)
  {
    if (strictXML)
    {
      rtn += " "
	+ (*i).first
	+ "="
	+ XMLifyAttVal((*i).second);
    }
    else
    {
      rtn += " " + (*i).first + "=\"" + (*i).second + "\"";
    }
  }

  rtn += "/>";
  return rtn;
}

std::string XMLObjectImplem::toString() const
{
  std::string rtn;
  if (content_.length()==0 && children_.length()==0)
  {
    rtn = terminatedHeader(true) + "\n";
  }
  else
  {
    rtn = header() + "\n";
    bool allBlankContent = true;
    for (int i=0; i<content_.length(); i++)
    {
      if (!StrUtils::isWhite(content_[i]))
      {
	allBlankContent=false;
	break;
      }
    }
    if (!allBlankContent)
    {
      for (int i=0; i<content_.length(); i++)
      {
	rtn += content_[i];
      }
      rtn += "\n";
    }
    for (int i=0; i<children_.length(); i++)
    {
      rtn += children_[i].toString();
    }
    rtn += "</" + tag_ + ">\n";
  }
  return rtn;
}

void XMLObjectImplem::print(std::ostream& os, int indent) const
{
  for (int i=0; i<indent; i++) os << " ";
  if (content_.length()==0 && children_.length()==0)
  {
    os << terminatedHeader(true) << std::endl;
    return;
  }
  else
  {
    os << header(true) << std::endl;
    printContent(os, indent+2);

    for (int i=0; i<children_.length(); i++)
    {
      children_[i].print(os, indent+2);
    }
    for (int i=0; i<indent; i++) os << " ";
    os << "</" << tag_ << ">\n";
  }
}

void XMLObjectImplem::printContent(std::ostream& os, int indent) const
{
  std::string space = "";
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

  if (!allBlankContent)
  {

    for (int i=0; i<content_.length(); i++)
    {
      // remove leading spaces, we will indent
      std::string s(content_[i]);
      s.erase(size_t(0), s.find_first_not_of(" \r\t"));
      if((s.length()>0) && (!StrUtils::isWhite(s)))
        os << space << s << '\n';
    }
  }
}

