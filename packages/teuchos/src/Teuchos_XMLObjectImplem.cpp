// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;


XMLObjectImplem::XMLObjectImplem(const std::string& tag)
  : tag_(tag), attributes_(), children_(0), content_(0)
{;}

XMLObjectImplem* XMLObjectImplem::deepCopy() const 
{
  XMLObjectImplem* rtn = new XMLObjectImplem(tag_);
  TEST_FOR_EXCEPTION(rtn==0, std::runtime_error, "XMLObjectImplem::deepCopy()");
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
    rtn = terminatedHeader() + "\n";
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
    os << space;
    for (int i=0; i<content_.length(); i++)
    {
      os << content_[i];
    }
    os << '\n';
  }
}

