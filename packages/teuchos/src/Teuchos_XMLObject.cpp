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


namespace Teuchos {


XMLObject::XMLObject(const std::string& tag)
  : ptr_(rcp(new XMLObjectImplem(tag)))
{}


XMLObject::XMLObject(XMLObjectImplem* ptr)
  : ptr_(rcp(ptr))
{}


XMLObject XMLObject::deepCopy() const
{
  if (is_null(ptr_))
  {
    return XMLObject();
  }
  return XMLObject(ptr_->deepCopy());
}


const std::string& XMLObject::getTag() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::getTag: XMLObject is empty");
  return ptr_->getTag();
}


bool XMLObject::hasAttribute(const std::string& name) const 
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::hasAttribute: XMLObject is empty");
  return ptr_->hasAttribute(name);
}


const std::string& XMLObject::getAttribute(const std::string& name) const 
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::getAttribute: XMLObject is empty");
  return ptr_->getAttribute(name);
}


const std::string& XMLObject::getRequired(const std::string& name) const 
{
  TEST_FOR_EXCEPTION(!hasAttribute(name), std::runtime_error,
                     "XMLObject::getRequired: key " 
                     << name << " not found");
  return getAttribute(name);
}


template<>
bool XMLObject::getRequired<bool>(const std::string& name) const
{
	return getRequiredBool(name);
}


template<>
int XMLObject::getRequired<int>(const std::string& name) const
{
	return getRequiredInt(name);
}


template<>
double XMLObject::getRequired<double>(const std::string& name) const
{
	return getRequiredDouble(name);
}


template<>
std::string XMLObject::getRequired<std::string>(const std::string& name) const
{
	return getRequired(name);
}


bool XMLObject::getRequiredBool(const std::string& name) const
{
  if (hasAttribute(name))
  {
    std::string val = StrUtils::allCaps(getRequired(name));
    if (val=="TRUE" || val=="YES" || val=="1")
    {
      return true;
    }
    else if (val=="FALSE" || val=="NO" || val=="0")
    {
      return false;
    }
    else
    {
      TEST_FOR_EXCEPTION(true, std::runtime_error, 
			 "XMLObject::getRequiredBool value [" << val 
			 << "] should have been {TRUE|FALSE|YES|NO|0|1}");
    }
  }
  return false; // -Wall
}


template<>
void XMLObject::addAttribute<const std::string&>(
  const std::string& name, const std::string& value) const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::addAttribute: XMLObject is empty");
  ptr_->addAttribute(name, value);
}


int XMLObject::numChildren() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::numChildren: XMLObject is empty");
  return ptr_->numChildren();
}


const XMLObject& XMLObject::getChild(int i) const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::getChild: XMLObject is empty");
  return ptr_->getChild(i);
}

int XMLObject::findFirstChild(std::string name) const{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::getChild: XMLObject is empty");
  for(int i = 0; i<numChildren(); ++i){
    if(getChild(i).getTag() == name){
      return i;
    }
  }
  return -1;
}

int XMLObject::numContentLines() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::numContentLines: XMLObject is empty");
  return ptr_->numContentLines();
}


const std::string& XMLObject::getContentLine(int i) const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::getContentLine: XMLObject is empty");
  return ptr_->getContentLine(i);
}


std::string XMLObject::toString() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::toString: XMLObject is empty");
  return ptr_->toString();
}


void XMLObject::print(std::ostream& os, int indent) const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::print: XMLObject is empty");
  ptr_->print(os, indent);
}


std::string XMLObject::header() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::header: XMLObject is empty");
  return ptr_->header();
}


std::string XMLObject::terminatedHeader() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::terminatedHeader: XMLObject is empty");
  return ptr_->terminatedHeader();
}


std::string XMLObject::footer() const
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::footer: XMLObject is empty");
  return ptr_->footer();
}


void XMLObject::checkTag(const std::string& expected) const 
{
  TEST_FOR_EXCEPTION(getTag() != expected, std::runtime_error,
                     "XMLObject::checkTag error: expected <"
                     << expected << ">, found <" 
                     << getTag() << ">");
}


void XMLObject::addChild(const XMLObject& child)
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::addChild: XMLObject is empty");
  ptr_->addChild(child);
}


void XMLObject::addContent(const std::string& contentLine)
{
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::addContent: XMLObject is empty");
  ptr_->addContent(contentLine);
}


} // namespace Teuchos
