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



XMLObject::XMLObject(const string& tag)
	: ptr_(rcp(new XMLObjectImplem(tag)))
{}

XMLObject::XMLObject(XMLObjectImplem* ptr)
	: ptr_(rcp(ptr))
{}

XMLObject XMLObject::deepCopy() const
{
	return XMLObject(ptr_->deepCopy());
}

const string& XMLObject::getRequired(const string& name) const 
{
	TEST_FOR_EXCEPTION(!hasAttribute(name), runtime_error,
                     "XMLObject::getRequired: key " 
                     << name << " not found");
  return getAttribute(name);
}

string XMLObject::getWithDefault(const string& name, 
																 const string& defaultValue) const
{
	if (hasAttribute(name)) return getRequired(name);
	else return defaultValue;
}

bool XMLObject::getRequiredBool(const string& name) const
{
	if (hasAttribute(name))
		{
			string val = StrUtils::allCaps(getRequired(name));
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
					TEST_FOR_EXCEPTION(true, runtime_error, 
                             "XMLObject::getRequiredBool value [" << val 
                             << "] should have been {TRUE|FALSE|YES|NO|0|1}");
				}
		}
	return false; // -Wall
}


void XMLObject::checkTag(const string& expected) const 
{
	TEST_FOR_EXCEPTION(getTag() != expected, std::runtime_error,
                     "XMLObject::checkTag error: expected <"
                     << expected << ">, found <" 
                     << getTag() << ">");
}



