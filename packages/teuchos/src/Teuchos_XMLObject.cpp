#include "Teuchos_XMLObject.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;



XMLObject::XMLObject(const string& tag)
	: ptr_(new XMLObjectImplem(tag))
{
	if (ptr_.get()==0) Error::raise("XMLObject ctor");
}

XMLObject::XMLObject(XMLObjectImplem* ptr)
	: ptr_(ptr)
{
	if (ptr_.get()==0) Error::raise("XMLObject ctor");
}

XMLObject XMLObject::deepCopy() const
{
	return XMLObject(ptr_->deepCopy());
}

const string& XMLObject::getRequired(const string& name) const 
{
	if (hasAttribute(name))
		{
			return getAttribute(name);
		}
	else
		{
			Error::raise("XMLObject::getRequired: key " 
													 + name + " not found");
			return getAttribute("-Wall"); // -Wall
		}
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
			if (val=="TRUE" || val=="YES")
				{
					return true;
				}
			else if (val=="FALSE" || val=="NO")
				{
					return false;
				}
			else
				{
					Error::raise("XMLObject::getRequiredBool value [" + val 
											 + "] should have been {TRUE|FALSE|YES|NO}");
				}
		}
	return false; // -Wall
}


void XMLObject::checkTag(const string& expected) const 
{
	if (getTag() != expected)
		{
			Error::raise("XMLObject::checkTag error: expected <"
													 + expected + ">, found <" + getTag() + ">");
		}
}



