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
					TEST_FOR_EXCEPTION(true, runtime_error, 
                             "XMLObject::getRequiredBool value [" << val 
                             << "] should have been {TRUE|FALSE|YES|NO}");
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



