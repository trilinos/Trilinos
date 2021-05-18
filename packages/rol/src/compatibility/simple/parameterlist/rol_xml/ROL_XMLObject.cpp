#include "ROL_XMLObject.hpp"
#include "ROL_StrUtils.hpp"


namespace ROL {


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
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::getTag: XMLObject is empty");
  return ptr_->getTag();
}


bool XMLObject::hasAttribute(const std::string& name) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::hasAttribute: XMLObject is empty");
  return ptr_->hasAttribute(name);
}


const std::string& XMLObject::getAttribute(const std::string& name) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::getAttribute: XMLObject is empty");
  return ptr_->getAttribute(name);
}


const std::string& XMLObject::getRequired(const std::string& name) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasAttribute(name), std::runtime_error,
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
  TEUCHOS_TEST_FOR_EXCEPTION(!hasAttribute(name), std::runtime_error,
                     "XMLObject::getRequired: key "
                     << name << " not found");
  std::string val = StrUtils::allCaps(getRequired(name));

  TEUCHOS_TEST_FOR_EXCEPTION( val!="TRUE" && val!="YES" && val!="1"
    && val!="FALSE" && val!="NO" && val!="0",
    std::runtime_error,
		"XMLObject::getRequiredBool value [" << val
		<< "] should have been {TRUE|FALSE|YES|NO|0|1}");

  if (val=="TRUE" || val=="YES" || val=="1")
  {
    return true;
  }
  else
  {
    return false;
  }
}


template<>
void XMLObject::addAttribute<const std::string&>(
  const std::string& name, const std::string& value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::addAttribute: XMLObject is empty");
  ptr_->addAttribute(name, value);
}


int XMLObject::numChildren() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::numChildren: XMLObject is empty");
  return ptr_->numChildren();
}


const XMLObject& XMLObject::getChild(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::getChild: XMLObject is empty");
  return ptr_->getChild(i);
}

int XMLObject::findFirstChild(std::string name) const{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
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
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::numContentLines: XMLObject is empty");
  return ptr_->numContentLines();
}


const std::string& XMLObject::getContentLine(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::getContentLine: XMLObject is empty");
  return ptr_->getContentLine(i);
}


std::string XMLObject::toString() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::toString: XMLObject is empty");
  return ptr_->toString();
}


void XMLObject::print(std::ostream& os, int indent) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::print: XMLObject is empty");
  ptr_->print(os, indent);
}


std::string XMLObject::header() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::header: XMLObject is empty");
  return ptr_->header();
}


std::string XMLObject::terminatedHeader() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::terminatedHeader: XMLObject is empty");
  return ptr_->terminatedHeader();
}


std::string XMLObject::footer() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::footer: XMLObject is empty");
  return ptr_->footer();
}


void XMLObject::checkTag(const std::string& expected) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(getTag() != expected, std::runtime_error,
                     "XMLObject::checkTag error: expected <"
                     << expected << ">, found <"
                     << getTag() << ">");
}


void XMLObject::addChild(const XMLObject& child)
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::addChild: XMLObject is empty");
  ptr_->addChild(child);
}


void XMLObject::addContent(const std::string& contentLine)
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), ROL::EmptyXMLError,
		     "XMLObject::addContent: XMLObject is empty");
  ptr_->addContent(contentLine);
}


} // namespace ROL
