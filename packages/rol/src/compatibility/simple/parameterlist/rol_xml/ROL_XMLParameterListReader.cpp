// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_XMLParameterListReader.hpp"
#include "ROL_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_Assert.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace ROL {

/**
 * \brief Thrown when an element inside a parameter list is bad.
 */
class BadParameterListElementException : public std::logic_error{
public:
  /**
   * \brief Constructs a BadParameterListElementException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadParameterListElementException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when the root xml tag for a parameter list is incorrect.
 */
class BadXMLParameterListRootElementException : public std::logic_error{
public:
  /**
   * \brief Constructs a BadXMLParameterListRootElementException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadXMLParameterListRootElementException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when a parameter entry tag is missing it's name attribute.
 */
class NoNameAttributeException : public std::logic_error{
public:
  /**
   * \brief Constructs a NoNameAttributeException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoNameAttributeException(const std::string& what_arg):std::logic_error(what_arg){}
};


XMLParameterListReader::XMLParameterListReader()
: _allowDuplicateSublists(true)
{;}

bool XMLParameterListReader::getAllowsDuplicateSublists() const
{ return _allowDuplicateSublists; }

void XMLParameterListReader::setAllowsDuplicateSublists(bool policy)
{ _allowDuplicateSublists = policy; }


ParameterList
XMLParameterListReader::toParameterList(const XMLObject& xml) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    xml.getTag()
    !=
    getParameterListTagName(),
    BadXMLParameterListRootElementException,
    "XMLParameterListReader expected tag " <<
    getParameterListTagName()
    <<", found " << xml.getTag());
  RCP<ParameterList> rtn = rcp(new ParameterList);
  convertParameterList(xml, rtn);
  ParameterList toReturn = ParameterList(*rtn);
  return toReturn;
}


void
XMLParameterListReader::convertParameterList(const XMLObject& xml,
  RCP<ParameterList> parentList) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    xml.getTag() != getParameterListTagName(),
    BadParameterListElementException,
    "XMLParameterListReader expected tag " <<
    getParameterListTagName()
    <<", found the tag "
    << xml.getTag());

  if(xml.hasAttribute(getNameAttributeName())){
    parentList->setName(
      xml.getRequired(getNameAttributeName()));
  }

  for (int i=0; i<xml.numChildren(); i++) {

      XMLObject child = xml.getChild(i);

      TEUCHOS_TEST_FOR_EXCEPTION(
        child.getTag() != getParameterListTagName()
        &&
        child.getTag() != ParameterEntry::getTagName(),
        BadParameterListElementException,
        "XMLParameterListReader expected tag "
        << getParameterListTagName() << " or "
        << ParameterEntry::getTagName() << ", but found "
        << child.getTag() << " tag.");


      if(
        child.getTag() == getParameterListTagName()
        ||
        child.getTag() == ParameterEntry::getTagName()
        )
      {

        std::string name;
        if (child.getTag()==getParameterListTagName()) {
          if ( child.hasAttribute(getNameAttributeName()) ) {
            name = child.getRequired(getNameAttributeName());
          }
          else {
            // the name needs to be unique: generate one
            std::ostringstream ss;
            ss << "child" << i;
            name = ss.str();
          }
          TEUCHOS_TEST_FOR_EXCEPTION(
            _allowDuplicateSublists == false
            &&
            parentList->isSublist(name) == true,
            DuplicateParameterSublist,
            "XMLParameterListReader encountered duplicate sublist \"" << name << "\", in violation"
            << " of the policy specified by XMLParameterListReader::setAllowsDuplicateSublists()." );
          RCP<ParameterList> newList = sublist(parentList, name);
          convertParameterList(child, newList);
        }
        else if (child.getTag() == ParameterEntry::getTagName()) {
          TEUCHOS_TEST_FOR_EXCEPTION(
              !child.hasAttribute(getNameAttributeName()),
              NoNameAttributeException,
              "All child nodes of a ParameterList must have a name attribute!" <<
              std::endl << std::endl);
          name = child.getRequired(getNameAttributeName());
          parentList->setEntry(
            name, ParameterEntryXMLConverterDB::convertXML(child));
      }
    }
  }
}

} // namespace ROL

