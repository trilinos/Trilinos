// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XMLOBJECTIMPLEM_H
#define TEUCHOS_XMLOBJECTIMPLEM_H

/*! \file Teuchos_XMLObjectImplem.hpp
  \brief Low level implementation of XMLObject
*/

#include "Teuchos_map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos
{

class XMLObject;

/**
 * \brief The XMLObjectImplem class takes care of the low-level implementation
 details of XMLObject
*/
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLObjectImplem
{
  typedef Teuchos::map<std::string, std::string> Map;

public:
  //! Construct with a 'tag'
  XMLObjectImplem(const std::string& tag);

  //! Deep copy
  XMLObjectImplem* deepCopy() const ;

  //! Add a [name, value] attribute
  void addAttribute(const std::string& name, const std::string& value);

  //! Add a child XMLObject
  void addChild(const XMLObject& child);

  //! Add a content line
  void addContent(const std::string& contentLine);

  //! Return the tag std::string
  const std::string& getTag() const {return tag_;}

  //! Determine whether an attribute exists
  bool hasAttribute(const std::string& name) const
    {return attributes_.find(name) != attributes_.end();}

  //! Look up an attribute by name
  const std::string& getAttribute(const std::string& name) const
    {return (*(attributes_.find(name))).second;}

  //! Return the number of children
  int numChildren() const ;

  //! Look up a child by its index
  const XMLObject& getChild(int i) const ;

  //! Get the number of content lines
  int numContentLines() const {return content_.length();}

  //! Get all attributes
  const Map& getAttributes() const {return attributes_;}

  //! Look up a content line by index
  const std::string& getContentLine(int i) const {return content_[i];}

  //! Add string at the the end of a content line
  void appendContentLine(const size_t& i, const std::string &str) {
    content_[i].append(str);
  }

  //! Remove content line by index
  void removeContentLine(const size_t& i);

  //!  Print to stream with the given indentation level. Output will be well-formed XML.
  void print(std::ostream& os, int indent) const ;

  //! Write as a std::string. Output may be ill-formed XML.
  std::string toString() const ;

  //! Write the header
  std::string header(bool strictXML = false) const ;

  //! Write the header terminated as <Header/>
  std::string terminatedHeader(bool strictXML = false) const ;

  //! Write the footer
  std::string footer() const {return "</" + getTag() + ">";}

private:

  //! Print content lines using the given indentation level
  void printContent(std::ostream& os, int indent) const ;

  //! Convert attribute value text into well-formed XML
  static std::string XMLifyAttVal(const std::string &attval);

  std::string tag_;
  Map attributes_;
  Array<XMLObject> children_;
  Array<std::string> content_;

};

}

#endif

