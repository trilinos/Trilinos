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
class TEUCHOS_LIB_DLL_EXPORT XMLObjectImplem
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

