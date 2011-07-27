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

#ifndef Teuchos_XMLOBJECT_H
#define Teuchos_XMLOBJECT_H

/*! \file Teuchos_XMLObject.hpp
    \brief An object representation of a subset of XML data
*/

#include "Teuchos_XMLObjectImplem.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos{

/** \brief Thrown when attempting to parse an empty XML std::string.*/
class EmptyXMLError : public std::runtime_error
{public: EmptyXMLError(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/** \ingroup XML 
 * \brief Representation of an XML data tree. XMLObject is a ref-counted
 * handle to a XMLObjectImplem object, allowing storage by reference.
 */
class TEUCHOS_LIB_DLL_EXPORT XMLObject{
public:

  //! @name Constructors 
  //@{

  //! Empty constructor
  XMLObject() : ptr_() {;}

  //! Construct using a node labeled by tag
  XMLObject(const std::string& tag);

  /** \brief Construct with a pointer to the low-level representation.
   *
   * This is used to allow construction of an XMLObject from the
   * XMLObjectImplem* return value of ExceptionBase::toXML().
   */
  XMLObject(XMLObjectImplem* ptr);
  //@}	

  //! @name Copy methods 
  //@{

  //! Make a deep copy of this object
  XMLObject deepCopy() const ;
  //@}

  //! @name Data Access methods 
  //@{

  //! Return the tag of the current node
  const std::string& getTag() const;

  //! Find out if the current node has an attribute of the specified name
  bool hasAttribute(const std::string& name) const;

  //! Return the value of the attribute with the specified name
  const std::string& getAttribute(const std::string& name) const;

  //! Get an attribute, throwing an std::exception if it is not found
  const std::string& getRequired(const std::string& name) const;

  //! Get a required attribute, returning it as a double
  double getRequiredDouble(const std::string& name) const 
    {return std::atof(getRequired(name).c_str());}
  
  //! Get a required attribute, returning it as an int
  int getRequiredInt(const std::string& name) const 
    {return std::atoi(getRequired(name).c_str());}

  //! Get a required attribute, returning it as T
  template<class T>
  T getRequired(const std::string& name) const{
    T toReturn;
    std::istringstream iss(getRequired(name));
    iss >> toReturn;
    return toReturn;
  }

  //! Get a required attribute, returning it as a bool
  bool getRequiredBool(const std::string& name) const ;

  /** \brief Get an attribute, assigning a default value if the requested
   * attribute does not exist */
  template<class T>
  T getWithDefault(const std::string& name, const T& defaultValue) const{
    if (hasAttribute(name)){
      return getRequired<T>(name);
    }
    else{
      return defaultValue;
    }
  }
  
  //! Return the number of child nodes owned by this node
  int numChildren() const;

  //! Return the i-th child node 
  const XMLObject& getChild(int i) const;

  /** \brief Returns the index of the first child found with the given tag name.
   * Returns -1 if no child is found.
   */
  int findFirstChild(std::string tagName) const;

  //! Return the number of lines of character content stored in this node 
  int numContentLines() const;

  //! Return the i-th line of character content stored in this node
  const std::string& getContentLine(int i) const;

  //! Represent this node and its children as a std::string
  std::string toString() const;

  //! Print this node and its children to stream with the given indentation
  void print(std::ostream& os, int indent) const;

  //! Write the header for this object to a std::string
  std::string header() const;

  //! Write the header for this object to a std::string
  std::string terminatedHeader() const;

  //! Write the footer for this object to a std::string
  std::string footer() const;

  //! Find out if a node is empty
  bool isEmpty() const { return ptr_.get()==0;}

  //! Check that a tag is equal to an expected std::string
  void checkTag(const std::string& expected) const ;
  //@}
	
  //! @name Tree-Assembly methods 
  //@{

  //! Add a double as an attribute
  void addDouble(const std::string& name, double val)
    {addAttribute(name, Teuchos::toString(val));}

  //! Add an int as an attribute
  void addInt(const std::string& name, int val)
    {addAttribute(name, Teuchos::toString(val));}

  //! Add a bool as an attribute
  void addBool(const std::string& name, bool val)
    {addAttribute(name, Teuchos::toString(val));}

  /** \brief Lookup whether or not Doubles are allowed.
   */
  template<class T>
  void addAttribute(const std::string& name, T value) {
  TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::addAttribute: XMLObject is empty");
  ptr_->addAttribute(name, Teuchos::toString(value));
  }

			
  //! Add a child node to the node
  void addChild(const XMLObject& child);

  //! Add a line of character content
  void addContent(const std::string& contentLine);
  //@}

private:

//use pragmas to disable some false-positive warnings for windows sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
  RCP<XMLObjectImplem> ptr_;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

};


template<>
TEUCHOS_LIB_DLL_EXPORT bool XMLObject::getRequired<bool>(const std::string& name) const;

template<>
TEUCHOS_LIB_DLL_EXPORT int XMLObject::getRequired<int>(const std::string& name) const;

template<>
TEUCHOS_LIB_DLL_EXPORT double XMLObject::getRequired<double>(const std::string& name) const;

template<>
TEUCHOS_LIB_DLL_EXPORT std::string XMLObject::getRequired<std::string>(const std::string& name) const;

template<>
TEUCHOS_LIB_DLL_EXPORT void XMLObject::addAttribute<const std::string&>(const std::string& name, const std::string& value);


/** \brief Write XMLObject to \c os stream.
 *
 * \relates XMLObject
 */
inline std::ostream& operator<<(std::ostream& os, const XMLObject& xml)
{
  xml.print(os, 0);
  return os;
}


/** \brief Write XMLObject to std::string.
 *
 * \relates XMLObject 
 */
inline std::string toString(const XMLObject& xml)
{
  return xml.toString();
}


} // namespace Teuchos


#endif
