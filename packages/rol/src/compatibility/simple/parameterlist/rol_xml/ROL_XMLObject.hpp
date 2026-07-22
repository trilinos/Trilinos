// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_XMLOBJECT_H
#define ROL_XMLOBJECT_H

/*! \file ROL_XMLObject.hpp
    \brief An object representation of a subset of XML data
*/

#include "ROL_XMLObjectImplem.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "ROL_toString.hpp"
using namespace Teuchos;

namespace ROL{

/** \brief Thrown when attempting to parse an empty XML std::string.*/
class EmptyXMLError : public std::runtime_error
{public: EmptyXMLError(const std::string& what_arg) : std::runtime_error(what_arg) {}};

class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLObject{
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
    {addAttribute(name, ROL::toString(val));}

  //! Add an int as an attribute
  void addInt(const std::string& name, int val)
    {addAttribute(name, ROL::toString(val));}

  //! Add a bool as an attribute
  void addBool(const std::string& name, bool val)
    {addAttribute(name, ROL::toString(val));}

  /** \brief Lookup whether or not Doubles are allowed.
   */
  template<class T>
  void addAttribute(const std::string& name, T value) {
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ptr_), Teuchos::EmptyXMLError,
		     "XMLObject::addAttribute: XMLObject is empty");
  ptr_->addAttribute(name, ROL::toString(value));
  }

			
  //! Add a child node to the node
  void addChild(const XMLObject& child);

  //! Add a line of character content
  void addContent(const std::string& contentLine);
  //@}

  void appendContentLine(const size_t& i, const std::string &str) {
    ptr_->appendContentLine(i, str);
  }

  void removeContentLine(const size_t& i) {
    ptr_->removeContentLine(i);
  }

protected:

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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT bool XMLObject::getRequired<bool>(const std::string& name) const;

template<>
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT int XMLObject::getRequired<int>(const std::string& name) const;

template<>
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT double XMLObject::getRequired<double>(const std::string& name) const;

template<>
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT std::string XMLObject::getRequired<std::string>(const std::string& name) const;

template<>
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void XMLObject::addAttribute<const std::string&>(const std::string& name, const std::string& value);


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


} // namespace ROL


#endif
