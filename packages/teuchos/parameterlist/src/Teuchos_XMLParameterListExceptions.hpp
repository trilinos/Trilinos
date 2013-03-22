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


#ifndef TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_
#define TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_

/*! \file Teuchos_XMLParameterListExceptions.hpp
 * \brief A collection of Exceptions that can be potentially
 * thrown when converting a ParameterList to and from XML
 */
#include <stdexcept>

namespace Teuchos {

/**
 * \brief Thrown when an appropriate ParameterEntryXMLConverter
 * can't be found.
 */
class CantFindParameterEntryConverterException : public std::logic_error{

public: 

  /**
   * \brief Constructs an CantFindParameterEntryConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindParameterEntryConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Constructs a CantFindParameterEntryConverterException */

/**
 * \brief Thrown when two validators in an XML file have the same ID.
 */
class DuplicateValidatorIDsException : public std::logic_error{

public: 

  /**
   * \brief Constructs an DuplicateValidatorIDsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  DuplicateValidatorIDsException(const std::string& what_arg):std::logic_error(what_arg){}

};

/**
 * \brief Thrown when two parameters in an XML file have the same ID.
 */
class DuplicateParameterIDsException : public std::logic_error{

public: 

  /**
   * \brief Constructs an DuplicateParameterIDsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  DuplicateParameterIDsException(const std::string& what_arg):std::logic_error(what_arg){}

};

/**
 * \brief Thrown when a bad validator xml converter is used.
 */
class BadValidatorXMLConverterException : public std::logic_error{

public: 
  /**
   * \brief Constructs an BadValidatorXMLConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadValidatorXMLConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};


/**
 * \brief Thrown when the ValidatorXMLConverterDB can't find an 
 * appropriate converter.
 */
class CantFindValidatorConverterException : public std::logic_error{

public: 
  /**
   * \brief Constructs a CantFindValidatorConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindValidatorConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};


/**
 * \brief Thrown when a converter is being used to convert either and XML tag or
 * ParameterEntry with an innappropriate type.
 */
class BadParameterEntryXMLConverterTypeException : public std::logic_error{

public: 

  /**
   * \brief Constructs a BadParmaeterEntryXMLConverterTypeException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadParameterEntryXMLConverterTypeException(const std::string& what_arg):std::logic_error(what_arg){}

};



/**
 * \brief Thrown when a parameter entry tag is missing it's value attribute.
 */
class NoValueAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoValueAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoValueAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};


/**
 * \brief Thrown when a parameter entry tag is missing it's type attribute.
 */
class NoTypeAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoTypeAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoTypeAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when a parameter entry tag is missing it's name attribute.
 */
class NoNameAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoNameAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoNameAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};

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
 * \brief Thrown when a referenced validator can't be found.
 */
class MissingValidatorDefinitionException : public std::logic_error{
public: 
  /**
   * \brief Constructs a MissingValidatorDefinitionException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValidatorDefinitionException(
    const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when a referenced ParameterEntry can't be found.
 */
class MissingParameterEntryDefinitionException : public std::logic_error{
public: 
  /**
   * \brief Constructs a MissingParameterEntryDefinitionException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingParameterEntryDefinitionException(
    const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when xml tag is encountered that is either unrecognized or 
 * inappropriate for a given context.
 */
class BadTagException : public std::logic_error{
public: 
  /**
   * \brief Constructs a MissingValidatorDefinitionException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadTagException(const std::string& what_arg):std::logic_error(what_arg){}
};


} // namespace Teuchos
#endif //TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_
