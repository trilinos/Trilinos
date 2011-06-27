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


#ifndef TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_
#define TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_

/*! \file Teuchos_XMLDependencyExceptions.hpp
 * \brief A collection of Exceptions thrown
 * when converting Dependencys to and from
 * XML.
 */
#include <stdexcept>

namespace Teuchos {

/** \brief Thrown when no dependes of a dependency can't be found 
 * when converting the dependency to or from XML.
 */
class MissingDependeeException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependeeException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependeeException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a dependent of a dependency cant be found 
 * when converting the dependency to or from XML.
 */
class MissingDependentException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependentException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependentException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when no dependess of a dependency are specified
 * when converting the dependency from XML.
 */
class MissingDependeesException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependeesException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependeesException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when no dependents of a dependency are specified
 * when converting the dependency from XML.
 */
class MissingDependentsException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependentsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependentsException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a Dependency has too many dependees specified in
 * its XML
 */
class TooManyDependeesException : public std::logic_error{

public:

  /**
   * \brief Constructs an TooManyDependeesException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  TooManyDependeesException( 
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a StringVisualDependency is being converted
 * from XML and no Values tag is found
 */
class ValuesTagMissingException : public std::logic_error{

public:

  /**
   * \brief Constructs an ValuesTagMissingException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  ValuesTagMissingException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when  the rangesAndValidators tag for
 * the RangeValidatorDepencyConverter can't be found.
 */
class MissingRangesAndValidatorsTagException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingRangesAndValidatorsTagException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingRangesAndValidatorsTagException(
    const std::string& what_arg):std::logic_error(what_arg){}

};


/** \brief Thrown when converting a StrinvValidatorDependcny from XML
 * and no valuesAndValidators tag is found.
 */
class MissingValuesAndValidatorsTagException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingValuesAndValidatorsTagException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValuesAndValidatorsTagException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when no condtion tag is found when converting a 
 * ConditionVisualDependency from XML.
 */
class MissingConditionTagException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingConditionTagException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingConditionTagException(
    const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when converting a dependency that has validaotrs to 
 * and from XML. This excetpion indicates that a specified validator
 * could not be found*/
class MissingValidatorException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingValidatorException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValidatorException(
    const std::string& what_arg):std::logic_error(what_arg){}

};


/** \brief Thrown when an appropriate Dependency Converter can't be found. */
class CantFindDependencyConverterException : public std::logic_error{

public:

  /**
   * \brief Constructs an CantFindDependencyConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindDependencyConverterException(
    const std::string& what_arg):std::logic_error(what_arg){}

};



} // namespace Teuchos
#endif //TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_
