// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_
#define TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_

/*! \file Teuchos_XMLConditionExceptions.hpp
 * \brief A collection of Exceptions thrown
 * when converting Conditions to and from
 * XML.
 */
#include <stdexcept>

namespace Teuchos {

/** \brief Thrown when a StringConditon is missing it's Value tag.
 */
class MissingValuesTagException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingValuesTagException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValuesTagException(const std::string& what_arg):
    std::logic_error(what_arg){}

};

/** \brief Thrown when an appropriate Condition Converter can't be found */
class CantFindConditionConverterException : public std::logic_error{

public:

  /**
   * \brief Constructs an CantFindConditionConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindConditionConverterException(const std::string& what_arg):
    std::logic_error(what_arg){}

};



} // namespace Teuchos
#endif //TEUCHOS_XMLCONDITIONEXCEPTIONS_HPP_

