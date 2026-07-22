// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XMLFUNCTIONOBJECTEXCEPTIONS_HPP_
#define TEUCHOS_XMLFUNCTIONOBJECTEXCEPTIONS_HPP_

/*! \file Teuchos_XMLFunctionObjectExceptions.hpp
 * \brief A collection of Exceptions thrown
 * when converting FunctionObjects to and from
 * XML.
 */
#include <stdexcept>

namespace Teuchos {

/** \brief Thrown when an appropriate FunctionObject Converter can't be found */
class CantFindFunctionObjectConverterException : public std::logic_error{

public:

  /**
   * \brief Constructs an CantFindFunctionObjectConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindFunctionObjectConverterException(const std::string& what_arg):
    std::logic_error(what_arg){}

};



} // namespace Teuchos
#endif //TEUCHOS_XMLFUNCTIONOBJECTEXCEPTIONS_HPP_

