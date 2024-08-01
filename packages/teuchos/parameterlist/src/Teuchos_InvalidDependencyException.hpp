// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_INVALIDDEPENDENCYEXCEPTION_HPP_
#define TEUCHOS_INVALIDDEPENDENCYEXCEPTION_HPP_
#include <stdexcept>

namespace Teuchos {

/**
 * Thrown when some aspect of a Dependency has been determined to be invalid.
 */
class InvalidDependencyException : public std::logic_error{
public:
  /**
   * Constructs an InvalidDependencyException
   *
   * @param what_arg The error message to be associated with this error.
   */
  InvalidDependencyException(const std::string& what_arg):std::logic_error(what_arg){}
};

}
#endif //TEUCHOS_INVALIDDEPENDENCYEXCEPTION_HPP_
