// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_INVALIDCONDITIONEXCEPTION_HPP_
#define TEUCHOS_INVALIDCONDITIONEXCEPTION_HPP_
#include <stdexcept>

namespace Teuchos {

/**
 * Thrown when some aspect of a Condition has been determined to be invalid.
 */
class InvalidConditionException : public std::logic_error{
public:
  /**
   * Constructs an InvalidConditionException
   *
   * @param what_arg The error message to be associated with this error.
   */
  InvalidConditionException(const std::string& what_arg):std::logic_error(what_arg){}
};

}
#endif //TEUCHOS_INVALIDCONDITIONEXCEPTION_HPP_
