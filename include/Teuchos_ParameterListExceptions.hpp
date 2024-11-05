// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETER_LIST_EXCEPTIONS_H
#define TEUCHOS_PARAMETER_LIST_EXCEPTIONS_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

namespace Exceptions {

/** \brief .
 * \relates ParameterList
 */
class InvalidArgument : public std::invalid_argument
{public: InvalidArgument(const std::string& what_arg) : std::invalid_argument(what_arg) {}};

/** \brief .
 * \relates ParameterList
 */
class InvalidParameter : public std::logic_error
{public: InvalidParameter(const std::string& what_arg) : std::logic_error(what_arg) {}};

/** \brief .
 * \relates ParameterList
 */
class InvalidParameterName : public InvalidParameter
{public: InvalidParameterName(const std::string& what_arg) : InvalidParameter(what_arg) {}};

/** \brief .
 * \relates ParameterList
 */
class InvalidParameterType : public InvalidParameter
{public: InvalidParameterType(const std::string& what_arg) : InvalidParameter(what_arg) {}};

/** \brief .
 * \relates ParameterList
 */
class InvalidParameterValue : public InvalidParameter
{public: InvalidParameterValue(const std::string& what_arg) : InvalidParameter(what_arg) {}};

} // namespace Exceptions

} // end of Teuchos namespace

#endif // TEUCHOS_PARAMETER_LIST_EXCEPTIONS_H


