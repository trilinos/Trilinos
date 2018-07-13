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


