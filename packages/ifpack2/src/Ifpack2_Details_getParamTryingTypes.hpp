/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DETAILS_GETPARAMTRYINGTYPES_HPP
#define IFPACK2_DETAILS_GETPARAMTRYINGTYPES_HPP

#include "Ifpack2_config.h"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Ifpack2 {
namespace Details {

template<class ... CandidateTypes>
struct GetParamTryingTypes {
  template<class ResultType>
  static void
  get (ResultType& result,
       const Teuchos::ParameterEntry& ent,
       const std::string& paramName,
       const char prefix[]);
};

template<>
struct GetParamTryingTypes<> {
  template<class ResultType>
  static void
  get (ResultType& /* result */,
       const Teuchos::ParameterEntry& /* ent */,
       const std::string& paramName,
       const char prefix[])
  {
    using Teuchos::TypeNameTraits;
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, prefix << "\"" << paramName
       << "\" parameter exists in input ParameterList, but does not "
       "have the right type.  The proper type is "
       << TypeNameTraits<ResultType>::name () << ".");
  }
};

template<class First, class ... Rest>
struct GetParamTryingTypes<First, Rest...> {
  template<class ResultType>
  static void
  get (ResultType& result,
       const Teuchos::ParameterEntry& ent,
       const std::string& paramName,
       const char prefix[])
  {
    if (ent.template isType<First> ()) {
      result = static_cast<ResultType> (Teuchos::getValue<First> (ent));
    }
    else {
      using rest_type = GetParamTryingTypes<Rest...>;
      rest_type::template get<ResultType> (result, ent, paramName, prefix);
    }
  }
};

template<class ResultType, class ... CandidateTypes>
void
getParamTryingTypes (ResultType& result,
                     const Teuchos::ParameterList& params,
                     const std::string& paramName,
                     const char prefix[])
{
  using Teuchos::ParameterEntry;
  const ParameterEntry* ent = params.getEntryPtr (paramName);
  if (ent != nullptr) {
    using impl_type = GetParamTryingTypes<CandidateTypes...>;
    impl_type::template get<ResultType> (result, *ent, paramName, prefix);
  }
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_GETPARAMTRYINGTYPES_HPP
