// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
