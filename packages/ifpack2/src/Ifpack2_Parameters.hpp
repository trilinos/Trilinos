// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_PARAMETERS_HPP
#define IFPACK2_PARAMETERS_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Ifpack2 {

//! Fills a list which contains all the parameters possibly used by Ifpack2.
void getValidParameters(Teuchos::ParameterList& params);

//! Set a value from a ParameterList if a parameter with the specified name exists.
/** If the specified name does not name a parameter in the list, then 'value' is
  not referenced.
*/
template<typename T>
void getParameter(const Teuchos::ParameterList& params, const std::string& name, T& value)
{
  if (params.isParameter(name)) {
    if (params.isType<T>(name)) {
      value = params.get<T>(name);
    }
  }
}

}//namespace Ifpack2

#endif
