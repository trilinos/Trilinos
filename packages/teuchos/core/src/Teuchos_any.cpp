// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_any.hpp"

template std::string& Teuchos::any_cast<std::string>(Teuchos::any&);
template int& Teuchos::any_cast<int>(Teuchos::any&);
template long long& Teuchos::any_cast<long long>(Teuchos::any&);
template bool& Teuchos::any_cast<bool>(Teuchos::any&);
template double& Teuchos::any_cast<double>(Teuchos::any&);
