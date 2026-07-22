// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_TEUCHOS_ETI
#define PYROL_TEUCHOS_ETI

#include <Teuchos_ParameterList.hpp>

namespace Teuchos {
    inline void initiate(ParameterList p) {}

#  define PARAMETERLIST_MF(T) \
	template T& ParameterList::get<T>(const std::string&); \
	template ParameterList& ParameterList::set<T>(std::string const&, T const&, std::string const&, RCP<const ParameterEntryValidator> const&);

PARAMETERLIST_MF(bool)
PARAMETERLIST_MF(int)
PARAMETERLIST_MF(double)
PARAMETERLIST_MF(std::string)
PARAMETERLIST_MF(Teuchos::ParameterList)
}

#endif // PYROL_TEUCHOS_ETI