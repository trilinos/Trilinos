// @HEADER
// *****************************************************************************
//          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
//
// Copyright 2022 NTESS and the PyTrilinos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYTRILINOS2_TEUCHOS_ETI
#define PYTRILINOS2_TEUCHOS_ETI

#include <Teuchos_ParameterList.hpp>

namespace Teuchos {
    inline void initiate(ParameterList p) {};

#  define PARAMETERLIST_MF(T) \
	template T& ParameterList::get<T>(const std::string&); \
	template ParameterList& ParameterList::set<T>(std::string const&, T const&, std::string const&, RCP<const ParameterEntryValidator> const&);

PARAMETERLIST_MF(int)
PARAMETERLIST_MF(double)
PARAMETERLIST_MF(std::string)
PARAMETERLIST_MF(Teuchos::ParameterList)
}

#endif // PYTRILINOS2_TEUCHOS_ETI
