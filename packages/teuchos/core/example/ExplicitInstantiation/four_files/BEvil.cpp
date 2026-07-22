// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BEvil_decl.hpp"

#ifdef DO_EXPLICIT_INSTANTIATION

#include "BEvil_def.hpp"

namespace EvilPack {

template class BEvil<double>;
template class BEvil<int>;

} // namespace EvilPack

#endif // DO_EXPLICIT_INSTANTIATION
