// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "AEvil_decl.hpp"

#ifdef DO_EXPLICIT_INSTANTIATION

#include "AEvil_def.hpp"

namespace EvilPack {

template class AEvil<double>;
template class AEvil<int>;

} // namespace EvilPack

#endif // DO_EXPLICIT_INSTANTIATION
