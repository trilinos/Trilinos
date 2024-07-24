// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "EvilBase_decl.hpp"

#ifdef DO_EXPLICIT_INSTANTIATION

#include "EvilBase_def.hpp"

namespace EvilPack {

template class EvilBase<double>;
template class EvilBase<int>;

} // namespace EvilPack

#endif // DO_EXPLICIT_INSTANTIATION
