// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_SLFAD_HPP
#define SACADO_CACHEFAD_SLFAD_HPP

#include "Sacado_CacheFad_GeneralFadExpr.hpp"
#include "Sacado_CacheFad_SLFadTraits.hpp"
#include "Sacado_Fad_StaticStorage.hpp"

#define FAD_NS CacheFad
#include "Sacado_Fad_SLFad_tmpl.hpp"
#undef FAD_NS

#include "Sacado_CacheFad_ViewFad.hpp"

#endif // SACADO_CACHEFAD_SLFAD_HPP
