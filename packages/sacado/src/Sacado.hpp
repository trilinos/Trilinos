// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_HPP
#define SACADO_HPP

// Top-level Sacado include file for Sacado classes that work with Kokkos.
// Users should use this file instead of Sacado_No_Kokkos.hpp when working
// with Kokkos.

// Ensure "Sacado.hpp" and "Sacado_No_Kokkos.hpp" are not both included
#ifdef SACADO_NO_KOKKOS_HPP
#error "Do not include Sacado.hpp and Sacado_No_Kokkos.hpp in the same file."
#endif

// Version string
#include "Sacado_Version.hpp"

// Declarations of all overloaded math functions
#include "Sacado_MathFunctions.hpp"

// Traits for all of the Sacado classes -- Include these first so they are all
// defined before any nesting of AD classes
#ifdef SACADO_ENABLE_NEW_DESIGN
#include "Sacado_Fad_Exp_ExpressionTraits.hpp"
#include "Sacado_Fad_Exp_GeneralFadTraits.hpp"
#endif
#include "Sacado_Fad_ExpressionTraits.hpp"
#include "Sacado_Fad_DFadTraits.hpp"
#include "Sacado_Fad_SFadTraits.hpp"
#include "Sacado_Fad_SLFadTraits.hpp"
#include "Sacado_ELRFad_ExpressionTraits.hpp"
#include "Sacado_ELRFad_DFadTraits.hpp"
#include "Sacado_ELRFad_SFadTraits.hpp"
#include "Sacado_ELRFad_SLFadTraits.hpp"
#include "Sacado_CacheFad_ExpressionTraits.hpp"
#include "Sacado_CacheFad_DFadTraits.hpp"
#include "Sacado_CacheFad_SFadTraits.hpp"
#include "Sacado_CacheFad_SLFadTraits.hpp"
#include "Sacado_ELRCacheFad_ExpressionTraits.hpp"
#include "Sacado_ELRCacheFad_DFadTraits.hpp"
#include "Sacado_ELRCacheFad_SFadTraits.hpp"
#include "Sacado_ELRCacheFad_SLFadTraits.hpp"
#if !defined(__CUDACC__) && !defined(__HIPCC__ )
#include "Sacado_Fad_DVFadTraits.hpp"
#include "Sacado_LFad_LogicalSparseTraits.hpp"
#include "Sacado_ScalarFlopCounterTraits.hpp"
#include "Sacado_Tay_TaylorTraits.hpp"
#include "Sacado_trad_Traits.hpp"
#include "Sacado_trad2_Traits.hpp"
#include "Sacado_tradvec_Traits.hpp"
#endif

// Standard forward AD classes
#ifdef SACADO_ENABLE_NEW_DESIGN
#include "Sacado_Fad_Exp_DFad.hpp"
#include "Sacado_Fad_Exp_SFad.hpp"
#include "Sacado_Fad_Exp_SLFad.hpp"
#include "Sacado_Fad_Exp_ViewFad.hpp"
#include "Sacado_Fad_Exp_Atomic.hpp"
#endif
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_Fad_SFad.hpp"
#include "Sacado_Fad_SLFad.hpp"

// Expression-level-reverse forward AD classes
#include "Sacado_ELRFad_DFad.hpp"
#include "Sacado_ELRFad_SFad.hpp"
#include "Sacado_ELRFad_SLFad.hpp"

// Caching forward AD classes
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_CacheFad_SFad.hpp"
#include "Sacado_CacheFad_SLFad.hpp"

// Caching expression-level reverse mode forward AD classes
#include "Sacado_ELRCacheFad_DFad.hpp"
#include "Sacado_ELRCacheFad_SFad.hpp"
#include "Sacado_ELRCacheFad_SLFad.hpp"

// Kokkos::View specialization for Sacado AD classes
#include "Kokkos_View_Fad.hpp"

//
// AD classes that won't or don't build under Cuda
//
#ifndef __CUDACC__

// Standard forward AD classes
#include "Sacado_Fad_DVFad.hpp"
#include "Sacado_Fad_Vector.hpp"
#include "Sacado_LFad_LogicalSparse.hpp"

// Reverse AD classes
#include "Sacado_trad.hpp"
#include "Sacado_trad2.hpp"
#include "Sacado_tradvec.hpp"

// Taylor polynomial AD classes
#include "Sacado_Tay_Taylor.hpp"

// Flop-counting classes
#include "Sacado_ScalarFlopCounter.hpp"

#endif

#endif // SACADO_KOKKOS_HPP
