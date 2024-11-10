// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#pragma once

//#define SACADO_DISABLE_FAD_VIEW_SPEC

#include "Kokkos_Timer.hpp"

struct Perf {
  double time;
  double flops;
  double throughput;
};

const int SFadSize  = 8;
const int SLFadSize = SFadSize;
#if defined (KOKKOS_ENABLE_HIP)
const int HierSFadSize  = 64;
#else
const int HierSFadSize  = 32;
#endif
const int HierSLFadSize = HierSFadSize;
