// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#pragma once

#include "common.hpp"

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad_hierarchical_dfad(const size_t m, const size_t n, const size_t p,
                              const size_t nloop, const bool check);

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad_hierarchical_dfad_scratch(
  const size_t m, const size_t n, const size_t p, const size_t nloop,
  const bool check);
