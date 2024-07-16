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

template <typename ... ViewArgs>
Perf
do_time_val(const size_t m, const size_t n, const size_t nloop,
            const bool check);

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad(const size_t m, const size_t n, const size_t p, const size_t nloop,
            const bool check);

template <typename FadType, typename ... ViewArgs>
Perf
do_time_scratch(const size_t m, const size_t n, const size_t p,
                const size_t nloop, const bool check);

template <typename ... ViewArgs>
Perf
do_time_analytic(const size_t m, const size_t n, const size_t p,
                 const size_t nloop, const bool check);

template <int MaxP, typename ... ViewArgs>
Perf
do_time_analytic_sl(const size_t m, const size_t n, const size_t p,
                    const size_t nloop, const bool check);

template <int p, typename ... ViewArgs>
Perf
do_time_analytic_s(const size_t m, const size_t n,
                   const size_t nloop, const bool check);
