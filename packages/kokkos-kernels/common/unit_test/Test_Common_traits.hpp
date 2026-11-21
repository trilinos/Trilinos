// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file Test_Common_traits.hpp
/// \brief Tests that traits are correctly defined

#ifndef TEST_COMMON_MODE_TRAITS_HPP
#define TEST_COMMON_MODE_TRAITS_HPP

#include <Kokkos_Core.hpp>
#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>

// Resource Mode
static_assert(KokkosBlas::is_mode_v<KokkosBlas::Mode::Serial>);
static_assert(KokkosBlas::is_mode_v<KokkosBlas::Mode::Team>);
static_assert(KokkosBlas::is_mode_v<KokkosBlas::Mode::TeamVector>);

static_assert(!KokkosBlas::is_mode_v<KokkosBatched::Diag::Unit>);
static_assert(!KokkosBlas::is_mode_v<KokkosBatched::Uplo::Upper>);

// Transpose Modes
static_assert(KokkosBlas::is_trans_v<KokkosBlas::Trans::NoTranspose>);
static_assert(KokkosBlas::is_trans_v<KokkosBlas::Trans::Transpose>);
static_assert(KokkosBlas::is_trans_v<KokkosBlas::Trans::ConjNoTranspose>);
static_assert(KokkosBlas::is_trans_v<KokkosBlas::Trans::ConjTranspose>);

static_assert(!KokkosBlas::is_trans_v<KokkosBlas::Mode::Team>);
static_assert(!KokkosBlas::is_trans_v<KokkosBatched::Side::Left>);

// Side Modes
static_assert(KokkosBatched::is_side_v<KokkosBatched::Side::Left>);
static_assert(KokkosBatched::is_side_v<KokkosBatched::Side::Right>);

static_assert(!KokkosBatched::is_side_v<double>);
static_assert(!KokkosBatched::is_side_v<KokkosBatched::Uplo::Lower>);

// Uplo Modes
static_assert(KokkosBatched::is_uplo_v<KokkosBatched::Uplo::Upper>);
static_assert(KokkosBatched::is_uplo_v<KokkosBatched::Uplo::Lower>);

static_assert(!KokkosBatched::is_uplo_v<KokkosBatched::Diag::NonUnit>);
static_assert(!KokkosBatched::is_uplo_v<KokkosBlas::Trans::Transpose>);

// Diag Modes
static_assert(KokkosBatched::is_diag_v<KokkosBatched::Diag::Unit>);
static_assert(KokkosBatched::is_diag_v<KokkosBatched::Diag::NonUnit>);

static_assert(!KokkosBatched::is_diag_v<KokkosBlas::Trans::ConjTranspose>);
static_assert(!KokkosBatched::is_diag_v<KokkosBatched::Side::Right>);

#endif  // TEST_COMMON_MODE_TRAITS_HPP
