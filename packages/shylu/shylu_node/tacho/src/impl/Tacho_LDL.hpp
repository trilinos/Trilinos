// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_LDL_HPP__
#define __TACHO_LDL_HPP__

/// \file Tacho_LDL.hpp
/// \brief Front interface for LDL^t dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// LDL:
///
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgAlgo> struct LDL;

struct LDL_Algorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
