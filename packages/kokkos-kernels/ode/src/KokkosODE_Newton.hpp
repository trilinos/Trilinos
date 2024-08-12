//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSODE_NEWTON_HPP
#define KOKKOSODE_NEWTON_HPP

/// \author Luc Berger-Vergiat (lberge@sandia.gov)
/// \file KokkosODE_Newton.hpp

#include "Kokkos_Core.hpp"

#include "KokkosODE_Types.hpp"
#include "KokkosODE_Newton_impl.hpp"

namespace KokkosODE {
namespace Experimental {

/// \brief Newton solver for non-linear system of equations
struct Newton {
  template <class system_type, class mat_type, class ini_vec_type, class rhs_vec_type, class update_type,
            class scale_type>
  KOKKOS_FUNCTION static newton_solver_status Solve(const system_type& sys, const Newton_params& params,
                                                    const mat_type& J, const mat_type& tmp, const ini_vec_type& y0,
                                                    const rhs_vec_type& rhs, const update_type& update,
                                                    const scale_type& scale) {
    return KokkosODE::Impl::NewtonSolve(sys, params, J, tmp, y0, rhs, update, scale);
  }
};

}  // namespace Experimental
}  // namespace KokkosODE

#endif  // KOKKOSODE_NEWTON_HPP
