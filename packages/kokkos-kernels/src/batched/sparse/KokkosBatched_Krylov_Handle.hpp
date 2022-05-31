/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>

#ifndef __KOKKOSBATCHED_KRYLOV_HANDLE_HPP__
#define __KOKKOSBATCHED_KRYLOV_HANDLE_HPP__
//#define VERBOSE

namespace KokkosBatched {

/// \brief KrylovHandle
///
/// \tparam scalar_type: Scalar type of the linear solver

template <class scalar_type>
class KrylovHandle {
 public:
  using norm_type =
      typename Kokkos::Details::ArithTraits<scalar_type>::mag_type;

 private:
  norm_type tolerance;
  int max_iteration;

 public:
  KOKKOS_INLINE_FUNCTION
  KrylovHandle() {
    tolerance     = Kokkos::Details::ArithTraits<norm_type>::epsilon();
    max_iteration = 200;
  }

  /// \brief set_tolerance
  ///   Set the tolerance of the batched Krylov solver
  ///
  /// \param _tolerance [in]: New tolerance

  KOKKOS_INLINE_FUNCTION
  void set_tolerance(norm_type _tolerance) { tolerance = _tolerance; }

  /// \brief get_tolerance
  ///   Get the tolerance of the batched Krylov solver

  KOKKOS_INLINE_FUNCTION
  norm_type get_tolerance() const { return tolerance; }

  /// \brief set_max_iteration
  ///   Set the maximum number of iterations of the batched Krylov solver
  ///
  /// \param _max_iteration [in]: New maximum number of iterations

  KOKKOS_INLINE_FUNCTION
  void set_max_iteration(norm_type _max_iteration) {
    max_iteration = _max_iteration;
  }

  /// \brief get_max_iteration
  ///   Get the maximum number of iterations of the batched Krylov solver

  KOKKOS_INLINE_FUNCTION
  int get_max_iteration() const { return max_iteration; }
};

}  // namespace KokkosBatched

#endif
