//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
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
#ifndef __KOKKOSBATCHED_GMRES_HPP__
#define __KOKKOSBATCHED_GMRES_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

/// \brief Batched GMRES: Selective Interface
///
/// \tparam OperatorType: The type of the operator of the system
/// \tparam VectorViewType: Input type for the right-hand side and the solution,
/// needs to be a 2D view
///
/// \param member [in]: TeamPolicy member
/// \param A [in]: batched operator (can be a batched matrix or a (left or right
/// or both) preconditioned batched matrix) \param B [in]: right-hand side, a
/// rank 2 view \param X [in/out]: initial guess and solution, a rank 2 view
/// \param handle [in]: a handle which provides different information such as
/// the tolerance or the maximal number of iterations of the solver.

#include "KokkosBatched_Krylov_Handle.hpp"
#include "KokkosBatched_GMRES_Team_Impl.hpp"
#include "KokkosBatched_GMRES_TeamVector_Impl.hpp"

namespace KokkosBatched {

template <typename MemberType, typename ArgMode>
struct GMRES {
  template <typename OperatorType, typename VectorViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const OperatorType &A, const VectorViewType &B,
      const VectorViewType &X,
      const KrylovHandle<typename VectorViewType::non_const_value_type>
          &handle) {
    int status = 0;
    if (std::is_same<ArgMode, Mode::Team>::value) {
      status =
          TeamGMRES<MemberType>::template invoke<OperatorType, VectorViewType>(
              member, A, B, X, handle);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      status = TeamVectorGMRES<MemberType>::template invoke<OperatorType,
                                                            VectorViewType>(
          member, A, B, X, handle);
    }
    return status;
  }
};

}  // namespace KokkosBatched
#endif
