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
// Questions? Contact Luc Berger-Vergiat (lberge@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ARMPL_SPMV_HPP_
#define ARMPL_SPMV_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
#include <armpl.h>

void armpl_matvec_wrapper(armpl_spmat_t A, float* x, float* y) {
  const float alpha = 1.0;
  const float beta  = 0.0;
  armpl_spmv_exec_s(ARMPL_SPARSE_OPERATION_NOTRANS, alpha, A, x, beta, y);
}

void armpl_matvec_wrapper(armpl_spmat_t A, double* x, double* y) {
  const double alpha = 1.0;
  const double beta  = 0.0;
  armpl_spmv_exec_d(ARMPL_SPARSE_OPERATION_NOTRANS, alpha, A, x, beta, y);
}

template <typename AType, typename XType, typename YType>
void armpl_matvec(AType /*A*/, XType x, YType y, spmv_additional_data* data) {
  // using Scalar = typename AType::non_const_value_type;
  // Run armpl spmv corresponding to scalar type
  armpl_matvec_wrapper(data->A, x.data(), y.data());
}

#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
#endif  // ARMPL_SPMV_HPP_
