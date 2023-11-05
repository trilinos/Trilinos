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
#ifndef __TACHO_BLAS_SERIAL_HPP__
#define __TACHO_BLAS_SERIAL_HPP__

/// \file  Tacho_Blas_Serial.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace Tacho {

template <typename T> struct BlasSerial {

  // GEMM
  inline static void gemm(const char transa, const char transb, int m, int n, int k,
                          const T alpha, const T *A, int lda,
                                         const T *B, int ldb,
                          const T beta,  /* */ T *C, int ldc) {

    const T one(1), zero(0);

    if (alpha == zero) {
      if (beta == zero) {
        for (int j = 0; j < n; j++) {
          for (int i = 0; i < m; i++) C[i + j*ldc] = zero;
        }
      } else if (beta != one) {
        for (int j = 0; j < n; j++) {
          for (int i = 0; i < m; i++) C[i + j*ldc] *= beta;
        }
      }
      return;
    }
    if (m <= 0 || n <= 0 || k <= 0)
      return;

    if (alpha != zero) {
      if (transa == 'N' || transa == 'n') {
        if (transb == 'N' || transb == 'n') {
          for (int j = 0; j < n; j++) {
            if (beta == zero) {
              for (int i = 0; i < m; i++) C[i + j*ldc] = zero;
            } else if (beta != one) {
              for (int i = 0; i < m; i++) C[i + j*ldc] *= beta;
            }
            for (int l = 0; l < k; l++) {
              T val = alpha * B[l + j*ldb] ;
              for (int i = 0; i < m; i++) {
                C[i + j*ldc] += (A[i + l*lda] * val);
              }
            }
          }
        } else {
          Kokkos::abort("transb is not valid");
        }
      } else {
        Kokkos::abort("transa is not valid");
      }
    }
  }

  // HERK
  inline static void herk(const char uplo, const char trans, int n, int k,
                          const T alpha, const T *A, int lda,
                          const T beta,  /* */ T *C, int ldc) {

    typedef ArithTraits<T> arith_traits;
    const T zero(0);

    if (n <= 0 || k <= 0)
      return;

    if (uplo == 'U' || uplo == 'u') {
      if (trans == 'N' || trans == 'n') {
        Kokkos::abort("trans is not valid");
      } else if (trans == 'T' || trans == 't') {
        // C = alpha * (A^T * A) + beta * C
        for (int j = 0; j < n; j++) {
          for (int i = 0; i <= j; i++) {
            // initialize
            if (beta == zero)
              C[i + j*ldc] = zero;
            else
              C[i + j*ldc] *= beta;

            // update
            if (alpha != zero) {
              for (int l = 0; l < k; ++l) {
                C[i + j*ldc] += alpha * (A[l + i*lda] * A[l * j*lda]);
              }
            }
          }
        }
      } else if (trans == 'C' || trans == 'c') {
        // C = alpha * (A^T * A) + beta * C
        for (int j = 0; j < n; j++) {
          for (int i = 0; i <= j; i++) {
            // initialize
            if (beta == zero)
              C[i + j*ldc] = zero;
            else
              C[i + j*ldc] *= beta;

            // update
            if (alpha != zero) {
              for (int l = 0; l < k; ++l) {
                C[i + j*ldc] += alpha * arith_traits::conj(A[l + i*lda]) * A[l + j*lda];
              }
            }
          }
        }
      }
    } else {
      Kokkos::abort("uplo is not valid");
    }
    return;
  }

  // TRSM
  inline static void trsm(const char side, const char uplo, const char transa, const char diag, int m, int n,
                          const T alpha, const T *A, int lda,
                                         /* */ T *B, int ldb) {

    typedef ArithTraits<T> arith_traits;
    const T one(1.0), zero(0.0);
    
    if (alpha == zero) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          B[i + j*ldb] = zero;
        }
      }
      return;
    } else if (alpha != one) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          B[i + j*ldb] *= alpha;
        }
      }
    }

    ///
    /// side left
    ///
    if (side == 'L' || side == 'l') {
      if (uplo == 'U' || uplo == 'u') {
        if (transa == 'T' || transa == 't') {
          // Transpose-solve with Upper-triangular matrix from Left, U^{-T} * B
          for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
              for (int l = 0; l < i; l++) {
                B[i + j*ldb] -= A[l + i*lda] * B[l + j*ldb];
              }
              if (diag != 'U' && diag != 'u') {
                B[i + j*ldb] /= A[i + i*lda];
              }
            }
          }
        } else if (transa == 'C' || transa == 'c') {
          // Conjugate-solve with Upper-triangular matrix from Left, U^{-T} * B
          for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
              for (int l = 0; l < i; l++) {
                B[i + j*ldb] -= arith_traits::conj(A[l + i*lda]) * B[l + j*ldb];
              }
              if (diag != 'U' && diag != 'u') {
                B[i + j*ldb] /= A[i + i*lda];
              }
            }
          }
        } else {
          Kokkos::abort("transa is not valid");
        }
      } else if (uplo == 'L' || uplo == 'l') {
        Kokkos::abort("uplo is not valid");
      }
    }

    ///
    /// side right
    ///
    else if (side == 'R' || side == 'r') {
      if (uplo == 'U' || uplo == 'u') {
        if (transa == 'N' || transa == 'n') {
          // NonTranspose-solve with Upper-triangular matrix from Right, B * U^{-1}
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
              for (int l = 0; l < j; l++) {
                B[i + j*ldb] -= A[l + j*lda] * B[i + l*ldb];
              }
              if (diag != 'U' && diag != 'u') {
                B[i + j*ldb] /= A[j + j*lda];
              }
            }
          }
        } else {
          Kokkos::abort("transa is not valid");
        }
      } else if (uplo == 'L' || uplo == 'l') {
        Kokkos::abort("uplo is not valid");
      }
    }
  }
};

} // namespace Tacho

#endif
