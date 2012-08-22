//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __Kokkos_Mkl_RawSparseKernels_hpp
#define __Kokkos_Mkl_RawSparseKernels_hpp

#include "Kokkos_ConfigDefs.hpp"

namespace Kokkos {
  /// \namespace Mkl
  /// \brief Access to Intel Math Kernel Library (MKL) functionality.
  ///
  /// \warning This class is _not_ meant for end users.  Kokkos
  ///   developers use this namespace for access to MKL functionality.
  ///   MKL support is optional for Kokkos, so don't include this
  ///   header file if you didn't build Trilinos with MKL enabled.
  namespace Mkl {
    /// \class RawSparseKernels
    /// \brief Thin wrapper around Intel MKL's sparse kernels.
    ///
    /// \note This class is not meant for end users.  We use it in the
    ///   implementation of MklSparseOps.
    ///
    /// Intel's Math Kernel Library (MKL) implements sparse
    /// matrix-(multi)vector multiply and sparse triangular solve.
    /// The MKL only defines sparse matrix kernels for a limited set
    /// of Scalar types, which correspond to LAPACK's four data types:
    /// float, double, std::complex<float>, and std::complex<double>.
    /// We only give you access to the complex types if Trilinos was
    /// build with complex arithmetic support (i.e., if the
    /// Teuchos_ENABLE_COMPLEX CMake configure option was set to ON).
    ///
    /// MKL only defines sparse kernels for one index (Ordinal) type.
    /// The Ordinal type depends on the MKL library with which you
    /// link.  You can test the Ordinal type by #including mkl.h and
    /// checking whether MKL_ILP64 is #defined.  If it is, then the
    /// only valid Ordinal type is a 64-bit signed integer type
    /// (int64_t in C99 terms; long long in C++ terms); otherwise,
    /// Ordinal is a 32-bit signed integer type (int, or int32_t in
    /// C99 terms).
    ///
    /// \note We have only tested linking these interfaces to MKL
    ///   version 10.3.  If you have trouble linking to future
    ///   versions of MKL, please let us Trilinos developers know and
    ///   we will fix the interfaces.
    ///
    /// \tparam Scalar The type of entries in the sparse matrix and
    ///   dense (multi)vectors.  MKL requires that the entries in the
    ///   matrix and (multi)vectors all have the same type.
    ///
    /// \tparam Ordinal The integer index type (MKL_INT) used by MKL.
    ///   MKL uses the same type for both row offsets and column
    ///   indices.
    ///
    /// \warning MKL only defines the sparse kernels for one Ordinal
    ///   type, which is either a 32-bit signed integer, or a 64-bit
    ///   signed integer.  Which one depends on link options and may
    ///   even be dependent on run-time choices made outside Trilinos.
    ///
    /// \note For zero-based indexing (matdescra[3] == 'C', rather than
    ///   'F'), the multivector inputs and outputs must be in row-major
    ///   order.  This means if you would normally use column-major
    ///   multivectors, then you have to convert ind, ptrBegin, and
    ///   ptrEnd to use one-based indices.
    template<class Scalar, class Ordinal>
    class RawSparseKernels {
    public:
      //! Sparse matrix-vector multiply (one vector): y = beta*y + alpha*A*x.
      static void
      csrmv (const char* const transa,
             const Ordinal m, // Number of rows in A
             const Ordinal k, // Number of columns in A
             const Scalar& alpha,
             const char* const matdescra,
             const Scalar* const val,
             const Ordinal* const ind,
             const Ordinal* const ptrBegin,
             const Ordinal* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
             const Scalar* const x,
             const Scalar& beta,
             Scalar* const y);

      //! Sparse matrix-multivector multiply: C = alpha*A*B + beta*C.
      static void
      csrmm (const char* const transa,
             const Ordinal m, // number of rows of A
             const Ordinal n, // number of columns of C
             const Ordinal k, // number of columns of A
             const Scalar& alpha,
             const char* const matdescra,
             const Scalar* const val,
             const Ordinal* const ind,
             const Ordinal* const ptrBegin,
             const Ordinal* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
             const Scalar* const B,
             const Ordinal LDB,
             const Scalar& beta,
             Scalar* const C,
             const Ordinal LDC);

      //! Sparse triangular solve (one vector): y = alpha*inv(op(A))*x.
      static void
      csrsv (const char* const transa,
             const Ordinal m,
             const Scalar& alpha,
             const char* const matdescra,
             const Scalar* const val,
             const Ordinal* const ind,
             const Ordinal* const ptrBegin,
             const Ordinal* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
             const Scalar* const x,
             Scalar* const y);

      //! Sparse triangular solve (multiple vectors): C = alpha*inv(op(A))*B.
      static void
      csrsm (const char* const transa,
             const Ordinal m, // Number of columns in A
             const Ordinal n, // Number of columns in C
             const Scalar& alpha,
             const char* const matdescra,
             const Scalar* const val,
             const Ordinal* const ind,
             const Ordinal* const ptrBegin,
             const Ordinal* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
             const Scalar* const B,
             const Ordinal LDB,
             Scalar* const C,
             const Ordinal LDC);
    };
  } // namespace Mkl
} // namespace Kokkos

#endif // __Kokkos_Mkl_RawSparseKernels_hpp

