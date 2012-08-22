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

#ifndef __TSQR_TsqrBlas_hpp
#define __TSQR_TsqrBlas_hpp

#include <Tsqr_ConfigDefs.hpp>

namespace TSQR {

  /// \class BLAS
  /// \brief Wrappers for BLAS routines used by the Tall Skinny QR factorization.
  ///
  template<class Ordinal, class Scalar>
  class BLAS {
  public:
    BLAS () {}

    void 
    GEMV (const char* const trans, 
	  const Ordinal m, 
	  const Ordinal n,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  const Scalar x[],
	  const Ordinal incx,
	  const Scalar beta,
	  Scalar y[],
	  const Ordinal incy);

    void
    GEMM (const char* const transa,
	  const char* const transb,
	  const Ordinal m,
	  const Ordinal n,
	  const Ordinal k,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  const Scalar B[],
	  const Ordinal ldb,
	  const Scalar beta,
	  Scalar C[],
	  const Ordinal ldc);

    ///
    /// If ScalarTraits< Scalar >::is_complex, calls _GERC.
    /// Otherwise, calls _GER.
    void
    GER (const Ordinal m,
	 const Ordinal n,
	 const Scalar alpha,
	 const Scalar x[],
	 const Ordinal incx,
	 const Scalar y[],
	 const Ordinal incy,
	 Scalar A[],
	 const Ordinal lda);

    void
    TRSM (const char* const side,
	  const char* const uplo,
	  const char* const transa,
	  const char* const diag,
	  const Ordinal m,
	  const Ordinal n,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  Scalar B[],
	  const Ordinal ldb);

  private:
    BLAS (const BLAS&);
    BLAS& operator= (const BLAS&);
  };

} // namespace TSQR

#endif // __TSQR_TsqrBlas_hpp
