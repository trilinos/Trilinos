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

#ifndef __TSQR_Tsqr_Lapack_hpp
#define __TSQR_Tsqr_Lapack_hpp

#include <Tsqr_ConfigDefs.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace TSQR {

  template<class Ordinal, class Scalar>
  class LAPACK {
  public:
    typedef Ordinal ordinal_type;
    typedef Scalar scalar_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    /// The type of the absolute value (or magnitude, if Scalar is
    /// complex) of a Scalar.
    typedef typename STS::magnitudeType magnitude_type;

    LAPACK () {}

    /// Whether or not the QR factorizations computed by LAPACK::GEQRF()
    /// and LAPACK::GEQR2() produce an R factor with all nonnegative
    /// diagonal entries.  This also corresponds to whether
    /// LAPACK::LARFP() always produces a nonnegative BETA output, and
    /// therefore whether the QR factorizations in the TSQR::Combine
    /// class produce R factors with all negative diagonal entries.
    static bool QR_produces_R_factor_with_nonnegative_diagonal();

    /// If the LAPACK library has _LARFGP, calls that.  Else, if the
    /// LAPACK library has _LARFP, calls that.  Otherwise, calls
    /// _LARF.  The last choice means that the alpha output may be
    /// negative if Scalar is real.
    void 
    LARFP (const Ordinal n, 
	   Scalar& alpha, 
	   Scalar x[], 
	   const Ordinal incx, 
	   Scalar& tau);

    /// If the LAPACK library has _GEQRFP, calls that.  Otherwise,
    /// calls _GEQRF.  _GEQRFP always computes an R factor with
    /// nonnegative diagonal entries.  _GEQRF does this in LAPACK 3.2
    /// and 3.2.1, but not in LAPACK <= 3.1.1 or LAPACK >= 3.2.2.
    void
    GEQRF  (const Ordinal m,
	    const Ordinal n, 
	    Scalar A[],
	    const Ordinal lda,
	    Scalar tau[],
	    Scalar work[],
	    const int lwork,
	    int* const INFO);

    /// If the LAPACK library has _GEQR2P, calls that.  Otherwise,
    /// calls _GEQR2.  _GEQR2P always computes an R factor with
    /// nonnegative diagonal entries.  _GEQR2 does this in LAPACK 3.2
    /// and 3.2.1, but not in LAPACK <= 3.1.1 or LAPACK >= 3.2.2.
    void 
    GEQR2 (const Ordinal m, 
	   const Ordinal n, 
	   Scalar A[],
	   const Ordinal lda, 
	   Scalar tau[],
	   Scalar work[],
	   int* const INFO);

    void
    ORM2R (const char* const side,
	   const char* const trans,
	   const Ordinal m,
	   const Ordinal n,
	   const Ordinal k,
	   const Scalar A[],
	   const Ordinal lda,
	   const Scalar tau[],
	   Scalar C[],
	   const Ordinal ldc,
	   Scalar work[],
	   int* const info);

    void
    ORMQR (const char* const side,
	   const char* const trans,
	   const Ordinal m,
	   const Ordinal n,
	   const Ordinal k,
	   const Scalar A[],
	   const Ordinal lda,
	   const Scalar tau[],
	   Scalar C[],
	   const Ordinal ldc,
	   Scalar work[],
	   const int lwork,
	   int* const INFO);

    void
    ORGQR (const Ordinal m,
	   const Ordinal n,
	   const Ordinal k,
	   Scalar A[],
	   const Ordinal lda,
	   Scalar tau[],
	   Scalar work[],
	   const int lwork,
	   int* const INFO);

    void
    POTRF (const char* const uplo,
	   const Ordinal n,
	   Scalar A[],
	   const Ordinal lda,
	   int* const INFO);

    void
    POTRS (const char* const uplo,
	   const Ordinal n,
	   const Ordinal nrhs,
	   const Scalar A[],
	   const Ordinal lda,
	   Scalar B[],
	   const Ordinal ldb,
	   int* const INFO);

    void
    POTRI (const char* const uplo, 
	   const Ordinal n, 
	   Scalar A[], 
	   const Ordinal lda, 
	   int* const INFO);

    void
    LARNV (const int idist, 
	   int iseed[],
	   const Ordinal n,
	   Scalar x[]);

    void 
    GESVD (const char* const jobu,
	   const char* const jobvt,
	   const Ordinal m,
	   const Ordinal n,
	   Scalar A[],
	   const Ordinal lda,
	   magnitude_type s[],
	   Scalar U[],
	   const Ordinal ldu,
	   Scalar VT[],
	   const Ordinal ldvt,
	   Scalar work[],
	   const Ordinal lwork,
	   magnitude_type rwork[],
	   int* const INFO);

  private:
    LAPACK (const LAPACK&);
    LAPACK& operator= (const LAPACK&);
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_Lapack_hpp
