//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
