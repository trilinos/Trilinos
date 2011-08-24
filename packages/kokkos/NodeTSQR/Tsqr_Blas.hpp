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
