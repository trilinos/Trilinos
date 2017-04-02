/*
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
*/

#include "Tpetra_Details_libGemm.hpp"
#include "KokkosKernels_config.h"
// For now, just use Teuchos, whatever.  Later, wrap it ourselves.
// This is way less trivial than you might think, given all the
// variations in link options and the platforms we must support.
#include "Teuchos_BLAS.hpp"

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Lib {
namespace Impl {

void
cgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<float>& alpha,
       const ::Kokkos::complex<float> A[],
       const int lda,
       const ::Kokkos::complex<float> B[],
       const int ldb,
       const ::Kokkos::complex<float>& beta,
       ::Kokkos::complex<float> C[],
       const int ldc)
{
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  typedef std::complex<float> IST; // impl_scalar_type;
  typedef Teuchos::BLAS<int, IST> blas_type;

  const ::Teuchos::ETransp transA_ =
    (transA == 'T' || transA == 't') ? ::Teuchos::TRANS :
    ((transA == 'C' || transA == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const ::Teuchos::ETransp transB_ =
    (transB == 'T' || transB == 't') ? ::Teuchos::TRANS :
    ((transB == 'C' || transB == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const IST alpha_ (alpha.real (), alpha.imag ());
  const IST beta_ (beta.real (), beta.imag ());

  blas_type blas;
  blas.GEMM (transA_, transB_, m, n, k,
             alpha_, reinterpret_cast<const IST*> (A), lda,
             reinterpret_cast<const IST*> (B), ldb,
             beta_, reinterpret_cast<IST*> (C), ldc);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::cgemm: "
                            "You must configure Teuchos and Tpetra with complex support.");
#endif
}

void
dgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const double alpha,
       const double A[],
       const int lda,
       const double B[],
       const int ldb,
       const double beta,
       double C[],
       const int ldc)
{
  typedef double IST; // impl_scalar_type;
  typedef Teuchos::BLAS<int, IST> blas_type;

  const ::Teuchos::ETransp transA_ =
    (transA == 'T' || transA == 't') ? ::Teuchos::TRANS :
    ((transA == 'C' || transA == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const ::Teuchos::ETransp transB_ =
    (transB == 'T' || transB == 't') ? ::Teuchos::TRANS :
    ((transB == 'C' || transB == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);

  blas_type blas;
  blas.GEMM (transA_, transB_, m, n, k,
             alpha, A, lda,
             B, ldb,
             beta, C, ldc);
}

void
sgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const float alpha,
       const float A[],
       const int lda,
       const float B[],
       const int ldb,
       const float beta,
       float C[],
       const int ldc)
{
  typedef float IST; // impl_scalar_type;
  typedef Teuchos::BLAS<int, IST> blas_type;

  const ::Teuchos::ETransp transA_ =
    (transA == 'T' || transA == 't') ? ::Teuchos::TRANS :
    ((transA == 'C' || transA == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const ::Teuchos::ETransp transB_ =
    (transB == 'T' || transB == 't') ? ::Teuchos::TRANS :
    ((transB == 'C' || transB == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);

  blas_type blas;
  blas.GEMM (transA_, transB_, m, n, k,
             alpha, A, lda,
             B, ldb,
             beta, C, ldc);
}

void
zgemm (const char transA,
       const char transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<double>& alpha,
       const ::Kokkos::complex<double> A[],
       const int lda,
       const ::Kokkos::complex<double> B[],
       const int ldb,
       const ::Kokkos::complex<double>& beta,
       ::Kokkos::complex<double> C[],
       const int ldc)
{
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  typedef std::complex<double> IST; // impl_scalar_type;
  typedef Teuchos::BLAS<int, IST> blas_type;

  const ::Teuchos::ETransp transA_ =
    (transA == 'T' || transA == 't') ? ::Teuchos::TRANS :
    ((transA == 'C' || transA == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const ::Teuchos::ETransp transB_ =
    (transB == 'T' || transB == 't') ? ::Teuchos::TRANS :
    ((transB == 'C' || transB == 'c') ? ::Teuchos::CONJ_TRANS :
     ::Teuchos::NO_TRANS);
  const IST alpha_ (alpha.real (), alpha.imag ());
  const IST beta_ (beta.real (), beta.imag ());

  blas_type blas;
  blas.GEMM (transA_, transB_, m, n, k,
             alpha_, reinterpret_cast<const IST*> (A), lda,
             reinterpret_cast<const IST*> (B), ldb,
             beta_, reinterpret_cast<IST*> (C), ldc);
#else
  throw std::runtime_error ("Tpetra::Details::Blas::Lib::Impl::cgemm: "
                            "You must configure Teuchos and Tpetra with complex support.");
#endif
}

} // namespace Impl
} // namespace Lib
} // namespace Blas
} // namespace Details
} // namespace Tpetra
