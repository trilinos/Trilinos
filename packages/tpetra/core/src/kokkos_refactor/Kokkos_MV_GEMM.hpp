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

#ifndef KOKKOS_MV_GEMM_HPP
#define KOKKOS_MV_GEMM_HPP

// Note this code lives only temporarily in TpetraCore.  As soon as
// GEMM kernels exist in the TpetraKernels subpackage, and thus a
// dependency on Teuchos can be eliminated, the code will move to
// TpetraKernels.

#include <Teuchos_BLAS.hpp>
#include <KokkosBlas2_gemv.hpp>
#include "Tpetra_Details_gemm.hpp"

namespace Teuchos {

  // mfh 11 Nov 2014: The DeviceGEMM specializations below need to be
  // able to use Teuchos::BLAS::{GEMM, GEMV}.  We provide just enough
  // of a specialization for Kokkos::complex<{float, double}> to make
  // DeviceGEMM work.  They just defer to BLAS<int,
  // std::complex<{float, double}> > via reinterpret_cast.  Please
  // feel free to expand these specializations if you need to.

  template<>
  class BLAS<int, ::Kokkos::complex<float> > {
  public:
    typedef float mag_type;
    typedef ::Kokkos::complex<float> val_type;
    typedef std::complex<float> impl_type;

    BLAS () {}
    BLAS (const BLAS<int, val_type>&) {}
    virtual ~BLAS () {}

    // void ROTG (val_type* da, val_type* db, mag_type* c, val_type* s) const;
    // void ROT (const int n, val_type* dx, const int incx, val_type* dy, const int incy, RealType* c, val_type* s) const;
    // RealType ASUM (const int n, const val_type* x, const int incx) const;
    //void AXPY (const int n, const val_type alpha, const val_type* x, const int incx, val_type* y, const int incy) const;
    //void COPY (const int n, const val_type* x, const int incx, val_type* y, const int incy) const;
    //val_type DOT(const int n, const val_type* x, const int incx, const val_type* y, const int incy) const;
    //RealType NRM2(const int n, const val_type* x, const int incx) const;
    //void SCAL(const int n, const val_type alpha, val_type* x, const int incx) const;
    //int IAMAX(const int n, const val_type* x, const int incx) const;

    void
    GEMV (ETransp trans, const int m, const int n, const val_type alpha,
          const val_type* A, const int lda, const val_type* x, const int incx,
          const val_type beta, val_type* y, const int incy) const
    {
      BLAS<int, impl_type> blas;
      blas.GEMV (trans, m, n, static_cast<impl_type> (alpha),
                 reinterpret_cast<const impl_type*> (A), lda,
                 reinterpret_cast<const impl_type*> (x), incx,
                 static_cast<impl_type> (beta),
                 reinterpret_cast<impl_type*> (y), incy);
    }

    //void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const val_type* A, const int lda, val_type* x, const int incx) const;
    //void GER(const int m, const int n, const val_type alpha, const val_type* x, const int incx, const val_type* y, const int incy, val_type* A, const int lda) const;

    void
    GEMM (ETransp transa, ETransp transb, const int m, const int n, const int k,
          const val_type alpha, const val_type* A, const int lda,
          const val_type* B, const int ldb, const val_type beta, val_type* C,
          const int ldc) const
    {
      BLAS<int, impl_type> blas;
      blas.GEMM (transa, transb, m, n, k,
                 static_cast<impl_type> (alpha),
                 reinterpret_cast<const impl_type*> (A), lda,
                 reinterpret_cast<const impl_type*> (B), ldb,
                 static_cast<impl_type> (beta),
                 reinterpret_cast<impl_type*> (C), ldc);
    }

    //void SYMM(ESide side, EUplo uplo, const int m, const int n, const val_type alpha, const val_type* A, const int lda, const val_type *B, const int ldb, const val_type beta, val_type *C, const int ldc) const;
    //void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const val_type alpha, const val_type* A, const int lda, const val_type beta, val_type* C, const int ldc) const;
    //void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const val_type alpha, const val_type* A, const int lda, val_type* B, const int ldb) const;
    //void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const val_type alpha, const val_type* A, const int lda, val_type* B, const int ldb) const;
  };

  template<>
  class BLAS<int, ::Kokkos::complex<double> > {
  public:
    typedef double mag_type;
    typedef ::Kokkos::complex<double> val_type;
    typedef std::complex<double> impl_type;

    BLAS () {}
    BLAS (const BLAS<int, val_type>&) {}
    virtual ~BLAS () {}

    // void ROTG (val_type* da, val_type* db, mag_type* c, val_type* s) const;
    // void ROT (const int n, val_type* dx, const int incx, val_type* dy, const int incy, RealType* c, val_type* s) const;
    // RealType ASUM (const int n, const val_type* x, const int incx) const;
    //void AXPY (const int n, const val_type alpha, const val_type* x, const int incx, val_type* y, const int incy) const;
    //void COPY (const int n, const val_type* x, const int incx, val_type* y, const int incy) const;
    //val_type DOT(const int n, const val_type* x, const int incx, const val_type* y, const int incy) const;
    //RealType NRM2(const int n, const val_type* x, const int incx) const;
    //void SCAL(const int n, const val_type alpha, val_type* x, const int incx) const;
    //int IAMAX(const int n, const val_type* x, const int incx) const;

    void
    GEMV (ETransp trans, const int m, const int n, const val_type alpha,
          const val_type* A, const int lda, const val_type* x, const int incx,
          const val_type beta, val_type* y, const int incy) const
    {
      BLAS<int, impl_type> blas;
      blas.GEMV (trans, m, n, static_cast<impl_type> (alpha),
                 reinterpret_cast<const impl_type*> (A), lda,
                 reinterpret_cast<const impl_type*> (x), incx,
                 static_cast<impl_type> (beta),
                 reinterpret_cast<impl_type*> (y), incy);
    }

    //void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const val_type* A, const int lda, val_type* x, const int incx) const;
    //void GER(const int m, const int n, const val_type alpha, const val_type* x, const int incx, const val_type* y, const int incy, val_type* A, const int lda) const;

    void
    GEMM (ETransp transa, ETransp transb, const int m, const int n, const int k,
          const val_type alpha, const val_type* A, const int lda,
          const val_type* B, const int ldb, const val_type beta, val_type* C,
          const int ldc) const
    {
      BLAS<int, impl_type> blas;
      blas.GEMM (transa, transb, m, n, k,
                 static_cast<impl_type> (alpha),
                 reinterpret_cast<const impl_type*> (A), lda,
                 reinterpret_cast<const impl_type*> (B), ldb,
                 static_cast<impl_type> (beta),
                 reinterpret_cast<impl_type*> (C), ldc);
    }

    //void SYMM(ESide side, EUplo uplo, const int m, const int n, const val_type alpha, const val_type* A, const int lda, const val_type *B, const int ldb, const val_type beta, val_type *C, const int ldc) const;
    //void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const val_type alpha, const val_type* A, const int lda, const val_type beta, val_type* C, const int ldc) const;
    //void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const val_type alpha, const val_type* A, const int lda, val_type* B, const int ldb) const;
    //void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const val_type alpha, const val_type* A, const int lda, val_type* B, const int ldb) const;
  };

} // namespace Teuchos


namespace Kokkos {
  namespace Impl {

    template<class ViewType>
    size_t getStride2DView (ViewType A) {
      size_t stride[8];
      A.stride (stride);
      return A.extent (1) > 1 ? stride[1] : A.extent (0);
    }
  }

  /// \struct DeviceGEMM
  /// \brief Class that provides GEMM for a particular Kokkos Device.
  /// \tparam Scalar The type of the entries in the matrices.
  /// \tparam DeviceType The Kokkos Device type.
  ///
  /// GEMM refers to the BLAS' dense matrix-matrix multiply routine.
  template <typename Scalar, typename DeviceType>
  struct DeviceGEMM {
  public:
    static void
    GEMM (const Teuchos::ETransp transA,
          const Teuchos::ETransp transB,
          const Scalar& alpha,
          const View<const Scalar**, LayoutLeft, DeviceType>& A,
          const View<const Scalar**, LayoutLeft, DeviceType>& B,
          const Scalar& beta,
          const View<Scalar**, LayoutLeft, DeviceType>& C)
    {
      const int n = static_cast<int> (C.extent (1));

      // For some BLAS implementations (e.g., MKL), GEMM when B has
      // one column may be signficantly less efficient than GEMV.
      if (n == 1 && transB == Teuchos::NO_TRANS) {
        const int lda = static_cast<int> (Impl::getStride2DView (A));
        Teuchos::BLAS<int,Scalar> blas;
        blas.GEMV (transA, A.extent (0), A.extent (1),
                   alpha, A.data (), lda,
                   B.data (), static_cast<int> (1),
                   beta, C.data (), static_cast<int> (1));
      }
      else {
        const char ctransA = (transA == Teuchos::CONJ_TRANS ? 'C' :
                              (transA == Teuchos::TRANS ? 'T' : 'N'));
        const char ctransB = (transB == Teuchos::CONJ_TRANS ? 'C' :
                              (transB == Teuchos::TRANS ? 'T' : 'N'));
        ::Tpetra::Details::Blas::gemm (ctransA, ctransB, alpha, A, B, beta, C);
      }
    }
  };

  // FIXME (mfh 10 May 2016) Temporary work-around for #243.
  // Don't call MKL for this case.
#ifdef HAVE_KOKKOSKERNELS_MKL
  template <typename DeviceType>
  struct DeviceGEMM<double, DeviceType> {
  public:
    static void
    GEMM (const Teuchos::ETransp transA,
          const Teuchos::ETransp transB,
          const double& alpha,
          const View<const double**, LayoutLeft, DeviceType>& A,
          const View<const double**, LayoutLeft, DeviceType>& B,
          const double& beta,
          const View<double**, LayoutLeft, DeviceType>& C)
    {
      const int n = static_cast<int> (C.extent (1));

      // For some BLAS implementations (e.g., MKL), GEMM when B has
      // one column may be signficantly less efficient than GEMV.
      if (n == 1 && transB == Teuchos::NO_TRANS) {
        char trans = 'N';
        if (transA == Teuchos::TRANS) {
          trans = 'T';
        }
        else if (transA == Teuchos::CONJ_TRANS) {
          trans = 'C';
        }
        auto B_0 = Kokkos::subview (B, Kokkos::ALL (), 0);
        auto C_0 = Kokkos::subview (C, Kokkos::ALL (), 0);
        KokkosBlas::gemv (&trans, alpha, A, B_0, beta, C_0);
      }
      else {
        const char ctransA = (transA == Teuchos::CONJ_TRANS ? 'C' :
                              (transA == Teuchos::TRANS ? 'T' : 'N'));
        const char ctransB = (transB == Teuchos::CONJ_TRANS ? 'C' :
                              (transB == Teuchos::TRANS ? 'T' : 'N'));
        ::Tpetra::Details::Blas::gemm (ctransA, ctransB,
                                       alpha, A, B, beta, C);
      }
    }
  };
#endif // HAVE_KOKKOSKERNELS_MKL

#ifdef KOKKOS_ENABLE_CUDA
  template <typename Scalar>
  struct DeviceGEMM<Scalar, Cuda> {
  public:
    static void
    GEMM (const Teuchos::ETransp transA,
          const Teuchos::ETransp transB,
          const Scalar& alpha,
          const View<const Scalar**, LayoutLeft, Cuda>& A,
          const View<const Scalar**,LayoutLeft, Cuda>& B,
          const Scalar& beta,
          const View<Scalar**,LayoutLeft,Cuda>& C)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "DeviceGEMM: Kokkos::Cuda has no support "
         "for GEMM operations over Scalar=" << Teuchos::typeName(alpha) << ".");
    }
  };

  template <>
  struct DeviceGEMM<float,Cuda> {
    public:
      static void
      GEMM (const Teuchos::ETransp transA,
            const Teuchos::ETransp transB,
            const float alpha,
            const View<const float**,LayoutLeft,Cuda>& A,
            const View<const float**,LayoutLeft,Cuda>& B,
            const float beta,
            const View<float**,LayoutLeft,Cuda>& C)
    {
      const char ctransA = (transA == Teuchos::NO_TRANS ? 'N' : 'T');
      const char ctransB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');

      ::Tpetra::Details::Blas::gemm (ctransA, ctransB,
                                     alpha, A, B, beta, C);
    }
  };

  template <>
  struct DeviceGEMM<double,Cuda> {
  public:
    static void
    GEMM (const Teuchos::ETransp transA,
          const Teuchos::ETransp transB,
          const double alpha,
          const View<const double**, LayoutLeft, Cuda>& A,
          const View<const double**, LayoutLeft, Cuda>& B,
          const double beta,
          const View<double**, LayoutLeft, Cuda>& C)
    {
      const char ctransA = (transA == Teuchos::NO_TRANS ? 'N' : 'T');
      const char ctransB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');

      ::Tpetra::Details::Blas::gemm (ctransA, ctransB,
                                     alpha, A, B, beta, C);
    }
  };
#endif // KOKKOS_ENABLE_CUDA

} // namespace Kokkos
#endif // KOKKOS_MV_GEMM_HPP

